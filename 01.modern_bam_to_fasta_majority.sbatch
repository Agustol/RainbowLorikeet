#!/bin/bash
#SBATCH --job-name=modern_bam2fa
#SBATCH --partition=gpu
#SBATCH --gres=none
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --output=modern_bam2fa_%j.out
#SBATCH --error=modern_bam2fa_%j.err

# -----------------------
# USER SETTINGS
# -----------------------
BAMLIST="psitacara.bamlist"
REF="GCA_045345405.1_bGuaGua1_haplotype_1_genomic.fna"

OUTDIR="modern_consensus_fasta"

# TOTAL CPUS AVAILABLE
TOTAL_CPUS=50

# NUMBER OF PARALLEL SAMPLES
N_JOBS=5              # 5 samples at once

# THREADS PER SAMPLE (50 / 5 = 10)
THREADS=$(( TOTAL_CPUS / N_JOBS ))

# FILTERS (modern samples)
MIN_MQ=30
MIN_BQ=30
MIN_DP=3

# -----------------------
# ENVIRONMENT
# -----------------------
source ~/.bashrc
# conda activate your_env_if_needed

mkdir -p ${OUTDIR}/{vcf,mask,fasta,logs}

# Index reference if needed
if [ ! -f ${REF}.fai ]; then
    samtools faidx ${REF}
fi

# -----------------------
# FUNCTION
# -----------------------
process_one() {

    BAM=$1
    SAMPLE=$(basename ${BAM} .bam)

    # Index BAM if missing
    if [ ! -f ${BAM}.bai ] && [ ! -f ${BAM%.bam}.bai ]; then
        samtools index ${BAM}
    fi

    # 1) BAM -> majority-rule VCF
    bcftools mpileup \
        -Ou \
        -f ${REF} \
        -q ${MIN_MQ} \
        -Q ${MIN_BQ} \
        -d 10000 \
        -a FORMAT/DP \
        --threads ${THREADS} \
        ${BAM} \
    | bcftools call -c -Ou \
    | bcftools view -Oz -o ${OUTDIR}/vcf/${SAMPLE}.raw.vcf.gz

    bcftools index -f ${OUTDIR}/vcf/${SAMPLE}.raw.vcf.gz

    # 2) Keep only biallelic SNPs
    bcftools view \
        -v snps -m2 -M2 \
        -Oz -o ${OUTDIR}/vcf/${SAMPLE}.snps.vcf.gz \
        ${OUTDIR}/vcf/${SAMPLE}.raw.vcf.gz

    bcftools index -f ${OUTDIR}/vcf/${SAMPLE}.snps.vcf.gz

    # 3) Zero coverage mask
    bedtools genomecov -ibam ${BAM} -bga \
    | awk '$4==0 {print $1"\t"$2"\t"$3}' \
    > ${OUTDIR}/mask/${SAMPLE}.zero.bed

    # 4) Low DP mask
    bcftools query -f '%CHROM\t%POS0\t%POS[\t%DP]\n' \
        ${OUTDIR}/vcf/${SAMPLE}.snps.vcf.gz \
    | awk -v mdp=${MIN_DP} '$4 < mdp {print $1"\t"$2"\t"$3}' \
    > ${OUTDIR}/mask/${SAMPLE}.lowdp.bed

    # 5) Heterozygous mask
    bcftools query -f '%CHROM\t%POS0\t%POS[\t%GT]\n' \
        ${OUTDIR}/vcf/${SAMPLE}.snps.vcf.gz \
    | awk '$4=="0/1" {print $1"\t"$2"\t"$3}' \
    > ${OUTDIR}/mask/${SAMPLE}.het.bed

    # 6) Merge masks
    cat \
        ${OUTDIR}/mask/${SAMPLE}.zero.bed \
        ${OUTDIR}/mask/${SAMPLE}.lowdp.bed \
        ${OUTDIR}/mask/${SAMPLE}.het.bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > ${OUTDIR}/mask/${SAMPLE}.final.bed

    # 7) FASTA consensus
    bcftools consensus \
        -f ${REF} \
        -m ${OUTDIR}/mask/${SAMPLE}.final.bed \
        ${OUTDIR}/vcf/${SAMPLE}.snps.vcf.gz \
    > ${OUTDIR}/fasta/${SAMPLE}.fa
}

export -f process_one
export REF OUTDIR THREADS MIN_MQ MIN_BQ MIN_DP

# -----------------------
# RUN PARALLEL
# -----------------------
grep -v '^#' ${BAMLIST} | grep -v '^\s*$' \
| parallel -j ${N_JOBS} process_one {}

# -----------------------
# QC SUMMARY
# -----------------------
seqkit stats ${OUTDIR}/fasta/*.fa > ${OUTDIR}/fasta_stats.tsv

echo "DONE"
echo "FASTAs: ${OUTDIR}/fasta/"
