#!/bin/bash
#SBATCH --job-name=modern_bam2fa_invariant
#SBATCH --partition=gpu
#SBATCH --gres=0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=45
#SBATCH --mem=250G
#SBATCH --time=72:00:00
#SBATCH --output=modern_bam2fa_%j.out
#SBATCH --error=modern_bam2fa_%j.err

#set -euo pipefail

# ============================================================
# USER SETTINGS
# ============================================================
BAMLIST="psitacara.bamlist"
REF="GCA_045345405.1_bGuaGua1_haplotype_1_genomic.fna"
OUTDIR="modern_consensus_fasta_not_patch"

TOTAL_CPUS=45
N_JOBS=7
THREADS=$(( TOTAL_CPUS / N_JOBS ))

MIN_MQ=30
MIN_BQ=30
MIN_DP=3

# ============================================================
# ENVIRONMENT
# ============================================================
source ~/.bashrc
module load BCFtools/1.21-GCC-13.3.0
# conda activate your_env_if_needed

mkdir -p ${OUTDIR}/{vcf,mask,fasta,logs}

# Index reference if needed
[ -f ${REF}.fai ] || samtools faidx ${REF}

# ============================================================
# FUNCTION
# ============================================================
process_one() {

    BAM=$1
    SAMPLE=$(basename ${BAM} .bam)

    echo "[${SAMPLE}] START $(date)" >&2

    # Index BAM if missing
    [ -f ${BAM}.bai ] || samtools index ${BAM}

    # ============================================================
    # 1) BAM → ALL-SITES VCF (variants + invariants)
    # ============================================================
    bcftools mpileup \
        -Ou \
        -f ${REF} \
        -q ${MIN_MQ} \
        -Q ${MIN_BQ} \
        -d 10000 \
        -a FORMAT/DP \
        --threads ${THREADS} \
        ${BAM} \
    | bcftools call \
        -c -A -Ou \
    | bcftools view \
        -Oz -o ${OUTDIR}/vcf/${SAMPLE}.allsites.raw.vcf.gz

    bcftools index -t ${OUTDIR}/vcf/${SAMPLE}.allsites.raw.vcf.gz

    # ============================================================
    # 2) REMOVE INDELS (keep SNPs + invariant sites)
    # ============================================================
    bcftools view \
        -V indels \
        -Oz -o ${OUTDIR}/vcf/${SAMPLE}.allsites.noindels.vcf.gz \
        ${OUTDIR}/vcf/${SAMPLE}.allsites.raw.vcf.gz

    bcftools index -t ${OUTDIR}/vcf/${SAMPLE}.allsites.noindels.vcf.gz

    # ============================================================
    # 3) MASK ONLY LOW-COVERAGE SITES (DP < MIN_DP → ./.)
    #    HETEROZYGOUS SITES ARE KEPT
    # ============================================================
    bcftools +setGT \
        ${OUTDIR}/vcf/${SAMPLE}.allsites.noindels.vcf.gz \
        -Oz -o ${OUTDIR}/vcf/${SAMPLE}.allsites.masked.vcf.gz -- \
        -t q -n . -i "FORMAT/DP<${MIN_DP}"

    bcftools index -t ${OUTDIR}/vcf/${SAMPLE}.allsites.masked.vcf.gz

    # ============================================================
    # 4) ZERO-COVERAGE BED MASK (true no-data regions)
    # ============================================================
    bedtools genomecov -ibam ${BAM} -bga \
    | awk '$4==0 {print $1"\t"$2"\t"$3}' \
    > ${OUTDIR}/mask/${SAMPLE}.zero.bed

    # ============================================================
    # 5) CONSENSUS (DIPLOID, INVARIANTS KEPT, NO PATCHING)
    # ============================================================
    bcftools consensus \
        -f ${REF} \
        --missing N \
        -m ${OUTDIR}/mask/${SAMPLE}.zero.bed \
        ${OUTDIR}/vcf/${SAMPLE}.allsites.masked.vcf.gz \
    > ${OUTDIR}/fasta/${SAMPLE}.fa

    echo "[${SAMPLE}] END   $(date)" >&2
}

export -f process_one
export REF OUTDIR THREADS MIN_MQ MIN_BQ MIN_DP

# ============================================================
# RUN PARALLEL
# ============================================================
grep -v '^#' ${BAMLIST} | grep -v '^\s*$' \
| parallel -j ${N_JOBS} process_one {}

# ============================================================
# QC SUMMARY
# ============================================================
seqkit stats ${OUTDIR}/fasta/*.fa > ${OUTDIR}/fasta_stats.tsv

echo "DONE"
echo "FASTAs written to: ${OUTDIR}/fasta/"


