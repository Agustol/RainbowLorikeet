
#!/bin/bash
#SBATCH --job-name=merge_snps_miss60_density
#SBATCH --output=merge_snps_miss60_density.%j.out
#SBATCH --error=merge_snps_miss60_density.%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=batch

set -euo pipefail

# ==========================
# USER VARIABLES
# ==========================
VCF_DIR="vcf"                          # per-sample *.snps.vcf.gz
REF="reference.fa"                     # reference FASTA
OUTDIR="merged_vcf"
PREFIX="merged.snps.miss60"
THREADS=8
MAX_MISS=0.60
WINDOW_BP=10000                        # 10 kb windows

mkdir -p "$OUTDIR"

# ==========================
# LOAD MODULES
# ==========================
module load bcftools/1.18 || module load bcftools
module load bedtools || true
module load samtools || true

# ==========================
# CHECK INPUTS
# ==========================
[[ -f "$REF" ]] || { echo "ERROR: reference not found"; exit 1; }

# ==========================
# INDEX REFERENCE (if needed)
# ==========================
if [[ ! -f "${REF}.fai" ]]; then
  samtools faidx "$REF"
fi

# ==========================
# LIST INPUT SNP VCFs
# ==========================
ls "$VCF_DIR"/*.snps.vcf.gz > vcf.list
NVCF=$(wc -l < vcf.list)

echo "Merging $NVCF SNP VCFs"

# ==========================
# STEP 1 — MERGE SNP VCFs
# ==========================
bcftools merge \
  --threads "$THREADS" \
  -m none \
  -l vcf.list \
  -Oz \
  -o "$OUTDIR/$PREFIX.raw.vcf.gz"

bcftools index -f "$OUTDIR/$PREFIX.raw.vcf.gz"

# ==========================
# STEP 2 — FILTER BY MISSINGNESS
# ==========================
bcftools view \
  --threads "$THREADS" \
  -i "F_MISSING <= $MAX_MISS" \
  -v snps \
  -Oz \
  -o "$OUTDIR/$PREFIX.vcf.gz" \
  "$OUTDIR/$PREFIX.raw.vcf.gz"

bcftools index -f "$OUTDIR/$PREFIX.vcf.gz"

# ==========================
# STEP 3 — BUILD 10 kb WINDOWS
# ==========================
bedtools makewindows \
  -g "${REF}.fai" \
  -w "$WINDOW_BP" \
> "$OUTDIR/genome_10kb.bed"

# ==========================
# STEP 4 — SNP DENSITY PER 10 kb
# ==========================
bcftools query -f '%CHROM\t%POS\n' \
  "$OUTDIR/$PREFIX.vcf.gz" \
| awk '{print $1"\t"$2-1"\t"$2}' \
| bedtools intersect \
    -a "$OUTDIR/genome_10kb.bed" \
    -b - \
    -c \
> "$OUTDIR/${PREFIX}.snp_density_10kb.bed"

# ==========================
# QC SUMMARY
# ==========================
echo "Final SNP count:"
bcftools view -H "$OUTDIR/$PREFIX.vcf.gz" | wc -l

echo "Example SNP density (first 5 windows):"
head -5 "$OUTDIR/${PREFIX}.snp_density_10kb.bed"

echo "DONE"

