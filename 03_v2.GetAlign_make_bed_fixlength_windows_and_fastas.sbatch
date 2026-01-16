#!/bin/bash
#SBATCH --job-name=windows_from_snps_to_fastas
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=windows_from_snps_to_fastas.%j.out
#SBATCH --error=windows_from_snps_to_fastas.%j.err

set -euo pipefail

# ==========================
# USER SETTINGS
# ==========================
VCF="merged.snps.miss60.vcf.gz"
REF="reference.fa"

WINDOW_SIZE=10000

FASTADIR="fasta"
OUT_FASTA_DIR="windows_per_sample"
WINDOWS="windows_${WINDOW_SIZE}bp.bed"

mkdir -p "$OUT_FASTA_DIR"

# ==========================
# CHECK INPUTS
# ==========================
[[ -f "$VCF" ]] || { echo "ERROR: VCF not found"; exit 1; }
[[ -f "$REF" ]] || { echo "ERROR: reference FASTA not found"; exit 1; }

# ==========================
# INDEX REFERENCE (FOR FASTA EXTRACTION ONLY)
# ==========================
if [[ ! -f "${REF}.fai" ]]; then
  samtools faidx "$REF"
fi

# ==========================
# STEP 1 — DERIVE CHR EXTENTS FROM VCF (NOT REF)
# ==========================
echo "Deriving chromosome extents from VCF"

bcftools view -H "$VCF" \
| awk '{print $1, $2}' \
| sort -k1,1 -k2,2n \
| awk '
  $1 != prev {
    if (NR > 1) print prev, max
    prev = $1
    max = $2
  }
  { max = $2 }
  END { print prev, max }
' > chr.lengths.from_vcf.txt

# ==========================
# STEP 2 — BUILD FIXED 10 kb WINDOWS (CHR-SAFE)
# ==========================
bedtools makewindows \
  -g chr.lengths.from_vcf.txt \
  -w "$WINDOW_SIZE" \
> "$WINDOWS"

echo "Generated fixed windows: $WINDOWS"

# ==========================
# STEP 3 — EXTRACT WINDOWS PER SAMPLE
# ==========================
for FA in "$FASTADIR"/*.fa; do
  SAMPLE=$(basename "$FA" .fa)

  echo "Extracting windows for $SAMPLE"

  bedtools getfasta \
    -fi "$FA" \
    -bed "$WINDOWS" \
    -name \
    -fo "$OUT_FASTA_DIR/${SAMPLE}.windows.fa"
done
