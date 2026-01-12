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
NSNP=50

FASTADIR="fasta"
OUT_FASTA_DIR="windows_per_sample"
WINDOWS="windows_${NSNP}snps.bed"

mkdir -p "$OUT_FASTA_DIR"

# ==========================
# CHECK INPUTS
# ==========================
[[ -f "$VCF" ]] || { echo "ERROR: VCF not found"; exit 1; }
[[ -f "$REF" ]] || { echo "ERROR: reference FASTA not found"; exit 1; }

# ==========================
# INDEX REFERENCE (IF NEEDED)
# ==========================
if [[ ! -f "${REF}.fai" ]]; then
  samtools faidx "$REF"
fi

# ==========================
# LOAD CHROM LENGTHS
# ==========================
awk '{print $1, $2}' "${REF}.fai" > chr.lengths

# ==========================
# STEP 1 — BUILD SNP-COUNT WINDOWS (BOUNDARY-SAFE)
# ==========================
bcftools query -f '%CHROM\t%POS\n' "$VCF" \
| sort -k1,1 -k2,2n \
| awk -v NSNP="$NSNP" '
BEGIN {
  while ((getline < "chr.lengths") > 0) {
    clen[$1] = $2
  }
}
{
  chr = $1
  pos = $2

  if (chr != prev_chr) {
    i = 0
    w = 1
  }

  i++
  if (i == 1) start = pos - 1
  end = pos

  if (i == NSNP) {
    if (end <= clen[chr]) {
      printf "%s\t%d\t%d\twindow_%d\n", chr, start, end, w
    }
    w++
    i = 0
  }

  prev_chr = chr
}
END {
  if (i > 0 && end <= clen[chr]) {
    printf "%s\t%d\t%d\twindow_%d\n", chr, start, end, w
  }
}
' > "$WINDOWS"

echo "Generated SNP-count windows: $WINDOWS"

# ==========================
# STEP 2 — EXTRACT WINDOWS PER SAMPLE
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

echo "DONE"

