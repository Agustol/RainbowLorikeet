#!/bin/bash
#SBATCH --job-name=windows_to_alignments
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=49
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --gres=0
#SBATCH --time=300:00:00
#SBATCH --output=windows_to_alignments.%j.out
#SBATCH --error=windows_to_alignments.%j.err

set -euo pipefail

INDIR="windows_per_sample"
TMPDIR="tmp_by_sample"
OUTDIR="alignments"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$TMPDIR" "$OUTDIR"
rm -rf "$TMPDIR"/*
rm -f "$OUTDIR"/*.fa

command -v gawk >/dev/null || { echo "ERROR: gawk not found"; exit 1; }
command -v parallel >/dev/null || { echo "ERROR: GNU parallel not found"; exit 1; }

# -------------------------------------------------
# WRITE THE AWK SCRIPT (CORRECT, SINGLE-LINE STRINGS)
# -------------------------------------------------
cat > extract_windows.awk <<'AWK'
function flush() {
  if (have && seq != "") {
    frag[chr][start] = frag[chr][start] seq
  }
  seq=""
}

BEGIN { have=0; seq="" }

/^>/ {
  flush()
  h = substr($0,2)
  pos = index(h,"::")
  if (pos==0) { have=0; next }

  rest = substr(h,pos+2)
  split(rest,a,":")
  chr=a[1]
  split(a[2],b,"-")
  start=b[1]+0
  end=b[2]+0
  coord[chr][start]=end
  have=1
  next
}

{
  if (have) {
    gsub(/[ \t\r\n]/,"",$0)
    seq=seq $0
  }
}

END {
  flush()
  for (c in frag) {
    n=asorti(frag[c],idx,"@ind_num_asc")
    for (i=1;i<=n;i++) {
      s=idx[i]
      e=coord[c][s]
      w=i
      fname=outdir "/window" w "_" c "-" s "-" e ".fa"
      print ">" sample "\n" frag[c][s] > fname
    }
  }
}
AWK

# ----------------------------
# PER-SAMPLE EXTRACTION
# ----------------------------
find "$INDIR" -name "*.windows.fa" | sort | \
parallel -j "$THREADS" '
  FA={}
  SAMPLE=$(basename "$FA" .windows.fa)
  SDIR="'"$TMPDIR"'/$SAMPLE"
  mkdir -p "$SDIR"
  gawk -v sample="$SAMPLE" -v outdir="$SDIR" -f extract_windows.awk "$FA"
'

# ----------------------------
# MERGE WINDOWS (SAFE, SERIAL)
# ----------------------------
find "$TMPDIR" -mindepth 2 -name "*.fa" -printf "%f\n" | sort -u > keys.txt
find "$TMPDIR" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > samples.txt

while read -r key; do
  out="$OUTDIR/$key"
  : > "$out"
  while read -r s; do
    f="$TMPDIR/$s/$key"
    [ -f "$f" ] && cat "$f" >> "$out"
  done < samples.txt
done < keys.txt

echo "DONE. Windows renumbered per scaffold."

