#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orig_ptr
#SBATCH --cpus-per-task=24
#SBATCH --mem 100G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-59%10
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/pytransrate_orig_orp/logs/ptr_%A_%a.log

# pyTransRate on each ORIGINAL assembly of the lineage_410 samples, called
# exactly as oyster.py's pytransrate step scores the ORP assembly, so the
# original and ORP sides of the comparison (and the components, scored by
# assembly_eval.py the same way) differ only in the assembly:
#
#   conda run -n orp pytransrate -o <out> -a <fasta> \
#       --left <cor1> --right <cor2> -t <cpu> --max-memory <mem>G
#
# * the orp env's pyTransRate (2.2.2, ORP's PYTRANSRATE_PINNED_VERSION);
# * the reads are that sample's ORP run's own trimmed, error-corrected pair,
#   lineage_410/<SRR>/rcorr/<runout>.TRIM_{1,2}P.cor.fq.gz.
#
# This replaces two earlier original-side scores that did not match ORP's:
# compare/old_results/ (ORP 3's Ruby TransRate, raw reads) and
# compare/pytransrate/results/ (pyTransRate from the transrate2026 env, mostly
# from before its 2.2.2.dev1 install on 2026-09-18, raw reads).
#
# Contig names are normalised the way run_pytransrate.sbatch did
# (ENA|ACC|ACC.V -> ACC, other '|' -> '_'), so contigs.csv rows keep the same
# ids as before; the task fails on duplicate names rather than score them.
#
#   mkdir -p /mnt/home/macmaneslab/macmanes/compare/pytransrate_orig_orp/logs
#   sbatch ~/Oyster_River_Protocol/experiments/compare108/orig_pytransrate.sbatch
#
# Writes compare/pytransrate_orig_orp/<pair>/assemblies.csv (+ contigs.csv),
# the layout collect_lineage410.py --pytransrate-original reads. A pair with
# an assemblies.csv already is skipped.

set -euo pipefail

MANIFEST="${MANIFEST:-/mnt/home/macmaneslab/macmanes/compare/lineage_410/manifest.tsv}"
RUNS="${RUNS:-/mnt/home/macmaneslab/macmanes/compare/lineage_410}"
OUT="${OUT:-/mnt/home/macmaneslab/macmanes/compare/pytransrate_orig_orp}"
ORPDIR="${ORPDIR:-$HOME/Oyster_River_Protocol}"
CPU="${SLURM_CPUS_PER_TASK:-24}"
MEM=$(( ${SLURM_MEM_PER_NODE:-102400} / 1024 ))

echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="

line=$(grep . "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$line" ] || { echo "manifest has no line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r PAIR ASM _ _ SRR _ <<< "$line"
[ -f "$ASM" ] || { echo "missing original assembly for $PAIR: $ASM" >&2; exit 1; }

# The run name, from the run's finished qualreport (v4_1_0 or v4_0_0).
shopt -s nullglob
done_files=("$RUNS/$SRR"/reports/qualreport."$SRR"_v*.done)
shopt -u nullglob
[ "${#done_files[@]}" -eq 1 ] || {
    echo "expected one finished run in $RUNS/$SRR/reports, found ${#done_files[@]}" >&2; exit 1; }
RUN=$(basename "${done_files[0]}" .done); RUN="${RUN#qualreport.}"
C1="$RUNS/$SRR/rcorr/$RUN.TRIM_1P.cor.fq.gz"
C2="$RUNS/$SRR/rcorr/$RUN.TRIM_2P.cor.fq.gz"
for f in "$C1" "$C2"; do
    [ -s "$f" ] || { echo "missing corrected reads: $f" >&2; exit 1; }
done

DEST="$OUT/$PAIR"
mkdir -p "$OUT/logs"
ln -sfn "ptr_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$OUT/logs/$PAIR.log"

echo "    pair     $PAIR ($SRR, run $RUN)"
echo "    assembly $ASM"
echo "    reads    $C1"
echo "             $C2"
echo "    out      $DEST"
echo "    ORP      $(cat "$ORPDIR/version.txt") ($(git -C "$ORPDIR" rev-parse --short HEAD 2>/dev/null || echo '?'))"

if [ -s "$DEST/assemblies.csv" ]; then
    echo "already complete, skipping"
    exit 0
fi

module load anaconda/colsa
echo "    $(conda run -n orp pytransrate --version 2>&1 | head -1), ${CPU} cpus, ${MEM}G"

WORK="${TMPDIR:-/tmp}/orig_ptr_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
rm -rf "$WORK"; mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT
FA="$WORK/$PAIR.fasta"
if [[ "$ASM" == *.gz ]]; then READER="gunzip -c"; else READER="cat"; fi
$READER "$ASM" | awk '
  /^>/ {
    split(substr($0, 2), a, " "); id = a[1]
    n = split(id, b, "|")
    if (n >= 3 && b[1] == "ENA") id = b[2]
    else if (n > 1) gsub(/\|/, "_", id)
    print ">" id; next
  }
  { print }' > "$FA"
total=$(grep -c '^>' "$FA")
uniq_n=$(grep '^>' "$FA" | sort -u | wc -l)
echo "    contigs  $total ($uniq_n unique names)"
[ "$total" -eq "$uniq_n" ] || { echo "duplicate contig names after normalising" >&2; exit 1; }

# A failed attempt's partial output would be taken as a cache by pyTransRate.
rm -rf "$DEST"
conda run --no-capture-output -n orp pytransrate \
    -o "$DEST" -a "$FA" \
    --left "$C1" --right "$C2" -t "$CPU" \
    --max-memory "${MEM}G"

[ -s "$DEST/assemblies.csv" ] || { echo "no assemblies.csv written" >&2; exit 1; }
# The BAMs are tens of GB and nothing downstream reads them.
find "$DEST" -name '*.bam' -delete
cat "$DEST/assemblies.csv"
echo "=== done $PAIR at $(date) ==="
