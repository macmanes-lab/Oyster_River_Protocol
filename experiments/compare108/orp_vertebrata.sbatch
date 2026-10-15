#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp_vert
#SBATCH --cpus-per-task=32
#SBATCH --mem 190G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-9%2
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/orp_vertebrata/logs/orp_%A_%a.log

# ORP 4.1.0 on the nine non-mammal vertebrate pairs, scored against
# vertebrata_odb12 instead of the default eukaryota_odb12.2. The companion to
# orp_mammalia.sbatch; only the pairs, the run tree and the lineage differ.
#
# Nine, not ten: tsa_GAQK (Hynobius chinensis, the one amphibian) has 62 bp R2
# reads, under make_manifest.sh's 76 bp minimum, so it is not in manifest.tsv.
# Anolis carolinensis appears twice (GAGA, GAGB), the only species that does.
#
# Same manifest format and per-task logic as orp_array.sbatch. What differs:
# the manifest, the run tree, --lineage, a versioned --runout, and filtering
# by --min-cov 2 alone (no --tpm-filt 1), and Trans-ABySS under MPI
# (--transabyss-mpi on). RUNS must not be orp_runs/:
# these SRRs already have finished 4.0 runs there, so the done-guard below
# would skip those tasks.
#
# Build the manifest and log directory once, before submitting:
#
#   cd /mnt/home/macmaneslab/macmanes/compare
#   mkdir -p orp_vertebrata/logs
#   grep -P '^(tsa_GARA|tsa_GAID|tsa_GAGA|tsa_GAJJ|tsa_GAAA|tsa_GAGB|tsa_GBAA|tsa_GAZY|orp_DRR046632)\t' \
#       manifest.tsv > orp_vertebrata/manifest.tsv
#   sbatch ~/Oyster_River_Protocol/experiments/compare108/orp_vertebrata.sbatch

set -euo pipefail

MANIFEST="${MANIFEST:-/mnt/home/macmaneslab/macmanes/compare/orp_vertebrata/manifest.tsv}"
RUNS="${RUNS:-/mnt/home/macmaneslab/macmanes/compare/orp_vertebrata}"
ORP="${ORP:-$HOME/Oyster_River_Protocol/oyster.py}"
LINEAGE="${LINEAGE:-vertebrata_odb12}"
LOGDIR="$RUNS/logs"

echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="

module load anaconda/colsa

line=$(grep . "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "$line" ]; then
    echo "manifest has no line ${SLURM_ARRAY_TASK_ID}; array range is wider than the manifest" >&2
    exit 1
fi
IFS=$'\037' read -r TSA TSA_FASTA R1 R2 SRR EXTRA <<< "$(printf '%s' "$line" | tr '\t' '\037')"
if [ -n "${EXTRA:-}" ]; then
    echo "manifest line $SLURM_ARRAY_TASK_ID has more than 5 fields" >&2
    exit 1
fi
[ -n "$SRR" ] || { echo "manifest line $SLURM_ARRAY_TASK_ID has no run name in column 5" >&2; exit 1; }

for f in "$R1" "$R2"; do
    [ -f "$f" ] || { echo "missing read file: $f" >&2; exit 1; }
done

mkdir -p "$RUNS/$SRR" "$LOGDIR"
ln -sfn "orp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$LOGDIR/${SRR}.log"
ln -sfn "orp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$LOGDIR/${TSA}.log"

cd "$RUNS/$SRR" || { echo "cannot cd to $RUNS/$SRR" >&2; exit 1; }

# The version goes into --runout so every output file says which ORP made it.
# Read from the checkout this task runs, not hardcoded, and once, here: a pull
# mid-array then labels each task with what it actually ran. Dots become
# underscores because the runout is the stem of every filename ORP writes.
VERSION=$(cat "$(dirname "$ORP")/version.txt")
RUNNAME="${SRR}_v${VERSION//./_}"

echo "=== sample: $SRR ($TSA) ==="
echo "    R1      $R1"
echo "    R2      $R2"
echo "    dir     $RUNS/$SRR"
echo "    lineage $LINEAGE"
echo "    ORP     $VERSION (runout $RUNNAME)"

if [ -f "$RUNS/$SRR/reports/qualreport.$RUNNAME.done" ]; then
    echo "already complete, skipping"
    exit 0
fi

CPU="${SLURM_CPUS_PER_TASK:-32}"
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    MEM=$(( SLURM_MEM_PER_NODE / 1024 ))
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    MEM=$(( SLURM_MEM_PER_CPU * CPU / 1024 ))
else
    MEM=190
fi
echo "    allocation: ${CPU} cpus, ${MEM}G"

"$ORP" \
    --normalize-reads \
    --min-cov 2 \
    --lineage "$LINEAGE" \
    --mem "$MEM" \
    --cpu "$CPU" \
    --read1 "$R1" \
    --read2 "$R2" \
    --max-parallel 2 \
    --transabyss-mpi on \
    --dir "$RUNS/$SRR" \
    --runout "$RUNNAME"

echo "=== done $RUNNAME at $(date) ==="
