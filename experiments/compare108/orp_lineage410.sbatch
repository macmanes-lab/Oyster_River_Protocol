#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp_lin410
#SBATCH --cpus-per-task=32
#SBATCH --mem 190G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-50%20
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/lineage_410/logs/orp_%A_%a.log

# ORP 4.1.0 rerun of the 50 lineage-specific samples (mammals 7, non-mammal
# vertebrates 9, insects 10, land plants 10, echinoderms 10, other 4), in one
# array instead of the six per-group ones (orp_mammalia.sbatch etc).
#
# What changed from those dev30 runs: --min-cov auto in place of --min-cov 2,
# and no --normalize-reads (Trinity gets --no_normalize_reads). --min-cov auto
# is ORP's default; it is passed explicitly so the timing log's command line
# says so. Everything else is the same: per-sample BUSCO lineage, no
# --tpm-filt, --transabyss-mpi on, --max-parallel 2, 32 cpus and 190G.
#
# Runs go to compare/lineage_410/<SRR>/, not the orp_<group>/ trees, so the
# dev30 results stay intact for the before/after comparison and the two runs
# of a sample never share a directory.
#
# The manifest has the five orp_array.sbatch columns plus two:
#   pair_id  original_assembly  R1  R2  SRR  group  lineage
# Build it once with build_lineage410_manifest.sh before submitting:
#
#   cd /mnt/home/macmaneslab/macmanes/compare
#   ~/Oyster_River_Protocol/experiments/compare108/build_lineage410_manifest.sh
#   sbatch ~/Oyster_River_Protocol/experiments/compare108/orp_lineage410.sbatch
#
# The task reads its line by number, so do not rewrite the manifest while the
# array has pending tasks.

set -euo pipefail

MANIFEST="${MANIFEST:-/mnt/home/macmaneslab/macmanes/compare/lineage_410/manifest.tsv}"
RUNS="${RUNS:-/mnt/home/macmaneslab/macmanes/compare/lineage_410}"
ORP="${ORP:-$HOME/Oyster_River_Protocol/oyster.py}"
LOGDIR="$RUNS/logs"

echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="

module load anaconda/colsa

line=$(grep . "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "$line" ]; then
    echo "manifest has no line ${SLURM_ARRAY_TASK_ID}; array range is wider than the manifest" >&2
    exit 1
fi
IFS=$'\037' read -r TSA TSA_FASTA R1 R2 SRR GROUP LINEAGE EXTRA <<< "$(printf '%s' "$line" | tr '\t' '\037')"
if [ -n "${EXTRA:-}" ]; then
    echo "manifest line $SLURM_ARRAY_TASK_ID has more than 7 fields" >&2
    exit 1
fi
[ -n "$SRR" ] || { echo "manifest line $SLURM_ARRAY_TASK_ID has no run name in column 5" >&2; exit 1; }
[ -n "${LINEAGE:-}" ] || { echo "manifest line $SLURM_ARRAY_TASK_ID has no lineage in column 7" >&2; exit 1; }

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

echo "=== sample: $SRR ($TSA, $GROUP) ==="
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
    --min-cov auto \
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
