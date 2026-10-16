#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp_echino
#SBATCH --cpus-per-task=32
#SBATCH --mem 190G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-10%2
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/orp_echinodermata/logs/orp_%A_%a.log

# ORP 4.1.0 on ten echinoderm pairs, one per species, scored against
# metazoa_odb12.2 instead of the default eukaryota_odb12.2. BUSCO has no
# echinoderm or deuterostome dataset, so metazoa is the closest one that
# includes them. The companion to the mammal, vertebrate, insect and plant
# arrays; only the pairs, the run tree and the lineage differ.
#
# The ten span 2.0-7.5 GB of reads (all 20 echinoderm pairs fall in that
# range) and all five classes: sea urchins (Eucidaris, Lytechinus,
# Sphaerechinus), a feather star (Oxycomanthus), a brittle star (Ophiocoma),
# sea cucumbers (Sclerodactyla, Apostichopus) and sea stars (Patiria miniata,
# Luidia, Asterias amurensis).
#
# metazoa_odb12.2 is OrthoDB 12.2, like embryophyta_odb12.2; the mammal,
# vertebrate and insect arrays use 12.1 datasets.
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
#   mkdir -p orp_echinodermata/logs
#   grep -P '^(tsa_GAZP|tsa_GAZO|tsa_GAUQ|tsa_GAUR|tsa_GAUT|tsa_GAWB|tsa_GAVQ|tsa_GAVS|tsa_GAVJ|tsa_GAVR)\t' \
#       manifest.tsv > orp_echinodermata/manifest.tsv
#   sbatch ~/Oyster_River_Protocol/experiments/compare108/orp_echinodermata.sbatch

set -euo pipefail

MANIFEST="${MANIFEST:-/mnt/home/macmaneslab/macmanes/compare/orp_echinodermata/manifest.tsv}"
RUNS="${RUNS:-/mnt/home/macmaneslab/macmanes/compare/orp_echinodermata}"
ORP="${ORP:-$HOME/Oyster_River_Protocol/oyster.py}"
LINEAGE="${LINEAGE:-metazoa_odb12.2}"
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
