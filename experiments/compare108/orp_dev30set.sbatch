#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp_dev30set
#SBATCH --cpus-per-task=32
#SBATCH --mem 190G
#SBATCH --exclude=node117,node118
#SBATCH --array=51-59%9
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/dev30_settings/logs/orp_%A_%a.log

# Section 4.4 for the nine orp_ libraries added on 2026-10-22 (lineage_410
# manifest lines 51-59). The first 50 have a dev30 run (orp_<group>.sbatch:
# --normalize-reads, --min-cov 2); these nine do not. This runs them with the
# dev30 settings so the settings comparison covers all 59:
#
#   --normalize-reads --min-cov 2      (the two settings dev30 differed in)
#   no --tpm-filt, --transabyss-mpi on, --max-parallel 2, 32 cpus, 190G,
#   the run's lineage                  (the same as every other run)
#
# Code: the ORP checkout this task runs (the 4.0.0 release by default), so for
# these nine the before/after differs in those two settings and nothing else.
# To use the dev30 build itself instead, point ORP at a worktree of the commit
# whose version.txt reads 4.1.0-dev30, e.g.
#   git -C ~/Oyster_River_Protocol worktree add ~/orp_dev30 <commit>
#   sbatch --export=ALL,ORP=$HOME/orp_dev30/oyster.py orp_dev30set.sbatch
#
# Runs go to compare/dev30_settings/<SRR>/ with --runout <SRR>_v<version>-dev30set,
# so nothing in lineage_410/ or orp_<group>/ is touched. Collect with
#   collect_lineage410.py --runs '$COMPARE/dev30_settings' \
#       --runout-tag v4_0_0-dev30set -o dev30set_metrics.csv
# (submit_dev30set_trace.sh collect does this and rebuilds the settings table).

set -euo pipefail

MANIFEST="${MANIFEST:-/mnt/home/macmaneslab/macmanes/compare/lineage_410/manifest.tsv}"
RUNS="${RUNS:-/mnt/home/macmaneslab/macmanes/compare/dev30_settings}"
ORP="${ORP:-$HOME/Oyster_River_Protocol/oyster.py}"
LOGDIR="$RUNS/logs"

echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="
module load anaconda/colsa

# Fail at once, not hours in, if this checkout has dropped the flag.
if ! "$ORP" --help 2>&1 | grep -q -- '--normalize-reads'; then
    echo "$ORP has no --normalize-reads option; cannot reproduce the dev30 settings" >&2
    exit 1
fi

line=$(grep . "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$line" ] || { echo "manifest has no line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\037' read -r PAIR _ORIG R1 R2 SRR GROUP LINEAGE EXTRA <<< "$(printf '%s' "$line" | tr '\t' '\037')"
[ -z "${EXTRA:-}" ] || { echo "manifest line $SLURM_ARRAY_TASK_ID has more than 7 fields" >&2; exit 1; }
for f in "$R1" "$R2"; do
    [ -f "$f" ] || { echo "missing read file: $f" >&2; exit 1; }
done

mkdir -p "$RUNS/$SRR" "$LOGDIR"
ln -sfn "orp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$LOGDIR/${SRR}.log"
ln -sfn "orp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$LOGDIR/${PAIR}.log"
cd "$RUNS/$SRR"

VERSION=$(cat "$(dirname "$ORP")/version.txt")
RUNNAME="${SRR}_v${VERSION//./_}-dev30set"
echo "=== sample: $SRR ($PAIR, $GROUP), lineage $LINEAGE, ORP $VERSION ($(git -C "$(dirname "$ORP")" rev-parse --short HEAD 2>/dev/null || echo '?')), runout $RUNNAME ==="

if [ -f "$RUNS/$SRR/reports/qualreport.$RUNNAME.done" ]; then
    echo "already complete, skipping"
    exit 0
fi

CPU="${SLURM_CPUS_PER_TASK:-32}"
MEM=$(( ${SLURM_MEM_PER_NODE:-194560} / 1024 ))

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
