#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J reads_lost
#SBATCH --cpus-per-task=8
#SBATCH --mem 32G
#SBATCH --exclude=node117,node118
#SBATCH --output=%x_%A_%a.log
#
# reads_lost.py for one sample per array task: which control contigs carried
# the reads a trial no longer maps.
#
# usage: sbatch --array=1-N reads_lost.sbatch SAMPLES_FILE [TRIAL]   (default two-track)
set -euo pipefail
SAMPLES=${1:?usage: sbatch reads_lost.sbatch SAMPLES_FILE [TRIAL]}
TRIAL=${2:-twotrack}
RUN=$(grep -v '^\s*$' "$SAMPLES" | sed -n "${SLURM_ARRAY_TASK_ID}p")
V=$HOME/redundancy_tests/validate_arms/$RUN
OUT=$HOME/redundancy_tests/reads_lost
mkdir -p "$OUT"
module load anaconda/colsa
python3 "$HOME/redundancy_tests/reads_lost.py" "$V/control" "$V/$TRIAL" "$OUT/$RUN.$TRIAL" \
    --threads "$SLURM_CPUS_PER_TASK"
