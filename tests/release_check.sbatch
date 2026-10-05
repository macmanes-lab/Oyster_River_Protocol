#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp_release_check
#SBATCH --cpus-per-task=32
#SBATCH --mem 140G
#SBATCH --output=orp_release_check.%j.out
# The release check (tests/release_check.py) on one node: every tier, four
# full runs at a time on 8 cores and 32 GB each. Submit from the directory
# the runs should go in; anything after the script name is passed on.
#
#     sbatch /path/to/Oyster_River_Protocol/tests/release_check.sbatch
#     sbatch .../release_check.sbatch --only chowder_default chowder_options
#
# Results: ./orp_release_check/summary.txt, logs in ./orp_release_check/logs/.

set -euo pipefail
REPO="$(cd "$(dirname "$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')")/.." && pwd)"
python3 "$REPO/tests/release_check.py" --jobs 4 --cpu 8 --mem 32 "$@"
