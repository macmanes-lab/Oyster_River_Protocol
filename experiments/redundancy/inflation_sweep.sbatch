#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J of_inflation
#SBATCH --cpus-per-task=4
#SBATCH --mem 60G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-12
#SBATCH --output=%x_%A_%a.log
#
# Re-cluster the saved all-vs-all searches from time_orthofinder_search.sbatch
# at several MCL inflations, one array task per (search program, -I) pair.
# OrthoFinder's -b reuses the searches, so each task is only the clustering.
#
# usage: sbatch inflation_sweep.sbatch [TIMING_DIR]
#   TIMING_DIR  default ~/redundancy_tests/search_timing/SRR1139197
#
# -og is not honoured by OrthoFinder 3.1.5 (NOTES.md): it goes on to
# MSA/trees after writing the orthogroups. This is a benchmark, not an ORP
# run, so the task stops OrthoFinder once its log says "Done orthogroups"
# rather than paying for work nothing reads.
set -euo pipefail

TIMING=${1:-$HOME/redundancy_tests/search_timing/SRR1139197}
PROGRAMS=(diamond blastn)
INFLATIONS=(1.5 2 3 5 8 12)

i=$((SLURM_ARRAY_TASK_ID - 1))
prog=${PROGRAMS[$((i / ${#INFLATIONS[@]}))]}
infl=${INFLATIONS[$((i % ${#INFLATIONS[@]}))]}
wd=$(ls -d "$TIMING/$prog"/in/OrthoFinder/Results_*/WorkingDirectory | head -1)
name="${prog}_I${infl}"
log=$TIMING/sweep/$name.log
mkdir -p "$TIMING/sweep"

module load anaconda/colsa
echo "$name on $(hostname), searches from $wd"
SECONDS=0
setsid conda run --no-capture-output -n orp_orthofinder \
    orthofinder -b "$wd" -og -I "$infl" -a "$SLURM_CPUS_PER_TASK" -n "$name" \
    > "$log" 2>&1 &
pid=$!
while kill -0 "$pid" 2>/dev/null; do
    if grep -q "Done orthogroups" "$log"; then
        sleep 20   # let the remaining orthogroup files land
        kill -- -"$pid" 2>/dev/null || true
        break
    fi
    sleep 10
done
wait "$pid" 2>/dev/null || true
grep -q "Done orthogroups" "$log" || { echo "FAILED: no orthogroups, see $log"; exit 1; }

# -b writes a fresh results tree under the reused working directory.
groups=$wd/OrthoFinder/Results_$name/Orthogroups/Orthogroups.txt
[[ -s $groups ]] || { echo "FAILED: cannot find Orthogroups.txt for $name"; exit 1; }
ln -sfn "$groups" "$TIMING/sweep/$name.Orthogroups.txt"
echo "$name: orthogroups in $((SECONDS / 60))m $((SECONDS % 60))s -> $groups"
