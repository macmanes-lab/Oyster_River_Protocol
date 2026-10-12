#!/bin/bash
#SBATCH -p macmanes
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J trin_full
#SBATCH -o %x.%j.out
# Full Phase 2 replay with the production (conda) Trinity, launched the way
# production launches it: conda run ... ParaFly -CPU 40 -v -shuffle.
# usage: sbatch full_replay.sbatch TAG ROOT SLOTS [extra flags]
set -euo pipefail
export LC_ALL=en_US.UTF-8
TAG=$1; ROOT=$2; SLOTS=$3; shift 3; EXTRA="$*"
ENG=~/trinity_eng; E=~/orp_envs/orp_trinity
WORK=$ROOT/$TAG; rm -rf $WORK; mkdir -p $WORK
echo "node=$(hostname) start=$(date) tag=$TAG slots=$SLOTS extra=$EXTRA trinity=${TRINITY:-$E/bin/Trinity}"
t0=$(date +%s)
python3 $ENG/replay.py stage --sample $ENG/sample_all.tsv --work $WORK --trinity ${TRINITY:-$E/bin/Trinity} --extra="$EXTRA"
echo "STAGE_WALL $(( $(date +%s) - t0 )) s"
s=$(date +%s)
conda run --no-capture-output -n orp_trinity ${PARAFLY:-$E/bin/trinity-plugins/BIN/ParaFly} -c $WORK/cmds -CPU $SLOTS -v -shuffle -failed_cmds $WORK/failed.txt 2>&1 \
  | tr "\r" "\n" | grep --line-buffered succeeded | awk -v s=$s "NR%500==1 {print systime()-s, \$0; fflush()}" > $WORK/progress.txt || true
echo "PARAFLY_WALL $(( $(date +%s) - s )) s"
python3 $ENG/replay.py compare --sample $ENG/sample_all.tsv --work $WORK
