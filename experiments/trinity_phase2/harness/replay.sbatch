#!/bin/bash
#SBATCH -p macmanes
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J trin_replay
#SBATCH -o %x.%j.out
# usage: sbatch replay.sbatch TAG WORKROOT SLOTS [extra trinity flags...]
set -euo pipefail
export LC_ALL=${REPLAY_LC:-en_US.UTF-8}
TAG=$1; ROOT=$2; SLOTS=$3; shift 3
EXTRA="$*"
E=~/orp_envs/orp_trinity
ENG=~/trinity_eng
SHIM=$ENG/shim; mkdir -p $SHIM
ln -sf $ENG/javashim $SHIM/java; ln -sf $ENG/pyshim $SHIM/python
export REAL_JAVA=$E/bin/java REAL_PYTHON=$E/bin/python
export PATH=$SHIM:$E/bin:$PATH
WORK=$ROOT/$TAG
rm -rf $WORK; mkdir -p $WORK
export TRINITY_PROF_JAVA=$WORK/java_py.tsv
echo "node=$(hostname) cpus=$(nproc) tag=$TAG work=$WORK slots=$SLOTS extra=$EXTRA"
python3 $ENG/replay.py stage --sample $ENG/${SAMPLE:-sample.tsv} --work $WORK --trinity ${TRINITY:-$ENG/src/Trinity} --extra="$EXTRA"
s=$(date +%s)
$ENG/src/trinity-plugins/BIN/ParaFly -c $WORK/cmds -CPU $SLOTS -shuffle -failed_cmds $WORK/failed.txt > $WORK/parafly.log 2>&1 || true
e=$(date +%s)
echo "PARAFLY_WALL $((e-s)) s"
python3 $ENG/replay.py compare --sample $ENG/${SAMPLE:-sample.tsv} --work $WORK
wc -l < $WORK/failed.txt 2>/dev/null || true
