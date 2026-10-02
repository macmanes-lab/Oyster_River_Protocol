#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp410
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
#SBATCH --output=%x_%A_%a.log
#
# ORP 4.1.0-dev0 end to end through chowder.py, production settings: the
# default merge (--merge-method twotrack), the normal end-of-run cleanup, and
# the settings the original 4.0 runs used (--tpm-filt 1, --cpu 24, --mem 120).
# Inputs are the sample's existing assemblies and corrected reads, so only the
# merge half runs. One array task per line of SAMPLES_FILE.
#
# usage: sbatch --array=1-N prod_410.sbatch SAMPLES_FILE CODE_DIR [OUT_BASE]
set -euo pipefail

SAMPLES=${1:?usage: sbatch prod_410.sbatch SAMPLES_FILE CODE_DIR [OUT_BASE]}
CODE=${2:?usage: sbatch prod_410.sbatch SAMPLES_FILE CODE_DIR [OUT_BASE]}
OUT_BASE=${3:-$HOME/redundancy_tests/prod_410}
COMPARE=${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}
RUN=$(grep -v '^\s*$' "$SAMPLES" | sed -n "${SLURM_ARRAY_TASK_ID}p")
SRC=$COMPARE/orp_runs/$RUN

labels=(spades55 spades75 transabyss trinity)
assemblies=()
for l in "${labels[@]}"; do
    f=$SRC/assemblies/$RUN.$l.fasta.gz
    [[ $l == trinity ]] && f=$SRC/assemblies/$RUN.trinity.Trinity.fasta.gz
    [[ -s $f ]] || { echo "missing $f" >&2; exit 1; }
    assemblies+=("$f")
done

module load anaconda/colsa
mkdir -p "$OUT_BASE/$RUN"
echo "$RUN 4.1.0 on $(hostname): $(cat "$CODE/version.txt"), $(git -C "$CODE" rev-parse --short HEAD)"
python3 "$CODE/chowder.py" \
    --assemblies "${assemblies[@]}" --labels "${labels[@]}" --assembly-order given \
    --read1 "$SRC/rcorr/$RUN.TRIM_1P.cor.fq.gz" --read2 "$SRC/rcorr/$RUN.TRIM_2P.cor.fq.gz" \
    --corrected-reads --cpu 24 --mem 120 --tpm-filt 1 \
    --runout "${RUN}_v410" --dir "$OUT_BASE/$RUN"
