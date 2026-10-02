#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J chowder_arms
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-3
#SBATCH --output=%x_%A_%a.log
#
# The merge half of ORP, end to end, three ways, on one sample's existing
# assemblies: the orthogroup step as ORP runs it now (diamond, -I 12) as the
# control, then blastn at -I 2 and -I 3. Everything else -- pytransrate,
# makeorthout, the diamond rescue, cd-hit-est, salmon, the TPM filter, BUSCO,
# transrate -- is the same code in all three, through chowder.py.
#
# usage: sbatch chowder_arms.sbatch RUN [CODE_DIR] [OUT_DIR]
#   RUN       run name under $COMPARE/orp_runs; its four assemblies and its
#             trimmed+corrected reads are reused, so no assembler reruns
#   CODE_DIR  ORP code to run (default ~/redundancy_tests/orp_code, a copy
#             carrying --orthofinder-program/--orthofinder-inflation, kept
#             apart from the checkout the 108-sample array runs from)
#
# Settings copied from the original run that matter after the assemblers:
# --tpm-filt 1, --cpu 24, --mem 120, --max-parallel 2 (the default). Each arm
# runs with --no-cleanup, so the orthogroups, the rescue lists and cd-hit's
# clusters survive for tracing where each duplicate copy came from.
set -euo pipefail

RUN=${1:?usage: sbatch chowder_arms.sbatch RUN [CODE_DIR] [OUT_DIR]}
CODE=${2:-$HOME/redundancy_tests/orp_code}
OUT=${3:-$HOME/redundancy_tests/chowder_arms/$RUN}
COMPARE=${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}
SRC=$COMPARE/orp_runs/$RUN

ARMS=(diamond_I12 blastn_I2 blastn_I3)
arm=${ARMS[$((SLURM_ARRAY_TASK_ID - 1))]}
program=${arm%_I*}
inflation=${arm#*_I}

# ORP's own concatenation order (ASSEMBLY_ORDER in oyster.py), kept with
# --assembly-order given so cd-hit-est's tie-breaks see what oyster.py would.
labels=(spades55 spades75 transabyss trinity)
assemblies=()
for l in "${labels[@]}"; do
    f=$SRC/assemblies/$RUN.$l.fasta.gz
    [[ $l == trinity ]] && f=$SRC/assemblies/$RUN.trinity.Trinity.fasta.gz
    [[ -s $f ]] || { echo "missing $f" >&2; exit 1; }
    assemblies+=("$f")
done

module load anaconda/colsa
mkdir -p "$OUT/$arm"
echo "$arm on $(hostname): orthofinder -S $program -I $inflation, code $CODE"
python3 "$CODE/chowder.py" \
    --assemblies "${assemblies[@]}" --labels "${labels[@]}" --assembly-order given \
    --read1 "$SRC/rcorr/$RUN.TRIM_1P.cor.fq.gz" --read2 "$SRC/rcorr/$RUN.TRIM_2P.cor.fq.gz" \
    --corrected-reads --cpu "$SLURM_CPUS_PER_TASK" --mem 120 --tpm-filt 1 --no-cleanup \
    --orthofinder-program "$program" --orthofinder-inflation "$inflation" \
    --runout "${RUN}_$arm" --dir "$OUT/$arm"
