#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp_ta_mpi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem 500G
#SBATCH --output=orp_ta_mpi_%j.log
# A full oyster.py run with Trans-ABySS under MPI. With --transabyss-mpi on,
# preflight stops the run in its first minute if mpirun cannot start ranks
# here; use auto to fall back to threaded instead.
#
#     sbatch --export=ALL,R1=<reads_1.fq.gz>,R2=<reads_2.fq.gz>,RUN=DRR031870_mpi \
#         experiments/transabyss_mpi/oyster_mpi.sbatch
#
# Afterwards: grep "=== Assemblers\|transabyss" in the log,
# reports/$RUN.transabyss.mode says "mpi 16", and reports/$RUN.timing.log has
# run_transabyss's time to set against a threaded run of the same reads.

set -euo pipefail
: "${R1:?set R1}" "${R2:?set R2}" "${RUN:?set RUN}"
ORP="${ORP:-$HOME/Oyster_River_Protocol/oyster.py}"
DIR="${DIR:-$PWD/$RUN}"

module load anaconda/colsa

# --cpu and --mem from what slurm granted, as compare108/orp_array.sbatch does.
CPU="${SLURM_CPUS_PER_TASK:-40}"
MEM=$(( ${SLURM_MEM_PER_NODE:-512000} / 1024 ))
mkdir -p "$DIR"
cd "$DIR"

"$ORP" \
    --normalize-reads \
    --tpm-filt 1 \
    --mem "$MEM" \
    --cpu "$CPU" \
    --read1 "$R1" \
    --read2 "$R2" \
    --max-parallel 2 \
    --transabyss-mpi on \
    --dir "$DIR" \
    --runout "$RUN"
