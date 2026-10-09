#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J ta_mpi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem 200G
#SBATCH --array=1-2
#SBATCH --output=ta_mpi_%A_%a.log
# Trans-ABySS as oyster.py runs it (task 1, "threads") against the same
# command with --mpi (task 2, "mpi"), on the same 16 cores: 40% of --cpu 40,
# what TRANSABYSS_SHARE gives it. See README.md for why and what to read off.
#
#     sbatch --export=ALL,R1=<cor_1.fq.gz>,R2=<cor_2.fq.gz>,NAME=DRR031870 \
#         experiments/transabyss_mpi/run.sbatch
#
# R1/R2: the corrected pair an ORP run assembled from,
#   <dir>/rcorr/<run>.TRIM_{1,2}P.cor.fq(.gz).
# Optional: STRAND=RF|FR (adds --SS, as oyster.py --strand does), OUT (default
#   ./ta_mpi_<NAME>), KMER (default 32, oyster.py's default).

set -euo pipefail

: "${R1:?set R1}" "${R2:?set R2}" "${NAME:?set NAME}"
OUT="${OUT:-$PWD/ta_mpi_$NAME}"
KMER="${KMER:-32}"
STRAND="${STRAND:-}"
REPO="$(cd "$(dirname "$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')")/../.." && pwd)"
CPU="${SLURM_CPUS_PER_TASK:-16}"
MEM=$(( ${SLURM_MEM_PER_NODE:-204800} / 1024 ))

case "${SLURM_ARRAY_TASK_ID:-1}" in
    1) ARM=threads; MPI=() ;;
    2) ARM=mpi;     MPI=(--mpi "$CPU") ;;
    *) echo "array task ${SLURM_ARRAY_TASK_ID} has no arm" >&2; exit 1 ;;
esac

echo "=== $NAME / $ARM on $(hostname) at $(date): ${CPU} cpus, ${MEM}G ==="
for f in "$R1" "$R2"; do [ -f "$f" ] || { echo "missing read file: $f" >&2; exit 1; }; done

module load anaconda/colsa

# One rank per core under a one-task allocation: OpenMPI counts slurm's single
# task as one slot and would refuse -np 16, and would pin every rank to it.
# Both spellings, for OpenMPI 4 (OMPI_) and 5 (PRTE_).
export OMPI_MCA_rmaps_base_oversubscribe=1
export OMPI_MCA_hwloc_base_binding_policy=none
export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
export PRTE_MCA_hwloc_default_binding_policy=none

run="$OUT/$ARM"
work="$run/transabyss"
mkdir -p "$run"
rm -rf "$work"
cd "$run"

# The command oyster.py's run_transabyss builds, plus --mpi in the mpi arm.
cmd=(transabyss)
if [ "$STRAND" = RF ] || [ "$STRAND" = FR ]; then cmd+=(--SS); fi
# ${a[@]+"${a[@]}"}: an empty array is "unbound" under set -u before bash 4.4.
cmd+=(--threads "$CPU" ${MPI[@]+"${MPI[@]}"} --outdir "$work" --kmer "$KMER" --length 250
      --name "$NAME.transabyss.fasta" --pe "$R1" "$R2")
echo "command: ${cmd[*]}"

timer=()
[ -x /usr/bin/time ] && timer=(/usr/bin/time -v -o "$run/time.txt")

start=$(date +%s)
${timer[@]+"${timer[@]}"} conda run --no-capture-output -n orp_transabyss "${cmd[@]}" &
pid=$!
python3 "$REPO/experiments/transabyss_mpi/cpu_sampler.py" "$pid" "$run/cpu.tsv" 60 &
sampler=$!
rc=0
wait "$pid" || rc=$?
wait "$sampler" || true
end=$(date +%s)
if [ "$rc" -ne 0 ]; then
    echo "transabyss exited $rc after $(( end - start ))s" >&2
    exit "$rc"
fi

# oyster.py keeps only the first word of each header (awk_first_field).
awk '/^>/{print $1; next} {print}' "$work/$NAME.transabyss.fasta-final.fa" > "$run/$NAME.$ARM.fasta"
contigs=$(grep -c '^>' "$run/$NAME.$ARM.fasta" || true)
bases=$(grep -v '^>' "$run/$NAME.$ARM.fasta" | tr -d '\n' | wc -c)
maxrss=$(awk -F': ' '/Maximum resident/{printf "%.1fG", $2/1048576}' "$run/time.txt" 2>/dev/null || true)
one_core=$(awk -F'\t' 'NR>1 && $2<1.5{n++} END{print n+0}' "$run/cpu.tsv")
samples=$(awk 'NR>1{n++} END{print n+0}' "$run/cpu.tsv")

{
    printf 'name\tarm\tcpu\twall_s\twall\tmax_rss\tminutes_under_1.5_cores\tminutes_sampled\tcontigs\tbases\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$NAME" "$ARM" "$CPU" "$(( end - start ))" \
        "$(date -u -d @$(( end - start )) +%H:%M:%S)" "${maxrss:-?}" "$one_core" "$samples" "$contigs" "$bases"
} > "$run/summary.tsv"
cat "$run/summary.tsv"
rm -rf "$work"

# Score it the way an ORP run would, on the same reads.
python3 "$REPO/assembly_eval.py" --assembly "$run/$NAME.$ARM.fasta" \
    --read1 "$R1" --read2 "$R2" --cpu "$CPU" --mem "$MEM" \
    --runout "$NAME.$ARM" --outdir "$run/eval"

echo "=== done $NAME / $ARM at $(date) ==="
