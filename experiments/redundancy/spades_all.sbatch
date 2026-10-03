#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J spades_all
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
#
# spades_branch.sbatch over every dataset with a Trans-ABySS and a Trinity
# assembly. One array task per line of SAMPLES_FILE (run, transabyss,
# trinity, R1, R2; tab-separated, written by spades_all.sh):
#   1. rnaSPAdes twice, spadesauto and spadeshigh (4.1.0 defaults), from the
#      run's existing corrected reads (spades_only.py; the "[kmer]" lines give
#      the k values used, copied to kmers.txt);
#   2. chowder.py's default twotrack merge of those two with the run's
#      existing Trans-ABySS and Trinity assemblies, again on the corrected
#      reads, production settings (--tpm-filt 1).
# Works in $OUT_BASE/<run>/ (GPFS home), runout <run>_spadesbr.
#
# Submit with spades_all.sh, which builds SAMPLES_FILE and the log directory.
# usage: sbatch --array=1-N%20 spades_all.sbatch SAMPLES_FILE CODE_DIR
set -euo pipefail

SAMPLES=${1:?usage: sbatch spades_all.sbatch SAMPLES_FILE CODE_DIR}
CODE=${2:?usage: sbatch spades_all.sbatch SAMPLES_FILE CODE_DIR}
OUT_BASE=${OUT_BASE:-$HOME/redundancy_tests/spades_all}
echo "=== task $SLURM_ARRAY_TASK_ID of $SLURM_ARRAY_JOB_ID on $(hostname) at $(date) ==="

line=$(grep . "$SAMPLES" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[[ -n $line ]] || { echo "no line $SLURM_ARRAY_TASK_ID in $SAMPLES" >&2; exit 1; }
IFS=$'\t' read -r RUN TA TR R1 R2 <<< "$line"
for f in "$TA" "$TR" "$R1" "$R2"; do
    [[ -s $f ]] || { echo "missing $f" >&2; exit 1; }
done
OUT=$OUT_BASE/$RUN
R=${RUN}_spadesbr

CPU=${SLURM_CPUS_PER_TASK:-24}
if [[ -n ${SLURM_MEM_PER_NODE:-} ]]; then MEM=$(( SLURM_MEM_PER_NODE / 1024 )); else MEM=120; fi

ln -sfn "spades_all_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$OUT_BASE/logs/$RUN.log"
module load anaconda/colsa
rm -rf "$OUT"
mkdir -p "$OUT"
cd "$OUT"
echo "$RUN spades_all: $(cat "$CODE/version.txt"), $(git -C "$CODE" rev-parse --short HEAD), ${CPU} cpus, ${MEM}G"
echo "    transabyss $TA"
echo "    trinity    $TR"
echo "    R1         $R1"
echo "    R2         $R2"
echo "writing to $OUT"
status=0
(
set -o pipefail
python3 "$CODE/experiments/redundancy/spades_only.py" --oyster "$CODE/oyster.py" \
    --read1 "$R1" --read2 "$R2" --trimmed-corrected-reads --cpu "$CPU" --mem "$MEM" \
    --runout "$R" --dir "$OUT" 2>&1 | tee "$OUT/spades_only.log"
python3 "$CODE/chowder.py" \
    --assemblies "$OUT/assemblies/$R.spadesauto.fasta" "$OUT/assemblies/$R.spadeshigh.fasta" "$TA" "$TR" \
    --labels spadesauto spadeshigh transabyss trinity --max-parallel 2 \
    --read1 "$R1" --read2 "$R2" \
    --corrected-reads --cpu "$CPU" --mem "$MEM" --tpm-filt 1 --runout "$R" --dir "$OUT"
) || status=$?
echo "exit status $status"

grep -h '^\[kmer\]' "$OUT/spades_only.log" > "$OUT/kmers.txt" || true

# Free the bulky files; the rest stays where it is.
find "$OUT/rcorr" -maxdepth 1 -type f -name '*.cor.fq' -delete 2>/dev/null || true
find "$OUT" -type f \( -name '*.bam' -o -name '*.bam.bai' -o -name GenomeIndex -o -name GenomeIndexHash \
     -o -name OverflowTable -o -name Genome \) -delete
echo "=== done $RUN at $(date), results in $OUT ==="
exit "$status"
