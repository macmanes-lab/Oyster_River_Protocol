#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J spades_branch
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
#SBATCH --output=%x_%A_%a.log
#
# The spades branch on one sample (one array task per line of SAMPLES_FILE):
#   1. rnaSPAdes twice, spadesauto and spadeshigh, from the sample's existing
#      corrected reads (spades_only.py; the "[kmer]" lines in the log give the
#      k values used);
#   2. chowder.py's default twotrack merge of those two with the sample's
#      existing Trans-ABySS and Trinity assemblies, production settings.
# Works directly in $HOME/redundancy_tests/spades_branch/<run>/ (GPFS home,
# no scratch step). Run 3 at a time: --array=1-10%3.
#
# usage: sbatch --array=1-N%3 spades_branch.sbatch SAMPLES_FILE CODE_DIR
set -euo pipefail

SAMPLES=${1:?usage: sbatch spades_branch.sbatch SAMPLES_FILE CODE_DIR}
CODE=${2:?usage: sbatch spades_branch.sbatch SAMPLES_FILE CODE_DIR}
COMPARE=${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}
RUN=$(grep -v '^\s*$' "$SAMPLES" | sed -n "${SLURM_ARRAY_TASK_ID}p")
SRC=$COMPARE/orp_runs/$RUN
OUT=$HOME/redundancy_tests/spades_branch/$RUN
WORK=$OUT
R=${RUN}_spadesbr
R1=$SRC/rcorr/$RUN.TRIM_1P.cor.fq.gz
R2=$SRC/rcorr/$RUN.TRIM_2P.cor.fq.gz

for f in "$R1" "$R2" "$SRC/assemblies/$RUN.transabyss.fasta.gz" "$SRC/assemblies/$RUN.trinity.Trinity.fasta.gz"; do
    [[ -s $f ]] || { echo "missing $f" >&2; exit 1; }
done

module load anaconda/colsa
rm -rf "$OUT"
mkdir -p "$OUT"
echo "$RUN spades branch on $(hostname): $(cat "$CODE/version.txt"), $(git -C "$CODE" rev-parse --short HEAD)"
echo "writing to $OUT"
status=0
(
python3 "$CODE/experiments/redundancy/spades_only.py" --oyster "$CODE/oyster.py" \
    --read1 "$R1" --read2 "$R2" --trimmed-corrected-reads --cpu 24 --mem 120 \
    --runout "$R" --dir "$WORK"
python3 "$CODE/chowder.py" \
    --assemblies "$WORK/assemblies/$R.spadesauto.fasta" "$WORK/assemblies/$R.spadeshigh.fasta" \
                 "$SRC/assemblies/$RUN.transabyss.fasta.gz" "$SRC/assemblies/$RUN.trinity.Trinity.fasta.gz" \
    --labels spadesauto spadeshigh transabyss trinity --assembly-order given \
    --read1 "$R1" --read2 "$R2" \
    --corrected-reads --cpu 24 --mem 120 --tpm-filt 1 --runout "$R" --dir "$WORK"
) || status=$?
echo "exit status $status"

# k values used, pulled out for easy collection
grep -h '^\[kmer\]' "$SLURM_SUBMIT_DIR/spades_branch_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" > "$OUT/kmers.txt" || true

# Free the bulky files; the rest stays where it is.
find "$WORK/rcorr" -maxdepth 1 -type f -name '*.cor.fq' -delete 2>/dev/null || true
find "$WORK" -type f \( -name '*.bam' -o -name '*.bam.bai' -o -name GenomeIndex -o -name GenomeIndexHash \
     -o -name OverflowTable -o -name Genome \) -delete
echo "results in $OUT"
exit "$status"
