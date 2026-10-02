#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J repick_arms
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-10
#SBATCH --output=%x_%A_%a.log
#
# Re-run the steps after makeorthout under an alternative pick rule
# (repick.py), starting from a finished chowder_arms.sbatch arm. Each task
# gets a fresh run directory holding copies of everything up to the pick --
# the ingested assemblies and their diamond hits, the filtered and search
# fastas, the orthofuser sentinel, merged.fasta and pytransrate's csvs, all
# copied with their timestamps so every step up to makeorthout reads as up
# to date -- plus a new good.<run>.list written after them. chowder.py then
# reruns orthofusing onwards: rescue, cd-hit, salmon, TPM filter, BUSCO,
# transrate. Nothing in the source arm is written to: the reads are the only
# symlinks, and they are only read.
#
# Task 1 is the sanity check: blastn_I2 under today's rule must reproduce
# that arm's results exactly, or the comparison means nothing.
#
# usage: sbatch repick_arms.sbatch [RUN] [ARMS_DIR] [CODE_DIR]
set -euo pipefail

RUN=${1:-SRR1139197}
ARMS=${2:-$HOME/redundancy_tests/chowder_arms/$RUN}
CODE=${3:-$HOME/redundancy_tests/orp_code}
HERE=$HOME/redundancy_tests
COMPARE=${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}
SRC=$COMPARE/orp_runs/$RUN

TASKS=(blastn_I2:score
       diamond_I12:score_len diamond_I12:score_orf diamond_I12:near_best
       blastn_I2:score_len   blastn_I2:score_orf   blastn_I2:near_best
       blastn_I3:score_len   blastn_I3:score_orf   blastn_I3:near_best)
task=${TASKS[$((SLURM_ARRAY_TASK_ID - 1))]}
arm=${task%%:*}
rule=${task#*:}
runout=${RUN}_$arm
from=$ARMS/$arm
to=$HERE/repick_arms/$RUN/${arm}.$rule
ofrom=$from/orthofuse/$runout
oto=$to/orthofuse/$runout

echo "$arm under rule $rule on $(hostname): $from -> $to"
rm -rf "$to"
mkdir -p "$to/assemblies/diamond" "$oto/working" "$oto/search" "$oto/merged" "$to/rcorr"

labels=(spades55 spades75 transabyss trinity)
for l in "${labels[@]}"; do
    cp -p "$from/assemblies/$runout.$l.fasta" "$to/assemblies/"
    cp -p "$from/assemblies/diamond/$runout.$l.diamond.txt" "$to/assemblies/diamond/"
    cp -p "$ofrom/working/$runout.$l.fasta.short.fasta" "$oto/working/"
    cp -p "$ofrom/search/$runout.$l.fasta.short.fasta" "$oto/search/"
done
cp -p "$from/assemblies/$runout.ingest.done" "$to/assemblies/"
cp -p "$ofrom/orthofuser.done" "$ofrom/merged.fasta" "$oto/"
cp -p "$ofrom/merged/assemblies.csv" "$ofrom/merged/contigs.csv" "$oto/merged/"
for m in 1 2; do
    ln -s "$from/rcorr/$runout.TRIM_${m}P.cor.fq" "$to/rcorr/$runout.TRIM_${m}P.cor.fq"
done

groups=$(ls "$ofrom"/search/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt | head -1)
module load anaconda/colsa
python3 "$HERE/repick.py" --pick-script "$CODE/scripts/pick_best_contigs.py" \
    --contigs-csv "$oto/merged/contigs.csv" --orthogroups "$groups" \
    --rule "$rule" --out "$oto/good.$runout.list"
echo "picks: $(wc -l < "$oto/good.$runout.list")"
if [[ $rule == score ]]; then
    cmp "$oto/good.$runout.list" "$ofrom/good.$runout.list" \
        && echo "sanity: pick list identical to the arm's"
fi

assemblies=()
for l in "${labels[@]}"; do
    f=$SRC/assemblies/$RUN.$l.fasta.gz
    [[ $l == trinity ]] && f=$SRC/assemblies/$RUN.trinity.Trinity.fasta.gz
    assemblies+=("$f")
done
program=${arm%_I*}
inflation=${arm#*_I}
python3 "$CODE/chowder.py" \
    --assemblies "${assemblies[@]}" --labels "${labels[@]}" --assembly-order given \
    --read1 "$SRC/rcorr/$RUN.TRIM_1P.cor.fq.gz" --read2 "$SRC/rcorr/$RUN.TRIM_2P.cor.fq.gz" \
    --corrected-reads --cpu "$SLURM_CPUS_PER_TASK" --mem 120 --tpm-filt 1 --no-cleanup \
    --orthofinder-program "$program" --orthofinder-inflation "$inflation" \
    --runout "$runout" --dir "$to"
