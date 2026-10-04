#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J repick_validate
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
#SBATCH --output=%x_%A_%a.log
#
# Re-run everything after makeorthout under a new pick rule, starting from a
# finished arm's orthogroups and pytransrate scores, one array task per line
# of TASKS_FILE:
#
#   sample  source_arm_dir  runout  program  inflation  new_arm  rule
#
# The new arm keeps the source's runout (every file name carries it) and is
# written to validate_arms/<sample>/<new_arm> on home. Like
# validate_arms.sbatch, the work happens in fast scratch and is moved off at
# the end: the analysis files copied out, the rest in one run.tar, scratch
# deleted once the archive checks out.
#
# What is copied in, with timestamps, so every step up to makeorthout reads as
# up to date: the ingested assemblies and their diamond blastx output, the
# filtered and search fastas, orthofuser.done, merged.fasta and pytransrate's
# csvs. Members missing from the source directory (validate_arms.sbatch
# archived them) are taken from its run.tar. The corrected reads were pruned
# from every source, so they are decompressed again and given the ingest
# marker's timestamp: newer than the source .gz, so chowder keeps them, and
# older than merged/assemblies.csv, so pytransrate's scoring of the pool still
# counts as current. The new good.<runout>.list is then written by the real
# picker and is newer than its inputs, so chowder starts at orthofusing.
#
# The fast fileset has a hidden space quota as well as its 5M-inode limit.
# Each task peaks around 60 GB (decompressed reads, BAMs), and ten at once
# exceeded the quota on 2026-10-02 -- some steps failed outright, and some
# BAMs were silently truncated. Run at most 3 at a time (%3), and check logs
# for "quota exceeded" before trusting a result.
#
# usage: sbatch --array=1-N%3 repick_validate.sbatch TASKS_FILE [CODE_DIR]
set -euo pipefail

TASKS=${1:?usage: sbatch repick_validate.sbatch TASKS_FILE [CODE_DIR]}
CODE=${2:-$HOME/redundancy_tests/orp_code}
COMPARE=${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}
SCRATCH_BASE=${SCRATCH_BASE:-/mnt/gpfs01/fast/macmaneslab/$USER/scratch}

read -r RUN SRCDIR R program inflation arm rule < <(grep -v '^\s*\(#\|$\)' "$TASKS" | sed -n "${SLURM_ARRAY_TASK_ID}p")
SRC=$COMPARE/orp_runs/$RUN
OUT=$HOME/redundancy_tests/validate_arms/$RUN/$arm
WORK=$SCRATCH_BASE/repick_validate/$RUN/$arm.${SLURM_JOB_ID}
labels=(spades55 spades75 transabyss trinity)

module load anaconda/colsa
rm -rf "$OUT" "$WORK"
mkdir -p "$OUT" "$WORK"
echo "$RUN $arm on $(hostname): rule $rule on $SRCDIR ($program -I $inflation, runout $R)"
echo "scratch $WORK -> $OUT"

# 1. Upstream files, from the source directory or its run.tar.
want=("assemblies/$R.ingest.done" "orthofuse/$R/orthofuser.done" "orthofuse/$R/merged.fasta"
      "orthofuse/$R/merged/assemblies.csv" "orthofuse/$R/merged/contigs.csv")
for l in "${labels[@]}"; do
    want+=("assemblies/$R.$l.fasta" "assemblies/diamond/$R.$l.diamond.txt"
           "orthofuse/$R/working/$R.$l.fasta.short.fasta" "orthofuse/$R/search/$R.$l.fasta.short.fasta")
done
from_tar=()
for f in "${want[@]}"; do
    if [[ -e $SRCDIR/$f ]]; then
        mkdir -p "$WORK/$(dirname "$f")"
        cp -p "$SRCDIR/$f" "$WORK/$f"
    else
        from_tar+=("./$f")
    fi
done
groups=$(ls "$SRCDIR"/orthofuse/"$R"/search/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt 2>/dev/null | head -1 || true)
if (( ${#from_tar[@]} )) || [[ -z $groups ]]; then
    [[ -s $SRCDIR/run.tar ]] || { echo "missing ${from_tar[*]} and no run.tar in $SRCDIR" >&2; exit 1; }
    printf '%s\n' "${from_tar[@]}" > "$WORK/.from_tar"
    tar -xf "$SRCDIR/run.tar" -C "$WORK" -T "$WORK/.from_tar"
    rm -f "$WORK/.from_tar"
fi
[[ -n $groups ]] || { echo "no Orthogroups.txt in $SRCDIR" >&2; exit 1; }
for f in "${want[@]}"; do [[ -e $WORK/$f ]] || { echo "could not stage $f" >&2; exit 1; }; done
echo "staged ${#want[@]} upstream files (${#from_tar[@]} from run.tar)"

# 2. Corrected reads, re-decompressed and dated to the ingest marker.
mkdir -p "$WORK/rcorr"
for m in 1 2; do
    gzip -dc "$SRC/rcorr/$RUN.TRIM_${m}P.cor.fq.gz" > "$WORK/rcorr/$R.TRIM_${m}P.cor.fq"
    touch -r "$WORK/assemblies/$R.ingest.done" "$WORK/rcorr/$R.TRIM_${m}P.cor.fq"
done

# 3. The new pick, by the picker itself -- or, for rule twotrack, by
#    twotrack_select.py, which replaces the orthogroups altogether (one contig
#    per swissprot gene, plus the no-hit contigs after cd-hit-est).
diamond=()
for l in "${labels[@]}"; do diamond+=("$WORK/assemblies/diamond/$R.$l.diamond.txt"); done
chowder_rule=$rule
if [[ $rule == twotrack ]]; then
    python3 "$HOME/redundancy_tests/twotrack_select.py" --merged "$WORK/orthofuse/$R/merged.fasta" \
        --contigs-csv "$WORK/orthofuse/$R/merged/contigs.csv" --diamond "${diamond[@]}" \
        --sprot "$CODE/software/diamond/uniprot_sprot.fasta" --out "$WORK/orthofuse/$R/good.$R.list" \
        --cdhit /mnt/home/macmaneslab/macmanes/orp_envs/orp/bin/cd-hit-est --threads "$SLURM_CPUS_PER_TASK"
    chowder_rule=score  # makeorthout is skipped either way; the list above is newer
else
    python3 "$CODE/experiments/redundancy/pick_best_contigs.py" "$WORK/orthofuse/$R/merged/contigs.csv" "$groups" \
        "$WORK/orthofuse/$R/good.$R.list" --rule "$rule" --diamond "${diamond[@]}"
fi
echo "picks: $(wc -l < "$WORK/orthofuse/$R/good.$R.list") (source arm: $(wc -l < "$SRCDIR/orthofuse/$R/good.$R.list"))"

# 4. Everything after the pick.
assemblies=()
for l in "${labels[@]}"; do
    f=$SRC/assemblies/$RUN.$l.fasta.gz
    [[ $l == trinity ]] && f=$SRC/assemblies/$RUN.trinity.Trinity.fasta.gz
    assemblies+=("$f")
done
status=0
python3 "$CODE/chowder.py" \
    --assemblies "${assemblies[@]}" --labels "${labels[@]}" --assembly-order given \
    --read1 "$SRC/rcorr/$RUN.TRIM_1P.cor.fq.gz" --read2 "$SRC/rcorr/$RUN.TRIM_2P.cor.fq.gz" \
    --corrected-reads --cpu 24 --mem 120 --tpm-filt 1 --no-cleanup \
    --orthofinder-program "$program" --orthofinder-inflation "$inflation" --pick-rule "$chowder_rule" \
    --runout "$R" --dir "$WORK" || status=$?
echo "chowder.py exited $status"

# 5. Move off scratch, as validate_arms.sbatch does -- unless chowder failed,
#    in which case scratch is left exactly as it is, BAMs and all, for
#    debugging; only the log points at it.
if [[ $status -ne 0 ]]; then
    echo "*** chowder.py failed; scratch kept untouched at $WORK ***"
    exit "$status"
fi
find "$WORK/rcorr" -maxdepth 1 -type f -name '*.cor.fq' -delete
find "$WORK" -type f \( -name '*.bam' -o -name '*.bam.bai' -o -name 'GenomeIndex' -o -name 'GenomeIndexHash' \
     -o -name 'OverflowTable' -o -name 'Genome' \) -delete
mkdir -p "$OUT/orthofuse/$R/search/OrthoFinder/Results_source/Orthogroups"
cp -p "$groups" "$OUT/orthofuse/$R/search/OrthoFinder/Results_source/Orthogroups/Orthogroups.txt"
KEEP=$(mktemp)
( cd "$WORK" && find . -type f \( \
      -path "./assemblies/$R.ORP.fasta" -o -path "./assemblies/$R.ORP.intermediate.fasta" \
      -o -path "./assemblies/$R.ORP.diamond.txt" -o -path "./assemblies/diamond/$R.list7" \
      -o -path "./orthofuse/$R/good.$R.list" -o -path "./orthofuse/$R/merged.fasta" \
      -o -path "./orthofuse/$R/merged/*.csv" \
      -o -path "./reports/run_$R.ORP/short_summary*" -o -path "./reports/run_$R.ORP/run_*/full_table.tsv" \
      -o -path "./reports/transrate_$R/*.csv" -o -path "./reports/$R.timing.log" \
      -o -path "./reports/qualreport.$R" \) ) > "$KEEP"
rsync -a --files-from="$KEEP" "$WORK"/ "$OUT"/
echo "copied $(wc -l < "$KEEP") analysis files to $OUT"
rm -f "$KEEP"
tar -C "$WORK" -cf "$OUT/run.tar" .
n_scratch=$(find "$WORK" -mindepth 1 | wc -l)
n_tar=$(tar -tf "$OUT/run.tar" | grep -vc '^\./$' || true)
if [[ $n_scratch -eq $n_tar ]]; then
    rm -rf "$WORK"
    echo "moved off scratch: $n_tar entries in $OUT/run.tar ($(du -h "$OUT/run.tar" | cut -f1))"
else
    echo "*** archive check failed: $n_scratch in scratch, $n_tar in run.tar; scratch kept at $WORK ***"
    [[ $status -eq 0 ]] && status=1
fi
exit "$status"
