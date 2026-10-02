#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J validate_arms
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-20
#SBATCH --output=%x_%A_%a.log
#
# The candidate merge (blastn, -I 3, pick rule score_orf) against a same-code
# control (diamond, -I 12, pick rule score: ORP as it is) on several samples,
# two array tasks per sample. Like chowder_arms.sbatch, each task runs the
# merge half through chowder.py on the sample's existing assemblies and
# corrected reads, so no assembler reruns.
#
# usage: sbatch [--array=...] [--mem=...] validate_arms.sbatch SAMPLES_FILE [CODE_DIR]
#   SAMPLES_FILE  one run name per line; task 2k-1 is sample k's control,
#                 task 2k its candidate
#
# ARM=<name> (environment) runs one named arm instead, one task per sample
# (task k is sample k), into validate_arms/<run>/<name>. Known arms:
#   blastn_I12   blastn, -I 12, pick rule score -- ORP's behaviour before
#                the OrthoFinder 3.1.5 upgrade made -d mean diamond
#
# chowder.py always gets --cpu 24 --mem 120 --tpm-filt 1, the settings of
# the original runs, whatever Slurm allocates: --mem steers the pipeline's
# own planning (diamond concurrency, pytransrate's memory cap), so it must
# not differ between the arms or from the original. The Slurm allocation is
# only headroom -- give the largest samples more with sbatch --mem -- and
# sacct's MaxRSS is then the memory each arm really needed.
#
# CHOWDER_EXTRA (environment) is appended to chowder.py's arguments, for a
# rerun that needs something the first attempt did not -- e.g.
#   sbatch --array=4 --mem=240G --export=ALL,CHOWDER_EXTRA="--orthofinder-analysis 4" ...
# when OrthoFinder's 200 s stall watchdog kills a blastn run at -a 1: with
# one worker, a single assembly's hits taking over 200 s to process looks
# like a hang. -a sets parallelism only, not the orthogroups.
#
# The run itself happens in scratch, SCRATCH_BASE (default the lab's fast
# GPFS fileset), not in home: OrthoFinder writes one file per orthogroup, tens
# to hundreds of thousands of them, and from 2026-10-01 home GPFS took up to
# seconds per file create while the fast fileset stayed under a millisecond.
# When chowder.py exits, success or not, the arm is moved off scratch:
#   1. prune what is large and rebuildable (decompressed reads, BAMs, snap
#      indexes);
#   2. copy the files the analysis scripts read (KEEP below) into
#      validate_arms/<run>/<arm>/ on home as ordinary files;
#   3. pack everything else, pruned, into one run.tar beside them;
#   4. check the archive lists as many files as scratch holds, and only then
#      delete the scratch directory.
# Home therefore sees ~20 file creates and one streamed archive per arm. A
# task killed outright (scancel, an OOM kill) never reaches step 1 and leaves
# its scratch directory behind: <SCRATCH_BASE>/validate_arms/<run>/<arm>.<job>.
set -euo pipefail

SAMPLES=${1:?usage: sbatch validate_arms.sbatch SAMPLES_FILE [CODE_DIR]}
CODE=${2:-$HOME/redundancy_tests/orp_code}
COMPARE=${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}

mapfile -t runs < <(grep -v '^\s*$' "$SAMPLES")
i=$((SLURM_ARRAY_TASK_ID - 1))
if [[ -n ${ARM:-} ]]; then
    RUN=${runs[$i]}
    case $ARM in
        blastn_I12) arm=blastn_I12 program=blastn inflation=12 rule=score ;;
        *) echo "unknown ARM=$ARM" >&2; exit 1 ;;
    esac
else
    RUN=${runs[$((i / 2))]}
    if (( i % 2 == 0 )); then
        arm=control program=diamond inflation=12 rule=score
    else
        arm=candidate program=blastn inflation=3 rule=score_orf
    fi
fi
SRC=$COMPARE/orp_runs/$RUN
OUT=$HOME/redundancy_tests/validate_arms/$RUN/$arm
SCRATCH_BASE=${SCRATCH_BASE:-/mnt/gpfs01/fast/macmaneslab/$USER/scratch}
WORK=$SCRATCH_BASE/validate_arms/$RUN/$arm.${SLURM_JOB_ID}
R=${RUN}_$arm

labels=(spades55 spades75 transabyss trinity)
assemblies=()
for l in "${labels[@]}"; do
    f=$SRC/assemblies/$RUN.$l.fasta.gz
    [[ $l == trinity ]] && f=$SRC/assemblies/$RUN.trinity.Trinity.fasta.gz
    [[ -s $f ]] || { echo "missing $f" >&2; exit 1; }
    assemblies+=("$f")
done

module load anaconda/colsa
rm -rf "$OUT" "$WORK"
mkdir -p "$OUT" "$WORK"
echo "$RUN $arm on $(hostname), Slurm mem ${SLURM_MEM_PER_NODE:-?}M:" \
     "orthofinder -S $program -I $inflation, pick rule $rule${CHOWDER_EXTRA:+, extra: $CHOWDER_EXTRA}"
echo "scratch $WORK -> $OUT"
status=0
python3 "$CODE/chowder.py" \
    --assemblies "${assemblies[@]}" --labels "${labels[@]}" --assembly-order given \
    --read1 "$SRC/rcorr/$RUN.TRIM_1P.cor.fq.gz" --read2 "$SRC/rcorr/$RUN.TRIM_2P.cor.fq.gz" \
    --corrected-reads --cpu 24 --mem 120 --tpm-filt 1 --no-cleanup \
    --orthofinder-program "$program" --orthofinder-inflation "$inflation" --pick-rule "$rule" \
    --runout "$R" --dir "$WORK" ${CHOWDER_EXTRA:-} || status=$?
echo "chowder.py exited $status"

# 1. Prune what is large and rebuildable, only inside this arm's scratch.
before=$(du -sh "$WORK" | cut -f1)
find "$WORK/rcorr" -maxdepth 1 -type f -name '*.cor.fq' -delete 2>/dev/null || true
find "$WORK" -type f \( -name '*.bam' -o -name '*.bam.bai' \) -delete
find "$WORK" -type f \( -name 'GenomeIndex' -o -name 'GenomeIndexHash' -o -name 'OverflowTable' \
     -o -name 'Genome' \) -delete
echo "pruned $WORK: $before -> $(du -sh "$WORK" | cut -f1)"

# 2. What compare_arms.py, trace_arms.py, dup_busco_pairs.py and
#    validate_summary.py read, copied out as ordinary files.
KEEP=$(mktemp)
( cd "$WORK" && find . -type f \( \
      -path "./assemblies/$R.ORP.fasta" -o -path "./assemblies/$R.ORP.intermediate.fasta" \
      -o -path "./assemblies/diamond/$R.list7" -o -path "./orthofuse/$R/good.$R.list" \
      -o -path "./orthofuse/$R/merged.fasta" -o -path "./orthofuse/$R/merged/*.csv" \
      -o -path "./orthofuse/$R/search/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt" \
      -o -path "./reports/run_$R.ORP/short_summary*" -o -path "./reports/run_$R.ORP/run_*/full_table.tsv" \
      -o -path "./reports/transrate_$R/*.csv" -o -path "./reports/$R.timing.log" \
      -o -path "./reports/qualreport.$R" \) ) > "$KEEP"
rsync -a --files-from="$KEEP" "$WORK"/ "$OUT"/
echo "copied $(wc -l < "$KEEP") analysis files to $OUT"
rm -f "$KEEP"

# 3. Everything, pruned, as one archive.
tar -C "$WORK" -cf "$OUT/run.tar" .

# 4. Delete scratch only once the archive holds what scratch does.
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
