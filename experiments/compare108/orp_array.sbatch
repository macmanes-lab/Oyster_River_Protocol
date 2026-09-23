#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp108
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-108%6
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/orp_runs/logs/orp_%A_%a.log

# One array task per line of the manifest, which is tab-separated:
#
#     tsa_GADU <TAB> assembly.fasta.gz <TAB> SRR527277_1.fastq.gz <TAB> SRR527277_2.fastq.gz <TAB> SRR527277
#
# Column 2 is the TSA assembly, which oyster.py has no use for -- that is
# chowder.py's input -- so it is read and ignored here.
#
# Each sample gets its own --dir, so the runs never touch each other's
# assemblies/, reports/ or orthofuse/ trees, and --runout is column 5, so every
# file inside carries the SRR rather than the tsa code.
#
# Submit with ./submit_all_orp.sh, not with sbatch directly. That script creates
# the log directory, checks the manifest and passes RUNS/MANIFEST/ORP in via
# --export.
# The --array and --output directives above are only a fallback for a bare
# `sbatch orp_array.sbatch`, and that path needs the log directory to exist
# already:
# slurm opens the --output file before this script runs, so a missing log
# directory kills the task before line one with nothing written anywhere.

set -euo pipefail

# The manifest layout this script expects. submit_all_orp.sh refuses to submit
# unless this matches what it is about to hand over, so a half-copied set of
# scripts is caught at submit time rather than by 108 tasks misparsing a line.
MANIFEST_COLUMNS=5

MANIFEST="${MANIFEST:-/mnt/home/macmaneslab/macmanes/compare/manifest.tsv}"
RUNS="${RUNS:-/mnt/home/macmaneslab/macmanes/compare/orp_runs}"
ORP="${ORP:-$HOME/Oyster_River_Protocol/oyster.py}"
LOGDIR="$RUNS/logs"

# Before anything that can fail, so a task that dies early still says which
# sample it was. If this line is missing from a log, the failure is upstream of
# the job script -- slurm could not launch it at all.
echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="

module load anaconda/colsa

# grep -v '^$' first, so a blank line in the manifest does not silently shift
# every task after it onto the wrong sample.
line=$(grep . "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "$line" ]; then
    echo "manifest has no line ${SLURM_ARRAY_TASK_ID}; array range is wider than the manifest" >&2
    exit 1
fi
# The trailing EXTRA is deliberate: read assigns everything left over to its
# last variable, separators and all, so parsing a 5-field line with 4 variables
# silently produces an R2 of "<path><tab><runname>". Catching the overflow here
# turns that into a stated error.
IFS=$'\t' read -r TSA TSA_FASTA R1 R2 SRR EXTRA <<< "$line"
if [ -n "${EXTRA:-}" ]; then
    echo "manifest line $SLURM_ARRAY_TASK_ID has more than $MANIFEST_COLUMNS fields;" \
         "leftover: [$EXTRA]" >&2
    exit 1
fi

# The run name comes from column 5 rather than being parsed out of $R1 here.
# A task sees only its own line, so it cannot tell that the accession in its
# filename is shared with another sample; make_manifest.sh settles that once,
# with all 108 in view. See its collision-resolution block.
[ -n "$SRR" ] || { echo "manifest line $SLURM_ARRAY_TASK_ID has no run name in column 5" >&2; exit 1; }

for f in "$R1" "$R2"; do
    [ -f "$f" ] || { echo "missing read file: $f" >&2; exit 1; }
done

mkdir -p "$RUNS/$SRR" "$LOGDIR"

# SLURM can only name its log by job/array id, so point an SRR-named symlink at
# it: `less logs/SRR527277.log` works without looking the task number up.
ln -sfn "orp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$LOGDIR/${SRR}.log"
ln -sfn "orp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$LOGDIR/${TSA}.log"

# Slurm starts a task in the directory it was submitted from, which is now
# inside the git checkout. Everything below addresses absolute paths, but a
# sub-tool that drops a scratch file relative to cwd would leave it in the
# repo; from here on cwd is the run's own directory instead.
cd "$RUNS/$SRR" || { echo "cannot cd to $RUNS/$SRR" >&2; exit 1; }

echo "=== sample: $SRR ($TSA) ==="
echo "    R1  $R1"
echo "    R2  $R2"
echo "    dir $RUNS/$SRR"

# ORP skips steps whose outputs are already current, so a resubmitted task
# resumes. This guard is the cheaper version of that for samples that finished:
# it costs nothing and keeps a rerun of the whole array short.
if [ -f "$RUNS/$SRR/reports/qualreport.$SRR.done" ]; then
    echo "already complete, skipping"
    exit 0
fi

# --cpu and --mem come from what slurm actually granted this task, not from
# numbers repeated here: the #SBATCH directives above are the single place the
# allocation is stated, and a job that asked for 24 cpus while telling ORP it
# had 16 would leave a third of the node idle for a day per sample. Editing a
# directive now changes both.
CPU="${SLURM_CPUS_PER_TASK:-24}"
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    MEM=$(( SLURM_MEM_PER_NODE / 1024 ))          # slurm reports MB
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    MEM=$(( SLURM_MEM_PER_CPU * CPU / 1024 ))
else
    MEM=120
fi
echo "    allocation: ${CPU} cpus, ${MEM}G"

"$ORP" \
    --normalize-reads \
    --tpm-filt 1 \
    --mem "$MEM" \
    --cpu "$CPU" \
    --read1 "$R1" \
    --read2 "$R2" \
    --max-parallel 2 \
    --dir "$RUNS/$SRR" \
    --runout "$SRR"

echo "=== done $SRR at $(date) ==="
