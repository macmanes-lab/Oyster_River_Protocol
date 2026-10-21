#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp_comp
#SBATCH --cpus-per-task=24
#SBATCH --mem 100G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-200%20
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/lineage_410/components/logs/comp_%A_%a.log

# Score each of the four component assemblies (Trinity, rnaSPAdes auto-k,
# rnaSPAdes high-k, Trans-ABySS) of the 50 ORP 4.1.0 lineage runs with ORP's
# own assembly_eval.py, so the paper can compare the merged ORP assembly with
# each assembler alone. 50 runs x 4 assemblers = 200 tasks, 20 at a time.
#
# assembly_eval.py runs the three evaluations an ORP run ends with --
# pyTransRate, BUSCO and strandeval -- in the same conda envs (orp, orp_busco)
# with the same BUSCO config and busco_dbs, and writes the same qualreport,
# so the component scores line up with the ORP ones in lineage410_metrics.csv.
# The reads are the run's own trimmed, error-corrected pairs from rcorr/,
# which is what oyster.py's pyTransRate step used on the ORP assembly.
#
# Outputs, per task, in components/<SRR>/<assembler>/ (--runout <SRR>_<assembler>):
#   reports/qualreport.<run>          BUSCO, pyTransRate scores, proper pairs
#   reports/run_<run>.ORP/            BUSCO
#   reports/pytransrate_<run>/        pyTransRate (assemblies.csv, contigs.csv)
#   reports/<run>.strandeval_summary.txt, reports/<run>.flagstat
# assembly_eval.py's .done sentinels make a resubmitted task resume.
#
# Build the manifest once, then submit:
#   ~/Oyster_River_Protocol/experiments/compare108/build_components_manifest.sh
# (later additions: rerun it and submit only the new lines, e.g. --array=201-236%20,
# as submit_orp_originals.sh does)
#   sbatch ~/Oyster_River_Protocol/experiments/compare108/components_score.sbatch
#
# The task reads its line by number, so do not rewrite the manifest while the
# array has pending tasks.

set -euo pipefail

RUNS="${RUNS:-/mnt/home/macmaneslab/macmanes/compare/lineage_410}"
MANIFEST="${MANIFEST:-$RUNS/components/manifest.tsv}"
EVAL="${EVAL:-$HOME/Oyster_River_Protocol/assembly_eval.py}"

echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="

module load anaconda/colsa

line=$(grep . "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$line" ] || { echo "manifest has no line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r SRR ASM FASTA R1 R2 LINEAGE EXTRA <<< "$line"
[ -z "${EXTRA:-}" ] || { echo "manifest line ${SLURM_ARRAY_TASK_ID} has more than 6 fields" >&2; exit 1; }
for f in "$FASTA" "$R1" "$R2"; do
    [ -s "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done

CPU="${SLURM_CPUS_PER_TASK:-24}"
MEM=$(( ${SLURM_MEM_PER_NODE:-102400} / 1024 ))
NAME="${SRR}_${ASM}"
OUT="$RUNS/components/$SRR/$ASM"
mkdir -p "$OUT"
ln -sfn "comp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$RUNS/components/logs/$NAME.log"

if [ -f "$OUT/reports/qualreport.$NAME" ]; then
    echo "already complete, skipping"
    exit 0
fi

# assembly_eval.py takes a plain fasta, so the gzipped component is unpacked
# into the output tree (not node scratch: its .done sentinels record the
# assembly by size and mtime, and a fresh copy each attempt would defeat the
# resume). Removed once the report is written.
FA="$OUT/$NAME.fasta"
if [ ! -s "$FA" ]; then
    gunzip -c "$FASTA" > "$FA.tmp" && mv "$FA.tmp" "$FA"
fi
echo "    assembler $ASM ($(grep -c '^>' "$FA") contigs), lineage $LINEAGE, ${CPU} cpus, ${MEM}G"

"$EVAL" \
    --assembly "$FA" \
    --read1 "$R1" \
    --read2 "$R2" \
    --runout "$NAME" \
    --outdir "$OUT" \
    --lineage "$LINEAGE" \
    --cpu "$CPU" \
    --mem "$MEM"

[ -f "$OUT/reports/qualreport.$NAME" ] && rm -f "$FA"
echo "=== done $NAME at $(date) ==="
