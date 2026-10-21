#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orig_uniq
#SBATCH --cpus-per-task=16
#SBATCH --mem 48G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-50%10
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/orig_unique_genes/logs/orig_uniq_%A_%a.log

# Unique swissprot genes in each ORIGINAL assembly of the 50 lineage samples,
# counted exactly as ORP counts its own (orp_diamond / orp_uniq in oyster.py):
#
#   diamond blastx --quiet -p <cpu> -e 1e-8 --top 0.1 -q <fasta> -d <swissprot> -o <out>
#
# in the orp conda env, against the swissprot database of the ORP checkout
# the lineage runs use, then oyster.py's own parse_unique_count over the hits
# (imported, not re-implemented): the distinct gene names, the part of the
# subject id's third |-field before the "_", e.g. sp|P12345|ACTB_HUMAN -> ACTB.
#
# The original is column 2 of the lineage manifest (pair, original assembly,
# R1, R2, SRR, group, lineage). diamond reads the gzipped TSA files directly;
# contig names don't enter the count, so the ENA headers need no rewriting.
#
#   mkdir -p /mnt/home/macmaneslab/macmanes/compare/orig_unique_genes/logs
#   sbatch orig_unique_genes.sbatch
#
# Writes, per pair, to compare/orig_unique_genes/:
#   <pair>.diamond.txt   the hits
#   <pair>.unique.txt    the count, one line, as ORP's unique.ORP.txt
# A pair with a .unique.txt already is skipped, so a rerun only fills gaps.

set -euo pipefail

MANIFEST="${MANIFEST:-/mnt/home/macmaneslab/macmanes/compare/lineage_410/manifest.tsv}"
OUT="${OUT:-/mnt/home/macmaneslab/macmanes/compare/orig_unique_genes}"
ORPDIR="${ORPDIR:-$HOME/Oyster_River_Protocol}"
DB="$ORPDIR/software/diamond/swissprot"
CPU="${SLURM_CPUS_PER_TASK:-16}"

echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="

line=$(grep . "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "$line" ]; then
    echo "manifest has no line ${SLURM_ARRAY_TASK_ID}; array range is wider than the manifest" >&2
    exit 1
fi
IFS=$'\t' read -r PAIR ASM _ <<< "$line"
[ -f "$ASM" ] || { echo "missing original assembly for $PAIR: $ASM" >&2; exit 1; }
[ -f "$DB.dmnd" ] || { echo "missing swissprot database: $DB.dmnd" >&2; exit 1; }

echo "    pair     $PAIR"
echo "    assembly $ASM"
echo "    db       $DB.dmnd"
echo "    ORP      $(cat "$ORPDIR/version.txt") ($(git -C "$ORPDIR" rev-parse --short HEAD 2>/dev/null || echo '?'))"

mkdir -p "$OUT"
if [ -s "$OUT/$PAIR.unique.txt" ]; then
    echo "already done: $(cat "$OUT/$PAIR.unique.txt") unique genes; skipping"
    exit 0
fi

module load anaconda/colsa

# .part then rename, so a killed task never leaves a short hit file that
# looks finished.
conda run --no-capture-output -n orp diamond blastx --quiet -p "$CPU" \
    -e 1e-8 --top 0.1 -q "$ASM" -d "$DB" -o "$OUT/$PAIR.diamond.txt.part"
mv "$OUT/$PAIR.diamond.txt.part" "$OUT/$PAIR.diamond.txt"

conda run --no-capture-output -n orp python - "$ORPDIR" "$OUT/$PAIR.diamond.txt" \
    "$OUT/$PAIR.unique.txt" <<'EOF'
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from oyster import parse_unique_count
count = parse_unique_count(Path(sys.argv[2]))
Path(sys.argv[3]).write_text(f"{count}\n")
EOF

echo "=== $PAIR: $(cat "$OUT/$PAIR.unique.txt") unique genes, $(wc -l < "$OUT/$PAIR.diamond.txt") hits; done at $(date) ==="
