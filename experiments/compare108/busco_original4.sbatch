#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J busco_orig4
#SBATCH --cpus-per-task=24
#SBATCH --mem 100G
#SBATCH --exclude=node117,node118
#SBATCH --array=51-59%5
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/busco_original/logs/busco_%A_%a.log

# Lineage-specific BUSCO on the original assemblies of the nine orp_ pairs
# that add_orp_originals.sh appended to lineage_410/manifest.tsv (lines
# 51-59): the 2018 ORP orthomerged assemblies, pairs/<pair>/assembly/. The
# companion to busco_original.sbatch, 2 and 3, which scored the other 50
# originals (six of them orp_); all four write to busco_original/<pair>/.
#
# Unlike those three it reads the pair (column 1) and lineage (column 7) from
# the manifest by array index, so the original and its ORP run are scored
# against the same dataset by construction. The BUSCO call is unchanged:
# ORP's own (the orp_busco env, busco_dbs, --offline, transcriptome mode,
# ORP's config.ini).
#
#   sbatch ~/Oyster_River_Protocol/experiments/compare108/busco_original4.sbatch

set -euo pipefail

MANIFEST="${MANIFEST:-/mnt/home/macmaneslab/macmanes/compare/lineage_410/manifest.tsv}"
PAIRS="${PAIRS:-/mnt/home/macmaneslab/macmanes/compare/pairs}"
OUT="${OUT:-/mnt/home/macmaneslab/macmanes/compare/busco_original}"
ORPDIR="${ORPDIR:-$HOME/Oyster_River_Protocol}"
CPU="${SLURM_CPUS_PER_TASK:-24}"

echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="

line=$(grep . "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$line" ] || { echo "manifest has no line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r PAIR FASTA _ _ _ _ LINEAGE <<< "$line"
[ -n "$LINEAGE" ] || { echo "manifest line ${SLURM_ARRAY_TASK_ID} has no lineage" >&2; exit 1; }
[ -f "$FASTA" ] || { echo "missing original assembly for $PAIR: $FASTA" >&2; exit 1; }
[ -d "$ORPDIR/busco_dbs/lineages/$LINEAGE" ] || {
    echo "$LINEAGE is not in $ORPDIR/busco_dbs; download it on the login node" >&2; exit 1; }

DEST="$OUT/$PAIR"
mkdir -p "$OUT/logs"
ln -sfn "busco_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$OUT/logs/${PAIR}.log"

echo "    pair     $PAIR"
echo "    assembly $FASTA"
echo "    lineage  $LINEAGE"
echo "    out      $DEST"

if compgen -G "$DEST/short_summary*.txt" >/dev/null; then
    echo "already complete, skipping"
    exit 0
fi

module load anaconda/colsa

WORK="${TMPDIR:-/tmp}/busco_orig4_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
rm -rf "$WORK"; mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

IN="$FASTA"
if [[ "$FASTA" == *.gz ]]; then
    IN="$WORK/$PAIR.fasta"
    gzip -dc "$FASTA" > "$IN"
fi
echo "    contigs  $(grep -c '^>' "$IN")"

export BUSCO_CONFIG_FILE="$ORPDIR/software/config.ini"
echo "    $(conda run -n orp_busco busco --version)"

conda run --no-capture-output -n orp_busco busco --offline \
    --lineage "$LINEAGE" \
    --download_path "$ORPDIR/busco_dbs" \
    -i "$IN" -m transcriptome --cpu "$CPU" \
    -o "$PAIR" --config "$BUSCO_CONFIG_FILE"

mkdir -p "$DEST"
cp "$WORK/$PAIR"/short_summary* "$DEST/"
find "$WORK/$PAIR" -name full_table.tsv -exec cp {} "$DEST/" \; -quit
find "$WORK/$PAIR" -name missing_busco_list.tsv -exec cp {} "$DEST/" \; -quit
cp "$WORK/$PAIR"/logs/busco.log "$DEST/" 2>/dev/null || true

grep -h "C:" "$DEST"/short_summary*.txt | head -1
echo "=== done $PAIR at $(date) ==="
