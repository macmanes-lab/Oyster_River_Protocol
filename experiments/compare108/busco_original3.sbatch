#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J busco_orig3
#SBATCH --cpus-per-task=24
#SBATCH --mem 100G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-4
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/busco_original/logs/busco_%A_%a.log

# BUSCO on the original assemblies of the four pairs in orp_other.sbatch
# (one crustacean, fungus, arachnid and flatworm), each with the lineage its
# ORP run uses, so the old and new assemblies are scored alike. The companion
# to busco_original.sbatch and busco_original2.sbatch; all three write to the
# same busco_original/ tree. "Original" is pairs/<pair>/assembly/: the
# published TSA assembly for tsa_*, and the 2018 ORP orthomerged assembly
# (Zenodo) for orp_*.
#
# The BUSCO call is ORP's own (oyster.py's busco step): the orp_busco env, the
# odb12 lineages in the ORP checkout's busco_dbs, --offline, transcriptome
# mode and ORP's config.ini. Only the input differs.
#
#   mkdir -p /mnt/home/macmaneslab/macmanes/compare/busco_original/logs
#   sbatch ~/Oyster_River_Protocol/experiments/compare108/busco_original3.sbatch

set -euo pipefail

# pair_id and lineage, one task per line. The order is the array order.
SAMPLES=(
    "tsa_GARO        crustacea_odb12.2"
    "tsa_GAGE        agaricales_odb12.2"
    "tsa_GAZR        arachnida_odb12.2"
    "orp_ERR1674585  lophotrochozoa_odb12.2"
)

PAIRS="${PAIRS:-/mnt/home/macmaneslab/macmanes/compare/pairs}"
OUT="${OUT:-/mnt/home/macmaneslab/macmanes/compare/busco_original}"
ORPDIR="${ORPDIR:-$HOME/Oyster_River_Protocol}"
CPU="${SLURM_CPUS_PER_TASK:-24}"

echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="

entry="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]:-}"
[ -n "$entry" ] || { echo "no sample for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
read -r PAIR LINEAGE <<< "$entry"

# Exactly one assembly per pair; more than one would leave it unclear which
# was scored.
shopt -s nullglob
fastas=("$PAIRS/$PAIR"/assembly/*.fasta "$PAIRS/$PAIR"/assembly/*.fasta.gz)
shopt -u nullglob
[ "${#fastas[@]}" -eq 1 ] || {
    echo "expected one assembly in $PAIRS/$PAIR/assembly, found ${#fastas[@]}" >&2; exit 1; }
FASTA="${fastas[0]}"

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

# Node-local work directory: BUSCO writes many small files, and only the
# summaries and tables are kept.
WORK="${TMPDIR:-/tmp}/busco_orig3_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
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
