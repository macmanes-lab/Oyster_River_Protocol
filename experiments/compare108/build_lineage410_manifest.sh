#!/bin/bash
# build_lineage410_manifest.sh - the 50-line manifest for orp_lineage410.sbatch.
#
# Concatenates the six per-group manifests the dev30 lineage runs used
# (compare/orp_<group>/manifest.tsv) and appends two columns, the group and
# its BUSCO lineage, so one array can carry all six groups. The lineages are
# the ones those runs used (LINEAGE_RUNS.md); the four "other" samples each
# get their own.
#
#   ./build_lineage410_manifest.sh            # writes $COMPARE/lineage_410/manifest.tsv
#
# Refuses to overwrite a manifest while an orp_lin410 array is queued or
# running: the tasks read it by line number.
set -euo pipefail

COMPARE="${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}"
OUT="${OUT:-$COMPARE/lineage_410/manifest.tsv}"

if [ -f "$OUT" ] && squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -q '^orp_lin410$'; then
    echo "refusing to rewrite $OUT while an orp_lin410 array is queued or running" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT")/logs"
tmp=$(mktemp "$OUT.XXXX")
trap 'rm -f "$tmp"' EXIT

for group in mammalia vertebrata insecta embryophyta echinodermata other; do
    case "$group" in
        mammalia)      lineage=mammalia_odb12 ;;
        vertebrata)    lineage=vertebrata_odb12 ;;
        insecta)       lineage=insecta_odb12 ;;
        embryophyta)   lineage=embryophyta_odb12.2 ;;
        echinodermata) lineage=metazoa_odb12.2 ;;
        other)         lineage="" ;;
    esac
    awk -F'\t' -v OFS='\t' -v g="$group" -v lin="$lineage" '
        NF == 0 { next }
        NF != 5 { printf "%s: line %d has %d fields, expected 5\n", FILENAME, FNR, NF > "/dev/stderr"; exit 1 }
        {
            l = lin
            if (g == "other") {
                if      ($1 == "tsa_GARO")       l = "crustacea_odb12.2"
                else if ($1 == "tsa_GAGE")       l = "agaricales_odb12.2"
                else if ($1 == "tsa_GAZR")       l = "arachnida_odb12.2"
                else if ($1 == "orp_ERR1674585") l = "lophotrochozoa_odb12.2"
                else { printf "no lineage for %s in group other\n", $1 > "/dev/stderr"; exit 1 }
            }
            print $1, $2, $3, $4, $5, g, l
        }' "$COMPARE/orp_$group/manifest.tsv" >> "$tmp"
done

n=$(grep -c . "$tmp")
dups=$(cut -f5 "$tmp" | sort | uniq -d)
if [ -n "$dups" ]; then
    echo "run name (column 5) appears more than once: $dups" >&2
    exit 1
fi
missing=0
while IFS=$'\t' read -r pair asm r1 r2 srr group lineage; do
    for f in "$r1" "$r2"; do
        [ -f "$f" ] || { echo "missing reads for $pair: $f" >&2; missing=1; }
    done
done < "$tmp"
[ "$missing" -eq 0 ] || exit 1

mv "$tmp" "$OUT"
trap - EXIT
echo "wrote $n samples to $OUT"
cut -f6,7 "$OUT" | sort | uniq -c
