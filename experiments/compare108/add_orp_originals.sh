#!/bin/bash
# add_orp_originals.sh - append the nine orp_ pairs the lineage runs left out
# to lineage_410/manifest.tsv, as lines 51-59.
#
# compare/pairs has 15 orp_ pairs: the samples of the 2018 ORP paper, whose
# "original" assembly is that paper's orthomerged assembly (Zenodo
# doi:10.5281/zenodo.1320141). Six were already in the 50 lineage samples
# (DRR036858, SRR1789336, DRR046632, ERR489297, SRR3499127, ERR1674585).
# This adds the other nine, in the 7-column format build_lineage410_manifest.sh
# writes (pair, original assembly, R1, R2, SRR, group, lineage), with the most
# specific odb12 dataset for each:
#
#   orp_DRR030368   Echinococcus multilocularis  other        lophotrochozoa_odb12.2
#   orp_DRR031870   Vigna angularis              embryophyta  embryophyta_odb12.2
#   orp_DRR053698   Cephalotus follicularis      embryophyta  embryophyta_odb12.2
#   orp_DRR069093   Hevea brasiliensis           embryophyta  embryophyta_odb12.2
#   orp_DRR082659   Aeginetia indica             embryophyta  embryophyta_odb12.2
#   orp_ERR058009   Entamoeba histolytica        other        amoebozoa_odb12.2
#   orp_ERR1016675  Heterorhabditis indica       other        nematoda_odb12.2
#   orp_SRR2016923  Phyllodoce medipapillata     other        lophotrochozoa_odb12.2
#   orp_SRR2086412  Mus musculus                 mammalia     mammalia_odb12
#
# The plants and the mouse take their group's dataset, as the lineage runs do
# (embryophyta, not fabales/malpighiales; mammalia_odb12 from OrthoDB 12.1, not
# glires). BUSCO has no flatworm or annelid dataset, so the tapeworm and the
# polychaete take lophotrochozoa, as Schistosoma does. amoebozoa_odb12.2 and
# nematoda_odb12.2 are new to busco_dbs; download them on the login node first
# (compute nodes have no internet) -- submit_orp_originals.sh does.
#
# Appending, not rebuilding, keeps lines 1-50 where the finished arrays read
# them, so every manifest-driven script takes the new samples with
# --array=51-59. Re-running is a no-op once the nine are present.
#
#   ./add_orp_originals.sh
set -euo pipefail

COMPARE="${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}"
MANIFEST="${MANIFEST:-$COMPARE/lineage_410/manifest.tsv}"
SOURCE="${SOURCE:-$COMPARE/manifest.tsv}"

SAMPLES=(
    "orp_DRR030368   other        lophotrochozoa_odb12.2"
    "orp_DRR031870   embryophyta  embryophyta_odb12.2"
    "orp_DRR053698   embryophyta  embryophyta_odb12.2"
    "orp_DRR069093   embryophyta  embryophyta_odb12.2"
    "orp_DRR082659   embryophyta  embryophyta_odb12.2"
    "orp_ERR058009   other        amoebozoa_odb12.2"
    "orp_ERR1016675  other        nematoda_odb12.2"
    "orp_SRR2016923  other        lophotrochozoa_odb12.2"
    "orp_SRR2086412  mammalia     mammalia_odb12"
)

if squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -qE '^(orp_lin410|orig_uniq|busco_orig4)$'; then
    echo "refusing to change $MANIFEST while an array that reads it is queued or running" >&2
    exit 1
fi

n=$(grep -c . "$MANIFEST")
present=0
for s in "${SAMPLES[@]}"; do
    read -r pair _ <<< "$s"
    cut -f1 "$MANIFEST" | grep -qx "$pair" && present=$((present + 1))
done
if [ "$present" -eq "${#SAMPLES[@]}" ]; then
    echo "all ${#SAMPLES[@]} already in $MANIFEST ($n lines); nothing to do"
    exit 0
elif [ "$present" -ne 0 ]; then
    echo "$present of ${#SAMPLES[@]} are already in $MANIFEST; fix it by hand" >&2
    exit 1
fi
[ "$n" -eq 50 ] || { echo "$MANIFEST has $n lines, expected the 50 lineage samples" >&2; exit 1; }

tmp=$(mktemp "$MANIFEST.XXXX")
trap 'rm -f "$tmp"' EXIT
cp "$MANIFEST" "$tmp"

missing=0
for s in "${SAMPLES[@]}"; do
    read -r pair group lineage <<< "$s"
    # Columns 1-5 come from the 108-sample manifest, as for the other 50.
    line=$(awk -F'\t' -v p="$pair" '$1 == p' "$SOURCE")
    [ "$(printf '%s\n' "$line" | grep -c .)" -eq 1 ] || {
        echo "expected one line for $pair in $SOURCE" >&2; missing=1; continue; }
    IFS=$'\t' read -r _ asm r1 r2 srr <<< "$line"
    for f in "$asm" "$r1" "$r2"; do
        [ -f "$f" ] || { echo "missing for $pair: $f" >&2; missing=1; }
    done
    [ -d "$HOME/Oyster_River_Protocol/busco_dbs/lineages/$lineage" ] || \
        echo "warning: $lineage is not in busco_dbs yet; download it before submitting" >&2
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$pair" "$asm" "$r1" "$r2" "$srr" "$group" "$lineage" >> "$tmp"
done
[ "$missing" -eq 0 ] || exit 1

dups=$(cut -f5 "$tmp" | sort | uniq -d)
[ -z "$dups" ] || { echo "run name (column 5) appears more than once: $dups" >&2; exit 1; }

mv "$tmp" "$MANIFEST"
trap - EXIT
echo "appended ${#SAMPLES[@]} samples to $MANIFEST (now $(grep -c . "$MANIFEST") lines)"
cut -f6,7 "$MANIFEST" | sort | uniq -c
