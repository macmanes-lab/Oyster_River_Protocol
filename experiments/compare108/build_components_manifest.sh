#!/bin/bash
# build_components_manifest.sh - the manifest for components_score.sbatch.
#
# One line per component assembly of the lineage_410 runs (4 per sample: 200
# for the first 50, 236 once add_orp_originals.sh has added nine): the four
# assemblies ORP merges (Trinity, rnaSPAdes auto-k, rnaSPAdes high-k,
# Trans-ABySS), each to be scored with BUSCO and pyTransRate exactly as ORP
# scores its own merged assembly. Columns:
#
#   SRR  assembler  assembly.fasta.gz  R1.cor.fq.gz  R2.cor.fq.gz  lineage
#
# The reads are the run's own trimmed, error-corrected pairs (rcorr/), the
# same ones ORP's pyTransRate step used on the ORP assembly; the lineage is
# column 7 of lineage_410/manifest.tsv.
#
# The run name is read from each run's finished qualreport, not assumed: the
# first 50 ran as <SRR>_v4_1_0, the nine orp_ additions as <SRR>_v4_0_0 (the
# same code, renumbered; only resume handling and the banner changed).
#
# Samples keep their manifest order, so appending samples to the lineage
# manifest only appends lines here, and the tasks already run keep their
# numbers. The script checks that: a rewrite that would move an existing
# line is refused.
#
#   ./build_components_manifest.sh     # writes $RUNS/components/manifest.tsv
set -euo pipefail

RUNS="${RUNS:-/mnt/home/macmaneslab/macmanes/compare/lineage_410}"
OUT="${OUT:-$RUNS/components/manifest.tsv}"

if [ -f "$OUT" ] && squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -q '^orp_comp$'; then
    echo "refusing to rewrite $OUT while an orp_comp array is queued or running" >&2
    exit 1
fi
mkdir -p "$(dirname "$OUT")/logs"
tmp=$(mktemp "$OUT.XXXX")
trap 'rm -f "$tmp"' EXIT

missing=0
while IFS=$'\t' read -r pair asm r1 r2 srr group lineage; do
    [ -n "$srr" ] || continue
    d="$RUNS/$srr"
    shopt -s nullglob
    done_files=("$d"/reports/qualreport."$srr"_v*.done)
    shopt -u nullglob
    if [ "${#done_files[@]}" -ne 1 ]; then
        echo "expected one finished run in $d/reports, found ${#done_files[@]}" >&2
        missing=1
        continue
    fi
    run=$(basename "${done_files[0]}" .done)
    run="${run#qualreport.}"
    c1="$d/rcorr/$run.TRIM_1P.cor.fq.gz"
    c2="$d/rcorr/$run.TRIM_2P.cor.fq.gz"
    for f in "$c1" "$c2"; do
        [ -s "$f" ] || { echo "missing reads: $f" >&2; missing=1; }
    done
    for a in trinity spadesauto spadeshigh transabyss; do
        case "$a" in
            trinity) fa="$d/assemblies/$run.trinity.Trinity.fasta.gz" ;;
            *)       fa="$d/assemblies/$run.$a.fasta.gz" ;;
        esac
        [ -s "$fa" ] || { echo "missing assembly: $fa" >&2; missing=1; }
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$srr" "$a" "$fa" "$c1" "$c2" "$lineage" >> "$tmp"
    done
done < "$RUNS/manifest.tsv"
[ "$missing" -eq 0 ] || exit 1

if [ -f "$OUT" ]; then
    old=$(grep -c . "$OUT")
    if ! head -n "$old" "$tmp" | cmp -s - <(grep . "$OUT"); then
        echo "the first $old lines differ from the existing $OUT; refusing to renumber tasks" >&2
        exit 1
    fi
    new=$(grep -c . "$tmp")
    if [ "$new" -gt "$old" ]; then
        echo "existing $old lines unchanged; new tasks are $((old + 1))-$new"
    else
        echo "existing $old lines unchanged; nothing new"
    fi
fi

mv "$tmp" "$OUT"
trap - EXIT
echo "wrote $(grep -c . "$OUT") component assemblies to $OUT"
cut -f2 "$OUT" | sort | uniq -c
