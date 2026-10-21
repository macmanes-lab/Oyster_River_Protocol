#!/bin/bash
# submit_orp_originals.sh - extend the lineage_410 comparison from 6 to all 15
# orp_ pairs (the 2018 ORP paper's samples), in three steps run by hand from
# the Premise login node:
#
#   ./submit_orp_originals.sh runs        # once
#   ./submit_orp_originals.sh components  # after the nine ORP runs finish
#   ./submit_orp_originals.sh collect     # after everything finishes
#
# runs:
#   1. downloads amoebozoa_odb12.2 and nematoda_odb12.2 into busco_dbs if they
#      are missing (login node only: compute nodes have no internet),
#   2. appends the nine pairs to lineage_410/manifest.tsv as lines 51-59
#      (add_orp_originals.sh),
#   3. submits, all on lines 51-59 of that manifest:
#        orp_lineage410.sbatch    ORP (the checkout's version, 4.0.0), the same
#                                 settings as the 50: --min-cov auto, no TPM
#                                 filter, no normalization, Trans-ABySS MPI
#        busco_original4.sbatch   lineage BUSCO on the nine 2018 originals
#        orig_unique_genes.sbatch unique swissprot genes in the nine originals
#      pyTransRate on the originals is already done for all 15
#      (compare/pytransrate/results/orp_*).
#
# components: rebuilds components/manifest.tsv (200 -> 236 lines, the first
#   200 unchanged) and submits components_score.sbatch on the 36 new lines:
#   assembly_eval.py on each run's Trinity, rnaSPAdes auto-k, rnaSPAdes high-k
#   and Trans-ABySS assembly, with the run's own corrected reads and lineage.
#
# collect: rewrites lineage410_metrics.csv (59 rows) and components_metrics.csv
#   (236 rows) in lineage_410/, keeping the previous copies as *.50.csv.
set -euo pipefail

COMPARE="${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}"
RUNS="$COMPARE/lineage_410"
ORPDIR="${ORPDIR:-$HOME/Oyster_River_Protocol}"
HERE="$ORPDIR/experiments/compare108"
NEW_LINES="51-59"

case "${1:-}" in
runs)
    module load anaconda/colsa
    for lin in amoebozoa_odb12.2 nematoda_odb12.2; do
        if [ ! -d "$ORPDIR/busco_dbs/lineages/$lin" ]; then
            echo "downloading $lin"
            conda run -n orp_busco busco --download "$lin" --download_path "$ORPDIR/busco_dbs"
        fi
    done
    "$HERE/add_orp_originals.sh"
    mkdir -p "$RUNS/logs" "$COMPARE/busco_original/logs" "$COMPARE/orig_unique_genes/logs"
    sbatch --array="$NEW_LINES%9" "$HERE/orp_lineage410.sbatch"
    sbatch --array="$NEW_LINES%5" "$HERE/busco_original4.sbatch"
    sbatch --array="$NEW_LINES%5" "$HERE/orig_unique_genes.sbatch"
    ;;
components)
    old=$(grep -c . "$RUNS/components/manifest.tsv")
    "$HERE/build_components_manifest.sh"
    new=$(grep -c . "$RUNS/components/manifest.tsv")
    [ "$new" -gt "$old" ] || { echo "no new component lines (still $old)" >&2; exit 1; }
    sbatch --array="$((old + 1))-$new%20" "$HERE/components_score.sbatch"
    ;;
collect)
    module load anaconda/colsa
    cd "$RUNS"
    for f in lineage410_metrics.csv components_metrics.csv; do
        [ -f "${f%.csv}.50.csv" ] || cp "$f" "${f%.csv}.50.csv"
    done
    python3 "$HERE/collect_lineage410.py" -o lineage410_metrics.csv
    python3 "$HERE/collect_components.py" -o components_metrics.csv
    ;;
*)
    sed -n '2,32p' "$0"
    exit 1
    ;;
esac
