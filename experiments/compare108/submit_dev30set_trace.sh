#!/bin/bash
# submit_dev30set_trace.sh - the two follow-ups to the 59-run handoff, run by
# hand from the Premise login node:
#
#   ./submit_dev30set_trace.sh submit    # once: both arrays
#   ./submit_dev30set_trace.sh collect   # after both finish
#
# submit:
#   orp_dev30set.sbatch        dev30 settings (--normalize-reads, --min-cov 2) on
#                              lineage manifest lines 51-59, the nine orp_
#                              libraries that have no dev30 run; ~10 h each
#   trace_merge_losses.sbatch  every BUSCO and Swiss-Prot gene a component
#                              assembly has and the ORP merge lost, all 59 runs
#
# collect:
#   merge_trace/merge_losses_summary.csv and merge_losses_detail.tsv
#   dev30_metrics.csv, dev30set_metrics.csv and settings_dev30_vs_400.csv in
#   lineage_410/, and the Section 4.4 numbers on stdout (also settings_summary.txt)
set -euo pipefail

COMPARE="${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}"
RUNS="$COMPARE/lineage_410"
HERE="${HERE:-$HOME/Oyster_River_Protocol/experiments/compare108}"

case "${1:-}" in
submit)
    [ "$(grep -c . "$RUNS/manifest.tsv")" -eq 59 ] || {
        echo "$RUNS/manifest.tsv does not have 59 lines; lines 51-59 must be the nine orp_ additions" >&2; exit 1; }
    mkdir -p "$COMPARE/dev30_settings/logs" "$RUNS/merge_trace/logs"
    sbatch "$HERE/orp_dev30set.sbatch"
    sbatch "$HERE/trace_merge_losses.sbatch"
    ;;
collect)
    module load anaconda/colsa
    cd "$RUNS"
    python3 "$HERE/trace_merge_losses.py" --collate
    # dev30 metrics for the first 50 (keep an existing table if there is one)
    [ -s dev30_metrics.csv ] || python3 "$HERE/collect_lineage410.py" \
        --runs '$COMPARE/orp_{group}' --runout-tag v4_1_0-dev30 -o dev30_metrics.csv
    python3 "$HERE/collect_lineage410.py" --runs '$COMPARE/dev30_settings' -o dev30set_metrics.csv
    python3 "$HERE/settings_compare.py" --release lineage410_metrics.csv \
        --dev30 dev30_metrics.csv --dev30set dev30set_metrics.csv \
        -o settings_dev30_vs_400.csv 2>&1 | tee settings_summary.txt
    ;;
*)
    sed -n '2,20p' "$0"
    exit 1
    ;;
esac
