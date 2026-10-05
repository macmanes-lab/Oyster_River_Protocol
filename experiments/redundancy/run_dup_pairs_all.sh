#!/usr/bin/env bash
# Run dup_busco_pairs.py over every manifest sample that has both an ORP
# re-assembly and an original BUSCO run, once for each, and stack the
# summaries into one table.
#
# usage: run_dup_pairs_all.sh [OUTDIR]     (run inside the orp env, or with
#                                           it on PATH -- needs blastn)
#
# Paths follow experiments/compare108: the manifest's column 1 is the sample
# directory under pairs/, column 5 the ORP run name under orp_runs/.
set -euo pipefail

COMPARE=${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}
OUT=${1:-$PWD/dup_pairs}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"

while IFS=$'\t' read -r code _asm _r1 _r2 run; do
    new_table=$(ls "$COMPARE"/orp_runs/"$run"/reports/run_"$run".ORP/run_*/full_table.tsv 2>/dev/null | head -1 || true)
    new_fasta=$COMPARE/orp_runs/$run/assemblies/$run.ORP.fasta
    new_csv=$COMPARE/orp_runs/$run/reports/pytransrate_$run/contigs.csv
    [[ -s $new_csv ]] || new_csv=$COMPARE/orp_runs/$run/reports/transrate_$run/contigs.csv  # before 4.1.0-dev9
    old_table=$COMPARE/busco/busco_out/$code/full_table.tsv
    old_fasta=$(ls "$COMPARE"/pairs/"$code"/assembly/*.fasta* 2>/dev/null | head -1 || true)
    if [[ -z $new_table || ! -s $new_fasta || ! -s $old_table || -z $old_fasta ]]; then
        echo "skip $code ($run): missing ORP or original BUSCO/assembly" >&2
        continue
    fi
    csv_arg=()
    [[ -s $new_csv ]] && csv_arg=(--contigs-csv "$new_csv")
    python "$HERE"/dup_busco_pairs.py --full-table "$new_table" --fasta "$new_fasta" \
        "${csv_arg[@]}" --label "$code:orp" --out-prefix "$OUT/$code.orp" >/dev/null \
        || echo "FAILED $code orp" >&2
    python "$HERE"/dup_busco_pairs.py --full-table "$old_table" --fasta "$old_fasta" \
        --label "$code:original" --out-prefix "$OUT/$code.original" >/dev/null \
        || echo "FAILED $code original" >&2
done < "$COMPARE/manifest.tsv"

# One header, then every summary row; same for the per-pair rows.
for kind in summary pairs; do
    files=("$OUT"/*."$kind".tsv)
    { head -1 "${files[0]}"; tail -q -n +2 "${files[@]}"; } > "$OUT/all.$kind.tsv"
done
echo "wrote $OUT/all.summary.tsv ($(($(wc -l < "$OUT/all.summary.tsv") - 1)) rows)"
