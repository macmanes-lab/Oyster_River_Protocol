#!/bin/bash
# Watch the dev26 MPI runs on Premise from a laptop (NOTES.md 2026-10-12).
#
#     experiments/transabyss_mpi/watch_dev26.sh          # status of every run
#     experiments/transabyss_mpi/watch_dev26.sh --report # + balance.py and qualreports
#
# Needs `ssh premise` to work non-interactively (a Host premise entry in
# ~/.ssh/config and its key loaded). Status: slurm state, elapsed, node, the
# banner's version, and the run's latest assembler-stage lines -- including
# whether phase 2 grew. --report (automatic once every run has finished)
# copies the finished runs' timing and slurm logs here, runs balance.py
# locally (Premise's python3 is 3.6) for TOTAL against dev24 and "if
# balanced", and puts each qualreport's metrics beside the dev24 run's.
#
# RUNS: slurm job id and --runout, as submitted 2026-10-12 06:4x from ~ on
# Premise with oyster_mpi.sbatch (so each run's dir is ~/<run>, its slurm
# log ~/orp_ta_mpi_<job>.log). Edit to watch other runs.

set -uo pipefail
HOST="${HOST:-premise}"
RUNS=(
    1321775:SRR802670_dev26
    1321776:SRR868662_dev26
    1321777:SRR1139191_dev26
    1321778:SRR1138707_dev26
    1321779:SRR796576_dev26
)
# WATCH_RUNS="job:run job:run" replaces the list; REF_FROM/REF_TO turn a
# run's name into its dev24 counterpart's (default _dev26 -> _dev24).
[ -n "${WATCH_RUNS:-}" ] && read -r -a RUNS <<< "$WATCH_RUNS"
REF_FROM="${REF_FROM:-_dev26}"
REF_TO="${REF_TO:-_dev24}"
HERE="$(cd "$(dirname "$0")" && pwd)"
LOCAL="${LOCAL:-${TMPDIR:-/tmp}/watch_dev26}"
REPORT=0
[ "${1:-}" = --report ] && REPORT=1

banner_filter() { grep -v '^\*\* \|UNAUTHORIZED\|^By accessing\|transmitted through\|and will be property\|consent to\|erasing, or removing\|status of your user\|^otherwise\|University of New Hampshire\|^$' || true; }

remote=$(cat <<'EOF'
for spec in RUNS_HERE; do
    job=${spec%%:*}; run=${spec##*:}
    log=~/orp_ta_mpi_$job.log
    state=$(sacct -j "$job" -X -n --format=State,Elapsed,NodeList 2>/dev/null | head -1 | xargs)
    timing=~/$run/reports/$run.timing.log
    total=$(grep -m1 '^TOTAL' "$timing" 2>/dev/null | awk '{print $2}')
    echo "== $run (job $job): ${state:-unknown}${total:+, TOTAL $total}"
    [ -f "$log" ] || { echo "   no slurm log yet"; continue; }
    grep -a -m1 -o 'ORP[^ ]* [Vv]ersion [^ ]*\|version [0-9][^ ]*' "$log" | head -1 | sed 's/^/   /'
    grep -a '^\[assemblers\]\|=== Assemblers\|=== Stage [AB]\|=== Assemblers done\|stopped at .* relaunching\|lane failed\|\*\*\* step failed' "$log" \
        | sed 's/^/   /' | cut -c1-170 | tail -7
    [ -n "$total" ] || echo "   last: $(grep -a '^=== ' "$log" | tail -1 | cut -c1-150)"
done
EOF
)
remote=${remote//RUNS_HERE/${RUNS[*]}}
ssh -o BatchMode=yes "$HOST" "$remote" 2>&1 | banner_filter

finished=()
pending=0
states=$(ssh -o BatchMode=yes "$HOST" "for s in ${RUNS[*]}; do r=\${s##*:}; grep -q '^TOTAL' ~/\$r/reports/\$r.timing.log 2>/dev/null && echo \"done \$s\" || echo \"open \$s\"; done" 2>/dev/null | banner_filter)
while read -r what spec; do
    case "$what" in done) finished+=("$spec") ;; open) pending=$((pending + 1)) ;; esac
done <<< "$states"
echo
echo "${#finished[@]} of ${#RUNS[@]} runs finished"
[ "$pending" = 0 ] && REPORT=1
[ "$REPORT" = 1 ] && [ "${#finished[@]}" -gt 0 ] || exit 0

mkdir -p "$LOCAL"
args=()
for spec in "${finished[@]}"; do
    job=${spec%%:*}; run=${spec##*:}
    scp -q "$HOST:$run/reports/$run.timing.log" "$HOST:orp_ta_mpi_$job.log" "$LOCAL/" 2>/dev/null
    scp -q "$HOST:$run/reports/qualreport.$run" "$LOCAL/" 2>/dev/null
    ref=${run%"$REF_FROM"}$REF_TO
    scp -q "$HOST:$ref/reports/qualreport.$ref" "$LOCAL/" 2>/dev/null
    args+=("$LOCAL/$run.timing.log:$LOCAL/orp_ta_mpi_$job.log")
done
echo
echo "---------------- balance.py (logs copied to $LOCAL) ----------------"
python3 "$HERE/balance.py" "${args[@]}" | sed '/^fit over/,$d'
echo "---------------- qualreports: dev24 -> dev26 ----------------"
for spec in "${finished[@]}"; do
    run=${spec##*:}; ref=${run%"$REF_FROM"}$REF_TO
    echo "== $run vs $ref"
    paste -d'|' \
        <(grep '~~~>' "$LOCAL/qualreport.$ref" 2>/dev/null | sed 's/\*\*\*\*\* *//; s/ *~*> */\t/') \
        <(grep '~~~>' "$LOCAL/qualreport.$run" 2>/dev/null | sed 's/.*~*> *//') \
        | awk -F'[\t|]' '{printf "   %-30s %-48s %s%s\n", $1, $2, $3, ($2 == $3 ? "" : "   *")}'
done
