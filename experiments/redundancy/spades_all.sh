#!/bin/bash
# Submit spades_all.sbatch over every ORP run that has both a Trans-ABySS and
# a Trinity assembly. Run this on the Premise login node, not sbatch directly.
#
#     ./spades_all.sh [-n|--dry-run] [-f|--force] [CODE_DIR] [THROTTLE]
#
#     CODE_DIR  ORP checkout to run (default: the one this script is in)
#     THROTTLE  array tasks at once (default 20)
#
# Scans $COMPARE/orp_runs/*/ (default COMPARE=/mnt/home/macmaneslab/macmanes/compare)
# for assemblies/<run>.transabyss.fasta[.gz], assemblies/<run>.trinity.Trinity.fasta[.gz]
# and the corrected reads rcorr/<run>.TRIM_{1,2}P.cor.fq[.gz]. A run missing
# any of them is reported and left out. A run whose
# $OUT_BASE/<run>/reports/qualreport.<run>_spadesbr.done already exists is
# left out too, unless --force.
#
# The runs to do go to $OUT_BASE/samples_<date>.tsv (run, transabyss,
# trinity, R1, R2), and the array is 1-N%THROTTLE over that file. Logs:
# $OUT_BASE/logs/spades_all_<job>_<task>.log, plus a <run>.log symlink.
# OUT_BASE defaults to $HOME/redundancy_tests/spades_all.
set -euo pipefail

DRYRUN="" FORCE="" pos=()
for a in "$@"; do
    case "$a" in
        -n|--dry-run) DRYRUN=1 ;;
        -f|--force)   FORCE=1 ;;
        -h|--help)    sed -n '2,21p' "$0"; exit 0 ;;
        -*)           echo "unknown option $a" >&2; exit 2 ;;
        *)            pos+=("$a") ;;
    esac
done
HERE=$(cd "$(dirname "$0")" && pwd)
CODE=$(cd "${pos[0]:-$HERE/../..}" && pwd)
THROTTLE=${pos[1]:-20}
COMPARE=${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}
OUT_BASE=${OUT_BASE:-$HOME/redundancy_tests/spades_all}

[[ -f $CODE/oyster.py && -f $CODE/chowder.py ]] || { echo "no oyster.py/chowder.py in $CODE" >&2; exit 1; }

# First existing file of the candidates, or nothing.
first() { local f; for f in "$@"; do [[ -s $f ]] && { echo "$f"; return 0; }; done; return 0; }

rows="" n=0 ndone=0 nskip=0
for d in "$COMPARE"/orp_runs/*/; do
    d=${d%/}
    run=$(basename "$d")
    [[ -d $d/assemblies ]] || continue
    ta=$(first "$d/assemblies/$run.transabyss.fasta.gz" "$d/assemblies/$run.transabyss.fasta")
    tr=$(first "$d/assemblies/$run.trinity.Trinity.fasta.gz" "$d/assemblies/$run.trinity.Trinity.fasta")
    [[ -n $ta && -n $tr ]] || continue      # not a dataset with both assemblies
    r1=$(first "$d/rcorr/$run.TRIM_1P.cor.fq.gz" "$d/rcorr/$run.TRIM_1P.cor.fq")
    r2=$(first "$d/rcorr/$run.TRIM_2P.cor.fq.gz" "$d/rcorr/$run.TRIM_2P.cor.fq")
    if [[ -z $r1 || -z $r2 ]]; then
        echo "SKIP $run: has both assemblies but no corrected reads in $d/rcorr" >&2
        nskip=$((nskip + 1)); continue
    fi
    if [[ -z $FORCE && -f $OUT_BASE/$run/reports/qualreport.${run}_spadesbr.done ]]; then
        ndone=$((ndone + 1)); continue
    fi
    rows+="$run"$'\t'"$ta"$'\t'"$tr"$'\t'"$r1"$'\t'"$r2"$'\n'
    n=$((n + 1))
done

echo "spades_all: $n to run, $ndone already done, $nskip without corrected reads"
echo "spades_all: ORP $(cat "$CODE/version.txt") at $(git -C "$CODE" rev-parse --short HEAD) in $CODE"
[[ $n -gt 0 ]] || exit 0

if [[ -n $DRYRUN ]]; then
    printf '%s' "$rows" | cut -f1 | paste - - - - - - | column -t
    echo "spades_all: dry run; would submit --array=1-$n%$THROTTLE to $OUT_BASE"
    exit 0
fi

mkdir -p "$OUT_BASE/logs"
SAMPLES=$OUT_BASE/samples_$(date +%Y%m%d_%H%M%S).tsv
printf '%s' "$rows" > "$SAMPLES"
echo "spades_all: sample list $SAMPLES"
sbatch --array="1-$n%$THROTTLE" \
       --output="$OUT_BASE/logs/spades_all_%A_%a.log" \
       --export="ALL,OUT_BASE=$OUT_BASE" \
       "$HERE/spades_all.sbatch" "$SAMPLES" "$CODE"
