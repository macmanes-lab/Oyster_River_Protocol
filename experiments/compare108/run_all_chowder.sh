#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J chowder
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
# Re-merge finished ORP runs with chowder.py, to test a merge-side setting
# (OrthoFinder, the pick rule, rescue, cd-hit-est) without reassembling.
#
#     ./run_all_chowder.sh [options] <manifest.tsv> <orp runs dir> <out dir> [throttle] [-- chowder args...]
#
#     # the baseline: chowder with its defaults
#     ./run_all_chowder.sh --only SRR807360,SRR651040 \
#         manifest.tsv orp_runs chowder/baseline
#     # the same samples with the setting under test
#     ./run_all_chowder.sh --only SRR807360,SRR651040 \
#         manifest.tsv orp_runs chowder/blastn -- --orthofinder-program blastn
#
# Each sample's inputs come from the ORP run that orp_array.sbatch left in
# <orp runs dir>/<SRR>/: the four raw assemblies under assemblies/ and the
# trimmed, corrected reads under rcorr/. Nothing there is written to. Output
# goes to <out dir>/<SRR>/, laid out the way oyster.py lays out a run, so
# collect_metrics.py --runs <out dir> reads it unchanged.
#
# Compare a setting against a chowder baseline, not against the ORP run. chowder
# uses one assembly order for both concatenation and the diamond rescue, while
# oyster.py uses ASSEMBLY_ORDER for the first and DIAMOND_PRIORITY for the
# second, and every contig carries a <label>_ prefix. So a default chowder merge
# is close to the ORP's own .ORP.fasta but not identical, and only two chowder
# runs differ by the setting alone.
#
# One out dir holds one setting. The chowder args are written to
# <out dir>/chowder.args on the first submit, and a later submit to the same
# directory with different args is refused, so a directory can never hold
# samples merged two different ways.
#
# This file is both the submitter and the array job: run by hand it checks
# everything and submits itself with sbatch; as an array task (CHOWDER_TASK=1,
# set by the submit) it runs one manifest line. The #SBATCH lines above apply to
# the task.

set -euo pipefail

# The four assemblies, in oyster.py's ASSEMBLY_ORDER, with the names oyster.py
# gave their files: <SRR>.<file>.fasta(.gz) under assemblies/.
LABELS=(spades55 spades75 transabyss trinity)
FILES=(spades55 spades75 transabyss trinity.Trinity)

# The first of <base>.gz and <base> that exists and is not empty, or nothing.
# The ORP gzips as it goes, but an interrupted run can leave either.
find_input() {
    if [ -s "$1.gz" ]; then echo "$1.gz"
    elif [ -s "$1" ]; then echo "$1"
    fi
}

# Every input of one sample, one per line: four assemblies, then R1 and R2.
# Prints "MISSING <path>" in place of any that is absent.
sample_inputs() {
    local src="$1/$2" f base path
    for f in "${FILES[@]}" rcorr1 rcorr2; do
        case "$f" in
            rcorr1) base="$src/rcorr/$2.TRIM_1P.cor.fq" ;;
            rcorr2) base="$src/rcorr/$2.TRIM_2P.cor.fq" ;;
            *)      base="$src/assemblies/$2.$f.fasta" ;;
        esac
        path=$(find_input "$base")
        echo "${path:-MISSING $base}"
    done
}

###############################################################################
# Array task: one manifest line.
###############################################################################
if [ "${CHOWDER_TASK:-}" = 1 ]; then
    : "${MANIFEST:?} ${ORP_RUNS:?} ${OUT:?} ${CHOWDER:?}"
    echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="

    module load anaconda/colsa

    line=$(grep . "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
    [ -n "$line" ] || { echo "manifest has no line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
    # \037 rather than tab in IFS: see orp_array.sbatch.
    IFS=$'\037' read -r TSA _ _ _ SRR EXTRA <<< "$(printf '%s' "$line" | tr '\t' '\037')"
    [ -z "${EXTRA:-}" ] || { echo "manifest line ${SLURM_ARRAY_TASK_ID} has more than 5 fields" >&2; exit 1; }
    [ -n "$SRR" ] || { echo "manifest line ${SLURM_ARRAY_TASK_ID} has no run name" >&2; exit 1; }

    mapfile -t inputs < <(sample_inputs "$ORP_RUNS" "$SRR")
    for f in "${inputs[@]}"; do
        case "$f" in MISSING*) echo "input ${f#MISSING }(.gz) not found" >&2; exit 1 ;; esac
    done

    # Read here rather than passed through --export, which splits on commas.
    mapfile -t extra < "$OUT/chowder.args"

    mkdir -p "$OUT/$SRR"
    ln -sfn "chowder_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$OUT/logs/${SRR}.log"
    ln -sfn "chowder_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" "$OUT/logs/${TSA}.log"
    cd "$OUT/$SRR"

    echo "=== sample: $SRR ($TSA) ==="
    echo "    from    $ORP_RUNS/$SRR"
    echo "    dir     $OUT/$SRR"
    echo "    chowder $CHOWDER ($(git -C "$(dirname "$CHOWDER")" log -1 --format='%h %s' 2>/dev/null || echo '?'))"
    echo "    args    $([ ${#extra[@]} -eq 0 ] || printf '%q ' "${extra[@]}")"

    if [ -f "$OUT/$SRR/reports/qualreport.$SRR.done" ]; then
        echo "already complete, skipping"
        exit 0
    fi

    # From the allocation, as in orp_array.sbatch.
    CPU="${SLURM_CPUS_PER_TASK:-24}"
    if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
        MEM=$(( SLURM_MEM_PER_NODE / 1024 ))
    elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
        MEM=$(( SLURM_MEM_PER_CPU * CPU / 1024 ))
    else
        MEM=120
    fi
    echo "    allocation: ${CPU} cpus, ${MEM}G"

    # --tpm-filt 1 and --max-parallel 2 as the ORP runs had them. The extra
    # args go last, so one of them can override anything set here.
    "$CHOWDER" \
        --assemblies "${inputs[@]:0:4}" \
        --labels "${LABELS[@]}" \
        --assembly-order given \
        --read1 "${inputs[4]}" \
        --read2 "${inputs[5]}" \
        --corrected-reads \
        --tpm-filt 1 \
        --mem "$MEM" \
        --cpu "$CPU" \
        --max-parallel 2 \
        --dir "$OUT/$SRR" \
        --runout "$SRR" \
        ${extra[@]+"${extra[@]}"}

    echo "=== done $SRR at $(date) ==="
    exit 0
fi

###############################################################################
# Submit.
###############################################################################
usage() {
    cat >&2 <<'USAGE'
usage: run_all_chowder.sh [options] <manifest.tsv> <orp runs dir> <out dir> [throttle] [-- chowder args...]

  manifest.tsv     5-column manifest from make_manifest.sh
  orp runs dir     finished ORP runs (orp_array.sbatch's <runs>), read only
  out dir          where each sample's chowder run goes; one setting per dir
  throttle         max concurrent array tasks (default 2)
  -- ARGS...       passed to every chowder.py run, after the script's own
                   arguments, so they can override them too

  -f, --force      submit samples already finished in <out dir> as well
  -n, --dry-run    run every check, list what would be submitted, then stop;
                   nothing is created or submitted
  --only NAMES     only these samples: comma-separated run names (column 5)
                   or TSA names (column 1); repeatable

  CHOWDER=/path/to/chowder.py   override the merge used
                                (default: $HOME/Oyster_River_Protocol/chowder.py)
USAGE
    exit 2
}

die() { echo "chowder-submit: $*" >&2; exit 1; }
say() { echo "chowder-submit: $*"; }

FORCE="" DRYRUN="" ONLY="" pos=() extra=()
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)    usage ;;
        -f|--force)   FORCE=1; shift ;;
        -n|--dry-run) DRYRUN=1; shift ;;
        --only)       [ $# -ge 2 ] && [ -n "$2" ] || die "--only needs a list of names"
                      ONLY="${ONLY:+$ONLY,}$2"; shift 2 ;;
        --only=*)     [ -n "${1#--only=}" ] || die "--only needs a list of names"
                      ONLY="${ONLY:+$ONLY,}${1#--only=}"; shift ;;
        --)           shift; extra=("$@"); break ;;
        -*)           die "unknown option: $1 (chowder args go after --)" ;;
        *)            pos+=("$1"); shift ;;
    esac
done
set -- ${pos[@]+"${pos[@]}"}
[ $# -ge 3 ] && [ $# -le 4 ] || usage

abspath() { case "$1" in /*) echo "$1" ;; *) echo "$PWD/$1" ;; esac; }
MANIFEST=$(abspath "$1")
ORP_RUNS=$(abspath "$2")
OUT=$(abspath "$3")
THROTTLE="${4:-2}"
CHOWDER="${CHOWDER:-$HOME/Oyster_River_Protocol/chowder.py}"
HERE=$(cd "$(dirname "$0")" && pwd)

case "$THROTTLE" in ''|*[!0-9]*|0) die "throttle must be a positive number, got '$THROTTLE'" ;; esac
[ -f "$MANIFEST" ] || die "no manifest at $MANIFEST"
[ -d "$ORP_RUNS" ] || die "no ORP runs directory at $ORP_RUNS"
ORP_RUNS=$(cd "$ORP_RUNS" && pwd)
[ -x "$CHOWDER" ] || die "$CHOWDER is not executable"
[ "$OUT" != "$ORP_RUNS" ] || die "out dir is the ORP runs dir; chowder would write into the ORP runs"

bad=$(awk -F'\t' 'NF && NF != 5 {print NR": "NF" fields"}' "$MANIFEST" | head -5)
[ -z "$bad" ] || die "expected 5 tab-separated fields per line:"$'\n'"$bad"
n=$(grep -c . "$MANIFEST")

# One setting per out dir. The args are compared one per line, so a quoted
# argument containing a space survives the round trip.
want=$(printf '%s\n' ${extra[@]+"${extra[@]}"})
# Shell-quoted, for messages: a quoted argument with a space in it shows as one.
shown=""
[ ${#extra[@]} -eq 0 ] || shown=$(printf '%q ' "${extra[@]}")
if [ -f "$OUT/chowder.args" ]; then
    have=$(cat "$OUT/chowder.args")
    [ "$have" = "$want" ] || die "$OUT was run with different chowder args:"$'\n'\
"    there:     ${have//$'\n'/ }"$'\n'"    requested: ${want//$'\n'/ }"$'\n'\
"Each setting needs its own out dir."
fi

checkout=$(dirname "$CHOWDER")
commit=$(git -C "$checkout" log -1 --format='%h %s' 2>/dev/null || echo "?")
branch=$(git -C "$checkout" rev-parse --abbrev-ref HEAD 2>/dev/null || echo "?")
# Tasks run whatever the checkout holds when each one starts, so editing it
# while the array runs would merge different samples with different code.
if [ -n "$(git -C "$checkout" status --porcelain --untracked-files=no 2>/dev/null)" ]; then
    echo "chowder-submit: WARNING $checkout has uncommitted changes; tasks will run them as they stand" >&2
fi

# Same index scheme as submit_all_orp.sh: count non-blank lines, as the task's
# `grep . | sed -n Np` does.
todo="" listing="" matched="," missing="" i=0 ndone=0 nskip=0 nmiss=0
while IFS=$'\037' read -r tsa _ _ _ run; do
    i=$((i + 1))
    [ -n "$run" ] || continue
    if [ -n "$ONLY" ]; then
        case ",$ONLY," in
            *",$run,"*|*",$tsa,"*) matched="$matched$run,$tsa," ;;
            *) nskip=$((nskip + 1)); continue ;;
        esac
    fi
    gaps=$(sample_inputs "$ORP_RUNS" "$run" | sed -n 's/^MISSING //p')
    if [ -n "$gaps" ]; then
        nmiss=$((nmiss + 1))
        missing="$missing    $run: $(echo "$gaps" | sed "s|$ORP_RUNS/$run/||" | tr '\n' ' ')"$'\n'
        continue
    fi
    if [ -z "$FORCE" ] && [ -f "$OUT/$run/reports/qualreport.$run.done" ]; then
        ndone=$((ndone + 1))
        [ -n "$ONLY" ] && echo "chowder-submit: $run is already complete in $OUT; --force to rerun" >&2
        continue
    fi
    todo="${todo:+$todo,}$i"
    listing="$listing$(printf '%5s  %-16s  %s' "$i" "$run" "$tsa")"$'\n'
done < <(grep . "$MANIFEST" | tr '\t' '\037')

if [ -n "$ONLY" ]; then
    unknown=$(printf '%s\n' "${ONLY//,/$'\n'}" | grep . | while read -r name; do
        case "$matched" in (*",$name,"*) ;; (*) echo "    $name" ;; esac
    done)
    [ -z "$unknown" ] || die "--only names not in column 5 or column 1 of the manifest:"$'\n'"$unknown"
fi

# A sample named by --only must be runnable; across the whole manifest, ones
# with no finished ORP run are expected and only reported.
if [ -n "$missing" ]; then
    [ -z "$ONLY" ] || die "inputs missing (.gz or plain) under $ORP_RUNS:"$'\n'"$missing"
    echo "chowder-submit: leaving out $nmiss with inputs missing under $ORP_RUNS:" >&2
    printf '%s' "$missing" >&2
fi

if [ -z "$todo" ]; then
    say "nothing to submit ($ndone already complete in $OUT; --force reruns them)"
    exit 0
fi

nrun=$(( n - ndone - nskip - nmiss ))
say "$nrun of $n samples, $THROTTLE at a time"
say "chowder  $CHOWDER on $branch ($commit)"
say "args     ${shown:-(none: chowder defaults)}"
say "inputs   $ORP_RUNS"
say "out      $OUT"
[ "$nskip" -gt 0 ] && say "$nskip not named by --only"
[ "$ndone" -gt 0 ] && say "skipping $ndone already complete (--force overrides)"

cmd=(sbatch --array="${todo}%${THROTTLE}"
            --job-name="chowder-$(basename "$OUT")"
            --export="ALL,CHOWDER_TASK=1,MANIFEST=$MANIFEST,ORP_RUNS=$ORP_RUNS,OUT=$OUT,CHOWDER=$CHOWDER"
            --output="$OUT/logs/chowder_%A_%a.log"
            "$HERE/run_all_chowder.sh")

if [ -n "$DRYRUN" ]; then
    echo
    printf '%5s  %-16s  %s\n' task run tsa
    printf '%s' "$listing"
    echo
    say "dry run, nothing submitted. Would run:"
    printf '    %s \\\n' "${cmd[@]:0:${#cmd[@]}-1}"
    printf '    %s\n' "${cmd[${#cmd[@]}-1]}"
    [ -f "$OUT/chowder.args" ] || say "and would first create $OUT and write its chowder.args"
    exit 0
fi

# Before sbatch: slurm opens the --output file before the task starts, and the
# task reads chowder.args.
mkdir -p "$OUT/logs" || die "cannot create $OUT/logs"
[ -f "$OUT/chowder.args" ] || printf '%s' "${want:+$want$'\n'}" > "$OUT/chowder.args"
echo "$(date '+%F %T')  $branch $commit  ${shown:-(defaults)}" >> "$OUT/chowder.submits"

"${cmd[@]}"
