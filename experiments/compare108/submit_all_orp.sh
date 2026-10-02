#!/bin/bash
# Submit the ORP array. Run this instead of calling sbatch by hand.
#
#     ./submit_all_orp.sh <manifest.tsv> <runs directory> [throttle]
#
#     ./submit_all_orp.sh /mnt/.../compare/manifest.tsv /mnt/.../compare/orp_runs
#     ./submit_all_orp.sh manifest.tsv /scratch/orp_runs 12
#
# To see what would be submitted and then submit only some of it:
#
#     ./submit_all_orp.sh --dry-run manifest.tsv orp_runs
#     ./submit_all_orp.sh --only SRR527277,SRR1234567 manifest.tsv orp_runs
#
# The runs directory is created if it does not exist; the manifest must already,
# and comes from make_manifest.sh.
#
# Both paths are passed into the job through --export, so there is one
# definition of them rather than one here and another in the job script.
#
# It also creates <runs>/logs BEFORE submitting. That is not a convenience:
# slurm opens the --output file when the task launches, before the job script
# runs, so a missing log directory kills every task instantly with no log to
# say so.

set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
usage: submit_all_orp.sh [options] <manifest.tsv> <runs directory> [throttle]

  manifest.tsv     5-column manifest from make_manifest.sh
  runs directory   where each sample's run tree goes; created if absent
  throttle         max concurrent array tasks (default 2)

  -f, --force      submit every sample, including ones already finished
                   (default: samples with reports/qualreport.<run>.done are
                   left out of the array)
  -n, --dry-run    run every check, list the samples that would be submitted
                   and the sbatch command, then stop; nothing is created or
                   submitted
  --only NAMES     submit only these samples: comma-separated run names
                   (column 5) or TSA names (column 1); repeatable. Combine
                   with --dry-run to check the selection first

  ORP=/path/to/oyster.py   override the pipeline used
                           (default: $HOME/Oyster_River_Protocol/oyster.py)
USAGE
    exit 2
}

die() { echo "submit: $*" >&2; exit 1; }

# Kept for the stale-script message below, which repeats the whole command.
ARGS="$*"

# Options may come before, between or after the paths; everything that is not
# one is collected, in order, as a positional argument.
FORCE="" DRYRUN="" ONLY="" pos=()
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)  usage ;;
        -f|--force) FORCE=1; shift ;;
        -n|--dry-run) DRYRUN=1; shift ;;
        --only)     [ $# -ge 2 ] && [ -n "$2" ] || die "--only needs a list of names"
                    ONLY="${ONLY:+$ONLY,}$2"; shift 2 ;;
        --only=*)   [ -n "${1#--only=}" ] || die "--only needs a list of names"
                    ONLY="${ONLY:+$ONLY,}${1#--only=}"; shift ;;
        --)         shift; pos+=("$@"); break ;;
        -*)         die "unknown option: $1" ;;
        *)          pos+=("$1"); shift ;;
    esac
done
# The ${pos[@]+...} form because bash before 4.4 calls an empty array unbound.
set -- ${pos[@]+"${pos[@]}"}

[ $# -ge 2 ] && [ $# -le 3 ] || usage

# Absolute, so that what is handed to the job and printed in messages does not
# depend on where this was run from -- a task resolves neither against its own
# cwd, which is its run directory.
MANIFEST=$(cd "$(dirname "$1")" 2>/dev/null && pwd)/$(basename "$1") \
    || die "no such directory for manifest: $1"
RUNS="$2"
case "$RUNS" in /*) ;; *) RUNS="$PWD/$RUNS" ;; esac
THROTTLE="${3:-2}"
ORP="${ORP:-$HOME/Oyster_River_Protocol/oyster.py}"

case "$THROTTLE" in
    ''|*[!0-9]*) die "throttle must be a number, got '$THROTTLE'" ;;
esac
[ "$THROTTLE" -gt 0 ] || die "throttle must be at least 1"

# Absolute, so a message naming this directory is one you can paste into rsync.
HERE=$(cd "$(dirname "$0")" && pwd)

[ -f "$MANIFEST" ] || die "no manifest at $MANIFEST"
[ -x "$ORP" ] || die "$ORP is not executable"

# Catch a half-copied set of scripts. This file and orp_array.sbatch are only
# ever updated together, so a mismatch here means one of them did not make it
# across -- which otherwise shows up as every task misparsing its manifest line.
jobcols=$(awk -F= '/^MANIFEST_COLUMNS=/ {print $2; exit}' "$HERE/orp_array.sbatch")
[ "${jobcols:-}" = "5" ] || die \
"$HERE/orp_array.sbatch is from an older version than this script"$'\n'\
"(it expects ${jobcols:-an older manifest layout}; this script hands over 5 columns)."$'\n'\
"The two must come from the same checkout. Rather than copying them, run this"$'\n'\
"script from the checkout -- it works from any directory, and relative paths are"$'\n'\
"taken from where you are:"$'\n'\
"    \$HOME/Oyster_River_Protocol/experiments/compare108/submit_all_orp.sh $ARGS"

n=$(grep -c . "$MANIFEST")
[ "$n" -gt 0 ] || die "manifest $MANIFEST has no non-blank lines"

# The columns the job relies on, checked once here rather than 108 times at 24
# cpus apiece: a manifest that is space- rather than tab-separated, or that has
# the assembly column missing, otherwise fails one task at a time.
bad=$(awk -F'\t' 'NF && NF != 5 {print NR": "NF" fields"}' "$MANIFEST" | head -5)
# A count other than 5 means the manifest and these scripts are different
# vintages, and it can be either one that is behind -- say so both ways rather
# than sending someone to regenerate a manifest that was already correct.
[ -z "$bad" ] || die "expected 5 tab-separated fields per line, column 5 being the"$'\n'\
"run name. Lines below have a different count, so manifest and scripts are out of"$'\n'\
"step: either regenerate it with make_manifest.sh, or re-copy all of"$'\n'\
"$HERE to this machine if the manifest is the newer of the two."$'\n'"$bad"

dupes=$(cut -f5 "$MANIFEST" | grep . | sort | uniq -d | head -5)
[ -z "$dupes" ] || die "run names in column 5 are not unique; these samples would"$'\n'"share a run directory:"$'\n'"$dupes"

missing=$(awk -F'\t' 'NF {print $3"\n"$4}' "$MANIFEST" | while read -r f; do
    [ -f "$f" ] || echo "$f"
done | head -5)
[ -z "$missing" ] || die "read files missing, e.g.:"$'\n'"$missing"

# A dry run creates nothing, so the runs directory is resolved only if it is
# already there; otherwise the path is shown as given, made absolute above.
if [ -z "$DRYRUN" ]; then
    mkdir -p "$RUNS/logs" || die "cannot create $RUNS/logs"
    [ -w "$RUNS/logs" ] || die "$RUNS/logs is not writable"
fi
# Now that it exists it can be resolved properly: ../orp_runs/ becomes a clean
# absolute path, rather than one carrying the .. and a doubled slash into every
# log path and every task's --dir.
[ -d "$RUNS" ] && RUNS=$(cd "$RUNS" && pwd)

branch=$(git -C "$(dirname "$ORP")" rev-parse --abbrev-ref HEAD 2>/dev/null || echo "?")
version=$(cat "$(dirname "$ORP")/version.txt" 2>/dev/null || echo "?")
# byo-assemblies merged to master at 4.0.0, so master is the expected branch
# now; the two feature branches stay accepted because a checkout sitting on
# either of them still has the pytransrate code these metrics are about.
case "$branch" in
    master|byo-assemblies|pytransrate) ;;
    *) echo "submit: WARNING ORP checkout is on '$branch' ($version); pytransrate" \
            "metrics need master (4.0.0+), byo-assemblies or pytransrate" >&2 ;;
esac

# Which manifest lines still need running. The array is given an explicit index
# list rather than 1-n, so a finished sample is never scheduled at all: the job
# script's own qualreport guard would exit in a second, but only after slurm had
# granted it 24 cpus and 120G, and with the throttle at 6 those are six slots
# that a sample still waiting to run could have had.
#
# Indices count non-blank lines, exactly as the job's own `grep . | sed -n Np`
# does, so the two agree about which line task N is.
#
# --only narrows this to the named samples first; a finished one that was named
# is still left out without --force, and said so, rather than silently rerun.
todo="" ndone=0 nskip=0 i=0 matched="," listing=""
# Tabs become \037 before read: a tab in IFS is whitespace, so read would fold
# the empty column 2 of a sample with no assembly into its neighbour and find
# column 5 empty -- dropping that sample from the array without a word.
while IFS=$'\037' read -r tsa _ _ _ run; do
    i=$((i + 1))
    [ -n "$run" ] || continue
    if [ -n "$ONLY" ]; then
        case ",$ONLY," in
            *",$run,"*|*",$tsa,"*) matched="$matched$run,$tsa," ;;
            *) nskip=$((nskip + 1)); continue ;;
        esac
    fi
    if [ -z "$FORCE" ] && [ -f "$RUNS/$run/reports/qualreport.$run.done" ]; then
        ndone=$((ndone + 1))
        [ -n "$ONLY" ] && echo "submit: $run is already complete; --force to run it again" >&2
        continue
    fi
    todo="${todo:+$todo,}$i"
    listing="$listing$(printf '%5s  %-16s  %s' "$i" "$run" "$tsa")"$'\n'
done < <(grep . "$MANIFEST" | tr '\t' '\037')

# A misspelt name would otherwise just shrink the array without a word.
if [ -n "$ONLY" ]; then
    unknown=$(printf '%s\n' "${ONLY//,/$'\n'}" | grep . | while read -r name; do
        case "$matched" in (*",$name,"*) ;; (*) echo "    $name" ;; esac
    done)
    [ -z "$unknown" ] || die "--only names not in column 5 or column 1 of the manifest:"$'\n'"$unknown"
fi

nrun=$(( n - ndone - nskip ))
if [ -z "$todo" ]; then
    if [ -n "$ONLY" ]; then
        echo "submit: every sample named by --only is already complete, nothing to submit"
    else
        echo "submit: all $n samples already complete in $RUNS, nothing to submit"
    fi
    echo "submit: pass --force to run them again"
    exit 0
fi

echo "submit: $nrun of $n samples to run, $THROTTLE at a time, ORP $version on $branch"
[ "$nskip" -gt 0 ] && echo "submit: $nskip not named by --only"
[ "$ndone" -gt 0 ] && echo "submit: skipping $ndone already complete (--force overrides)"
echo "submit: manifest $MANIFEST"
echo "submit: runs -> $RUNS"

cmd=(sbatch --array="${todo}%${THROTTLE}"
            --export="ALL,RUNS=$RUNS,MANIFEST=$MANIFEST,ORP=$ORP"
            --output="$RUNS/logs/orp_%A_%a.log"
            "$HERE/orp_array.sbatch")

if [ -n "$DRYRUN" ]; then
    # The task number is what the log will be named by: orp_<jobid>_<task>.log.
    echo
    printf '%5s  %-16s  %s\n' task run tsa
    printf '%s' "$listing"
    echo
    echo "submit: dry run, nothing submitted. Would run:"
    printf '    %s \\\n' "${cmd[@]:0:${#cmd[@]}-1}"
    printf '    %s\n' "${cmd[${#cmd[@]}-1]}"
    [ -d "$RUNS/logs" ] || echo "submit: and would first create $RUNS/logs"
    echo "submit: to submit some of these, rerun without --dry-run and with"
    echo "submit:     --only <run>,<run>,..."
    exit 0
fi

"${cmd[@]}"
