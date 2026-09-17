#!/bin/bash
# Submit the ORP array. Run this instead of calling sbatch by hand.
#
# Everything the job needs is defined here and passed through --export, so there
# is one definition of RUNS/MANIFEST/ORP rather than one in the submitting shell
# and another in the job. Nothing needs to be exported before running this.
#
# It also creates $RUNS/logs BEFORE submitting. That is not a convenience: slurm
# opens the --output file when the task launches, before the job script runs, so
# a missing log directory kills every task instantly with no log to say so.

set -euo pipefail

# The manifest sits beside the run tree, not inside it, so its default is keyed
# on the compare directory rather than on RUNS -- orp_runs/ does not exist until
# this script creates it.
COMPARE="${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}"
RUNS="${RUNS:-$COMPARE/orp_runs}"
MANIFEST="${MANIFEST:-$COMPARE/manifest.tsv}"
ORP="${ORP:-$HOME/Oyster_River_Protocol/oyster.py}"
THROTTLE="${THROTTLE:-6}"

die() { echo "submit: $*" >&2; exit 1; }

# Absolute, so a message naming this directory is one you can paste into rsync.
HERE=$(cd "$(dirname "$0")" && pwd)

[ -f "$MANIFEST" ] || die "no manifest at $MANIFEST"
[ -x "$ORP" ] || die "$ORP is not executable"

# Catch a half-copied set of scripts. This file and orp_array.sbatch are only
# ever updated together, so a mismatch here means one of them did not make it
# across -- which otherwise shows up as every task misparsing its manifest line.
jobcols=$(awk -F= '/^MANIFEST_COLUMNS=/ {print $2; exit}' "$HERE/orp_array.sbatch")
[ "${jobcols:-}" = "5" ] || die \
    "$HERE/orp_array.sbatch expects ${jobcols:-an older layout}, this script writes 5."$'\n'\
"The two are out of step: re-copy the whole of $HERE to this machine."

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

mkdir -p "$RUNS/logs" || die "cannot create $RUNS/logs"
[ -w "$RUNS/logs" ] || die "$RUNS/logs is not writable"

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

echo "submit: $n samples, $THROTTLE at a time, ORP $version on $branch"
echo "submit: manifest $MANIFEST"
echo "submit: runs -> $RUNS"

sbatch --array="1-${n}%${THROTTLE}" \
       --export="ALL,RUNS=$RUNS,MANIFEST=$MANIFEST,ORP=$ORP" \
       --output="$RUNS/logs/orp_%A_%a.log" \
       "$HERE/orp_array.sbatch"
