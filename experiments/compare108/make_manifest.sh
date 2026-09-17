#!/bin/bash
# Build the manifest from the pairs/ tree.
#
# Walks $COMPARE/pairs/tsa_XXXX/ and writes one headerless tab-separated line
# per sample, in the layout the rest of these scripts read:
#
#     tsa_GADU <TAB> assembly.fasta.gz <TAB> SRR527277_1.fastq.gz <TAB> SRR527277_2.fastq.gz <TAB> SRR527277
#
#     ./make_manifest.sh                    # -> $COMPARE/manifest.tsv
#     OUT=/tmp/subset.tsv ./make_manifest.sh
#
# Column 5 is the run name: what --runout gets and what the run directory is
# called. It is settled here, once, rather than re-derived inside each array
# task, because only this script sees all the samples at once -- a task cannot
# tell that the SRR it just parsed out of its filename belongs to two samples.
#
# A sample with no usable read pair is left out and reported, because it cannot
# be run; a sample with no assembly is kept with column 2 empty, because
# oyster.py never reads that column -- only chowder.py does. Both go to stderr,
# so the file on stdout stays clean and the line count is what will actually run.

set -uo pipefail

COMPARE="${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}"
BASE="${BASE:-$COMPARE/pairs}"
OUT="${OUT:-$COMPARE/manifest.tsv}"

[ -d "$BASE" ] || { echo "make_manifest: no pairs directory at $BASE" >&2; exit 1; }

tmp=$(mktemp "${TMPDIR:-/tmp}/manifest.XXXXXX") || exit 1
trap 'rm -f "$tmp"' EXIT

found=0 skipped=0 noasm=0

for sampledir in "$BASE"/tsa_*; do
    [ -d "$sampledir" ] || continue
    tsa=$(basename "$sampledir")

    # R1 first: without a pair there is nothing to run, whatever else is there.
    # Several naming conventions are in play across these samples, so try each
    # and derive the mate from whichever matched, so the pair cannot disagree.
    r1=""
    for pat in '*_1.fastq.gz' '*_1.fq.gz' '*_R1.fastq.gz' '*_R1_*.fastq.gz' '*.1.fastq.gz'; do
        for c in "$sampledir"/reads/$pat; do
            [ -f "$c" ] && { r1="$c"; break 2; }
        done
    done
    if [ -z "$r1" ]; then
        # Say what is there, not just what is missing: these are usually either
        # an empty directory or a naming convention the patterns above miss, and
        # the two need opposite fixes.
        if [ ! -d "$sampledir/reads" ]; then
            what="no reads/ directory"
        else
            what=$(ls -A "$sampledir/reads" 2>/dev/null | head -6 | tr '\n' ' ')
            [ -n "$what" ] || what="(directory is empty)"
        fi
        echo "SKIP $tsa: no R1 in $sampledir/reads -- contains: $what" >&2
        skipped=$((skipped+1)); continue
    fi

    case "$r1" in
        *_1.fastq.gz)  r2="${r1%_1.fastq.gz}_2.fastq.gz" ;;
        *_1.fq.gz)     r2="${r1%_1.fq.gz}_2.fq.gz" ;;
        *_R1.fastq.gz) r2="${r1%_R1.fastq.gz}_R2.fastq.gz" ;;
        *.1.fastq.gz)  r2="${r1%.1.fastq.gz}.2.fastq.gz" ;;
        *_R1_*)        r2=$(echo "$r1" | sed 's/_R1_/_R2_/') ;;
        *)             r2="" ;;
    esac
    if [ -z "$r2" ] || [ ! -f "$r2" ]; then
        echo "SKIP $tsa: no mate for $(basename "$r1")" >&2
        skipped=$((skipped+1)); continue
    fi

    asm=""
    for c in "$sampledir"/assembly/*.fasta.gz "$sampledir"/assembly/*.fsa_nt.gz \
             "$sampledir"/assembly/*.fasta; do
        [ -f "$c" ] && { asm="$c"; break; }
    done
    if [ -z "$asm" ]; then
        echo "WARN $tsa: no assembly in $sampledir/assembly; column 2 left empty" \
             "(fine for oyster.py, not for chowder.py)" >&2
        noasm=$((noasm+1))
    fi

    srr=$(basename "$r1" | grep -o -m1 '[SED]RR[0-9]\{4,\}' | head -1)
    [ -n "$srr" ] || srr="$tsa"

    # Two signatures for the collision resolver, stripped back off below.
    # canon: the files the paths actually resolve to, so reads catalogued under
    # two accessions as symlinks or hardlinks are recognised as one set.
    # sizes: a weaker hint, for two independent copies of the same data.
    canon="$(readlink -f "$r1")|$(readlink -f "$r2")"
    sizes="$(wc -c < "$r1" | tr -d ' ')|$(wc -c < "$r2" | tr -d ' ')"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$tsa" "$asm" "$r1" "$r2" "$srr" "$canon" "$sizes" >> "$tmp"
    found=$((found+1))
done

# Resolve run-name collisions. Two samples can land on one SRR for two quite
# different reasons, and they want opposite treatment:
#
#   same SRR, same read files  -- one set of reads catalogued under two TSA
#       accessions. A full ORP run of each would read the same fastqs and
#       produce byte-identical output, so the second is dropped and the pairing
#       recorded in <manifest>.folded.tsv, so nothing about which accessions
#       share reads is lost.
#
#   same SRR, different read files -- distinct samples whose filenames happen to
#       carry the same accession. Both must run, so both get a run name
#       disambiguated with the tsa code: SRR516821_tsa_GBBZ.
#
folded="${OUT%.tsv}.folded.tsv"
: > "$folded"
awk -F'\t' -v OFS='\t' -v folded="$folded" '
    FNR == NR {
        key = $5 SUBSEP $6
        if (key in owner) {
            fold[FNR] = owner[key]          # same files as an earlier sample
        } else {
            owner[key] = $1
            distinct[$5]++                  # how many real samples share this name
        }
        if (!(($5 SUBSEP $7) in sizeseen)) { sizeseen[$5 SUBSEP $7] = 1; sizes[$5]++ }
        next
    }
    {
        if (FNR in fold) { print $1, $5, fold[FNR] >> folded; next }
        if (distinct[$5] > 1) {
            if (sizes[$5] == 1 && !(warned[$5]++))
                printf "WARN %s: %d samples share this accession with different\n" \
                       "     paths but identical file sizes -- possibly the same reads\n" \
                       "     copied twice. Both will run; check before spending the time.\n", \
                       $5, distinct[$5] > "/dev/stderr"
            $5 = $5 "_" $1
        }
        print $1, $2, $3, $4, $5
    }
' "$tmp" "$tmp" > "$tmp.resolved"
mv "$tmp.resolved" "$tmp"

nfolded=$(grep -c . "$folded" 2>/dev/null) || true
nfolded=${nfolded:-0}
if [ "$nfolded" -gt 0 ]; then
    echo "NOTE folded $nfolded sample(s) whose reads are the same files as another's;" \
         "see $folded (columns: dropped tsa, run name, tsa kept):" >&2
    sed 's/^/    /' "$folded" >&2
    found=$((found - nfolded))
else
    rm -f "$folded"
fi

renamed=$(awk -F'\t' 'index($5, "_tsa_") {print "    " $1 "  ->  " $5}' "$tmp")
if [ -n "$renamed" ]; then
    echo "NOTE these share an accession with another sample but have different reads," >&2
    echo "     so each runs under its own name:" >&2
    echo "$renamed" >&2
fi

# Belt and braces: after all that, the run names must be unique, because each
# one becomes a directory that a job writes into.
dupes=$(cut -f5 "$tmp" | sort | uniq -d)
if [ -n "$dupes" ]; then
    echo "ERROR run names still collide after resolution:" >&2
    echo "$dupes" | sed 's/^/    /' >&2
    exit 1
fi

# Only now overwrite whatever was there: a failed scan leaves the old manifest.
mv "$tmp" "$OUT"
trap - EXIT

echo "make_manifest: wrote $found samples to $OUT ($skipped skipped, $noasm without an assembly)" >&2
