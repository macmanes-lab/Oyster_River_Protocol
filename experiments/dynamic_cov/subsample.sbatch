#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J tuco_sub
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 64G
#SBATCH --array=1-5%3
#SBATCH --output=tuco_sub_%A_%a.log
#SBATCH --exclude=node117,node118
# How --min-cov should move with depth: salmon over the tuco merge's own
# pre-filter index, on random subsamples of its trimmed + corrected reads (the
# reads ORP's filter quant uses; never normalized ones). Nothing is
# re-assembled, so every depth scores the same contigs, and depth_curve.py
# can read off where the no-hit contigs stop looking like noise.
#
#     cd ~/dyncov && sbatch ~/Oyster_River_Protocol/experiments/dynamic_cov/subsample.sbatch
#
# Task i quantifies DEPTHS[i] million read pairs, drawn with seqtk sample at
# one seed. One seed means R1 and R2 keep the same pairs, and each depth's
# pairs are a subset of the next one's. The full library (480,093,281 pairs)
# is the chowder run's own quant, already in $ROOT/tuco_cov. Each task
# writes $ROOT/tuco_d<NNN>_cov laid out like tuco_cov, so preview.py and
# depth_curve.py read it as is:
#   assemblies/   links to the merge's intermediate and its diamond hits
#   nohit_orfs.tsv  link to tuco_cov's (ORFs depend only on the contigs)
#   reports/tuco_d<NNN>.quant.sf.gz, reports/tuco_d<NNN>.filter_candidates.tsv
#                 (just its '# frag_length_mean' line)
#   depth.txt     salmon's num_processed: the read pairs it was given
# Then: conda run -n orp python experiments/dynamic_cov/depth_curve.py ~/dyncov \
#           tuco_d015 tuco_d030 tuco_d060 tuco_d120 tuco_d240 tuco

set -euo pipefail
DEPTHS=(15 30 60 120 240)
FULL_PAIRS=480093281
SEED="${SEED:-11}"
ROOT="${ROOT:-$HOME/dyncov}"
TUCO="${TUCO:-$HOME/tuco/assembly}"
RUN=380C_0C5D_001F_chowder
BIN="${BIN:-$HOME/orp_envs/orp/bin}"

d=${DEPTHS[$(( SLURM_ARRAY_TASK_ID - 1 ))]}
name=$(printf "tuco_d%03d" "$d")
frac=$(awk -v d="$d" -v n="$FULL_PAIRS" 'BEGIN { printf "%.6f", d * 1e6 / n }')
CPU="${SLURM_CPUS_PER_TASK:-24}"
DIR="$ROOT/${name}_cov"
echo "$name: ${d}M of $FULL_PAIRS pairs (fraction $frac, seed $SEED) -> $DIR -- $(date '+%F %T')"

mkdir -p "$DIR/assemblies" "$DIR/reports"
ln -sfn "$TUCO/assemblies/$RUN.ORP.intermediate.fasta" "$DIR/assemblies/$name.ORP.intermediate.fasta"
ln -sfn "$TUCO/assemblies/$RUN.ORP.diamond.txt" "$DIR/assemblies/$name.ORP.diamond.txt"
ln -sfn "$ROOT/tuco_cov/nohit_orfs.tsv" "$DIR/nohit_orfs.tsv"

fifo=$(mktemp -d "${TMPDIR:-/tmp}/tuco_sub.XXXXXX")
trap 'rm -rf "$fifo"' EXIT
mkfifo "$fifo/r1.fq" "$fifo/r2.fq"
"$BIN/seqtk" sample -s "$SEED" "$TUCO/rcorr/$RUN.TRIM_1P.cor.fq" "$frac" > "$fifo/r1.fq" &
pid1=$!
"$BIN/seqtk" sample -s "$SEED" "$TUCO/rcorr/$RUN.TRIM_2P.cor.fq" "$frac" > "$fifo/r2.fq" &
pid2=$!

# ORP's own salmon quant (oyster.py, Oyster.salmon), on the subsample.
rm -rf "$DIR/quant"
"$BIN/salmon" quant --no-version-check \
    -p "$CPU" -i "$TUCO/quants/$RUN.shucked.idx" --seqBias --gcBias --libType A \
    -1 "$fifo/r1.fq" -2 "$fifo/r2.fq" -o "$DIR/quant"
# A seqtk that died early would hand salmon a short read stream and a
# plausible-looking quant; its exit status is the only sign.
wait "$pid1"
wait "$pid2"

meta="$DIR/quant/aux_info/meta_info.json"
python3 - "$meta" "$DIR" "$name" <<'EOF'
import json, sys
meta, out, name = sys.argv[1:]
m = json.load(open(meta))
open(f"{out}/depth.txt", "w").write(f"{m['num_processed']}\n")
with open(f"{out}/reports/{name}.filter_candidates.tsv", "w") as f:
    f.write(f"# frag_length_mean\t{m['frag_length_mean']}\n")
    f.write(f"# num_processed\t{m['num_processed']}\n# percent_mapped\t{m['percent_mapped']}\n")
print(f"{name}: {m['num_processed']} pairs, {m['percent_mapped']:.2f}% mapped, "
      f"fragment {m['frag_length_mean']:.1f}")
EOF
gzip -c "$DIR/quant/quant.sf" > "$DIR/reports/$name.quant.sf.gz"
echo "$name done -- $(date '+%F %T')"
