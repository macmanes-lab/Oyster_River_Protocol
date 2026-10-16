#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J dyncov_sub
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 64G
#SBATCH --array=1-5%3
#SBATCH --output=%x_%A_%a.log
#SBATCH --exclude=node117,node118
# How --min-cov should move with depth: salmon over a finished run's own
# pre-filter index, on random subsamples of its trimmed + corrected reads (the
# reads ORP's filter quant uses; never normalized ones). Nothing is
# re-assembled, so every depth scores the same contigs, and depth_curve.py
# can read off where the no-hit contigs stop looking like noise.
#
# tuco (the defaults; 15-240M of 480M pairs):
#     cd ~/dyncov && sbatch -J tuco_sub ~/Oyster_River_Protocol/experiments/dynamic_cov/subsample.sbatch
# another run: SRC holds its assemblies/, quants/ and rcorr/, RUN is its
# --runout, NAME prefixes the output dirs; one array task per depth:
#     cd ~/dyncov && sbatch -J ERR489297_sub --array=1-4%3 \
#         --export=ALL,SRC=$HOME/dyncov/keep_ERR489297,RUN=ERR489297_v4_1_0-dev30,NAME=ERR489297,DEPTHS="15 30 60 120" \
#         ~/Oyster_River_Protocol/experiments/dynamic_cov/subsample.sbatch
#
# Task i quantifies DEPTHS[i] million read pairs, drawn with seqtk sample at
# one seed. One seed means R1 and R2 keep the same pairs, and each depth's
# pairs are a subset of the next one's. The full library is the run's own
# quant: $ROOT/<NAME>_cov, which must exist (tuco_cov does; set one up as
# the README says), and gives the read count the fractions are taken of
# (its depth.txt). salmon 2.7.0 can't read a pipe (it seeks: "Illegal
# seek"), so the subsample is written to $SCRATCH and deleted after. seqtk
# exits 0 on a failed write (task 5 of 1323242 lost the end of r1 on fast
# scratch), so both files are counted before salmon runs. Each task writes
# $ROOT/<NAME>_d<NNN>_cov, laid out like <NAME>_cov, so preview.py and
# depth_curve.py read it as is:
#   assemblies/   links to the run's intermediate and its diamond hits
#   nohit_orfs.tsv  link to <NAME>_cov's (ORFs depend only on the contigs)
#   reports/<NAME>_d<NNN>.quant.sf.gz, reports/<NAME>_d<NNN>.filter_candidates.tsv
#                 (its '# frag_length_mean' line, and num_processed)
#   depth.txt     salmon's num_processed: the read pairs it was given
# Then: depth_curve.py ~/dyncov <NAME>_d015 ... <NAME>

set -euo pipefail
SRC="${SRC:-$HOME/tuco/assembly}"
RUN="${RUN:-380C_0C5D_001F_chowder}"
NAME="${NAME:-tuco}"
read -r -a DEPTHS <<< "${DEPTHS:-15 30 60 120 240}"
SEED="${SEED:-11}"
ROOT="${ROOT:-$HOME/dyncov}"
BIN="${BIN:-$HOME/orp_envs/orp/bin}"
SCRATCH="${SCRATCH:-$HOME/scratch_dyncov_sub}"

read1="$SRC/rcorr/$RUN.TRIM_1P.cor.fq"
read2="$SRC/rcorr/$RUN.TRIM_2P.cor.fq"
[ -e "$read1" ] || { read1="$read1.gz"; read2="$read2.gz"; }
BASE="$ROOT/${NAME}_cov"
FULL_PAIRS=$(cat "$BASE/depth.txt")

d=${DEPTHS[$(( SLURM_ARRAY_TASK_ID - 1 ))]}
name=$(printf "%s_d%03d" "$NAME" "$d")
frac=$(awk -v d="$d" -v n="$FULL_PAIRS" 'BEGIN { printf "%.6f", d * 1e6 / n }')
awk -v f="$frac" 'BEGIN { exit !(f < 1) }' || { echo "$d M pairs is not below the full $FULL_PAIRS" >&2; exit 1; }
CPU="${SLURM_CPUS_PER_TASK:-24}"
DIR="$ROOT/${name}_cov"
echo "$name: ${d}M of $FULL_PAIRS pairs (fraction $frac, seed $SEED) from $read1 -> $DIR -- $(date '+%F %T')"

mkdir -p "$DIR/assemblies" "$DIR/reports"
ln -sfn "$SRC/assemblies/$RUN.ORP.intermediate.fasta" "$DIR/assemblies/$name.ORP.intermediate.fasta"
ln -sfn "$SRC/assemblies/$RUN.ORP.diamond.txt" "$DIR/assemblies/$name.ORP.diamond.txt"
ln -sfn "$BASE/nohit_orfs.tsv" "$DIR/nohit_orfs.tsv"

sub="$SCRATCH/$name"
mkdir -p "$sub"
trap 'rm -rf "$sub"' EXIT
"$BIN/seqtk" sample -s "$SEED" "$read1" "$frac" > "$sub/r1.fq" &
pid1=$!
"$BIN/seqtk" sample -s "$SEED" "$read2" "$frac" > "$sub/r2.fq" &
pid2=$!
wait "$pid1"
wait "$pid2"
n1=$(wc -l < "$sub/r1.fq")
n2=$(wc -l < "$sub/r2.fq")
echo "subsample written: $(du -sh --apparent-size "$sub" | cut -f1), $((n1 / 4)) and $((n2 / 4)) reads -- $(date '+%F %T')"
if [ "$n1" -ne "$n2" ] || [ $((n1 % 4)) -ne 0 ]; then
    echo "r1 has $n1 lines and r2 $n2: a write failed" >&2
    exit 1
fi

# ORP's own salmon quant (oyster.py, Oyster.salmon), on the subsample.
rm -rf "$DIR/quant"
"$BIN/salmon" quant --no-version-check \
    -p "$CPU" -i "$SRC/quants/$RUN.shucked.idx" --seqBias --gcBias --libType A \
    -1 "$sub/r1.fq" -2 "$sub/r2.fq" -o "$DIR/quant"
rm -rf "$sub"

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
