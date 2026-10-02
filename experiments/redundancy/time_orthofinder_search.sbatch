#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J of_search_time
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --exclude=node117,node118
#SBATCH --output=%x_%j.log
#
# Time ORP's orthogroup search as it runs now (diamond blastp over DNA, one
# strand) against OrthoFinder's blastn entry (both strands) on one sample's
# four real assemblies. Same node, same allocation, run one after the other
# as separate job steps, so `sacct -j <job>` gives each one's Elapsed and
# MaxRSS on the same hardware.
#
# usage: sbatch time_orthofinder_search.sbatch RUN [OUTDIR]
#   RUN     run name under $COMPARE/orp_runs (its four .fasta.gz are used)
#
# Inputs match ORP's: contigs over 200 bp (scripts/long.seq.py), and for
# diamond only, N -> X as in write_search_inputs(). blastn gets the Ns.
#
# Threads: 24 cpus / 120G through ORP's own planner gives diamond 10
# searches x 2 threads (-t 10, diamond_orp_2, 12G floor per search). The
# stock blastn entry has no thread flag, so blastn gets -t 16: all 16
# searches at once, one core each. That is the most it can use without
# adding a threaded entry to the shared config.json.
set -euo pipefail

RUN=${1:?usage: sbatch time_orthofinder_search.sbatch RUN [OUTDIR]}
COMPARE=${COMPARE:-/mnt/home/macmaneslab/macmanes/compare}
OUT=${2:-$HOME/redundancy_tests/search_timing/$RUN}
ASM=$COMPARE/orp_runs/$RUN/assemblies

module load anaconda/colsa
rm -rf "$OUT"
mkdir -p "$OUT/diamond/in" "$OUT/blastn/in"

python3 - "$ASM" "$RUN" "$OUT" <<'EOF'
import gzip, sys
asm, run, out = sys.argv[1:]
mask = bytes.maketrans(b"Nn", b"XX")
for name in ("transabyss", "spades75", "spades55", "trinity.Trinity"):
    src = f"{asm}/{run}.{name}.fasta.gz"
    recs, head, seq = 0, None, []
    with gzip.open(src, "rb") as f, \
         open(f"{out}/diamond/in/{run}.{name}.fasta.short.fasta", "wb") as d, \
         open(f"{out}/blastn/in/{run}.{name}.fasta.short.fasta", "wb") as b:
        def flush():
            global recs
            s = b"".join(seq)
            if head is not None and len(s) > 200:
                b.write(head + s + b"\n")
                d.write(head + s.translate(mask) + b"\n")
                recs += 1
        for line in f:
            if line.startswith(b">"):
                flush()
                head, seq = line, []
            else:
                seq.append(line.strip())
        flush()
    print(f"{name}: {recs} contigs over 200 bp")
EOF

echo "node: $(hostname)"
for arm in diamond blastn; do
    if [[ $arm == diamond ]]; then flags=(-t 10 -S diamond_orp_2); else flags=(-t 16 -S blastn); fi
    echo "== $arm: orthofinder -d -I 12 -og -a 1 ${flags[*]}  start $(date '+%F %T')"
    SECONDS=0
    srun --ntasks=1 --cpus-per-task=24 --job-name="$arm" \
        conda run --no-capture-output -n orp_orthofinder \
        orthofinder -d -I 12 -f "$OUT/$arm/in" -og -a 1 "${flags[@]}" \
        > "$OUT/$arm/orthofinder.log" 2>&1
    echo "== $arm: done $(date '+%F %T'), $((SECONDS / 3600))h $((SECONDS % 3600 / 60))m $((SECONDS % 60))s"
done

sacct -j "$SLURM_JOB_ID" --format=JobName%12,Elapsed,MaxRSS,TotalCPU,State
