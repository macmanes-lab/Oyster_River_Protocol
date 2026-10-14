#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J dyncov
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 120G
#SBATCH --output=dyncov_%x_%j.log
#SBATCH --exclude=node117,node118
# One dataset's filter sweep: chowder.py over a finished run's four
# assemblies and corrected reads, once per filter level, in one run
# directory. The first level does the whole merge; each later one changes
# only the threshold, which (<run>.filter.settings) re-runs filter,
# secondfilter, BUSCO, pyTransRate, strandeval and reportgen and nothing
# before them. After each level its assembly and reports are copied to
# $DIR/levels/<level>/.
#
#     sbatch -J SRR1138704 --export=ALL,SRC=$HOME/SRR1138704_dev29 \
#         experiments/dynamic_cov/sweep.sbatch
#
# A level is none, tpm<X> (--tpm-filt X) or cov<X> (--min-cov X).
# Then: python experiments/dynamic_cov/collect.py $HOME/dyncov

set -euo pipefail
: "${SRC:?set SRC to a finished oyster.py run directory}"
LEVELS="${LEVELS:-tpm1 none cov0.5 cov1 cov2 cov5}"
ORP="${ORP:-$HOME/orp_dynamic_cov}"
SRCRUN=$(basename "$SRC")
RUN="${RUN:-${SRCRUN%_dev*}_cov}"
DIR="${DIR:-$HOME/dyncov/$RUN}"

module load anaconda/colsa

CPU="${SLURM_CPUS_PER_TASK:-24}"
MEM=$(( ${SLURM_MEM_PER_NODE:-122880} / 1024 - 10 ))
A="$SRC/assemblies/$SRCRUN"
mkdir -p "$DIR/levels"
cd "$DIR"
echo "chowder sweep: $SRCRUN -> $DIR ($RUN), levels: $LEVELS, ORP $(cat "$ORP/version.txt")"

for level in $LEVELS; do
    case "$level" in
        none)  filt=(--tpm-filt 0 --min-cov 0) ;;
        tpm*)  filt=(--tpm-filt "${level#tpm}" --min-cov 0) ;;
        cov*)  filt=(--tpm-filt 0 --min-cov "${level#cov}") ;;
        *)     echo "unknown level $level" >&2; exit 1 ;;
    esac
    echo "=== level $level: ${filt[*]} -- $(date '+%F %T') ==="
    "$ORP/chowder.py" \
        --assemblies "$A.trinity.Trinity.fasta.gz" "$A.transabyss.fasta.gz" \
                     "$A.spadesauto.fasta.gz" "$A.spadeshigh.fasta.gz" \
        --labels trinity transabyss spadesauto spadeshigh \
        --read1 "$SRC/rcorr/$SRCRUN.TRIM_1P.cor.fq.gz" \
        --read2 "$SRC/rcorr/$SRCRUN.TRIM_2P.cor.fq.gz" \
        --corrected-reads \
        "${filt[@]}" \
        --no-cleanup \
        --cpu "$CPU" --mem "$MEM" \
        --dir "$DIR" --runout "$RUN"
    out="$DIR/levels/$level"
    rm -rf "$out"
    mkdir -p "$out"
    cp "$DIR/assemblies/$RUN.ORP.fasta" "$DIR/reports/qualreport.$RUN" \
       "$DIR/reports/$RUN.filter_candidates.tsv" "$out/"
    cp "$DIR/reports/pytransrate_$RUN/assemblies.csv" "$out/pytransrate.assemblies.csv"
    cp "$DIR"/reports/run_"$RUN".ORP/short_summary*.txt "$out/" 2>/dev/null || true
    echo "=== level $level done -- $(date '+%F %T') ==="
done
