#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J busco_pool
#SBATCH --cpus-per-task=24
#SBATCH --mem 60G
#SBATCH --exclude=node117,node118
#SBATCH --output=%x_%j.log
#
# BUSCO over the pooled search input (all four assemblies, contigs over
# 200 bp) for one sample. It tells compare_orthogroups.py which contigs
# carry which BUSCO gene, so it can count orthogroups that merge two genes
# (makeorthout would then keep only one of them) and genes split across
# several orthogroups. Contig names do not clash between assemblers
# (TRINITY_*, NODE_*, R*/S*), so the pool needs no renaming.
#
# usage: sbatch busco_pool.sbatch [TIMING_DIR]
set -euo pipefail

TIMING=${1:-$HOME/redundancy_tests/search_timing/SRR1139197}
ORP=${ORP:-$HOME/Oyster_River_Protocol}
out=$TIMING/busco_pool

module load anaconda/colsa
mkdir -p "$out"
cat "$TIMING"/blastn/in/*.short.fasta > "$out/pool.fasta"
cd "$out"
conda run --no-capture-output -n orp_busco busco --offline \
    --lineage eukaryota_odb12.2 --download_path "$ORP/busco_dbs" \
    -i pool.fasta -m transcriptome --cpu "$SLURM_CPUS_PER_TASK" \
    -o run_pool --config "$ORP/software/config.ini" -f
ls run_pool/run_*/full_table.tsv
