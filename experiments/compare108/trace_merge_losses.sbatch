#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J orp_trace
#SBATCH --cpus-per-task=32
#SBATCH --mem 128G
#SBATCH --exclude=node117,node118
#SBATCH --array=1-59%20
#SBATCH --output=/mnt/home/macmaneslab/macmanes/compare/lineage_410/merge_trace/logs/trace_%A_%a.log

# What the ORP merge loses from its four component assemblies, one task per
# line of lineage_410/manifest.tsv (59 runs). trace_merge_losses.py lists every
# BUSCO complete in a component but fragmented or missing in ORP, and every
# Swiss-Prot gene hit in a component but not in ORP, with the evidence for why
# (best DIAMOND hit, ORP contigs kept for that gene, best match in the ORP
# assembly, merge files that mention the contig). DIAMOND is rerun on all five
# assemblies as ORP runs it (blastx -e 1e-8 --top 0.1) and cached in
# merge_trace/<SRR>/, so a resubmitted task resumes.
#
#   mkdir -p /mnt/home/macmaneslab/macmanes/compare/lineage_410/merge_trace/logs
#   sbatch ~/Oyster_River_Protocol/experiments/compare108/trace_merge_losses.sbatch
#   # then, from the login node:
#   python3 ~/Oyster_River_Protocol/experiments/compare108/trace_merge_losses.py --collate

set -euo pipefail
HERE="${HERE:-$HOME/Oyster_River_Protocol/experiments/compare108}"
echo "=== task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} on $(hostname) at $(date) ==="
module load anaconda/colsa
# the orp env has diamond; minimap2 or blastn is used if present, else the
# sequence match column is left empty
conda run --no-capture-output -n orp python3 "$HERE/trace_merge_losses.py" \
    --task "$SLURM_ARRAY_TASK_ID" --cpu "${SLURM_CPUS_PER_TASK:-32}"
echo "=== done at $(date) ==="
