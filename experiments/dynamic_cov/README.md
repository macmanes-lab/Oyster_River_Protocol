# `--min-cov` against `--tpm-filt`

Does a coverage filter (`--min-cov`: NumReads x mean fragment length /
contig length) keep a better set of contigs without a swissprot hit than
the TPM filter (`--tpm-filt 1`)? TPM is a share of a million: across the
119 final assemblies of the 108-sample set, 1 TPM is anywhere from 4 to 183
reads per kb, mostly following depth, so the same cutoff means different
things in different samples.

## Datasets

Five finished 4.1.0-dev29 runs on Premise (`~/<run>`), which kept their four
assemblies and corrected reads as `.gz`, spanning final size and depth:

| run | dev29 final contigs | pooled contigs | corrected R1 |
|---|---|---|---|
| DRR030368_dev29 | 26,782 | 227,879 | 6.2 GB |
| DRR053698_dev29 | 41,079 | 318,577 | 11.4 GB |
| SRR866209_dev29 | 45,836 | 246,715 | 2.6 GB |
| SRR1138704_dev29 | 108,256 | 364,748 | 1.1 GB |
| SRR868662_dev29 | 139,119 | 814,056 | 7.9 GB |

## Running

Premise checkout of this branch at `~/orp_dynamic_cov`, with `software` and
`busco_dbs` symlinked from `~/Oyster_River_Protocol`. One job per dataset:

    for r in DRR030368 DRR053698 SRR866209 SRR1138704 SRR868662; do
        sbatch -J dyncov_$r --export=ALL,SRC=$HOME/${r}_dev29 \
            ~/orp_dynamic_cov/experiments/dynamic_cov/sweep.sbatch
    done

Each job runs `chowder.py` in `~/dyncov/<run>_cov` once per level (`tpm1
none cov0.5 cov1 cov2 cov5`). Only the first does the merge; later levels
re-run from the filter on. Then:

    python ~/orp_dynamic_cov/experiments/dynamic_cov/collect.py ~/dyncov

Every level also writes `filter_candidates.tsv`: how many contigs each of
a range of TPM, coverage and read-count thresholds would keep, from the
same quantification.
