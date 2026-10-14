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

## Results (2026-10-14)

`results/sweep_table.tsv` is collect.py's table, `results/removed_by_rule.tsv`
preview.py's, and `results/depth.tsv` each dataset's filter_candidates
header. The first sweep (jobs 1322434-8) scored every level after tpm1
against tpm1's snap index (fixed in `clear_pytransrate_outdir`); every
level but tpm1 was rerun (1322478-82), and these are the reruns.

BUSCO and unique genes are identical at every level: the filter only ever
removes contigs without a swissprot hit. pyTransRate score, with no-hit
contigs removed and those among them with an ORF > 149 codons:

| dataset | 1 TPM = | tpm1 | cov2 | cov5 |
|---|---|---|---|---|
| SRR866209 | 1.8x | 0.516 (1,624 / 165) | 0.517 (3,763 / 140) | 0.524 (15,213 / 278) |
| SRR1138704 | 2.6x | 0.362 (6,361 / 237) | 0.366 (6,775 / 110) | 0.369 (33,766 / 435) |
| DRR053698 | 6.1x | 0.450 (34,399 / 2,737) | 0.472 (6,844 / 651) | 0.460 (26,769 / 1,276) |
| SRR868662 | 6.9x | 0.502 (40,739 / 3,425) | 0.475 (9,256 / 637) | 0.506 (32,129 / 1,273) |
| DRR030368 | 7.3x | 0.384 (27,439 / 6,175) | 0.375 (4,271 / 1,116) | 0.383 (16,502 / 2,392) |

- 1 TPM is 2-7x coverage depending on depth, so tpm1 removes 4-6x as many
  no-hit contigs as cov2 on deep libraries, including 2-5x as many with a
  long ORF.
- On the deep libraries every coverage level maps more fragments than tpm1
  (+0.1 to +0.45%), and keeps pyTransRate's optimal score as high or
  higher (tpm1 lowers it, e.g. DRR053698 0.490 vs 0.510).
- cov5 matches or beats tpm1's score everywhere but DRR030368 (-0.001) and
  removes about half the long-ORF contigs tpm1 does there, but on the
  shallow libraries it costs mapping (SRR1138704 -0.5%) and removes the
  most long-ORF contigs of any rule.
