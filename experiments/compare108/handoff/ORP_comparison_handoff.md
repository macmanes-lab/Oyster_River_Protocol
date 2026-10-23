# ORP vs original assemblies and single assemblers: hand-off for the manuscript

Status as of 2026-10-23: all 59 samples are complete. Each sample has an ORP assembly, the
original assembly and the four component assemblies inside the ORP run, all scored the same
way.

## Files in this folder

| File | What it is |
|---|---|
| `orp_vs_original.csv` | One row per sample (59). The ORP assembly and the original side by side; columns end in `_orp` or `_original`. |
| `orp_vs_components.csv` | Five rows per sample (295): the ORP assembly and its four component assemblies (Trinity, rnaSPAdes auto-k, rnaSPAdes high-k, Trans-ABySS). |
| `comparison_figure.png` / `.pdf` | Six scatter panels, original (x) vs ORP (y): BUSCO complete, duplicated and missing; pyTransRate score and optimal score; good mappings. |
| `orp_v5_comparison_tables.xlsx` | One sheet per group with both assemblies' metrics, plus a Summary sheet with paired tests. |
| `summary_stats.csv` | The Summary sheet as csv. |

The figure and workbook still say "ORP version 5" / "ORP v5". That label is out of date: the
release is numbered **4.0.0** (see Methods). Rename it before the figure goes into the paper.

## The samples

There are 59 paired-end Illumina RNA-seq libraries from the SRA. Each one has an original
assembly made by someone else:

- **44 TSA assemblies:** the assembly the original authors deposited in NCBI's Transcriptome
  Shotgun Assembly archive.
- **15 ORP 2018 assemblies:** the orthomerged assemblies from the 2018 Oyster River Protocol
  paper (Zenodo doi:10.5281/zenodo.1320141). These are all 15 samples from that paper, so the
  new ORP can be compared with the old one on the same reads.

| Group | Samples | BUSCO dataset |
|---|---:|---|
| Land plants | 14 | embryophyta_odb12.2 |
| Insects | 10 | insecta_odb12 |
| Echinoderms | 10 | metazoa_odb12.2 |
| Non-mammal vertebrates | 9 | vertebrata_odb12 |
| Mammals | 8 | mammalia_odb12 |
| Other | 8 | one dataset per sample (below) |

"Other" is made up of:

| Sample | Organism | BUSCO dataset |
|---|---|---|
| tsa_GARO | *Echinogammarus veneris* (amphipod) | crustacea_odb12.2 |
| tsa_GAZR | *Stegodyphus mimosarum* (spider) | arachnida_odb12.2 |
| tsa_GAGE | *Apiospora malaysiana* (fungus; see Known issues) | agaricales_odb12.2 |
| orp_ERR1674585 | *Schistosoma mansoni* (fluke) | lophotrochozoa_odb12.2 |
| orp_DRR030368 | *Echinococcus multilocularis* (tapeworm) | lophotrochozoa_odb12.2 |
| orp_SRR2016923 | *Phyllodoce medipapillata* (polychaete) | lophotrochozoa_odb12.2 |
| orp_ERR1016675 | *Heterorhabditis indica* (nematode) | nematoda_odb12.2 |
| orp_ERR058009 | *Entamoeba histolytica* (amoeba) | amoebozoa_odb12.2 |

Each sample is scored against the most specific BUSCO dataset that covers its whole group. BUSCO
has no echinoderm, flatworm or annelid dataset, so echinoderms use metazoa, and the two
flatworms and the polychaete use lophotrochozoa. The lophotrochozoa dataset includes
*S. mansoni* and *E. multilocularis*, and the amoebozoa dataset includes *E. histolytica*.
The nematoda dataset has no *Heterorhabditis*, but it does include its closest relatives, the
strongylids (*Haemonchus*, *Ancylostoma*, *Necator*).

The datasets come from two OrthoDB releases. `mammalia_odb12`, `vertebrata_odb12` and
`insecta_odb12` are OrthoDB 12.1; all the others are 12.2.

## Methods

### ORP runs

Every sample was assembled from its raw reads with the same command:

```
oyster.py --min-cov auto --lineage <dataset> --mem 190 --cpu 32 \
    --read1 <R1> --read2 <R2> --max-parallel 2 --transabyss-mpi on \
    --dir <run dir> --runout <run>
```

What the settings mean:

- **Defaults otherwise:** no read normalization, and no TPM filter.
- **`--min-cov auto`:** contigs are filtered by a coverage threshold that scales with library
  depth, max(1, 2 × (read pairs / 60M)^0.85). The read count comes from salmon on the trimmed,
  error-corrected reads.
- **Resources:** 32 CPUs and 190 GB per sample.
- **Run time:** median 7.1 h per sample (range 2.1–33.5 h).

Software versions:

| Tool | Version |
|---|---|
| Trimmomatic | 0.41 |
| Rcorrector | 1.0.7 |
| Trinity | 2.15.2 (ORP's bundled build) |
| rnaSPAdes | 4.3.0 |
| Trans-ABySS | 2.0.1 (under MPI) |
| DIAMOND | 2.2.5 |
| CD-HIT | 4.8.1 |
| BUSCO | 6.1.0 |
| pyTransRate | 2.2.2 |

**ORP version.** 50 samples ran on the 4.1.0 release candidate (commit e88d037). The 9 ORP 2018
samples added later ran on the tagged release 4.0.0 (b9459aa). These are the same pipeline:
between the two commits, the release was renumbered from 4.1.0 to 4.0.0, and the only code
changes were resume handling and a startup banner. Nothing that affects a fresh run's output
changed. The `orp_version` column records which one each sample used.

### Scoring

The ORP assembly, the original and the four components are all scored the same way.

- **BUSCO completeness:** BUSCO 6.1.0 in transcriptome mode, offline, with ORP's `config.ini`,
  against the dataset for that sample's group.
- **pyTransRate:** pyTransRate 2.2.2 with the same command ORP runs on its own assembly
  (`pytransrate -a … --left … --right … -t … --max-memory …`). It uses the trimmed,
  error-corrected reads from that sample's ORP run, so all six assemblies of a sample are scored
  against the same reads.
  - Before 2026-10-23, the original assemblies had been scored on raw reads with an older build,
    and the 47-pair draft figure used ORP 3's Ruby TransRate. Their scores were redone to match
    ORP's. Most changed by less than 0.02.
- **Unique genes:** the number of distinct SwissProt genes hit by `diamond blastx -e 1e-8
  --top 0.1`, counted the way ORP counts them. For each component, the count is the one ORP
  itself reports for that assembler.
- **Proper pairs:** the percentage of reads mapping as proper pairs in ORP's strand check.
  This is reported for the ORP assembly and the components only.

All paired tests are two-sided Wilcoxon signed-rank tests on the 59 pairs.

## Results

### ORP vs the original assembly (59 pairs)

| Metric | Median: original | Median: ORP | Median difference | ORP better | ORP worse | Wilcoxon p |
|---|---:|---:|---:|---:|---:|---:|
| BUSCO complete (%) | 61.3 | 80.5 | +13.0 | **59** | 0 | 2.4e-11 |
| BUSCO single-copy (%) | 39.8 | 67.4 | +26.4 | 57 | 2 | 3.8e-11 |
| BUSCO duplicated (%) | 15.7 | 5.5 | −9.7 | 44 (lower) | 15 | 1.4e-5 |
| BUSCO fragmented (%) | 13.3 | 7.2 | −5.6 | 53 (lower) | 6 | 4.2e-10 |
| BUSCO missing (%) | 23.0 | 10.8 | −6.8 | 58 (lower) | 0 | 3.5e-11 |
| pyTransRate score | 0.250 | 0.438 | +0.153 | 56 | 3 | 7.7e-11 |
| pyTransRate optimal score | 0.409 | 0.495 | +0.080 | 47 | 12 | 4.3e-8 |
| Good mappings (fraction) | 0.772 | 0.888 | +0.089 | 53 | 6 | 4.8e-10 |
| Unique SwissProt genes | 10,217 | 12,243 | +1,501 | 57 | 2 | 2.3e-10 |
| Contigs | 76,838 | 105,273 | +15,189 | — | — | 0.069 |
| N50 (bp) | 1,142 | 1,240 | +192 | 40 | 19 | 0.0071 |

Notes on reading the table:

- The tie for BUSCO missing (58 + 0 of 59) is one sample with the same value on both sides.
- The figure and workbook round good mappings to 0.1%, which turns a few small differences into
  ties: they show 51 better, where the unrounded csv gives 53.

**By group.** Mean BUSCO complete and duplicated, and median pyTransRate score:

| Group | n | Complete: ORP | Complete: original | Duplicated: ORP | Duplicated: original | pyTransRate: ORP | pyTransRate: original |
|---|---:|---:|---:|---:|---:|---:|---:|
| Land plants | 14 | 86.7 | 63.5 | 19.3 | 18.5 | 0.452 | 0.290 |
| Echinoderms | 10 | 78.5 | 66.5 | 4.4 | 22.3 | 0.389 | 0.308 |
| Insects | 10 | 75.0 | 54.7 | 9.0 | 16.9 | 0.462 | 0.276 |
| Non-mammal vertebrates | 9 | 71.5 | 51.8 | 4.6 | 15.4 | 0.413 | 0.185 |
| Other | 8 | 69.5 | 54.8 | 11.5 | 19.8 | 0.389 | 0.214 |
| Mammals | 8 | 49.0 | 41.5 | 1.2 | 14.3 | 0.486 | 0.224 |

Land plants are the one group where ORP does not reduce duplication; the mean stays at about
19%. The mammal scores are low on both sides because the mammal dataset is large (12,277 genes)
and several of the libraries are small: the two camels score 13% and 29%.

**New ORP vs ORP 2018 (the 15 samples from the 2018 paper):**

| Metric | Result |
|---|---|
| BUSCO complete | Higher in 15 of 15; mean 80.2% vs 67.0% |
| BUSCO duplicated | Mean halved, 13.9% vs 27.8% |
| pyTransRate score | Higher in 14 of 15 |
| Unique genes | Higher in 13 of 15 |
| pyTransRate optimal score | Higher in only 5 of 15 |
| Good mappings | Higher in only 9 of 15 |

The 2018 assemblies were also built from error-corrected reads, so they are the closest
competitors in the set. All 6 samples where the original has more good mappings are ORP 2018
assemblies, and so are 10 of the 12 where its optimal score is higher. The other two are
tsa_GAKE and tsa_GAVJ.

### ORP vs its own component assemblies (59 samples × 4 assemblers)

These compare the merged ORP assembly with the single-assembler builds it was made from, on the
same reads.

| Metric | ORP beats the best single assembler | Which single assembler is usually best |
|---|---:|---|
| BUSCO complete | 53 of 59 (2 ties) | rnaSPAdes auto-k 34, Trinity 24, rnaSPAdes high-k 1 |
| BUSCO duplicated (lower) | 38 of 59 (2 ties) | Trans-ABySS 32, rnaSPAdes high-k 26, rnaSPAdes auto-k 1 |
| Unique genes | 57 of 59 | rnaSPAdes auto-k 31, Trans-ABySS 20, rnaSPAdes high-k 5, Trinity 3 |
| pyTransRate score | **20 of 59** | rnaSPAdes high-k 43 |
| pyTransRate optimal score | **7 of 59** | rnaSPAdes auto-k 37, rnaSPAdes high-k 19 |
| Good mappings | **11 of 59** (1 tie) | rnaSPAdes auto-k 34, rnaSPAdes high-k 24 |

ORP against each assembler, as the number of samples where ORP is better, with the median of
ORP minus that assembler in brackets:

| Metric | vs Trinity | vs rnaSPAdes auto-k | vs rnaSPAdes high-k | vs Trans-ABySS |
|---|---:|---:|---:|---:|
| BUSCO complete | 54 (+2.9) | 57 (+1.4) | 59 (+7.2) | 58 (+5.9) |
| BUSCO duplicated (lower) | 59 (−28.3) | 58 (−13.7) | 54 (−7.1) | 42 (−5.0) |
| Unique genes | 59 (+911) | 59 (+639) | 59 (+1,147) | 57 (+998) |
| pyTransRate score | 54 (+0.056) | 46 (+0.079) | **24 (−0.011)** | 55 (+0.087) |
| pyTransRate optimal | 47 (+0.030) | **11 (−0.014)** | **23 (−0.012)** | 55 (+0.067) |
| Good mappings | 45 (+0.021) | **12 (−0.013)** | **18 (−0.006)** | 59 (+0.036) |

**What this means for the paper.** The merged ORP assembly is the best choice on gene content:
it has more complete BUSCOs and more unique genes than any single assembler, and less
duplication than Trinity or rnaSPAdes auto-k. On pyTransRate's read-mapping measures, though,
the rnaSPAdes assemblies, especially high-k, usually match or slightly beat ORP. The differences
are small (medians of about 0.01) but consistent, and the paper needs to address them rather
than claim ORP wins on every metric.

## Known issues and caveats

- **tsa_GAGE is the wrong organism for its dataset.** Its `pair.json` lists *Termitomyces
  clypeatus* (a mushroom), but ENA lists the reads (SRR746797) as *Apiospora malaysiana*, a sac
  fungus. Both assemblies were scored against `agaricales_odb12.2`, the mushroom dataset, so both
  are low (ORP 32%, original 30%). It should be rescored with `sordariomycetes_odb12.2`, or
  dropped.
- **tsa_GAGA reads had placeholder quality scores.** ENA's copy had SRA-lite qualities (every
  quality character was `?`). The reads were rebuilt from NCBI's original-quality data with
  `fasterq-dump`; the sequences are identical.
- **Duplicate contig names.** The two *Anolis* TSA assemblies (GAGA, GAGB) had duplicate names,
  and were de-duplicated before scoring. GAGA went from 1,135,207 to 735,899 sequences, and GAGB
  from 2,056,067 to 1,195,995.
- **Repeated species.** *Mus musculus* appears twice (orp_DRR036858, orp_SRR2086412) and
  *Anolis carolinensis* twice (tsa_GAGA, tsa_GAGB): different libraries of the same species.
- **Parasites with reduced gene sets** will score lower whatever the assembly: *Aeginetia
  indica* (a holoparasitic plant), *Entamoeba* and the flatworms. *Entamoeba*'s low unique-gene
  count (about 2,400) mostly reflects how few of its genes have close SwissProt matches.
- **Libraries that limit the score:** *Stegodyphus* (about 11,000 genes expressed),
  *Eucidaris* (an ovary-only GAIIx library) and the two camels (2–4 GB of reads).
- **tsa_GAKE (*Thlaspi*)** is the main exception on duplication: 24.4% for ORP against 2.8% for
  the original. It is also one of only two samples (with tsa_GAZR) where ORP has fewer
  single-copy BUSCOs than the original.
- **Excluded:** tsa_GAQK (*Hynobius*, R2 reads of only 62 bp).

## Column dictionary

`orp_vs_original.csv`:

| Column | Meaning |
|---|---|
| `pair_id` | Sample id: `tsa_XXXX` (original is the TSA assembly) or `orp_<run>` (original is the ORP 2018 assembly) |
| `run` | SRA/ENA/DDBJ run accession |
| `organism`, `group`, `busco_dataset` | as above |
| `original_source` | `TSA` or `ORP 2018 (orthomerged)` |
| `orp_version` | `4.1.0` (release candidate) or `4.0.0` (tagged); same code |
| `reads_gb`, `read_pairs` | raw read volume (gzipped GB) and read pairs before trimming |
| `min_cov_applied` | the coverage cutoff `--min-cov auto` chose |
| `busco_{complete,single,duplicated,fragmented,missing}_{orp,original}` | BUSCO percentages |
| `busco_n` | number of BUSCOs in the dataset |
| `pytransrate_score_*`, `pytransrate_optimal_*` | pyTransRate assembly score and optimal score |
| `p_good_mapping_*`, `p_fragments_mapped_*` | fraction of read pairs mapping well / at all (0–1) |
| `contigs_*`, `n50_*`, `mean_len_*`, `total_bases_*` | assembly statistics from pyTransRate |
| `unique_genes_orp`, `unique_genes_original` | distinct SwissProt genes hit |
| `proper_pairs_pct_orp` | % proper pairs mapping to the ORP assembly |
| `orp_hours`, `orp_peak_mem_gb` | ORP wall-clock hours and peak memory |

`orp_vs_components.csv` uses the same metric names without the side suffix, with `assembly` set
to one of `ORP`, `trinity`, `spadesauto`, `spadeshigh` or `transabyss`.

## Where everything lives

- **Raw outputs on Premise:**
  - ORP runs: `/mnt/home/macmaneslab/macmanes/compare/lineage_410/<run>/`
  - Component evaluations: `lineage_410/components/<run>/<assembler>/`
  - Original-assembly BUSCO: `compare/busco_original/<pair>/`
  - Original-assembly pyTransRate: `compare/pytransrate_orig_orp/<pair>/`
- **Collected tables:** `lineage_410/lineage410_metrics.csv` and `components_metrics.csv`.
  These csvs are curated subsets of those two.
- **Scripts:** `experiments/compare108/` in the Oyster_River_Protocol repository (commits
  88f245a and 48f0757): `orp_lineage410.sbatch`, `busco_original*.sbatch`,
  `orig_unique_genes.sbatch`, `components_score.sbatch`, `orig_pytransrate.sbatch`,
  `collect_lineage410.py`, `collect_components.py`, `submit_orp_originals.sh`.
- **Slurm jobs:**
  - ORP runs: 1323838 (50 samples) and 1324733 (9)
  - Components: 1324263 and 1324923
  - BUSCO on originals: 1322827, 1323174, 1323375, 1323209 and 1324734
  - Unique genes in originals: 1324735, plus its earlier 50-sample run
  - Original pyTransRate: 1325177
- **Figure and workbook:** built by `~/orp/compare_analysis/make_figure.py` from
  `lineage410_metrics.csv`.
