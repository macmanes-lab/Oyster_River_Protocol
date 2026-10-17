# ORP 4.1.0-dev30 on 50 samples, scored with lineage-specific BUSCO

Status as of 2026-10-17: 48 of 50 runs finished. Still running are *Brassica juncea*
(tsa_GBEP, plant task 10) and the *Anolis* GAGA rerun (job 1323199).

The ORP runs were assembled on Premise from pairs in `compare/pairs`. Each one is
BUSCO-scored against the most specific odb12 dataset for its group. Each pair's
original assembly is scored the same way. That's the published TSA assembly for
`tsa_*` pairs, and the 2018 ORP orthomerged assembly (Zenodo) for `orp_*` pairs.

## Summary

Every finished ORP assembly scores higher on BUSCO completeness than its original.

| Group | Samples | Finished | ORP higher | Mean complete: ORP | Mean complete: original |
|---|---:|---:|---:|---:|---:|
| Mammals | 7 | 7 | 7 of 7 | 48.8% | 41.2% |
| Non-mammal vertebrates | 9 | 8 | 8 of 8 | 73.8% | 54.8% |
| Insects | 10 | 10 | 10 of 10 | 74.9% | 54.7% |
| Land plants | 10 | 9 | 9 of 9 | 85.7% | 63.4% |
| Echinoderms | 10 | 10 | 10 of 10 | 78.5% | 66.5% |
| Other groups | 4 | 4 | 4 of 4 | 50.6% | 41.6% |
| **All** | **50** | **48** | **48 of 48** | | |

The gains usually come with much less duplication. Originals often carry 15–50%
duplicated BUSCOs, while ORP assemblies mostly have 1–15%. The one exception is
*Thlaspi* (see below).

## How the runs were set up

All runs use ORP 4.1.0-dev30 from `~/Oyster_River_Protocol` on Premise (master)
with the same settings:

- `--normalize-reads`, `--min-cov 2` and no `--tpm-filt`
- `--transabyss-mpi on`, `--max-parallel 2`
- 32 CPUs and 190G per sample
- `--runout <SRR>_v4_1_0-dev30`, so every output file carries the ORP version

Each group writes to its own `compare/orp_<group>/<SRR>/` tree, so the finished 4.0
runs in `compare/orp_runs/` are untouched.

BUSCO on the originals uses ORP's own BUSCO call: BUSCO 6.1.0 from the `orp_busco`
env, offline, transcriptome mode, ORP's `config.ini`, and the same odb12 dataset as
that pair's ORP run.

The BUSCO datasets come from two OrthoDB releases:

- **OrthoDB 12.1:** `mammalia_odb12`, `vertebrata_odb12` and `insecta_odb12`
- **OrthoDB 12.2:** `embryophyta_odb12.2`, `metazoa_odb12.2`, `crustacea_odb12.2`,
  `agaricales_odb12.2`, `arachnida_odb12.2` and `lophotrochozoa_odb12.2`

BUSCO has no echinoderm or flatworm dataset. Echinoderms use metazoa.
*Schistosoma* uses lophotrochozoa, whose species list includes *S. mansoni*.

### Scripts and jobs

| Script (`experiments/compare108/`) | Job | Samples | Output |
|---|---|---|---|
| `orp_mammalia.sbatch` | 1322800 | 7 mammals | `compare/orp_mammalia/` |
| `orp_vertebrata.sbatch` | 1322803, 1323199 | 9 non-mammal vertebrates | `compare/orp_vertebrata/` |
| `orp_insecta.sbatch` | 1322859 | 10 insects | `compare/orp_insecta/` |
| `orp_embryophyta.sbatch` | 1323168 | 10 land plants | `compare/orp_embryophyta/` |
| `orp_echinodermata.sbatch` | 1323171 | 10 echinoderms | `compare/orp_echinodermata/` |
| `orp_other.sbatch` | 1323371 | crustacean, fungus, arachnid, flatworm | `compare/orp_other/` |
| `busco_original.sbatch` | 1322827, 1323209 | originals of mammals and vertebrates | `compare/busco_original/` |
| `busco_original2.sbatch` | 1323174 | originals of insects, plants and echinoderms | `compare/busco_original/` |
| `busco_original3.sbatch` | 1323375 | originals of the four other groups | `compare/busco_original/` |

## Results by group

"Reads mapped" is ORP's percentage of reads mapping back to the ORP assembly as
proper pairs.

### Mammals

Job 1322800, `orp_mammalia/`. BUSCO dataset: `mammalia_odb12`.

| Pair | Organism | Reads (GB) | ORP contigs | Complete: ORP | Complete: original | Change | Duplicated: ORP / original | Reads mapped |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| tsa_GAMP | *Callithrix jacchus* | 4.7 | 47,203 | **70.4%** | 65.6% | +4.8 | 0.7% / 18.2% | 97.18% |
| tsa_GABC | *Pan troglodytes* | 9.0 | 24,003 | **65.6%** | 61.3% | +4.3 | 0.8% / 13.6% | 98.81% |
| tsa_GAOE | *Capra hircus* | 6.3 | 109,536 | **54.9%** | 44.5% | +10.4 | 1.2% / 33.1% | 95.50% |
| orp_DRR036858 | *Mus musculus* | 20.7 | 292,257 | **54.6%** | 41.2% | +13.4 | 2.5% / 19.4% | 86.16% |
| orp_SRR1789336 | *Oryctolagus cuniculus* | 5.8 | 72,083 | **54.1%** | 43.7% | +10.4 | 1.3% / 9.0% | 96.55% |
| tsa_GADW | *Camelus dromedarius* | 3.8 | 75,181 | **29.1%** | 24.1% | +5.0 | 0.3% / 3.4% | 88.78% |
| tsa_GAFA | *Camelus bactrianus* | 2.2 | 63,156 | **13.2%** | 7.8% | +5.4 | 0.1% / 0.6% | 81.95% |

### Non-mammal vertebrates

Job 1322803 (+1323199 for GAGA), `orp_vertebrata/`. BUSCO dataset: `vertebrata_odb12`.

| Pair | Organism | Reads (GB) | ORP contigs | Complete: ORP | Complete: original | Change | Duplicated: ORP / original | Reads mapped |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| tsa_GAZY | *Morone chrysops* | 10.9 | 168,425 | **97.8%** | 91.8% | +6.0 | 8.8% / 34.8% | 94.19% |
| tsa_GBAA | *Morone saxatilis* | 10.4 | 166,489 | **97.6%** | 91.3% | +6.3 | 8.8% / 35.5% | 93.93% |
| tsa_GAJJ | *Anoplopoma fimbria* | 7.9 | 99,690 | **85.9%** | 68.8% | +17.1 | 5.1% / 19.1% | 94.29% |
| tsa_GAID | *Oreochromis niloticus* | 3.5 | 117,764 | **75.0%** | 40.9% | +34.1 | 3.7% / 0.3% | 91.50% |
| tsa_GAAA | *Latimeria chalumnae* | 8.6 | 350,540 | **64.4%** | 40.2% | +24.2 | 2.7% / 6.3% | 81.06% |
| tsa_GARA | *Miichthys miiuy* | 1.9 | 65,928 | **59.9%** | 37.8% | +22.1 | 2.1% / 0.6% | 92.60% |
| orp_DRR046632 | *Oncorhynchus mykiss* | 11.3 | 120,929 | **59.7%** | 38.1% | +21.6 | 6.7% / 10.7% | 93.93% |
| tsa_GAGB | *Anolis carolinensis* | 9.5 | 1,100,803 | **49.8%** | 29.4% | +20.4 | 2.3% / 15.8% | 72.43% |
| tsa_GAGA | *Anolis carolinensis* | 8.9 | | still running | 27.9% | | | |

### Insects

Job 1322859, `orp_insecta/`. BUSCO dataset: `insecta_odb12`.

| Pair | Organism | Reads (GB) | ORP contigs | Complete: ORP | Complete: original | Change | Duplicated: ORP / original | Reads mapped |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| tsa_GAJC | *Leptopilina heterotoma* | 7.1 | 252,045 | **93.2%** | 88.2% | +5.0 | 15.9% / 21.0% | 96.12% |
| orp_ERR489297 | *Anopheles gambiae* | 35.7 | 76,075 | **92.2%** | 86.6% | +5.6 | 15.0% / 40.0% | 93.52% |
| tsa_GBDO | *Kerria lacca* | 4.2 | 99,579 | **88.2%** | 84.2% | +4.0 | 14.7% / 52.0% | 94.49% |
| tsa_GAKG | *Cotesia vestalis* | 5.3 | 42,266 | **84.9%** | 63.6% | +21.3 | 12.3% / 14.9% | 94.15% |
| tsa_GAGF | *Chrysopa pallens* | 3.9 | 86,592 | **82.6%** | 68.4% | +14.2 | 3.3% / 9.4% | 94.26% |
| tsa_GAWF | *Sipyloidea sipylus* | 8.5 | 126,156 | **73.0%** | 64.2% | +8.8 | 8.7% / 27.2% | 87.74% |
| tsa_GAJS | *Chilo suppressalis* | 2.1 | 74,570 | **71.7%** | 30.3% | +41.4 | 3.5% / 0.3% | 92.31% |
| tsa_GADR | *Drosophila pseudoananassae nigrens* | 3.0 | 45,995 | **62.4%** | 4.9% | +57.5 | 5.7% / 0.0% | 94.69% |
| tsa_GABB | *Archaeopsylla erinacei* | 9.7 | 48,926 | **53.1%** | 24.4% | +28.7 | 9.9% / 0.3% | 81.92% |
| tsa_GAPQ | *Bemisia tabaci* | 1.8 | 42,686 | **47.6%** | 32.1% | +15.5 | 2.1% / 3.6% | 84.17% |

### Land plants

Job 1323168, `orp_embryophyta/`. BUSCO dataset: `embryophyta_odb12.2`.

| Pair | Organism | Reads (GB) | ORP contigs | Complete: ORP | Complete: original | Change | Duplicated: ORP / original | Reads mapped |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| tsa_GAKE | *Thlaspi arvense* | 33.1 | 82,203 | **98.7%** | 92.1% | +6.6 | 26.5% / 2.8% | 96.59% |
| orp_SRR3499127 | *Nicotiana tabacum* | 5.6 | 149,904 | **97.3%** | 94.3% | +3.0 | 29.9% / 54.5% | 91.05% |
| tsa_GAXU | *Capsicum annuum* | 8.8 | 99,567 | **95.1%** | 87.9% | +7.2 | 12.7% / 9.0% | 95.79% |
| tsa_GABH | *Musa ABB Group* | 6.4 | 146,761 | **94.0%** | 73.2% | +20.8 | 21.4% / 22.2% | 91.45% |
| tsa_GAJE | *Avena sativa* | 10.7 | 160,500 | **91.6%** | 30.2% | +61.4 | 16.2% / 23.9% | 84.83% |
| tsa_GADC | *Plukenetia volubilis* | 3.5 | 65,576 | **89.1%** | 61.6% | +27.5 | 12.8% / 4.4% | 93.47% |
| tsa_GABT | *Cocos nucifera* | 4.4 | 81,943 | **83.8%** | 64.1% | +19.7 | 12.2% / 11.0% | 93.15% |
| tsa_GACB | *Lotus corniculatus* | 1.9 | 53,444 | **68.3%** | 42.8% | +25.5 | 6.0% / 2.8% | 91.90% |
| tsa_GAIG | *Arachis hypogaea* | 0.9 | 59,199 | **53.6%** | 24.6% | +29.0 | 3.6% / 0.2% | 87.00% |
| tsa_GBEP | *Brassica juncea* | 21.4 | | still running | 27.3% | | | |

### Echinoderms

Job 1323171, `orp_echinodermata/`. BUSCO dataset: `metazoa_odb12.2`.

| Pair | Organism | Reads (GB) | ORP contigs | Complete: ORP | Complete: original | Change | Duplicated: ORP / original | Reads mapped |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| tsa_GAUR | *Lytechinus variegatus* | 4.0 | 114,652 | **96.1%** | 83.2% | +12.9 | 7.5% / 35.1% | 90.33% |
| tsa_GAUT | *Sclerodactyla briareus* | 4.2 | 71,714 | **95.1%** | 87.7% | +7.4 | 5.3% / 28.0% | 77.09% |
| tsa_GAWB | *Patiria miniata* | 4.4 | 119,710 | **94.0%** | 85.9% | +8.1 | 8.7% / 31.3% | 89.44% |
| tsa_GAZO | *Oxycomanthus japonicus* | 2.8 | 40,959 | **91.7%** | 84.9% | +6.8 | 2.6% / 23.0% | 85.99% |
| tsa_GAVR | *Sphaerechinus granularis* | 7.5 | 163,021 | **88.0%** | 70.7% | +17.3 | 5.7% / 25.5% | 95.58% |
| tsa_GAUQ | *Ophiocoma echinata* | 3.7 | 189,492 | **80.9%** | 60.8% | +20.1 | 3.4% / 14.1% | 88.65% |
| tsa_GAVS | *Apostichopus japonicus* | 6.4 | 149,713 | **80.5%** | 62.6% | +17.9 | 3.5% / 27.3% | 95.38% |
| tsa_GAVJ | *Asterias amurensis* | 7.1 | 122,235 | **67.9%** | 54.6% | +13.3 | 1.6% / 14.6% | 97.20% |
| tsa_GAVQ | *Luidia clathrata* | 5.8 | 154,305 | **67.3%** | 56.1% | +11.2 | 4.6% / 19.0% | 95.19% |
| tsa_GAZP | *Eucidaris tribuloides* | 2.0 | 108,103 | **23.6%** | 18.3% | +5.3 | 1.3% / 4.9% | 86.56% |

### Other groups

Job 1323371, `orp_other/`. BUSCO dataset: `agaricales_odb12.2`, `arachnida_odb12.2`, `crustacea_odb12.2`, `lophotrochozoa_odb12.2`.

| Pair | Organism | Reads (GB) | ORP contigs | Complete: ORP | Complete: original | Change | Duplicated: ORP / original | Reads mapped |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| orp_ERR1674585 | *Schistosoma mansoni* | 2.4 | 134,978 | **97.1%** | 95.0% | +2.1 | 12.0% / 44.5% | 87.81% |
| tsa_GARO | *Echinogammarus veneris* | 6.1 | 130,204 | **53.7%** | 22.1% | +31.6 | 5.6% / 1.4% | 93.37% |
| tsa_GAGE | *Termitomyces clypeatus* in `pair.json`; really *Apiospora malaysiana* (see below) | 4.5 | 24,579 | **32.0%** | 30.3% | +1.7 | 3.0% / 9.2% | 98.31% |
| tsa_GAZR | *Stegodyphus mimosarum* | 6.1 | 63,886 | **19.8%** | 19.2% | +0.6 | 1.9% / 1.0% | 90.51% |

## Problems found and how they were handled

**tsa_GAGA (*Anolis*, SRR389086): reads had no real quality scores.** Every quality
character in the files was `?`, which is SRA-lite's placeholder. Trimmomatic stopped
with "Unable to detect quality encoding". ENA serves the same files.

- The reads were rebuilt from NCBI's original-quality "SRA Normalized" object with
  `fasterq-dump` (new `sratools` conda env). Names were rewritten to the ENA style.
- The sequences are identical to the old files; only the quality scores differ.
- The old files are kept in `pairs/tsa_GAGA/reads_sralite/`.
- The ORP run is rerunning as job 1323199. Its original-assembly BUSCO (27.9%) is
  listed above.

**tsa_GAGA and tsa_GAGB original assemblies: duplicate sequence names.** BUSCO
refused both TSA assemblies (job 1322827, tasks 10 and 13).

- Both files were rewritten on 2026-10-16 and rescored in job 1323209, outside this
  set of scripts.
- GAGA dropped from 1,135,207 to 735,899 sequences, and GAGB from 2,056,067 to
  1,195,995. Their original scores are for the rewritten files.

**tsa_GAGE is not *Termitomyces*.** `pair.json` lists *Termitomyces clypeatus*
MTCC 5091, a mushroom (Basidiomycota). ENA lists the reads (SRR746797) as
***Apiospora malaysiana*** (formerly *Arthrinium malaysianum*) strain MTCC 5091, a sac
fungus (Ascomycota).

- Both assemblies were scored against the mushroom dataset `agaricales_odb12.2`.
  That explains why both are low (32.0% and 30.3%).
- They should be rescored against `sordariomycetes_odb12.2`, or `xylariales_odb12.2`
  if it includes *Apiospora*. That hasn't been done yet.
- tsa_GAFV is listed as the same strain and is probably mislabelled too. It isn't in
  these runs.

**Low scores that come from the libraries, not the assembly:**

- ***Stegodyphus* (tsa_GAZR, 19.8%):** 90.5% of reads map back, but the assembly has
  only about 11,000 distinct genes.
- ***Echinogammarus* (tsa_GARO, 53.7%):** about 9,000 distinct genes, with 93.4% of
  reads mapping back.
- ***Eucidaris* (tsa_GAZP, 23.6%):** an ovary-only library of 17.5 million read
  pairs from a GAIIx.
- **The two camels (13.2% and 29.1%):** the smallest mammal read sets, against the
  12,277-gene mammal dataset.

**tsa_GAKE (*Thlaspi arvense*): duplication went up.** The ORP assembly is 98.7%
complete, but 26.5% of BUSCOs are duplicated. The original is 92.1% complete with
only 2.8% duplicated, so single-copy BUSCOs fall from 89.3% to 72.2%. That's the
reverse of every other run.

- The published assembly may have been clustered to one transcript per gene, or the
  ORP assembly may keep redundant copies from the largest read set in the study
  (33 GB). Not yet checked.

**tsa_GAGB (*Anolis*): lower read mapping.** 72.4% of reads map back as proper pairs,
the lowest in the set. Trinity phase 2 had 734,181 components, and the run took
35.8 hours. Not investigated further.

**Samples left out:**

- tsa_GAQK (*Hynobius*, the only amphibian) has 62 bp R2 reads, under the
  manifest's 76 bp minimum.
- *Spraguea lophii* (a microsporidian parasite with a greatly reduced gene set) was
  not chosen for the fungus.
- *Aeginetia indica* (a parasitic plant) was not chosen for the plants.
