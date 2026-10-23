# ORP 4.0.0 lineage comparison: results package

Rebuilt 2026-10-23 from the compare108 handoff (experiments/compare108/handoff/). Supersedes the 21 October package (50 datasets). What changed: nine more ORP 2018 libraries (59 datasets, 15 of them paired with the 2018 assemblies), every previous assembly rescored with pyTransRate 2.2.2 on ORP's own reads, and the release commit question settled. The companion workbook `ORP_4.0.0_lineage_results.xlsx` (Supplementary Data S1) holds every number below.

## Headline results

| Result | Value |
|---|---|
| Datasets | 59 (44 TSA-paired, 15 paired with the ORP 2018 assembly); 57 species; six groups |
| ORP vs TSA: BUSCO complete | median 58.5% → 74.0%; higher in 44/44; median gain 15.0 points |
| ORP vs TSA: unique genes | median 8,942 → 12,140; higher in 44/44 |
| ORP vs TSA: BUSCO duplicated | median 13.8% → 4.8% |
| ORP vs TSA: pyTransRate score | median 0.25 → 0.44; higher in 42/44 |
| ORP vs TSA: good mappings | median 71.4% → 88.5%; higher in 44/44 |
| ORP vs ORP 2018: BUSCO complete | median 75.7% → 86.7%; higher in 15/15 |
| ORP vs ORP 2018: BUSCO duplicated | median 23.8% → 13.0%; lower in 14/15 |
| ORP vs ORP 2018: contigs | median 147,000 → 89,000; fewer in 15/15 |
| ORP vs ORP 2018: read support | good mappings 89.6% vs 88.9% (higher in 9/15, q = 0.85); optimal score lower in 10/15 (q = 0.36) |
| ORP vs best single assembler (59 runs) | more complete in 53 (2 ties), more genes in 57; pyTransRate score higher in only 20 |
| Run time (ORP 4.0.0) | median 7.1 h wall (IQR 4.7–11.3; range 2.1–33.5), 46 GB median peak; 32 cores, 190 GB allocated |
| dev30 vs 4.0.0 (first 50 datasets) | BUSCO complete within 0.4 points; total wall time 409 → 442 h (+8%) |

## What changed since 21 October

- **Nine new datasets**, all from the 2018 ORP paper: *Vigna angularis*, *Cephalotus follicularis*, *Hevea brasiliensis*, *Aeginetia indica* (land plants), *Mus musculus* SRR2086412 (mammals), and *Echinococcus multilocularis*, *Phyllodoce medipapillata*, *Heterorhabditis indica* and *Entamoeba histolytica* (other). All 15 libraries from the 2018 paper are now included, so Table 1b has 15 pairs instead of 6.
- **New BUSCO datasets** for two of them: nematoda_odb12.2 (1,477 genes) for *Heterorhabditis* and amoebozoa_odb12.2 (309) for *Entamoeba*. The two flatworms and the polychaete use lophotrochozoa_odb12.2.
- **Previous assemblies rescored.** Every previous assembly (TSA and ORP 2018) was rescored with pyTransRate 2.2.2 using the trimmed, error-corrected reads of that library's ORP run, the same command and reads ORP uses on its own assembly. Most scores moved by less than 0.02. This closes the old open item and removes the "different pyTransRate release and raw reads" caveat. TSA medians that moved: good mappings 69.1 → 71.4%, pyTransRate score 0.24 → 0.25, optimal score 0.33 → 0.35.
- **Bases uncovered and low-coverage contigs dropped from Table 1.** The handoff tables do not carry them for the rescored previous assemblies. They are replaced by reads mapped (%) and pyTransRate optimal score, keeping 11 metrics. The full rescored outputs are on Premise (compare/pytransrate_orig_orp/) if you want them back.
- **ORP build.** 50 runs used the release candidate e88d037 (then labelled 4.1.0); the 9 new runs used the tagged v4.0.0 release, b9459aa. Between the two commits only resume handling and a startup banner changed.
- **Run time** now covers all 59 runs: median 7.1 h (was 5.9 h on 50). CPU-hours and per-step times exist only for the first 50.

## Design

- **Datasets.** 59 public paired-end RNA-seq runs in six groups (mammals 8, non-mammal vertebrates 9, insects 10, land plants 14, echinoderms 10, other 8). 44 paired with the deposited TSA assembly, 15 with the 2018 ORP orthomerged assembly (Zenodo doi:10.5281/zenodo.1320141). Two species appear twice as different libraries: *Anolis carolinensis* (GAGA, GAGB) and *Mus musculus* (DRR036858, SRR2086412). 6.6 to 206 million read pairs (median 40.0 million). Excluded: tsa_GAQK (*Hynobius*, R2 reads 62 bp).
- **Assembly.** `oyster.py --min-cov auto --lineage <dataset> --mem 190 --cpu 32 --max-parallel 2 --transabyss-mpi on`; no read normalization, no TPM filter. Trimmomatic 0.41, Rcorrector 1.0.7, Trinity 2.15.2 (ORP's bundled build), rnaSPAdes 4.3.0, Trans-ABySS 2.0.1 (MPI), DIAMOND 2.2.5, CD-HIT 4.8.1. Slurm jobs 1323838 (50 runs) and 1324733 (9).
- **Coverage threshold.** `--min-cov auto` gave 1 for 15 libraries (all ≤35 million pairs), above 2 for 15 of the 16 above 60 million (*Vigna*, 60.5 million, got 1.96), up to 5.58 (*Anopheles*, 206 million); median 1.40.
- **BUSCO.** v6.1.0, transcriptome mode, offline, ORP's config.ini. mammalia_odb12 (12,277), vertebrata_odb12 (3,390), insecta_odb12 (3,114) from OrthoDB 12.1; embryophyta (1,658), metazoa (932, for echinoderms), crustacea (2,255), arachnida (1,046), lophotrochozoa (1,329), nematoda (1,477), amoebozoa (309) and agaricales (3,186) odb12.2. Same dataset for all six assemblies of a library.
- **pyTransRate.** 2.2.2 for every assembly, on that library's trimmed, error-corrected reads.
- **Unique genes.** Distinct Swiss-Prot gene names among DIAMOND blastx hits (E ≤ 1e-8, --top 0.1), counted as ORP counts them; component counts are the ones each ORP run reports.
- **Component assemblies.** Trinity, rnaSPAdes auto-k, rnaSPAdes high-k and Trans-ABySS from each ORP run, scored with assembly_eval.py (Slurm jobs 1324263 and 1324923).
- **Statistics.** Two-sided Wilcoxon signed-rank tests on paired values; Benjamini-Hochberg within each Table 1 panel (11 metrics) or within each assembler (11 metrics). Medians with IQR.
- **Known caveats.** GAGE reads are *Apiospora malaysiana* (ascomycete) but both assemblies were scored against agaricales_odb12.2. TSA assemblies are curated (no contigs < 200 bp). The two *Anolis* TSA assemblies had duplicate sequence names and were de-duplicated before scoring (GAGA 1,135,207 → 735,899 sequences; GAGB 2,056,067 → 1,195,995). Parasites with reduced gene sets (*Aeginetia*, *Entamoeba*, the flatworms) score low whatever the assembler; *Entamoeba*'s low gene count (about 2,400) mostly reflects few close Swiss-Prot matches.

## Table 1a. TSA deposited vs ORP 4.0.0 (44 datasets)

Medians (IQR); q = Benjamini-Hochberg-adjusted Wilcoxon signed-rank P across the 11 metrics.

| Metric | TSA deposited | ORP 4.0.0 | q | ORP higher in |
|---|---|---|---|---|
| Contigs (thousands) | 56.7 (39.0–85.3) | 106.9 (63.3–153.4) | <0.001 | 40/44 |
| Reads mapped (%) | 85.2 (72.7–90.8) | 96.1 (92.7–97.0) | <0.001 | 44/44 |
| Good mappings (%) | 71.4 (45.3–81.7) | 88.5 (81.5–93.0) | <0.001 | 44/44 |
| pyTransRate score | 0.25 (0.12–0.34) | 0.44 (0.31–0.50) | <0.001 | 42/44 |
| pyTransRate optimal score | 0.35 (0.22–0.43) | 0.50 (0.43–0.56) | <0.001 | 42/44 |
| Unique Swiss-Prot genes | 8,942 (7,285–11,690) | 12,140 (9,023–14,109) | <0.001 | 44/44 |
| BUSCO complete (%) | 58.5 (30.0–71.3) | 74.0 (54.5–91.6) | <0.001 | 44/44 |
| BUSCO single-copy (%) | 38.6 (21.0–51.5) | 68.9 (50.2–77.2) | <0.001 | 42/44 |
| BUSCO duplicated (%) | 13.8 (3.2–23.2) | 4.8 (2.5–9.1) | 0.002 | 14/44 |
| BUSCO fragmented (%) | 15.2 (9.5–25.3) | 9.5 (3.5–13.8) | <0.001 | 5/44 |
| BUSCO missing (%) | 27.4 (10.1–41.6) | 10.9 (4.7–27.4) | <0.001 | 0/44 |

## Table 1b. ORP 2018 vs ORP 4.0.0 (15 datasets)

| Metric | ORP 2018 | ORP 4.0.0 | q | ORP higher in |
|---|---|---|---|---|
| Contigs (thousands) | 147.4 (95.6–303.0) | 88.8 (56.4–165.3) | <0.001 | 0/15 |
| Reads mapped (%) | 97.0 (95.5–97.3) | 96.8 (96.0–97.4) | 0.131 | 11/15 |
| Good mappings (%) | 89.6 (84.4–92.3) | 88.9 (84.7–91.3) | 0.847 | 9/15 |
| pyTransRate score | 0.28 (0.22–0.33) | 0.42 (0.37–0.46) | <0.001 | 14/15 |
| pyTransRate optimal score | 0.49 (0.43–0.52) | 0.46 (0.45–0.50) | 0.363 | 5/15 |
| Unique Swiss-Prot genes | 11,843 (9,822–12,360) | 12,785 (10,014–13,664) | 0.049 | 13/15 |
| BUSCO complete (%) | 75.7 (43.7–87.3) | 86.7 (67.5–92.9) | <0.001 | 15/15 |
| BUSCO single-copy (%) | 39.8 (30.5–48.3) | 64.7 (58.5–75.8) | <0.001 | 15/15 |
| BUSCO duplicated (%) | 23.8 (15.9–41.7) | 13.0 (7.0–20.4) | 0.003 | 1/15 |
| BUSCO fragmented (%) | 8.4 (4.8–14.5) | 3.8 (1.6–8.4) | <0.001 | 1/15 |
| BUSCO missing (%) | 22.4 (7.1–37.3) | 6.9 (4.5–21.0) | 0.002 | 0/15 (1 tie) |

Reading: the 2018 assemblies were also built from error-corrected reads and are the closest competitors in the set. All six libraries where the previous assembly has more good mappings are ORP 2018 assemblies, as are 10 of the 12 where its optimal score is higher (the other two are TSA *Thlaspi arvense*, GAKE, and *Asterias amurensis*, GAVJ). Duplication rose only in *Phyllodoce medipapillata* (3.0% → 19.0%). Genes fell in *Nicotiana tabacum* (16,423 → 14,786) and *Echinococcus multilocularis* (8,512 → 7,726). The pyTransRate score fell only in *Hevea brasiliensis* (0.30 → 0.22).

## By group (all 59 pairs)

| Group | Datasets (TSA / ORP 2018) | Median BUSCO complete, previous (%) | Median BUSCO complete, ORP (%) | Median gain (points) | Median BUSCO duplicated, previous (%) | Median BUSCO duplicated, ORP (%) | Duplication fell in | Median unique genes gain |
|---|---|---|---|---|---|---|---|---|
| mammalia | 8 (5 / 3) | 43.7 | 54.4 | 6.2 | 15.4 | 1.0 | 8 | 1,601 |
| vertebrata | 9 (8 / 1) | 40.2 | 64.5 | 21.7 | 15.7 | 3.8 | 7 | 1,972 |
| insecta | 10 (9 / 1) | 63.9 | 77.8 | 15.0 | 12.2 | 8.9 | 7 | 1,574 |
| embryophyta | 14 (9 / 5) | 62.8 | 92.2 | 20.2 | 12.8 | 14.3 | 7 | 1,710 |
| echinodermata | 10 (10 / 0) | 66.7 | 84.2 | 12.1 | 24.2 | 4.4 | 10 | 1,082 |
| other | 8 (3 / 5) | 53.0 | 85.2 | 8.0 | 13.3 | 9.7 | 5 | 485 |

Land plants are the one group where ORP does not reduce duplication: it rose in all seven of the plant datasets where it rose at all, every one a TSA pair.

## ORP vs each component assembler (59 runs)

Cells: runs in which ORP is better (median difference, ORP minus assembler; percentage points for %).

| Metric | Trinity | rnaSPAdes auto-k | rnaSPAdes high-k | Trans-ABySS | Best single |
|---|---|---|---|---|---|
| BUSCO complete (%) | 54 (+2.9) | 57 (+1.4) | 59 (+7.2) | 58 (+5.9) | 53, 2 ties (+1.1) |
| BUSCO duplicated (%) | 59 (−28.3) | 58 (−13.7) | 54 (−7.1) | 42 (−5.0) | 38, 2 ties (−2.2) |
| Unique Swiss-Prot genes | 59 (+911) | 59 (+639) | 59 (+1,147) | 57 (+998) | 57 (+434) |
| pyTransRate score | 54 (+0.056) | 46 (+0.079) | 24 (−0.011) | 55 (+0.087) | 20 (−0.014) |
| pyTransRate optimal score | 47 (+0.030) | 11 (−0.014) | 23 (−0.012) | 55 (+0.067) | 7 (−0.019) |
| Good mappings (%) | 45 (+2.1) | 12 (−1.3) | 18 (−0.6) | 59 (+3.6) | 11, 1 tie (−1.3) |
| Proper pairs (%) | 52 (+2.0) | 13 (−1.1) | 17 (−0.6) | 54 (+2.8) | 10 (−1.2) |

Best single assembler, by metric: completeness rnaSPAdes auto-k 34, Trinity 24, rnaSPAdes high-k 1; genes rnaSPAdes auto-k 31, Trans-ABySS 20, rnaSPAdes high-k 5, Trinity 3; pyTransRate score rnaSPAdes high-k 43; good mappings rnaSPAdes auto-k 34, rnaSPAdes high-k 24. Against rnaSPAdes high-k the pyTransRate score difference is not significant (q = 0.11); every other component comparison above has q < 0.02.

Exceptions on completeness: ORP is beaten in four runs, each by 0.4 points or less and each by an assembly two to five times more duplicated: *Nicotiana tabacum* (Trinity 97.6% vs 97.4%; D 88.2% vs 30.0%), *Avena sativa* (Trinity 91.8% vs 91.6%; D 78.3% vs 15.7%), *Echinogammarus veneris* (rnaSPAdes auto-k 53.9% vs 53.7%; D 12.2% vs 5.1%) and GAGE (Trinity 32.4% vs 32.0%; D 13.0% vs 2.9%). It ties Trinity on *Kerria lacca* (88.3%; D 57.1% vs 15.0%) and *Heterorhabditis indica* (93.0%; D 74.5% vs 29.6%). Trans-ABySS recovers more genes than ORP in *Patiria miniata* (+134) and *Anopheles gambiae* (+77).

### Genes missing from the ORP count (merge trace, 59 runs; S1 sheet Gene_masking)

Genes among some component's DIAMOND hits but not ORP's: median 798 per run (155–1,821), 45,699 in all. Decision (2026-10-23): report as a property of the gene count, not a merge bug.

| Cause | Genes | Share |
|---|---|---|
| Masked: sequence inside an ORP contig whose best hit beats the gene's by >10% (--top 0.1) | 35,426 | 78% |
| Inside an ORP contig with no hit (marginal hit fails -e 1e-8 on the longer query) | 1,718 | 4% |
| Inside an ORP contig, gene within 10% but not reported (unresolved) | 3,231 | 7% |
| No ORP match at ≥95% identity over ≥90% of the contig (upper bound on sequence lost) | 5,324 | 12% |

No match per run: median 77 (9–303); land plants median 149, other groups 16–80; *Avena sativa* 267 (second of 59, after ERR1016675 at 303). *P. miniata*: 1,094 of 1,329 masked, 91 no match. *A. gambiae*: 963 of 1,136 masked, 88 no match. Two-track selection keeps every best-hit gene; the removal is the final cd-hit-est -c .98, which clusters across genes. BUSCO side of the trace (1,061 BUSCOs complete in a component, F/M in ORP: 79% copy_selection, 20% seq_present, 12 gene_dropped, 3 nohit_path) is not yet written into the draft.

### Land plants and polyploidy

In all 14 land plants ORP recovered more distinct genes than any component, and it was the most complete assembly in 12 and within 0.2 points in the other two (*A. sativa*, *N. tabacum*). Concordant mapping trailed the best single assembler by 15.5 points in *Avena sativa* (79.2% vs 94.7%, rnaSPAdes high-k), 8.7 in *Brassica juncea*, 5.0 in *Musa* ABB and 4.8 in *Nicotiana tabacum*, against a median of 3.0 across plants and 0.8 outside them. Outside plants the largest gap is *Stegodyphus mimosarum* (7.3). Trinity and rnaSPAdes are 51–91% duplicated in *B. juncea*, *N. tabacum* and *A. sativa*; Trans-ABySS stays at 3.7–11.0% with the weakest read support (60.4% good mappings in *B. juncea*). The diploid *Thlaspi arvense* is 24.4% duplicated (Trinity and rnaSPAdes 56–73%). The four plants added in this round (*Vigna*, *Cephalotus*, *Hevea*, *Aeginetia*) fit the same pattern at lower amplitude (gaps 1.6–3.4 points; ORP duplication 10.5–23.9%).

## Settings: dev30 vs 4.0.0 (first 50 datasets)

Unchanged; dev30 was run only on the original 50. Trinity normalization on and fixed `--min-cov 2` in dev30. BUSCO completeness moved at most 0.4 points; unique genes −24 to +72; pyTransRate score median −0.005. Trinity time 291 → 343 h and total wall time 409 → 442 h (+8%; longer in 39 of 50). On the nine of those libraries above 60 million read pairs the depth-scaled filter (2.1–5.6) removed up to 35% of contigs with BUSCO within 0.2 points.

## Run time (ORP 4.0.0, 32 cores, 59 runs)

Median 7.1 h wall (IQR 4.7–11.3; range 2.1–33.5), median peak memory 46 GB (max 180 GB); wall time vs read pairs r = 0.83 (log-log). TSA libraries median 5.5 h; ORP 2018 libraries, which are deeper, median 9.8 h. Median CPU-hours (128) is for the first 50 runs only.

## Open items

- Rescore GAGE against sordariomycetes_odb12.2 (reads are an ascomycete), or drop it.
- Restore bases-uncovered and low-coverage-contig columns for the rescored previous assemblies if the paper should keep those two metrics.
- CPU-hours and per-step times for the 9 new runs (sacct for job 1324733) if the paper reports CPU-hours.
- Rename "ORP version 5" / "ORP v5" in the handoff figure and workbook; the Figure 1 here already reads 4.0.0.
- Decide whether Trinity normalization should stay off by default (+8% wall time, no quality gain).
- Optional: rerun the merge on the polyploid plants with a higher per-gene copy cap (twotrack_select.py --max-per-gene, default 3).
- Resolved 2026-10-23: genes missing from the ORP count are reported as a property of the count (Section 4.3, Table 2 caption, polyploid paragraph, S1 sheet Gene_masking).
- Resolved: previous assemblies rescored with pyTransRate 2.2.2; v4.0.0 tag is b9459aa, functionally identical to e88d037.
