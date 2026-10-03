# Where ORP's extra duplicated BUSCOs come from

Compared with the original assemblies, the ORP re-assemblies have a higher
median duplicated-BUSCO rate (17.6% to 28.8%) and ~30,000 more contigs. These
tests ask whether that comes from resolving isoforms or from keeping the same
transcript more than once.

## `dup_busco_pairs.py`: what are the copies behind each duplicated BUSCO?

For each duplicated BUSCO, it aligns every pair of its contigs with blastn
(both strands) and puts each pair in one class: `isoform` (near-identical
with an internal indel of 20 bp or more), `redundant` (near-identical across
90% or more of the shorter contig), `ends` (near-identical but staggered),
`variant` (90-98% identity) or `divergent`. The docstring gives the exact
rules. `run_dup_pairs_all.sh` runs it on the ORP and original assembly of
every manifest sample and stacks the results into `all.summary.tsv` and
`all.pairs.tsv`.

```
conda run --no-capture-output -n orp bash run_dup_pairs_all.sh ~/redundancy_tests/dup_pairs
```

Takes about 8 minutes for all 104 pairs on a login node.

First results (2026-09-30, 104 pairs). Each duplicated BUSCO is counted once,
under the strongest class among its pairs:

| Why the BUSCO is duplicated | Original | ORP | Change |
| --- | --- | --- | --- |
| same sequence twice (only redundant or staggered pairs) | 1,286 | 1,217 | -69 |
| isoform-like (at least one isoform pair) | 862 | 1,891 | +1,029 |
| variant (90-98% identity, no isoform pair) | 300 | 792 | +492 |
| divergent | 56 | 84 | +28 |

Pure redundancy is about as common in ORP as in the originals. ORP drops
some of the originals' redundant copies and adds others: 36% of the BUSCOs
that are duplicated only in ORP have nothing but redundant copies. The net
increase is in isoform-like and variant pairs. Two caveats:

- `isoform` is a generous call. 786 of ORP's 1,891 isoform-like BUSCOs rest
  only on pairs from two *different* assemblers, which can simply disagree
  about structure. Real isoforms should show up in at least two assemblers.
- The TPM share of the minor copy (`--contigs-csv`) does not separate the
  classes. Salmon splits reads evenly between redundant copies, so they look
  as well supported as isoforms do.

## `orthofinder_strand_test.py`: can the orthogroup step see both strands?

No. ORP runs `orthofinder -d ... -S diamond_orp_<threads>`. In the installed
OrthoFinder 3.1.5, `-d` only switches off the "looks like nucleotides" input
check (`utils/fasta_processor.py`). It does not choose the search program:
that comes from `-S` and defaults to `diamond` (`run/process_args.py:107`).
So the search is diamond **blastp** over the DNA, which covers one strand
only, and the run's `Log.txt` says `Search program: diamond`. OrthoFinder's
config does include a `blastn` entry (`makeblastdb -dbtype nucl`, then
`blastn`). That is the command a maintainer quotes for `-d` in
OrthoFinder/OrthoFinder#54, but here it only runs when you pass `-S blastn`. The test
takes N real contigs, makes forward, reverse-complement and 1%-mutated
copies of each, and runs OrthoFinder with ORP's flags:

```
conda run --no-capture-output -n orp_orthofinder python orthofinder_strand_test.py \
    --fasta <run>.trinity.Trinity.fasta.gz --workdir strand_test [--orient] [--search blastn]
```

| Copies that should share an orthogroup | `-I 12` (ORP) | `-I 1.5` | `-I 12`, oriented by longest ORF | `-I 12 -S blastn` |
| --- | --- | --- | --- | --- |
| same strand, 1% substitutions | 1.00 | 1.00 | 0.85 | 0.97 |
| opposite strand, exact copy | **0.01** | **0.01** | 0.98 | 0.96 |
| all four copies together | 0.01 | 0.01 | 0.83 | 0.95 |
| orthogroups (should be 200) | 403 | 403 | 235 | 218 |

Copies of one transcript in opposite orientations never compete in
`makeorthout`, so both are kept. The one thing that can still merge them is
cd-hit-est, which checks both strands but only merges a contig that matches
across its whole length. In `dup_busco_pairs.py`, 55% of ORP's
redundant/staggered pairs align on opposite strands. Only 240 of its 1,939
redundant pairs would pass cd-hit's global-identity test at `-c .98`.

`-S blastn` is the simplest fix, and it does better here than orienting by
longest ORF. That approach fixes the exact copies but is fragile: 1%
substitutions flip the call for about 15% of contigs. Its cost at full ORP
scale is still untested. The stock `blastn` entry has no thread flag, so each
of the 16 searches runs on one thread. The N-to-X masking in
`write_search_inputs` is only there for diamond, so blastn should get the
original Ns back.

## `time_orthofinder_search.sbatch`: is blastn affordable?

One median sample (SRR1139197, whose `run_orthofuser` took 49 minutes in the
real run). Its four real assemblies are filtered to contigs over 200 bp
(425,015 in all), and each search runs on the same node with 24 cpus and
120G: diamond as ORP's planner would set it (`-t 10 -S diamond_orp_2`), then
`-t 16 -S blastn`. `compare_orthogroups.py` compares the results.

| | diamond (ORP now) | blastn |
| --- | --- | --- |
| all-vs-all search | 5m 58s | 4m 16s |
| start to `Done orthogroups` | 16m 31s | 12m 30s |
| whole OrthoFinder run (see below) | 35m 12s | 35m 53s |
| CPU time | 3h 35m | 3h 48m |
| peak RSS (sacct MaxRSS) | 9.6 GB | 36.3 GB |
| orthogroups | 267,315 | 178,799 |
| singletons | 152,785 | 36,202 |

blastn is no slower, and it reaches orthogroups four minutes sooner. It
leaves a third fewer orthogroups and a quarter as many singletons, and
makeorthout keeps at most one contig per orthogroup. The cost is memory:
3.8x the peak. sacct doesn't say which phase that peak comes from, so check
how it scales on the large assemblies before switching. The chowder
all-vs-all already peaks near 550 GB with diamond.

This sample has 29 duplicated-BUSCO contig pairs in its final ORP assembly.
None of them share a diamond orthogroup. blastn would put 31% of them
together: 2 of 5 redundant pairs, 5 of 11 isoform-like pairs, 1 of 4
variant pairs and 1 of 9 staggered pairs. The rest stay split even with
blastn, presumably because `-I 12` splits them. That can be tested without
repeating the search, by re-clustering from the saved blastn hits at a
lower `-I`. That sample is small, so treat the 31% as a direction, not an
estimate.

Both runs spend more than half their time after `Done orthogroups`, on
MSA/trees, STRIDE and hierarchical orthogroups, which ORP never uses. That
is the "`-og` is not honoured" problem already in NOTES.md, and it applies
equally to both search programs.

## `inflation_sweep.sbatch`: does a lower `-I` help?

This re-clusters the two saved searches above at `-I` 1.5, 2, 3, 5, 8 and 12
with `orthofinder -b`. Each run takes 6-11 minutes on 4 cpus and stops once
the orthogroups are written. `busco_pool.sbatch` runs BUSCO over the pooled
425,015 contigs, finding 106 genes in 401 contigs. `compare_orthogroups.py
--busco` then counts **split genes** (one BUSCO gene in several orthogroups,
so potential extra copies) and **merged groups** (two BUSCO genes in one
orthogroup, so `makeorthout` drops one).

| | dia `-I 12` (ORP now) | dia `-I 1.5` | bn `-I 12` | bn `-I 5` | bn `-I 3` | bn `-I 2` | bn `-I 1.5` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| orthogroups | 266,927 | 252,571 | 177,572 | 162,955 | 151,776 | 143,199 | 134,987 |
| singletons | 152,741 | 148,082 | 36,036 | 24,550 | 18,736 | 17,737 | 17,706 |
| BUSCO genes split | 57 | 56 | 37 | 28 | 23 | 20 | 12 |
| extra orthogroups for them | 114 | 82 | 71 | 47 | 36 | 28 | 15 |
| merged BUSCO groups | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| duplicate pairs grouped (n=29) | 0.00 | 0.21 | 0.38 | 0.48 | 0.59 | 0.62 | 0.79 |
| groups of more than 20 contigs | 6 | n/a | 6 | n/a | 27 | 70 | 166 |
| largest group | 33 | n/a | 33 | n/a | 74 | 112 | 197 |

(dia = diamond, bn = blastn. "n/a" means not measured. The re-run at
diamond `-I 12` matches the original run, 266,927 groups against 267,315.
The re-run at blastn `-I 12` has 177,572 groups against the original run's
178,799.)

- **For diamond, inflation hardly matters.** A missed opposite-strand hit
  cannot be fixed by clustering, and 57 BUSCO genes stay split at every
  setting. Lowering `-I` alone is not a fix.
- **For blastn, it matters a lot.** Split BUSCO genes fall from 37 to 12,
  and their extra orthogroups from 71 to 15. Combined with the strand fix,
  that is 114 extra orthogroups down to 15 against ORP now.
- **BUSCO shows no over-merging at any setting, but that is a weak test.**
  BUSCO genes are conserved single-copy genes with few close paralogs, and
  close gene families are where over-merging would happen. The growth of
  very large groups at low `-I` (166 groups of more than 20 contigs at 1.5,
  against 6 at 12) is the warning sign. `makeorthout` keeps one contig from
  each, so any real gene family in there loses all but one member.
- **Candidate setting: blastn with `-I 2` or `-I 3`.** That keeps most of
  the gain and only a few dozen large groups. The final choice needs the
  end-to-end test below, because it depends on what BUSCO completeness and
  duplication do in the finished assembly, not on orthogroup counts.

## `chowder_arms.sbatch`: end to end on SRR1139197

This runs the whole merge half three ways through `chowder.py`, using the
same four assemblies, the same corrected reads, `--tpm-filt 1` and 24 cpus
throughout. The OrthoFinder settings come from two hidden options,
`--orthofinder-program` and `--orthofinder-inflation`, run from a separate
copy of the code. Each run takes about 1h10m.

The control reproduces the original ORP run of this sample exactly: same
BUSCO to the decimal, and the same 29 duplicate pairs.

| | diamond `-I 12` (control) | blastn `-I 2` | blastn `-I 3` |
| --- | --- | --- | --- |
| orthogroups (= makeorthout picks) | 267,313 | 149,705 | 155,656 |
| rescued (list7) | 1,068 | 1,658 | 1,600 |
| after cd-hit | 184,518 | 138,413 | 142,501 |
| **final contigs** | **167,672** | **136,012** (-19%) | **139,598** (-17%) |
| total bases | 66.2 Mb | 51.5 Mb | 53.4 Mb |
| BUSCO C / S / D | 51.2 / 36.8 / **14.4** | 49.6 / 40.0 / **9.6** | 49.6 / 38.4 / **11.2** |
| BUSCO F / M | 28.0 / **20.8** | 25.6 / **24.8** | 27.2 / **23.2** |
| transrate score | 0.410 | 0.458 | 0.460 |
| transrate optimal score | 0.582 | 0.570 | 0.582 |
| good mappings | 0.930 | 0.905 | 0.910 |
| low-covered contigs | 4.6% | 1.7% | 1.9% |

In BUSCO genes, `-I 2` turns 6 duplicated BUSCOs into 4 single and 2
missing; 3 fragmented go missing as well, so 5 more are missing in all.
`trace_arms.py` follows each change back through the merge:

- **Every remaining duplicate is two picks from separate orthogroups.**
  The diamond rescue adds none.
- **The lost BUSCOs are not over-merging.** In each of the 5, the orthogroup
  is right: one gene, assembled by several assemblers. The difference is
  which member `makeorthout` keeps. Looking contig by contig
  (`repick.py --compare --watch`):

  | Lost BUSCO (in control) | Contig that carried it | Picked instead |
  | --- | --- | --- |
  | 5000870 (fragmented) | 1,574 bp, score 0.853 | 222 bp, score 0.879 |
  | 5003022 (fragmented) | 355 bp, score 0.140 | 211 bp, score 0.822 |
  | 5005196 (complete) | 1,396 bp, score 0.061 | 2,371 bp, score 0.243, shorter ORF |
  | 5002573 (fragmented) | 387 bp, score 0.058 | 459 bp, score 0.855 |
  | 319730 (fragmented) | 3,287 bp isoform, score 0.01 (no reads) | the other isoform, 2,572 bp |

  Two are short fragments beating a longer contig. Two are longer or
  similar-length picks that do not carry the BUSCO match. One is a choice
  between isoforms. Under diamond, each of these contigs sat alone in an
  orthogroup built from a one-strand search, so it survived. Four of the
  five were only fragmented BUSCOs in the control.

## `repick.py` / `repick_arms.sbatch`: a length-aware pick rule

`repick.py` re-picks from an arm's saved `Orthogroups.txt` and `contigs.csv`
under `score` (today's rule, which reproduces each arm's `good.<run>.list`
byte for byte), `score_len` (score x length), `score_orf` (score x ORF
length) or `near_best` (longest member scoring at least 0.8 of the group's
best). Each rule changes 6,500-12,500 picks per arm, about 7-8% of picks
under blastn, and raises the total length picked by 3-7%. Of the 5 lost
contigs above, all three new rules recover only 5000870.
`repick_arms.sbatch` reruns everything after the pick under each rule, for
all three arms. Task 1 (blastn_I2 under `score`) must reproduce that arm
exactly.

### Results (job 1319548, 2026-10-01)

The sanity task reproduced the blastn `-I 2` arm exactly, so the re-pick
runs are directly comparable to the arms. BUSCO is in genes out of 125
(0.8% = 1 gene).

| Search | Pick rule | Final contigs | BUSCO C | D | M | Transrate score | Optimal | Good mappings |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| diamond `-I 12` | score (ORP now) | 167,672 | 64 | 18 | 26 | 0.410 | 0.582 | 0.930 |
| diamond `-I 12` | score_len | 165,841 | 64 | 19 | 24 | 0.407 | 0.573 | 0.932 |
| diamond `-I 12` | score_orf | 166,868 | 64 | 19 | 24 | 0.408 | 0.563 | 0.932 |
| diamond `-I 12` | near_best | 166,871 | 64 | 18 | 24 | 0.409 | 0.583 | 0.931 |
| blastn `-I 2` | score | 136,012 | 62 | 12 | 31 | 0.458 | 0.570 | 0.905 |
| blastn `-I 2` | score_len | 133,844 | 63 | 14 | 30 | 0.471 | 0.575 | 0.918 |
| blastn `-I 2` | score_orf | 134,806 | 63 | 13 | 28 | 0.464 | 0.563 | 0.914 |
| blastn `-I 2` | near_best | 135,141 | 62 | 13 | 30 | 0.468 | 0.573 | 0.913 |
| blastn `-I 3` | score | 139,598 | 62 | 14 | 29 | 0.460 | 0.582 | 0.910 |
| blastn `-I 3` | score_len | 137,538 | 64 | 15 | 28 | 0.470 | 0.573 | 0.921 |
| **blastn `-I 3`** | **score_orf** | **138,506** | **64** | **14** | **26** | **0.465** | **0.590** | 0.917 |
| blastn `-I 3` | near_best | 138,724 | 62 | 15 | 28 | 0.468 | 0.571 | 0.917 |

- **Under blastn, every length-aware rule helps.** Each recovers 1-3 missing
  BUSCOs, raises the transrate score by 0.006-0.013 and good mappings by
  about 0.01, and gives up at most 2 of the duplicates blastn removed.
- **blastn `-I 3` with `score_orf` is the best combination here.** Against
  ORP now, it is equal on completeness (64) and missing (26), has 4 fewer
  duplicated BUSCOs (14 against 18), 17% fewer contigs, a transrate score
  0.055 higher and a higher optimal score. The one metric that falls is good
  mappings, 0.917 against 0.930, as expected for an assembly 18% smaller.
- **Under diamond, the rules barely matter**: 2 fewer missing, 0-1 more
  duplicated, and transrate unchanged.
- **This is one sample, and the differences are 1-5 genes.** The pipeline is
  deterministic (the sanity task matched to the contig), so these are real
  differences for this sample. But whether they hold across samples is
  exactly what a single sample cannot show.

## `validate_arms.sbatch`: blastn `-I 3` + `score_orf` on 10 samples

`--pick-rule` is now a third hidden option, next to the two OrthoFinder ones
(`add_merge_experiment_args` in `oyster.py`, shared with `chowder.py`). It
passes `--rule` to `scripts/pick_best_contigs.py`, which now implements
`score_len`, `score_orf` and `near_best` itself. Its default output is
unchanged, and each rule matches `repick.py` byte for byte on all three
SRR1139197 arms.

The 10 samples are roughly one per decile of total assembly size among the
113 runs that kept their assemblies and corrected reads, with the largest
standing in for the top decile. Each was originally run with
`--tpm-filt 1 --cpu 24 --mem 120` (`validate_samples.txt`):

| Sample | Assemblies (gz) | Reads (gz) | Original orthogroup step |
| --- | --- | --- | --- |
| SRR954929 | 24 MB | 4.2 GB | 0:08 |
| SRR544889 | 41 MB | 3.3 GB | 0:37 |
| SRR747027 | 51 MB | 7.9 GB | 0:45 |
| SRR1139198 | 59 MB | 7.1 GB | 0:33 |
| SRR1060332 | 65 MB | 2.1 GB | 0:42 |
| SRR1176880 | 72 MB | 4.8 GB | 0:29 |
| SRR1174805 | 89 MB | 4.7 GB | 1:12 |
| SRR1123893 | 104 MB | 9.9 GB | 1:07 |
| SRR951913 | 138 MB | 8.9 GB | 0:56 |
| DRR036858 | 211 MB | 21.3 GB | 5:07 (largest) |

Each sample gets a control (diamond, `-I 12`, `score`) and a candidate
(blastn, `-I 3`, `score_orf`), both through `chowder.py` with the original
settings. DRR036858's tasks get a 400G Slurm allocation, so blastn's memory
on the largest sample is measured rather than cut short. chowder's own
`--mem` stays 120 in both arms. `validate_summary.py` makes the table:

```
python3 validate_summary.py validate_arms validate_samples.txt 1319575 1319576
```

### Results: 9 samples plus SRR1139197 (2026-10-01)

DRR036858, the largest, was still running when this was written (see
below). BUSCO is in genes out of 125. "Δ" is candidate minus control.

| Sample | Contigs Δ | Complete | Duplicated | Missing | Transrate | Peak GB (control → candidate) |
| --- | --- | --- | --- | --- | --- | --- |
| SRR954929 | -9.6% | 90 → 91 | 8 → 5 | 14 → 14 | 0.560 → 0.590 | 12 → 33 |
| SRR544889 | -8.2% | 97 → 98 | 18 → 16 | 14 → 14 | 0.518 → 0.522 | 13 → 33 |
| SRR747027 | -0.5% | 123 → 124 | 47 → 41 | 0 → 0 | 0.535 → 0.577 | 22 → 72 |
| SRR1139198 | -15.8% | 51 → **49** | 9 → **12** | 28 → **34** | 0.401 → 0.441 | 15 → 44 |
| SRR1060332 | -10.9% | 111 → **108** | 10 → **11** | 4 → **6** | 0.410 → 0.483 | 19 → 72 |
| SRR1176880 | -4.9% | 124 → 125 | 65 → **73** | 1 → 0 | 0.510 → 0.552 | 28 → 101 |
| SRR1174805 | -8.4% | 125 → 125 | 83 → 82 | 0 → 0 | 0.467 → 0.513 | 22 → 102 |
| SRR1123893 | -6.9% | 124 → 124 | 49 → 49 | 0 → 0 | 0.495 → 0.540 | 42 → 117 |
| SRR951913 | -4.0% | 111 → **108** | 35 → 20 | 3 → **4** | 0.345 → 0.377 | 40 → 112 |
| SRR1139197 | -17.4% | 64 → 64 | 18 → 14 | 26 → 26 | 0.410 → 0.465 | (other runs) |
| **Total / median** | **median -8.3%** | **-4 in total** | **342 → 323** | **90 → 98** | **10/10 up, median +0.044** | **2.6-5.4x** |

- **Transrate goes up in every sample** (median +0.044) and the assembly
  shrinks (median -8% contigs), with good mappings unchanged (median
  0.000). By transrate's measure the candidate is consistently better.
- **The duplicated-BUSCO reduction is small, and it is not consistent.**
  Over all 10 samples, duplicated BUSCOs fall 6% (342 to 323). Six samples
  fall, three rise and one is unchanged, and two samples provide most of
  the drop (SRR951913 -15, SRR747027 -6). Samples whose duplication is
  mostly biological, like SRR1176880 (half its duplicate pairs are
  divergent: paralogs, or a polyploid genome), cannot and should not lose
  it.
- **BUSCO completeness costs a little in some samples.** Missing BUSCOs
  rise in 3 samples (+6, +2, +1) and fall in 1; complete BUSCOs total -4
  over 10 samples. All 10 lost BUSCOs trace to the same step
  (`trace_arms.py`): the BUSCO-carrying contig sat in a correct
  multi-assembler orthogroup and was not the one `makeorthout` kept. Seven
  were only fragmented in the control. Three complete BUSCOs were lost,
  and 5005196 is lost in two samples.
- **blastn is expensive.** Peak memory is 2.6-5.4x the control's (33-117
  GB against 12-42 GB). Wall time is mixed on mid-sized samples (each arm
  faster in some samples), but on DRR036858 (1.01 M contigs) the blastn
  search alone runs about 8 hours, against under 2 for diamond.
- **Two failure modes appear only under blastn, and both have fixes.**
  OrthoFinder's 200 s stall watchdog (SRR544889; fixed by `-a` = number of
  assemblies, see `validate_arms.sbatch`), and its unused species-tree
  stage running past 120 GB (SRR1123893, where FastTree was OOM-killed
  after the orthogroups were written). The second is harmless here but
  needs the `-og` problem fixed before blastn can be a default.

### Where this leaves the decision

- **The strand bug is a regression, not a design choice.** OrthoFinder
  2.5.2's `-d` ran blastn (`options.search_program = "blast_nucl"` unless
  `-S` was given), and ORP's `-d -I 12` dates from 2021, under 2.5.2. The
  2026-08-14 move to 3.1.5 silently made the search one-strand diamond;
  NOTES.md recorded the side effect (`run_orthofuser` 5:15 -> 30:57 on
  SRR1789336) without the cause. Fixing it is justified on correctness
  alone.
- **blastn `-I 3` + `score_orf` is not ready to be the default.** It wins
  on transrate and size, but its BUSCO effect is mixed, and the loss
  mechanism (the pick) is still active.
- **Next, blastn `-I 12` with today's pick rule** on these same samples:
  literally the behaviour ORP had before the upgrade, with fewer pick
  decisions forced on it than at `-I 3`. If it keeps the transrate gain
  without the BUSCO losses, it is the default, shipped as a regression fix.
- **Two cheaper or better-targeted options are worth a test:**
  - Orient contigs before a diamond search, which gets the strand fix at
    diamond's cost. Orienting by swissprot blastx hit where there is one
    is more robust than by longest ORF, which flips about 15% of contigs
    under 1% substitutions (`orthofinder_strand_test.py`).
  - A pick rule that prefers protein evidence. Every BUSCO loss was the
    pick keeping a contig without the gene's protein match. ORP already
    has each assembly's swissprot diamond hits, so a member with a longer
    hit could win. `repick.py` can test that without new searches.

## Synthesis: search, inflation and pick rule (10 samples, 2026-10-02)

Three arms per sample, all through `chowder.py` with the original settings:

- **control:** diamond, `-I 12`, pick rule `score` (ORP now)
- **blastn_I12:** blastn, `-I 12`, `score` (ORP before the OrthoFinder
  3.1.5 upgrade)
- **candidate:** blastn, `-I 3`, `score_orf`

The samples are the 9 in `validate_samples.txt` plus SRR1139197. DRR036858,
the largest, has no finished arm: both blastn arms fail at OrthoFinder's
200 s stall watchdog even at `-a 4`, and its diamond control was cancelled
an hour from the end. BUSCO is in genes; "uniq" is qualreport's UNIQUE
GENES ORP.

| Totals over 10 samples | control | blastn_I12 | candidate |
| --- | --- | --- | --- |
| BUSCO complete | 1,020 | 998 (-22) | 1,016 (-4) |
| BUSCO duplicated | 342 | 352 (+10) | 323 (-19) |
| BUSCO missing | 90 | 110 (+20) | 98 (+8) |
| samples with more missing / fewer | - | 8 / 1 | 3 / 1 |
| unique genes | 96,168 | +425 (10/10 up) | +336 (10/10 up) |
| transrate score, median change | - | +0.024 (9/10 up) | +0.042 (10/10 up) |
| contigs, median change | - | -3,102 | -6,874 |
| peak memory | 12-42 GB | 33-109 GB | 33-117 GB |

**The search (control vs blastn_I12).** Fixing the strand regression on its
own, with today's pick rule, costs BUSCOs: complete -22, missing +20 (8 of
10 samples), duplicated +10. The search is not the problem in itself.
blastn puts opposite-strand copies of a transcript into one orthogroup at
last, and `score` then picks among them badly. Every lost BUSCO traces to
that pick (`trace_arms.py`): a correct multi-assembler orthogroup, and a
different member kept. The worst case is SRR1123893, where a 214 bp
fragment was kept over the 3,429 bp contig carrying BUSCO 530014.

**Inflation and pick rule (blastn_I12 vs candidate).** Moving to `-I 3`
with `score_orf` recovers most of it: complete +18, missing -12, duplicated
-29. The two changes do different jobs. On SRR1139197, the one sample with
every combination: `-I 12` to `-I 3` under `score` cut duplicated 20 to 14
(missing 28 to 29), and `score` to `score_orf` at `-I 3` cut missing 29 to
26 (duplicated unchanged). **Inflation drives the duplicate reduction; the
pick rule drives completeness.**

**The combination (control vs candidate).** 6% fewer duplicated BUSCOs, -4
complete, +8 missing, more unique genes in all 10 samples, transrate up in
all 10. A net improvement, but it still loses some BUSCOs, through the same
pick mechanism.

**Gene loss (unique genes).** Both blastn arms raise the count in every
sample, but the sets churn underneath (`gene_sets.py`). There are about
100 genes lost and 135-140 gained per sample, almost identically in both
arms, so the churn comes from the search, not from inflation or the pick
rule. Most of it is weak hits; at bitscore >= 200 it is about 15 lost
against 24 gained per sample. A trivial change (re-picking diamond's
orthogroups) loses 0-4 strong genes on SRR1139197, so the effect is real
but small. On SRR1123893 the control contig behind each strongly hit "lost"
gene was searched for in the arm: candidate 16 near-identical, 11 partial,
0 absent; blastn_I12 18, 6 and 1. These "losses" are mostly bookkeeping.
diamond's `--top 0.1` keeps every hit within 10% of a contig's best, so a
slightly different copy of the same contig can drop secondary names out of
that window. **BUSCO, not the unique-gene count, is where genuine gene loss
shows, and it traces to the pick rule.**

**Cost.** blastn's peak memory is 2.5-5x diamond's. OrthoFinder needs `-a`
of at least the number of assemblies to survive its stall watchdog, and
even that is not enough for the largest sample. Its unused MSA/tree stages
are heavier with blastn's groups (FastTree was OOM-killed past 120 GB on
SRR1123893).

### Recommendation

1. **Don't ship blastn at `-I 12` with `score`** ("restore the pre-upgrade
   behaviour"). With today's pick rule it loses BUSCOs in 8 of 10 samples.
2. **Any strand fix needs a better pick rule.** The pick is where genes are
   lost, and `score_orf` helps without closing the gap. Next test: a rule
   that prefers protein evidence, such as the member with the longest
   swissprot alignment, which ORP already computes. Every lost BUSCO was a
   kept contig lacking the gene's protein match. `repick.py` plus
   `repick_arms.sbatch` can test it on the saved arms without new searches.
3. **Lower inflation is the lever for duplication,** once picks are safe.
4. **Consider orient-then-diamond** instead of blastn. It gives the same
   strand fix without blastn's memory, run time on large samples, or
   watchdog failure.

## The protein-evidence pick rule (`protein`)

`scripts/pick_best_contigs.py --rule protein --diamond <per-assembly blastx>`
keeps, in each orthogroup, the member with the strongest swissprot hit (best
bitscore in ORP's per-assembly diamond blastx output). Ties, and groups with
no hit, fall back to the transrate score, and the `score > 0` floor stays.
`protein_len` ranks by aligned length instead. In `oyster.py`,
`--pick-rule protein` runs the four per-assembly diamond passes before
`makeorthout` (they normally run just after) and passes them to the picker.
The default rule's step order and output are unchanged.

**Preview** (`preview_pick.py`, offline re-pick of the saved orthogroups): of
the BUSCOs each blastn arm lost against the control, how many have the
control's BUSCO-carrying contig picked under each rule:

| Arm | today's rule (score) | score_orf | protein | protein_len |
| --- | --- | --- | --- | --- |
| candidate (blastn `-I 3`) | 0/10 | 0/10 | **9/10** | 8/10 |
| blastn_I12 | 0/22 | 9/22 | **18/22** | 17/22 |

`protein` changes 3,000-17,000 picks per arm and picks slightly more
sequence. Its end-to-end test is `repick_validate.sbatch` with
`protein_tasks.tsv`: the rule applied to all three arms' orthogroups on all
10 samples (job 1319839). Each run stages everything up to the pick from
the saved arm (from home or its run.tar), re-decompresses the reads dated
to the ingest marker, writes the new pick with the real picker and lets
chowder run from `orthofusing`.

## Two-track redundancy removal (`twotrack_select.py`)

This replaces OrthoFinder and the per-orthogroup pick, on Matt's proposal of
2026-10-02, after the search/inflation/pick tests kept trading duplication
against gene loss. OrthoFinder's orthogroups do not line up with genes: some
split one gene (the duplicates) and some mix several (the gene loss at the
pick).

- **Track 1, contigs with a swissprot hit:** grouped by the gene of their
  best hit, the name qualreport's UNIQUE GENES counts. One contig is kept
  per gene: the one whose best hit covers most of the protein, ties to
  transrate score, then TPM. Every gene in the pool keeps a representative.
- **Track 2, contigs with no hit:** `cd-hit-est -c 0.95 -G 0 -aS 0.9 -r 1`
  (both strands, local identity, shorter contig at least 90% covered).
- **Expression rescue:** ORP's existing TPM filter (`--tpm-filt 1`) after
  cd-hit and salmon on the reduced assembly. Contigs above 1 TPM stay, and
  so do all contigs with a swissprot hit, so a no-hit contig survives only
  above 1 TPM.

On SRR954929 the 117,613 pooled contigs reduce to 36,813 (8,104 genes plus
28,709 no-hit) going forward, against 77,417 orthogroup picks in the
control. The no-hit cd-hit-est took 15 s. The test is
`repick_validate.sbatch` with `twotrack_tasks.tsv`: rule `twotrack` on each
sample's control files (job 1319913, 10 samples, after the protein reruns).

## Results across every trial, and the move to production (2026-10-02)

`redundancy_trials.py` collects every finished run into
`results/redundancy_trials.csv` (80 rows), with a summary in
`results/redundancy_trials_summary.txt`. Totals over the 9 samples that
have every trial (BUSCO in genes; transrate and mapping rates are medians):

| Trial | Contigs | C | D | M | Unique genes | Transrate | Good mappings |
| --- | --- | --- | --- | --- | --- | --- | --- |
| control (diamond `-I 12`, score) | 761,160 | 909 | 307 | 87 | 85,805 | 0.495 | 0.953 |
| blastn `-I 12`, score | 708,324 | 894 | 329 | 104 | 86,204 | 0.493 | 0.945 |
| candidate (blastn `-I 3`, score_orf) | 676,289 | 908 | 303 | 94 | 86,127 | 0.522 | 0.953 |
| diamond `-I 12` + protein | 751,205 | 928 | 325 | 78 | 85,128 | 0.489 | 0.954 |
| blastn `-I 12` + protein | 702,566 | 924 | 351 | 79 | 85,104 | 0.494 | 0.952 |
| blastn `-I 3` + protein | 673,014 | 924 | 320 | 79 | 85,220 | 0.515 | 0.950 |
| two-track | 629,957 | 918 | 36 | 97 | 86,603 | 0.425 | 0.837 |

Two-track is the only trial that removes duplication (307 to 36) while
holding gene content. Its cost is read mapping: good mappings fall about
11 points. `mapping_breakdown.py` splits that into about 7 points of
fragments that no longer align anywhere and about 3 of pairs that align but
break. `reads_lost.py` (megablast of control contigs against the two-track
assembly, weighted by the control's salmon NumReads) puts the reads on
other versions of kept genes, not on lost genes or no-hit contigs: 9-33% of
control reads sit on contigs two-track keeps only partly (shorter kept
version, so UTRs and ends are missing) and 7-26% on contigs it drops
entirely (other isoforms, same-name paralogs). No-hit contigs and lost
genes each carry under 1%. These are upper bounds; reads on the covered
part of a partial copy still map.

**In production as ORP 4.1.0-dev0** (`scripts/twotrack_select.py`,
`--merge-method twotrack`, the default), with two rescue steps added:
step 1, the representative is the *longest* contig with near-best protein
coverage (restores ends, no added contigs); step 2, up to two distinct
(under 50% covered by the representative), expressed (>= 1 TPM) copies per
gene come back. On SRR954929 that keeps 8,104 representatives plus 2,620
rescued copies. The no-hit contigs are untouched by the rescue: they are
kept after cd-hit-est and left to `--tpm-filt`, which keeps every one above
1 TPM. **`--tpm-filt` defaults to 0 (no filter); the test runs all pass
1.** The end-to-end validation of 4.1.0-dev0 through `chowder.py` on 5 of
the 10 samples is described in NOTES.md.

## SPAdes with auto k (`spades_branch.sbatch`, 2026-10-03)

The same two-track merge, but with the two rnaSPAdes runs changed: the
first leaves `-k` off (rnaSPAdes' auto pair, k=33,49 on these 100 bp reads),
the second uses 60% and 75% of read length (k=59,75) instead of the fixed
k=55 and k=75. `spades_branch.sbatch` re-runs both from each sample's
existing corrected reads, then runs `chowder.py` (default merge,
`--tpm-filt 1`) with the sample's existing Trans-ABySS and Trinity
assemblies. Results land in `~/redundancy_tests/spades_branch/<run>/`. The
numbers below are from the ten `qualreport.*_spadesbr` files (BUSCO counts
are the report's percentages x 125; no contig count or good-mapping rate
is in a qualreport, so proper pairs are shown instead).

Same 5 samples as `v410` (BUSCO genes out of 625):

| arm | C | S | D | F | M | unique genes | transrate | proper pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| control (4.0) | 487 | 360 | 127 | 88 | 50 | 49,670 | 0.445 | 0.941 |
| two-track, no rescue | 499 | 475 | 24 | 63 | 63 | 49,884 | 0.404 | 0.839 |
| v410 (two-track + rescue) | 503 | 428 | 75 | 69 | 53 | 49,852 | 0.455 | 0.924 |
| **spades branch** | 512 | 433 | 79 | 70 | 43 | 50,556 | 0.469 | 0.914 |

Nine samples (everything but DRR036858, which has no earlier arm), against
the control: complete 983 vs 956, duplicated 148 vs 324, missing 57 vs 64,
unique genes 89,501 vs 87,203, transrate 0.481 vs 0.471, proper pairs 0.933
vs 0.953. DRR036858 alone: 124/125 complete, 16 duplicated, 16,183 unique
genes, transrate 0.266.

What this does and does not show:

- Complete BUSCOs were equal or higher than v410's on all 5 shared samples,
  unique genes higher on 4 (SRR1139198 was the exception, 8,553 vs 8,696).
  SRR954929 lost proper pairs (0.924 to 0.872).
- Two things changed together (the auto first run and the 60%/75% second
  run), so the gain is not yet credited to either. In every qualreport the
  spadesauto assembly has more unique genes than spadeshigh.
- Samples are 100 bp reads only; one BUSCO gene is 0.8% of 125.
- Still to do: compare spadesauto with the old spades55 per sample (the
  `v410` run directories hold both), and run auto with the old k=75 to
  isolate the second run's change.

## Not yet done

- **Check isoform calls across assemblers.** For each isoform-like pair, see
  whether each structure has a near-identical match in at least two of the
  four raw assemblies (still kept in each run directory).
- **Trace copies through the pipeline.** Run one small sample with
  `--no-cleanup` and follow every duplicated-BUSCO contig through
  `Orthogroups.txt`, the posthack rescue (`newbies.fasta`) and cd-hit's
  `.clstr`.
- **Test candidate fixes on real data.** Re-deduplicate `ORP.fasta` (for
  example, orient and then re-cluster, or `cd-hit-est -G 0 -aS 0.9`) and
  rerun BUSCO on about 10 samples. Keep a fix only if duplicates fall and
  completeness holds.
