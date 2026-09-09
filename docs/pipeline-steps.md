# oyster.py Step Reference

What every `self.step()` in `Pipeline.main()` actually reads and writes, and what it does functionally. For execution order and where concurrency kicks in, see [pipeline-schedule.html](pipeline-schedule.html) — this document is the companion piece: same steps, but focused on inputs/outputs/purpose rather than scheduling, so a reader can tell what each stage is *for* and where the pipeline's CPU/time actually goes.

Reflects `oyster.py` as of commit [53e436f](https://github.com/macmanes-lab/Oyster_River_Protocol/commit/53e436feeaa277ae23468f0c5e2c1decf8202118) (2026-08-18) · ORP 3.1.0.

All paths below are relative to the run directory (`--dir`) and use `<run>` for `--runout`. "Env" is the conda environment the step's tool runs in.

## Read prep

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `run_trimmomatic` | `orp` trimmomatic (or a bare jar on `bridges.psc.edu`) | `--read1`/`--read2` | `rcorr/<run>.TRIM_{1,2}P.fastq` | Adapter/quality trims raw reads (`LEADING:3 TRAILING:3 ILLUMINACLIP MINLEN:25`), paired output only. |
| `run_rcorrector` | `orp` `run_rcorrector.pl` | the two `TRIM_*P.fastq` files | `rcorr/<run>.TRIM_{1,2}P.cor.fq` | K-mer-based read error correction (k=31). This corrected pair (`c1`/`c2`) is what every assembler and downstream alignment step reads from here on — trimmed-but-uncorrected reads are never touched again. |

## Assembly lanes

Two sequential stage-pairings, each a fixed `ThreadPoolExecutor(max_workers=2)`, independent of `--max-parallel` — see [pipeline-schedule.html](pipeline-schedule.html) for why. All assemblers consume `c1`/`c2`.

**Stage A** (`TRINITY_PHASE1_SHARE`, 25/75 cpu/mem split):

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `run_trinity_phase1` | `orp_trinity` Trinity, `--no_distributed_trinity_exec` | `c1`, `c2` | `assemblies/<run>.trinity/recursive_trinity.cmds.ok` | Inchworm (de Bruijn contig graph) + Chrysalis (read partitioning per gene component), then stops. Doesn't assemble anything yet — just prepares the per-component job list Phase 2 dispatches. Largely insensitive to its CPU share above Inchworm's own `--inchworm_cpu=10` cap. |
| `run_spades55` / `run_spades75` | `orp_spades` rnaspades.py, `--only-assembler` | `c1`, `c2` | `<run>.spades{55,75}.fasta` | rnaSPAdes at k=55 and k=75 (`--spades1-kmer`/`--spades2-kmer`). Sequential, slowest first. |
| `diamond_spades55` / `diamond_spades75` | `orp` diamond blastx | the matching spades fasta | `diamond/<run>.spades{55,75}.diamond.txt` | Blastx against swissprot; fires immediately after each assembly since it only needs its own fasta, not the merge stage below. |

**Stage B** (`TRINITY_PHASE2_SHARE`, 95/5 cpu split; Trans-ABySS's mem share is not cut the same way — see `oyster.py`):

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `run_trinity_phase2` | `orp_trinity` Trinity (no stop flag, `--full_cleanup`) | Phase 1's checkpoint files | `<run>.trinity.Trinity.fasta` | Resumes from Phase 1 straight into the actual per-gene-component assembly — thousands of small independent jobs dispatched via ParaFly, and Trinity's dominant cost by far (~34h on the SRR1789336 benchmark). Runs alongside Trans-ABySS rather than waiting for Stage A's short lane to fully clear, since Trans-ABySS's own dominant cost is single-threaded regardless of CPU count (see NOTES.md 2026-08-19). |
| `run_transabyss` | `orp_transabyss` transabyss | `c1`, `c2` | `assemblies/<run>.transabyss.fasta` | De novo assembly at `--transabyss-kmer` (default 32). Paired with Phase 2 instead of Stage A's SPAdes pair, since its dominant cost (initial FASTQ read + De Bruijn graph build) can't use extra cores anyway. |
| `diamond_transabyss` | `orp` diamond blastx | `<run>.transabyss.fasta` | `diamond/<run>.transabyss.diamond.txt` | Blastx against swissprot; fires immediately after the assembly. |

## Merging into one assembly (OrthoFuser)

`run_filtershort` fans out to both branches; `orthofuser_branch` and `merge_branch` run concurrently when `--max-parallel ≥ 2`, then join at `makeorthout`.

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `run_filtershort` | `orp` `scripts/long.seq.py` (×4, one per assembly) | the 4 raw assembly fastas | `orthofuse/<run>/working/<run>.<name>.short.fasta` ×4 | Drops contigs ≤200bp from each of the four assemblies. Short junk contigs would otherwise dominate the orthogroup clustering below. |
| `run_orthofuser` | `orp_orthofinder` orthofinder, `-d` (DNA mode), `-og` (orthogroups only) | the 4 `*.short.fasta` | `Orthogroups.txt` under the OrthoFinder results tree; `orthofuser.done` sentinel | Clusters contigs from all four assemblers into orthogroups — the pipeline's "these are probably the same transcript, assembled four different ways" grouping. |
| `makelist` (untimed) | pure Python | `Orthogroups.txt` | `orthofuse/<run>/<run>.list` | Just an index `1..n` (n = orthogroup count) — feeds the range `makegroups` iterates. |
| `makegroups` (untimed) | pure Python, thread pool | `Orthogroups.txt`, the index list | `orthofuse/<run>/<i>.groups` ×n; `groups.done` | Splits `Orthogroups.txt` into one file per orthogroup, each listing its member contig IDs. |
| `merge` (untimed) | pure Python (file concat) | the 4 `*.short.fasta` | `orthofuse/<run>/merged.fasta` | Concatenates the four short-filtered assemblies into one pool — no dedup yet, just union. |
| `orthotransrate` | `orp` orp-transrate | `merged.fasta`, `c1`, `c2` | `orthofuse/<run>/merged/contigs.csv` | Scores every contig in the pooled fasta for assembly quality (read-support based), by aligning `c1`/`c2` back to it. |
| `makeorthout` | `orp` `scripts/pick_best_contigs.py` | `contigs.csv`, the `*.groups` files | `orthofuse/<run>/good.<run>.list` | For each orthogroup, keeps the single highest-transrate-scoring member contig (score must be > 0); this is the actual "best isoform per gene, chosen across all four assemblers" selection. |
| `orthofusing` | `orp` `scripts/filter.py` | `merged.fasta`, `good.<run>.list` | `assemblies/<run>.orthomerged.fasta` | Filters the pooled fasta down to just the winning contigs — the first cut of the merged, deduplicated assembly. |

## Gene-uniqueness accounting

This is the least obvious part of the pipeline: a set-algebra pass that finds genes OrthoFuser's clustering *dropped* (because their contig lost the transrate vote, or wasn't grouped at all) and rescues them back into the assembly.

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `diamond_orthomerged` | `orp` diamond blastx | `<run>.orthomerged.fasta` | `diamond/<run>.orthomerged.diamond.txt` | Blastx of the merged assembly. |
| `diamond_trinity` | `orp` diamond blastx | `<run>.trinity.Trinity.fasta` | `diamond/<run>.trinity.diamond.txt` | Blastx of the raw (un-filtered, un-merged) Trinity assembly — `diamond_{transabyss,spades55,spades75}` already ran earlier, in Stage A/Stage B above. |
| `diamond_uniq` (untimed) | pure Python | all 5 diamond outputs | `diamond/<run>.unique.{trinity,sp55,sp75,transabyss}.txt` | Counts distinct swissprot gene IDs hit by each individual assembler — reporting metrics only, doesn't gate anything downstream. |
| `make_list1` (untimed) | pure Python | `diamond_orthomerged` | `diamond/<run>.list1` | Gene IDs hit by the *merged* assembly. |
| `make_list2` (untimed) | pure Python | the 4 individual-assembler diamond outputs | `diamond/<run>.list2` | Union of gene IDs hit by *any* of the four raw assemblies. |
| `make_list3` (untimed) | pure Python | `list1`, `list2` | `diamond/<run>.list3` | `list2 − list1`: genes some individual assembler found that the merged assembly does **not** represent — i.e., genes OrthoFuser's selection accidentally dropped. |
| `make_list5` | `orp` `scripts/build_list5.py` | `list3`, the 4 individual diamond outputs | `diamond/<run>.list5` | For each dropped gene in `list3`, picks a rescue contig ID — the first individual-assembly contig (checked transabyss → spades75 → spades55 → trinity) that hit it. |
| `make_list6` (untimed) | pure Python | `<run>.orthomerged.fasta` | `diamond/<run>.list6` | Every sequence ID currently in the merged assembly. |
| `make_list7` (untimed) | pure Python | `list5`, `list6` | `diamond/<run>.list7` | `list5 − list6`: rescue contig IDs not already present in the merged assembly (belt-and-suspenders — should already be disjoint, but confirms it). |
| `posthack` | `orp` `scripts/filter.py` via a `bash -c` process substitution | the 4 **raw** (not short-filtered) assembly fastas, `list7` | `diamond/<run>.newbies.fasta`; `assemblies/working/<run>.orthomerged.fasta` | Pulls the `list7` rescue contigs back out of the original, un-filtered assemblies (not the ≤200bp-trimmed ones from `run_filtershort` — a rescued contig may be short) as `newbies.fasta`, then appends them onto `orthomerged.fasta` to produce the true working assembly used from here on. |

## Dedup & quantify

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `cdhit` | `orp` cd-hit-est | the rescued working assembly | `assemblies/<run>.ORP.intermediate.fasta` | Collapses near-duplicate contigs at 98% identity (`-c .98`) — the last dedup pass. |
| `orp_diamond` | `orp` diamond blastx | `ORP.intermediate.fasta` | `assemblies/<run>.ORP.diamond.txt` | Blastx of the (nearly) final assembly — used both for the unique-gene count below and for the low-TPM rescue logic in `secondfilter`. |
| `orp_uniq` (untimed) | pure Python | `ORP.diamond.txt` | `assemblies/working/<run>.unique.ORP.txt` | Counts distinct genes hit — the headline "unique genes (ORP)" metric in the final report. |
| `salmon_index` | `orp` salmon | `ORP.intermediate.fasta` | `quants/<run>.ortho.idx` | Builds a salmon index (k=31) over the intermediate assembly. |
| `salmon` | `orp` salmon quant | the index, `c1`, `c2` | `quants/salmon_orthomerged_<run>/quant.sf` | Quantifies expression (TPM) per contig by pseudo-aligning the corrected reads. |
| `filter` (untimed) | pure Python | `quant.sf` | `assemblies/working/<run>.{HIGH,LOW}EXP.txt` | Splits contigs into above/below `--tpm-filt` TPM lists. |
| `secondfilter` | `orp` `scripts/filter.py` (×2) + pure Python | `ORP.intermediate.fasta`, `LOWEXP.txt`, `HIGHEXP.txt`, `ORP.diamond.txt` | `assemblies/<run>.ORP.fasta`; a `*_BEFORE_TPM_FILT.fasta` backup copy | If any contigs fell below the TPM threshold: keeps all high-TPM contigs outright, but rescues a low-TPM contig anyway if it's the *only* one with a diamond hit to its gene (`donotremove.list`) — so a real-but-lowly-expressed transcript with no redundant coverage isn't thrown away just for being quiet. If nothing was below threshold, `ORP.intermediate.fasta` is simply copied through unchanged. This is the file every later step (BUSCO, transrate, strandeval, `reportgen`) treats as "the assembly." |

## QC / report

`transrate` and `strandeval` run concurrently when `--max-parallel ≥ 2`; `busco` runs alone just before them (see [pipeline-schedule.html](pipeline-schedule.html) for why it isn't folded into that pair).

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `busco` | `orp_busco` busco, `--offline`, `-m transcriptome` | `ORP.fasta` | `reports/run_<run>.ORP/` | Scores completeness against the `--lineage` ortholog set (default `eukaryota_odb12.2`). |
| `transrate` | `orp` orp-transrate | `ORP.fasta`, `c1`, `c2` | `reports/transrate_<run>/assemblies.csv` | Same read-support quality scoring as `orthotransrate` earlier, now on the final assembly rather than the mid-pipeline pool. |
| `strandeval` | `orp_trinity` bwa + `orp` samtools + `scripts/examine_strand.pl` | `ORP.fasta`, a 400k-read subsample of `c1`/`c2` | `reports/<run>.strandeval_summary.txt` | Aligns a read subsample back to the assembly and checks read-orientation-vs-transcript-strand agreement — a sanity check on whether `--strand` was set correctly. |
| `reportgen` (untimed) | pure Python | BUSCO/transrate/diamond/salmon/strandeval outputs above | `reports/qualreport.<run>` | Pulls one headline number from each prior report into a single human-readable summary (BUSCO score, transrate scores, unique-gene counts per assembler, proper-pair mapping rate, strand histogram). |

## Observations: where the time goes and where to look for further gains

Notes from reading the pipeline end to end with input/output in hand — some are candidate efficiencies, others are just useful context for anyone tuning this further.

- **`run_filtershort` is four sequential `conda run` subprocesses for a trivial length filter.** Each of the four `scripts/long.seq.py` calls pays full `conda run -n orp ...` activation overhead for what's otherwise a single-pass Biopython length check, and none of the four depend on each other. They're not part of any lane or `run_parallel` group, so they always run one after another regardless of `--max-parallel`. Candidate fix: either fold the four filters into one script invocation (one `conda run`, loop over the four fastas in Python), or give them the same `run_parallel`/thread-pool treatment as `orthofuser_branch`/`merge_branch`.
- **DIAMOND blastx runs six times against the same swissprot database over the course of a run** — once each for `transabyss`, `spades55`, `spades75`, `trinity`, `orthomerged`, and finally `ORP.intermediate` — and is the largest CPU cost after the assemblers themselves. This is inherent to how unique-gene accounting works (each assembly's sequence IDs are distinct, so hits can't be shared across runs), and each call already gets the full `--cpu` budget since none of these six are inside a `run_parallel` group — so there's little scheduling headroom left here. If there's a future win it's more likely algorithmic (e.g. only re-searching the *new* sequences introduced by cd-hit/`posthack` rather than the whole intermediate fasta again) than a concurrency change.
- **`secondfilter` writes a `*_BEFORE_TPM_FILT.fasta` snapshot that nothing downstream reads.** It looks like a manual-inspection safety copy rather than dead code, but worth confirming that's the intent — it's a full copy of the intermediate assembly written on every run that has any low-TPM contigs.
- **The bookkeeping chain (`make_list1`–`make_list7`, `diamond_uniq`, `orp_uniq`, `filter`, `makelist`, `makegroups`, `merge`) is already about as cheap as it can get** — pure single-pass Python, no subprocess overhead, sub-second on real data (see the timing tables in [benchmarks.md](../sampledata/benchmarks.md)). Not a place to spend further optimization effort.
- **`diamond_orthomerged`/`diamond_trinity` run sequentially** even though, by the time they're reached, both of their dependencies (the merge chain and the Trinity lane) are already satisfied and they don't depend on each other. This matches the pipeline's stated policy that CPU-bound DIAMOND/BUSCO/salmon stages don't benefit from splitting cores to overlap (see `orp_diamond`'s comment in `oyster.py` and the same claim in the README), so it's likely intentional rather than an oversight — but unlike `orp_diamond`, there's no comment on this pair confirming that reasoning locally.
