# oyster.py Step Reference

What every `self.step()` in `Pipeline.main()` actually reads and writes, and what it does functionally. (`main()` is three methods -- `prepare_reads`, `run_assemblers`, `merge_and_report` -- matching the three sections below.) For execution order and where concurrency kicks in, see [pipeline-schedule.html](pipeline-schedule.html) — this document is the companion piece: same steps, but focused on inputs/outputs/purpose rather than scheduling, so a reader can tell what each stage is *for* and where the pipeline's CPU/time actually goes.

Reflects `oyster.py` as of ORP 4.1.0-dev9 (two-track merge only; the OrthoFinder merge was removed). Every step below is timed and appears in `reports/<run>.timing.log`: the bookkeeping steps used to be exempt, which was fair when each was seconds against an assembler's tens of hours, but a `chowder.py` merge has no assemblers and the merge half is the whole run.

All paths below are relative to the run directory (`--dir`) and use `<run>` for `--runout`. "Env" is the conda environment the step's tool runs in.

A step is skipped when its outputs exist and are newer than its inputs. Each step also touches `reports/.<run>.<step>.running` before it starts and removes it on success, so a step that was killed or failed partway is re-run on the next invocation even if it left output files behind.

`chowder.py` runs this same reference from `run_filtershort` down, over assemblies it was handed rather than ones it built: it skips the whole Assembly lanes section and puts one `ingest` step in front of it (pure Python -- copies each input to `assemblies/ingested/<run>.<label>.fasta`, prefixing every contig name with `<label>_`, and declares those copies as its outputs so a re-invocation resumes; the inputs themselves are never written, and one inside `assemblies/ingested/` is refused). Everywhere a row below says "the 4 assemblies", read "the N assemblies given to `--assemblies`". See the README for what else differs.

## Read prep

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `run_trimmomatic` | `orp` trimmomatic | `--read1`/`--read2` | `rcorr/<run>.TRIM_{1,2}P.fastq` | Adapter/quality trims raw reads (`LEADING:3 TRAILING:3 ILLUMINACLIP MINLEN:25`), paired output only. Declares `rcorr/<run>.trim.done` as its output instead of the TRIM files themselves once the corrected pair below is present and current, since by then `reclaim_trimmed_reads` has deleted them. |
| `run_rcorrector` | `orp` `run_rcorrector.pl` | the two `TRIM_*P.fastq` files | `rcorr/<run>.TRIM_{1,2}P.cor.fq` | K-mer-based read error correction (k=31). This corrected pair (`c1`/`c2`) is what every assembler and downstream alignment step reads from here on — trimmed-but-uncorrected reads are never touched again. |
| `reclaim_trimmed_reads` (not a step) | pure Python | — | deletes all four `rcorr/<run>.TRIM_*.fastq`; touches `rcorr/<run>.trim.done` | Runs inline right after `run_rcorrector`, not at the end: those four files are the size of the raw input and, per the row above, nothing ever opens them again. The corrected pair is then queued for background gzipping (`compress_async`), which runs through the whole assembly phase below. Skipped under `--no-cleanup`. |
| `use_corrected_reads` (not a step) | pure Python, or `pigz`/`gzip -d` | `--read1`/`--read2` | `rcorr/<run>.TRIM_{1,2}P.cor.fq` | Replaces the three rows above under `--trimmed-corrected-reads` (`--corrected-reads` in `chowder.py`): the user's pair stands in for rcorrector's output. A plain fastq is symlinked; a gzipped one is decompressed, since Trinity and SPAdes trust the `.fq` name. Neither is queued for compression, and `cleanup` keeps the symlinks and deletes the decompressed copies. |

## Assembly lanes

Two lanes side by side for the whole assembly stage, a fixed `ThreadPoolExecutor(max_workers=2)` independent of `--max-parallel` — see [pipeline-schedule.html](pipeline-schedule.html) for why. The Trans-ABySS lane gets `TRANSABYSS_SHARE` (25%) of `--cpu`, or under `--transabyss-mpi` a share fitted to the uncompressed corrected R1's size (0.2–0.45, `transabyss_share()`; a fixed 33% without `--normalize-reads`), and `TRANSABYSS_MEM_SHARE` (25%) of `--mem`, from the start; the Trinity lane gets the rest and runs Stage A, then Stage B (all of `--cpu` if Trans-ABySS is done by then). All assemblers consume `c1`/`c2`.

**Stage A**, Trinity lane (`TRINITY_PHASE1_SHARE`, cores 67/33; `TRINITY_PHASE1_MEM_SHARE`, mem 25/75):

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `run_trinity_phase1` | `orp_trinity` Trinity, `--no_distributed_trinity_exec` | `c1`, `c2` | `assemblies/<run>.trinity.phase1.done` (a sentinel deliberately *outside* the working directory Phase 2's `--full_cleanup` deletes — see the observation at the bottom) | Inchworm (de Bruijn contig graph) + Chrysalis (read partitioning per gene component), then stops. Doesn't assemble anything yet — just prepares the per-component job list Phase 2 dispatches. Largely insensitive to its CPU share above Inchworm's own `--inchworm_cpu=10` cap. |
| `run_spadesauto` / `run_spadeshigh` | `orp_spades` rnaspades.py, `--only-assembler` | `c1`, `c2` | `<run>.spades{auto,high}.fasta` | Two rnaSPAdes runs covering different k bands (`--spades1-kmer`/`--spades2-kmer`). `spadesauto` omits `-k` so rnaSPAdes picks its documented default pair (~1/3 and ~1/2 of maximum read length). `spadeshigh` sits above that band at 60%/75% of maximum read length, resolved per dataset. Sequential, slowest first. Each deletes its `<run>.spades_k*/` working directory when done, unless `--no-cleanup`. |
| `diamond_spadesauto` / `diamond_spadeshigh` | `orp` diamond blastx | the matching spades fasta | `diamond/<run>.spades{auto,high}.diamond.txt` | Blastx against swissprot; fires immediately after each assembly since it only needs its own fasta, not the merge stage below. |

**Stage B**, Trinity lane (the whole lane), and the **Trans-ABySS lane** (from the start, beside both stages):

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `run_trinity_phase2` | `orp_trinity` Trinity (no stop flag, `--full_cleanup` unless `--no-cleanup`) | Phase 1's sentinel (it still resumes from the on-disk checkpoints inside `<run>.trinity/`; the sentinel is only what `needs_run()` compares) | `<run>.trinity.Trinity.fasta` | Resumes from Phase 1 straight into the actual per-gene-component assembly — thousands of small independent jobs dispatched via ParaFly, and Trinity's dominant cost by far (~2h on 38 cores on the SRR1789336 benchmark with an OpenMP ParaFly; ~34h without). Gets the whole Trinity lane once Stage A is done. |
| `run_transabyss` | `orp_transabyss` transabyss | `c1`, `c2` | `assemblies/<run>.transabyss.fasta` | De novo assembly at `--transabyss-kmer` (default 32). Starts at once in its own lane: its dominant cost (initial FASTQ read + De Bruijn graph build) is single-threaded and can't use extra cores, which makes it the longest assembler. `--transabyss-mpi auto|on` runs that stage under MPI instead (`--mpi <its cores>`, i.e. `mpirun -np` `ABYSS-P`) and writes `reports/<run>.transabyss.mode`. Deletes its `<run>.transabyss/` working directory when done, unless `--no-cleanup`. |
| `diamond_transabyss` | `orp` diamond blastx | `<run>.transabyss.fasta` | `diamond/<run>.transabyss.diamond.txt` | Blastx against swissprot; fires immediately after the assembly. |

## Merging into one assembly (two-track)

`build_pool` → `score_pool` (as `pool_branch`) and any per-assembly `diamond_<label>` pass not already done (as `assembly_diamonds`) are independent until `twotrack_select`, so at `--max-parallel ≥ 2` they run together: `score_pool` keeps all of `--cpu`, and the diamond passes run beside it on a quarter of it, with up to 16 GB of `--mem` taken from pyTransRate's budget (`run_beside`).

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `run_filtershort` | `orp` `scripts/long.seq.py` (×4, one per assembly, in parallel) | the 4 raw assembly fastas | `shuck/<run>/working/<run>.<name>.short.fasta` ×4 | Drops contigs ≤200bp from each of the four assemblies. |
| `build_pool` | pure Python (file concat) | the 4 `*.short.fasta` | `shuck/<run>/pool.fasta` | Concatenates the four short-filtered assemblies into one pool — no dedup yet, just union. |
| `score_pool` | `orp` pyTransRate | `pool.fasta`, `c1`, `c2` | `shuck/<run>/pool/contigs.csv` | Scores every contig in the pooled fasta for assembly quality (read-support based), by aligning `c1`/`c2` back to it. |
| `diamond_<label>` | `orp` diamond blastx | each assembly | `diamond/<run>.<label>.diamond.txt` | Any per-assembly swissprot pass not already run in an assembler lane (normally only Trinity's under oyster.py; all of them under chowder.py). |
| `twotrack_select` | `orp` `scripts/twotrack_select.py` (blastn, cd-hit-est inside) | `pool.fasta`, `pool/contigs.csv`, the 4 per-assembly diamond outputs, `software/diamond/uniprot_sprot.fasta` | `shuck/<run>/good.<run>.list`; `shuck/<run>/twotrack.<run>.tsv` (the fate of every pooled contig) | Contigs with a Swiss-Prot hit are grouped by best-hit gene; each gene keeps the longest contig with near-best protein coverage, plus up to two distinct (under 50% covered by it at >=95% identity), expressed (>=1 TPM) copies. Contigs without a hit go through cd-hit-est (`-c 0.95 -G 0 -aS 0.9 -r 1`). |
| `shuck` | `orp` `scripts/filter.py` | `pool.fasta`, `good.<run>.list` | `assemblies/<run>.shucked.fasta` | Filters the pooled fasta down to just the winning contigs — the first cut of the merged, deduplicated assembly. |

## Gene-uniqueness accounting

This is the least obvious part of the pipeline: a set-algebra pass that finds genes the pick *dropped* (because their contig lost the pyTransRate vote, or wasn't grouped at all) and rescues them back into the assembly.

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `diamond_shucked` | pure Python (`hits_for`) | `<run>.shucked.fasta`, the per-assembly diamond outputs | `diamond/<run>.shucked.diamond.txt` | The merged assembly's swissprot hits, looked up from the per-assembly blastx outputs rather than searched for again: every contig in it came whole out of one of the assemblies. |
| `diamond_trinity` | `orp` diamond blastx | `<run>.trinity.Trinity.fasta` | `diamond/<run>.trinity.diamond.txt` | Blastx of the raw (un-filtered, un-merged) Trinity assembly — `diamond_{transabyss,spadesauto,spadeshigh}` already ran earlier, in the assembler lanes above. |
| `diamond_uniq` | pure Python | all 5 diamond outputs | `diamond/<run>.unique.{trinity,spauto,sphigh,transabyss}.txt` | Counts distinct swissprot gene IDs hit by each individual assembler — reporting metrics only, doesn't gate anything downstream. |
| `make_list1` | pure Python | `diamond_shucked` | `diamond/<run>.list1` | Gene IDs hit by the *merged* assembly. |
| `make_list2` | pure Python | the 4 individual-assembler diamond outputs | `diamond/<run>.list2` | Union of gene IDs hit by *any* of the four raw assemblies. |
| `make_list3` | pure Python | `list1`, `list2` | `diamond/<run>.list3` | `list2 − list1`: genes some individual assembler found that the merged assembly does **not** represent — i.e., genes the pick accidentally dropped. |
| `make_list5` | `orp` `scripts/build_list5.py` | `list3`, the 4 individual diamond outputs | `diamond/<run>.list5` | For each dropped gene in `list3`, picks a rescue contig ID — the first individual-assembly contig (checked transabyss → spadeshigh → spadesauto → trinity) that hit it. |
| `make_list6` | pure Python | `<run>.shucked.fasta` | `diamond/<run>.list6` | Every sequence ID currently in the merged assembly. |
| `make_list7` | pure Python | `list5`, `list6` | `diamond/<run>.list7` | `list5 − list6`: rescue contig IDs not already present in the merged assembly (belt-and-suspenders — should already be disjoint, but confirms it). |
| `posthack` | `orp` `scripts/filter.py` via a `bash -c` process substitution | the 4 **raw** (not short-filtered) assembly fastas, `list7` | `diamond/<run>.newbies.fasta`; `assemblies/working/<run>.shucked.fasta` | Pulls the `list7` rescue contigs back out of the original, un-filtered assemblies (not the ≤200bp-trimmed ones from `run_filtershort` — a rescued contig may be short) as `newbies.fasta`, then appends them onto `shucked.fasta` to produce the true working assembly used from here on. |

## Dedup & quantify

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `cdhit` | `orp` cd-hit-est | the rescued working assembly | `assemblies/<run>.ORP.intermediate.fasta` (plus cd-hit-est's own `.clstr` beside it) | Collapses near-duplicate contigs at 98% identity (`-c .98`) — the last dedup pass. |
| `orp_diamond` | pure Python (`hits_for`) | `ORP.intermediate.fasta`, the per-assembly diamond outputs | `assemblies/<run>.ORP.diamond.txt` | The (nearly) final assembly's swissprot hits, looked up the same way as `diamond_shucked` — used both for the unique-gene count below and for the low-TPM rescue logic in `secondfilter`. |
| `orp_uniq` | pure Python | `ORP.diamond.txt` | `assemblies/working/<run>.unique.ORP.txt` | Counts distinct genes hit — the headline "unique genes (ORP)" metric in the final report. |
| `salmon_index` | `orp` salmon | `ORP.intermediate.fasta` | `quants/<run>.shucked.idx` | Builds a salmon index (k=31) over the intermediate assembly. |
| `salmon` | `orp` salmon quant | the index, `c1`, `c2` | `quants/salmon_shucked_<run>/quant.sf` | Quantifies expression (TPM) per contig by pseudo-aligning the corrected reads. |
| `filter` | pure Python | `quant.sf` | `assemblies/working/<run>.{HIGH,LOW}EXP.txt` | Splits contigs into at-or-above / below `--tpm-filt` TPM lists (a contig exactly at the threshold counts as high). |
| `secondfilter` | `orp` `scripts/filter.py` (×2) + pure Python | `ORP.intermediate.fasta`, `LOWEXP.txt`, `HIGHEXP.txt`, `ORP.diamond.txt` | `assemblies/<run>.ORP.fasta`; a `*_BEFORE_TPM_FILT.fasta` backup copy | If any contigs fell below the TPM threshold: keeps all high-TPM contigs outright, but rescues a low-TPM contig anyway if it's the *only* one with a diamond hit to its gene (`donotremove.list`) — so a real-but-lowly-expressed transcript with no redundant coverage isn't thrown away just for being quiet. If nothing was below threshold, `ORP.intermediate.fasta` is simply copied through unchanged. This is the file every later step (BUSCO, pyTransRate, strandeval, `reportgen`) treats as "the assembly." |

## QC / report

At `--max-parallel ≥ 2`, `strandeval` runs beside `pytransrate` on a quarter of `--cpu` (at most 8 threads) while `pytransrate` keeps all of it (`run_beside`); `busco` runs alone just before them with the whole of `--cpu`/`--busco-threads`.

| Step | Env / tool | Inputs | Outputs | What it does |
|---|---|---|---|---|
| `busco` | `orp_busco` busco, `--offline`, `-m transcriptome` | `ORP.fasta` | `reports/run_<run>.ORP/` | Scores completeness against the `--lineage` ortholog set (default `eukaryota_odb12.2`). A re-run replaces the previous report. |
| `pytransrate` | `orp` pyTransRate | `ORP.fasta`, `c1`, `c2` | `reports/pytransrate_<run>/assemblies.csv` (a `reports/transrate_<run>/` left by an older run is renamed to this) | Same read-support quality scoring as `score_pool` earlier, now on the final assembly rather than the mid-pipeline pool. |
| `strandeval` | `orp_trinity` bwa + `orp` samtools + `scripts/examine_strand.pl` | `ORP.fasta`, a 400k-read subsample of `c1`/`c2` | `reports/<run>.strandeval_summary.txt` | Aligns a read subsample back to the assembly and checks read-orientation-vs-transcript-strand agreement — a sanity check on whether `--strand` was set correctly. Deletes its sorted BAM and bwa index when done, unless `--no-cleanup`. |
| `reportgen` | pure Python | BUSCO/pyTransRate/diamond/salmon/strandeval outputs above | `reports/qualreport.<run>` | Pulls one headline number from each prior report into a single human-readable summary (BUSCO score, pyTransRate scores, unique-gene counts per assembler, proper-pair mapping rate, strand histogram). |
| `cleanup` | pure Python | `qualreport.<run>` | `reports/<run>.cleanup.done` (a manifest of what was kept and removed, with sizes) | Last, because it deletes files earlier steps declare as inputs. Keeps `reports/`, `.ORP.fasta`, and the four individual assemblies plus the corrected reads as the `.gz` that `compress_async` has been building in the background since each was written — so this only unlinks, and the compression cost was already paid in parallel with an assembler. Removes `shuck/`, `quants/`, `assemblies/diamond/`, `assemblies/working/`, and the working assemblies between `shuck` and `.ORP.fasta` (including cd-hit-est's `.ORP.intermediate.fasta.clstr`); all of it is reproducible from what's kept, and every number it fed is already in `qualreport.<run>`. Also sweeps up what the steps above normally delete themselves — the rnaSPAdes and Trans-ABySS working directories, Trinity's gene_trans_map, strandeval's BAM and bwa index — which are only still there after a `--no-cleanup` or interrupted run. A file whose background gzip didn't finish is kept uncompressed instead of deleted. Skipped under `--no-cleanup`, which also leaves `cleanup.done` unwritten, so re-running without the flag cleans up then. |

## Observations: where the time goes and where to look for further gains

Notes from reading the pipeline end to end with input/output in hand — some are candidate efficiencies, others are just useful context for anyone tuning this further.

- **DIAMOND blastx runs four times, once per assembly**, and is the largest CPU cost after the assemblers themselves. `shucked` and `ORP.intermediate` used to get searches of their own; since every contig in them came whole out of one of the assemblies, their hits are now looked up from the four per-assembly outputs (`hits_for`).
- **`secondfilter` writes a `*_BEFORE_TPM_FILT.fasta` snapshot that nothing downstream reads.** It looks like a manual-inspection safety copy rather than dead code, but worth confirming that's the intent — it's a full copy of the intermediate assembly written on every run that has any low-TPM contigs. It lands in `assemblies/working/`, so `cleanup` removes it at the end of the run along with the rest of that directory.
- **Trinity erasing its own evidence that it ran was worth ~35h on a resume.** `run_trinity_phase2` passes `--full_cleanup`, which deletes the whole `<run>.trinity/` working directory — including `recursive_trinity.cmds.ok`, which used to be both Phase 1's declared output and Phase 2's declared input. So on a resumed run `needs_run()` found Phase 1's output missing and re-ran it (~90min), which rewrote `cmds.ok` *newer* than `<run>.trinity.Trinity.fasta`, which made Phase 2 look out of date and re-run against an assembly that was already complete (~34h on the SRR1789336 benchmark). Both phases now hang off `assemblies/<run>.trinity.phase1.done` instead, which lives outside the directory `--full_cleanup` removes. The window was narrow but common: Stage B succeeded and the run failed or was killed later — a walltime kill near the end of a long run, which is exactly when a job gets resubmitted. `seed_trinity_phase1_sentinel()` back-fills the sentinel for run directories created before it existed, copying the mtime of whatever already proves Phase 1 ran (`cmds.ok`, else the finished assembly) so the ordering `needs_run()` compares is preserved rather than stamped `now`.
- **`cleanup` deleting most steps' outputs is what `already_complete()` exists for.** `needs_run()` decides each step purely from whether its outputs are present and current, so on a cleaned run directory nearly every stage reads as out of date. Re-invoking `oyster.py` there — a resubmitted cluster job, say — would therefore reassemble from scratch over a finished run, where before cleanup existed it no-opped. `already_complete()` short-circuits `main()` on `reports/<run>.cleanup.done` plus an `.ORP.fasta` still newer than the raw reads, before `timing_init()` gets a chance to truncate the finished run's timing log.
- **The bookkeeping chain (`make_list1`–`make_list7`, `diamond_uniq`, `orp_uniq`, `filter`, `build_pool`) is about as cheap as it can get** — pure single-pass Python, no subprocess overhead, sub-second on real data (see the timing tables in [benchmarks.md](../sampledata/benchmarks.md)). Not a place to spend further optimization effort.
