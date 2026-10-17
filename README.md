# Oyster River Protocol

Official Repository of the Oyster River Protocol for Transcriptome Assembly

## Installation

See [INSTALL.md](INSTALL.md) for step-by-step install directions, covering both the one-command `make` installer and the manual, step-by-step alternative.

## Usage

`--read1` and `--read2` are the only required flags; everything else has a default.

```bash
python3 oyster.py --read1 R1.fq.gz --read2 R2.fq.gz --mem 110 --cpu 24 --runout runname --strand RF
```

`python3 oyster.py --help` prints the same reference from the command line, and `python3 oyster.py --version` prints the installed ORP version.

### Trinity read normalization

By default, oyster.py runs Trinity with `--no_normalize_reads`, i.e. read normalization is disabled. To let Trinity normalize reads instead, pass `--normalize-reads` on the command line:

```bash
python3 oyster.py --read1 R1.fq.gz --read2 R2.fq.gz --mem 110 --cpu 24 --runout runname --strand RF --normalize-reads
```

### Reads you have already corrected

If your pair has already been through trimmomatic and rcorrector, pass `--trimmed-corrected-reads` and oyster.py assembles it as it is, skipping both steps (and not requiring either tool at preflight):

```bash
python3 oyster.py --read1 R1.cor.fq.gz --read2 R2.cor.fq.gz --mem 110 --cpu 24 --runout runname --trimmed-corrected-reads
```

A plain fastq is symlinked into `rcorr/` as the corrected pair. A gzipped one is decompressed there instead, because Trinity and SPAdes decide whether to gunzip from the file extension; that copy is deleted at cleanup rather than re-gzipped, since your original is still where you left it.

### Merging assemblies you already have (`chowder.py`)

`chowder.py` runs the second half of the ORP -- everything after the four assemblers -- over assemblies you bring yourself, so two or more existing transcriptomes can be merged into one the same way ORP merges its own:

```bash
python3 chowder.py --assemblies best.fasta other.fasta third.fasta.gz \
        --read1 R1.fq.gz --read2 R2.fq.gz --mem 110 --cpu 24 --runout runname
```

It is the same code rather than a copy of it -- `chowder.py` subclasses `oyster.py`'s pipeline and reuses the merge stages wholesale -- so the two-track selection rule, the contig ordering that reaches cd-hit-est's tie-breaks, and the pyTransRate scoring are identical to a full ORP run's by construction. The output is a `<run>.ORP.fasta` and a `reports/qualreport.<run>` in the usual layout, with one `UNIQUE GENES` line per input assembly.

Three things differ from `oyster.py`, and all three are worth knowing before you run it:

- **The reads are not optional.** The merge scores every contig against them (pyTransRate), quantifies the survivors (salmon) and strand-checks the result, so it needs the library the assemblies were built from. It trims and error-corrects that library the way ORP always does; pass `--corrected-reads` if your pair has already been through trimmomatic and rcorrector, and it will use it as-is.
- **Contig names are prefixed with the label of the assembly they came from** (`trinity_TRINITY_DN0_c0_g1_i1`). Two assemblies of one library routinely share contig names -- two Trinity runs both start at `TRINITY_DN0_c0_g1_i1` -- and every stage downstream joins on that name, so without the prefix two different contigs would silently be treated as one. It doubles as provenance: every contig in the final assembly says which input it survived from. Labels come from the filenames unless you pass `--labels`.
- **Assembly order matters, and by default it is not yours to get wrong.** It sets contig order in the pooled fasta, which reaches cd-hit-est, where input order breaks length ties between near-identical contigs; and it is the order the diamond rescue searches, so for genes the selection missed, earlier assemblies are preferred. Rather than let the sequence you happened to type decide either of those, chowder sorts the assemblies by label and permutes them with a fixed seed, so the order depends on the *set* of assemblies and not on how you listed them. It is a seeded shuffle rather than a random one on purpose: randomising outright would mean the same command gave a different assembly on a different day, and a resumed run disagreeing with the run it resumed. The order used and the seed that produced it are printed when the run starts and written to `assemblies/<run>.ingest.done`.

  This makes the order arbitrary and reproducible; it does not make the pipeline order-*independent*, which would mean breaking cd-hit-est's ties and the rescue ranking on merit rather than on position. To find out what the order is worth on your own data, run it twice with different `--seed` values and diff the assemblies. To rank the assemblies yourself, best first, pass `--assembly-order given`.

`python3 chowder.py --help` prints the full flag reference, and [All flags (`chowder.py`)](#all-flags-chowderpy) below is the same list. There are no assembler flags -- no k-mers, no `--strand`, no `--normalize-reads` -- and preflight does not require SPAdes, Trinity or Trans-ABySS, since a merge never runs them. It does check for `bwa`, which `strandeval` runs.

### Scoring an assembly you already have (`assembly_eval.py`)

`assembly_eval.py` scores one transcriptome assembly the way an ORP run scores its own, without assembling or merging anything. Use it to put an assembly from another pipeline, a published transcriptome, or an older ORP run on the same scale as a new one:

```bash
python3 assembly_eval.py --assembly my.fasta --read1 R1.fq.gz --read2 R2.fq.gz \
        --cpu 24 --mem 110 --runout myrun
```

It runs the three evaluations an ORP run ends with:

- **pyTransRate**, the contig-level read-mapping score (`PYTRANSRATE SCORE` and `PYTRANSRATE OPTIMAL SCORE`)
- **BUSCO** completeness against `--lineage` (default `eukaryota_odb12.2`)
- **strandeval**, the strand-specificity histogram, plus the proper-pair rate from the same mapping

These are the same commands, conda envs and seeds `oyster.py` uses, and the result is a `reports/qualreport.<run>` in the same layout, so a report from each diffs cleanly. The only lines missing are `UNIQUE GENES`, because those compare the assemblers inside one ORP run and a single assembly has nothing to compare against.

Some things to know before you run it:

- **The reads are used as given.** Nothing is trimmed or error-corrected, so pass the pair you want the assembly scored against, normally the one it was built from. The assembly must be uncompressed fasta.
- **The three steps run at the same time**, with `--cpu` and `--mem` split between them: pyTransRate gets half, and BUSCO and strandeval a quarter each. pyTransRate is the slow step, and BUSCO gains little past a few threads.
- **It prints as it goes.** The run opens with a banner, the command line, and the inputs and settings. pyTransRate and BUSCO output then streams to the screen, each line marked `[pytransrate]` or `[busco]` because the two run together. Every tool's full output is also kept in `logs/<step>.log`.
- **Interrupted runs resume.** Each step that finishes leaves `reports/<run>.<step>.done`, and a re-run skips steps that have one. `--force` redoes all of them, and `--skip pytransrate busco strandeval` leaves out any you don't need. If a step fails, the report is still written, with that step's lines blank.

Output goes to `--outdir` (default `./<runout>_eval`):

| Path | What it holds |
|------|---------------|
| `reports/qualreport.<run>` | the quality report |
| `reports/pytransrate_<run>/` | pyTransRate's output |
| `reports/run_<run>.ORP/` | BUSCO's output |
| `reports/<run>.strandeval_summary.txt` | the strand histogram |
| `reports/<run>.flagstat` | `samtools flagstat` of the sampled mapping |
| `logs/<step>.log` | each tool's full output |

[All flags (`assembly_eval.py`)](#all-flags-assembly_evalpy) below lists every option, and `python3 assembly_eval.py --help` prints the same list.

### How the assemblies are merged

From 4.1.0 the four assemblies are merged by Swiss-Prot gene, not by OrthoFinder orthogroup:

- **Contigs with a Swiss-Prot hit** are grouped by the gene of their best hit -- the same gene names the qualreport's UNIQUE GENES counts. Each gene keeps one representative (the longest of the contigs whose hit covers close to the best fraction of the protein), plus up to two more copies that are both distinct from it (it covers under half of them) and expressed (at least 1 TPM): other isoforms, paralogs under the same name.
- **Contigs without a hit** are deduplicated with cd-hit-est on both strands, and the expression filter (`--tpm-filt`, `--min-cov`) later drops those below threshold.

Everything after that -- the diamond rescue, cd-hit-est, salmon, BUSCO, pyTransRate -- works as before; the expression filter is now a coverage filter, described below. On 10 test samples this cut duplicated BUSCOs by about 90% against ORP 4.0 while keeping complete BUSCOs and raising unique genes; the experiments behind it are in `experiments/redundancy/`. The OrthoFinder merge ORP used through 4.0 was removed in 4.1, so a 4.1 assembly is not directly comparable with a 4.0 one.

### Parallel task management

See [docs/pipeline-schedule.html](docs/pipeline-schedule.html) for a full DAG of execution order and concurrency (download and open locally, or view via [htmlpreview](https://htmlpreview.github.io/?https://github.com/macmanes-lab/Oyster_River_Protocol/blob/master/docs/pipeline-schedule.html)), and [docs/pipeline-steps.md](docs/pipeline-steps.md) for what each step actually reads and writes.

The assemblers run as two lanes side by side, and Trinity itself runs in two stages within its lane, using its documented [multi-stage execution support](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity#running-trinity-in-multiple-sequential-stages):

- **Trans-ABySS lane**: Trans-ABySS starts at once and keeps its cores for the whole assembly stage, with a quarter of `--mem` (`TRANSABYSS_MEM_SHARE`). How many cores depends on `--transabyss-mpi`. Threaded (the default), its dominant cost (the unitig assembly) runs on one thread whatever it is given, so it gets a quarter of `--cpu` (`TRANSABYSS_SHARE`): more did not make it faster, and it is the longest assembler either way. Under MPI that stage uses its cores and both lanes scale, but no up-front split keeps them level: Trans-ABySS's work follows the reads, Trinity's the transcriptome. So Trans-ABySS gets half of `--cpu` (`TRANSABYSS_MPI_SHARE`), enough that it has not been the long pole on any dataset so far, and whatever it doesn't need comes back to Trinity when it finishes (Stage B below). The log prints the share and the corrected R1's size.
- **Trinity lane**, on the rest of the machine:
  - **Stage A**: Trinity Phase 1 (Inchworm + Chrysalis, building the whole-transcriptome graph and partitioning reads per gene component) runs alongside rnaSPAdes (the `spadesauto` and `spadeshigh` runs), splitting the lane's cores two thirds to Phase 1 (`TRINITY_PHASE1_SHARE`) and its memory 25/75 (`TRINITY_PHASE1_MEM_SHARE`) -- Phase 1 has been the slower of the two on every dataset, while SPAdes treats `--memory` as a hard cap.
  - **Stage B**: once Stage A finishes, Trinity Phase 2 -- the actual per-gene-component assembly, thousands of small independent jobs -- gets the whole lane, or all of `--cpu` if Trans-ABySS has already finished. If Trans-ABySS finishes while Phase 2 is running and at least 15% of its components are left (`TRINITY_GROW_MIN_LEFT`), Phase 2 is stopped and relaunched on all of `--cpu`; ParaFly skips the components already done, so only the few in flight are redone. The log says `run_trinity_phase2 grows to N cpu` (or `stays at`) with the count left.

Each assembly is immediately followed by its own diamond search rather than waiting for the merge stage below, since that search only ever needed its own assembly. On a resumed run whose Trans-ABySS is already done, the Trinity lane gets the whole machine.

This split is fixed and not affected by `--max-parallel`.

Two later pairs of steps are independent of each other, and each pairs a short job with a long one. By default (`--max-parallel 2`) the short job runs beside the long one on a few threads of its own (a quarter of `--cpu`, at most 8 for strandeval), while the long job keeps all of `--cpu`:

- the diamond passes not already run in an assembler lane (Trinity's under oyster.py, every assembly's under chowder.py) beside `score_pool`
- strandeval beside pyTransRate

The short job's memory (a quarter of `--mem`, at most 16 GB) comes out of the long job's budget, so memory is never overcommitted. This used to split the cores evenly, which left the short job's half of the machine idle once it finished while the long job carried on at half speed.

Everything else (`orp_diamond`, salmon, BUSCO and the rest) runs on its own with the full `--cpu`. `--max-parallel 1` runs each pair one after the other, every step with the whole machine.

### What a finished run leaves behind

A completed run keeps five things and reclaims the rest:

| Kept | |
|---|---|
| `assemblies/<run>.ORP.fasta` | the assembly — the point of the run, left uncompressed |
| `assemblies/<run>.{spadesauto,spadeshigh,transabyss,trinity.Trinity}.fasta.gz` | the four individual assemblies, gzipped |
| `rcorr/<run>.TRIM_{1,2}P.cor.fq.gz` | the trimmed **and error-corrected** reads, gzipped — the pair every assembler actually read |
| `reports/` | BUSCO, pyTransRate, strand evaluation, `qualreport.<run>`, timings |
| `reports/<run>.cleanup.done` | what was reclaimed and what was kept, with sizes |

Everything else goes: the trimmed-but-uncorrected reads (deleted as soon as read correction is done with them — nothing downstream ever reads them again), the `shuck/` tree (the pooled fasta and its pyTransRate scoring, normally the largest directory in a run), `quants/`, `assemblies/diamond/`, `assemblies/working/`, and the chain of working assemblies between `shuck` and `.ORP.fasta`. Every number any of those contributed is already in `reports/qualreport.<run>`.

The gzipping runs in the background, starting the moment each file is finished being written rather than at the end of the run — the corrected reads compress alongside the assemblers, and each assembly compresses while the next stage runs — so cleanup itself is just an unlink and adds no wall time. Pass `--no-cleanup` to switch all of this off and keep a run exactly as it was, which is what you want when debugging a run rather than shipping its results. It also keeps what individual steps would otherwise delete as they finish: Trinity's working directory (Phase 2 runs without `--full_cleanup`), the rnaSPAdes and Trans-ABySS working directories, and strandeval's BAM and bwa index. Running the same command again later without the flag does the cleanup then. `--keep-intermediates` is an older name for the same flag and still works.

Re-running `oyster.py` on a directory whose run already finished and was cleaned up is a no-op: it reports where the assembly and reports are and exits, rather than treating the reclaimed intermediates as work to redo. To assemble the same reads again, use a different `--runout`/`--dir`, or delete `reports/<run>.cleanup.done` to force a full re-run in place.

Re-running on an unfinished run resumes it: a step whose outputs exist and are newer than its inputs is skipped. A step that was killed or failed partway is re-run even if it left output files behind, since tools like rcorrector and diamond write as they go and a half-written file would otherwise look finished. Each step leaves `reports/.<run>.<step>.running` while it runs and removes it when it succeeds. Timestamps alone can't tell that a setting changed, or that a read file was swapped for a different, older one, so the settings each assembler and chowder's `ingest` depend on (`--strand`, the k-mers, `--normalize-reads`, `--seed`/`--assembly-order`) and the size and mtime of the reads and of chowder's assemblies are recorded in `assemblies/working/<run>.*.settings` and `rcorr/<run>.reads.source`. Change one and the steps built on it are redone; moving the files without changing them redoes nothing.

### The expression filter (`--min-cov`, `--tpm-filt`)

After salmon quantifies the merged assembly, contigs **without a Swiss-Prot hit** are dropped if they look like noise; a contig with a hit is never dropped by this step. Two thresholds decide it, and a contig has to pass both:

- **`--min-cov`** (default `auto`): coverage, salmon's NumReads x mean fragment length / contig length -- the average number of read-pair fragments over each base of the contig. `auto` picks the threshold from sequencing depth (below); a number is used as given.
- **`--tpm-filt`** (default `0`, off): TPM.

Before 4.1.0 nothing was filtered unless you gave `--tpm-filt` (1 was the usual value); now coverage is filtered by default. Coverage rather than TPM because TPM is a share of a million, so a fixed TPM cutoff means something different in every library: across the 119 final assemblies of the 108-sample test set, 1 TPM was anywhere from 4 to 183 reads per kb, mostly following sequencing depth. On a deep library `--tpm-filt 1` removes far more than intended (on a 451M-read tuco library, 1 TPM is 41x coverage, and it removed 91% of the contigs without a hit); on a shallow one it removes almost nothing. TPM also divides by salmon's effective length, which lets a short contig pass on a couple of reads. Coverage means the same thing in every library.

**`--min-cov auto` scales the threshold with sequencing depth.** The same share of noise fragments becomes more coverage the deeper you sequence: on the 451M-read tuco library a fixed `--min-cov 2` removed only 10% of the no-hit contigs. So `auto` (the default) counts the read pairs salmon quantified -- the trimmed and corrected reads, never Trinity's normalized ones -- and applies

    min-cov = max(1, 2 x (read pairs / 60M) ^ 0.85)

| read pairs | `--min-cov auto` resolves to |
|---|---|
| 30M or fewer | `1` |
| 60M | `2` |
| 120M | `3.6` |
| 240M | `6.5` |
| 480M | `11.7` |

The rule comes from re-quantifying two finished assemblies (tuco, 480M pairs; *Anopheles* ERR489297, 200M) on nested subsamples of their own reads, 15M pairs up. At every depth, the contigs without a hit were binned by coverage and compared with their own shuffled sequence: noise has no excess of long open reading frames, real coding sequence does. The coverage where that excess climbs out of the noise rose with depth at an exponent of 0.8-0.9 in both libraries, and the two sat within a factor of about 1.5 of each other at the same depth. The value 2 at 60M pairs is where the earlier five-library sweep (16-123M pairs) found `--min-cov 2` removed fewer no-hit contigs than `--tpm-filt 1` while keeping BUSCO and unique genes unchanged; the floor of 1 is there because at 15M pairs a threshold of 2 already cuts into the real signal. The log says what `auto` chose (`[filter] --min-cov auto: N read pairs -> --min-cov X`), and so does the `# min_cov` line of `reports/<run>.filter_candidates.tsv`. Two libraries is still a small basis: for an unusual library, or many tissues pooled, check the choice against `filter_candidates.tsv`.

To check or override the choice: every run writes `reports/<run>.filter_candidates.tsv`, which lists how many contigs, and how many without a hit, each of a range of coverage, TPM and read-count thresholds would keep, from the same quantification, and its header gives the library's `cov_at_1tpm` (how much coverage 1 TPM is here). If the applied threshold looks too loose or too strict, re-run the same command with a new `--min-cov`: the threshold is recorded in `assemblies/working/<run>.filter.settings`, so a resumed run re-filters and redoes only the steps from the filter on (BUSCO, pyTransRate, the report), not the assemblers or the merge. That needs the run's intermediates, which cleanup removes, so on a deep library you expect to tune, run with `--no-cleanup` and then re-run without it once you've settled on a value. Pass a number such as `--min-cov 2` for a fixed threshold, `--min-cov 0 --tpm-filt 1` for the old TPM behaviour, or `--min-cov 0` to switch filtering off.

The experiments behind this are in [experiments/dynamic_cov/](experiments/dynamic_cov/README.md).

### Tuning pyTransRate

Both pyTransRate runs -- the one that scores the pooled fasta for the two-track selection, and the one that scores the finished assembly for the report -- take a fixed argument list. `--pytransrate-args` appends to both, as one quoted string, `shlex`-split so the arguments arrive separately:

```bash
python3 chowder.py --assemblies a.fasta b.fasta --read1 R1.fq.gz --read2 R2.fq.gz \
        --pytransrate-args '--location-size 5'
```

Because it appends, a value given here overrides the same flag ORP passes above it. `pytransrate --help` (inside the `orp` env) lists the full set. Two are worth calling out for a large merge:

- **`--location-size`** is the one that matters. snap encodes genome locations in a fixed number of bytes and sweeps upward until they fit; passing `5` skips the sweep when you already know four bytes will not hold the merge. The sweep is not cheap to get wrong -- a failed four-byte attempt is a full index build thrown away.
- **`--padding` is not the memory lever it looks like.** snap writes padding between contigs as `N` and skips any seed containing `N`, so it grows the 1 byte/base genome array and nothing else. Measured on a 5.4M-contig, 5.7 Gbp merge: dropping it entirely saved 1% of the branch, and real sequence alone still exceeded the four-byte ceiling.

`oyster.py` takes the same flag. It is on both entry points because they carry separate parsers, but a merge is where an assembly large enough to need it normally shows up.

### All flags (`oyster.py`)

| Flag | Default | Description |
|------|---------|-------------|
| `--read1` | *(required)* | Path to R1 fastq(.gz) |
| `--read2` | *(required)* | Path to R2 fastq(.gz) |
| `--mem` | `110` | Memory in GB |
| `--cpu` | `16` | CPU threads |
| `--busco-threads` | same as `--cpu` | BUSCO threads |
| `--runout` | `USER_RUN` | Run name prefix: letters, digits, `_` and `-` only |
| `--strand` | unstranded (`""`) | Strand-specificity: `RF`, `FR`, or unset |
| `--lineage` | `eukaryota_odb12.2` | BUSCO lineage |
| `--normalize-reads` | off | Let Trinity normalize reads (default is `--no_normalize_reads`) |
| `--tpm-filt` | `0` | TPM filter threshold |
| `--min-cov` | `auto` | coverage filter threshold (NumReads x mean fragment length / contig length); a contig must also pass `--tpm-filt`. `auto` scales it with read pairs, `max(1, 2 x (pairs / 60M) ^ 0.85)`; a number is used as given -- see [The expression filter](#the-expression-filter---min-cov---tpm-filt) |
| `--trimmed-corrected-reads` | off | Treat `--read1`/`--read2` as already trimmed and error-corrected, and skip trimmomatic and rcorrector (see [Reads you have already corrected](#reads-you-have-already-corrected) above) |
| `--spades1-kmer` | `auto` | rnaSPAdes k-mer(s) for the spadesauto assembly — `auto` lets rnaSPAdes pick its documented default pair (~1/3 and ~1/2 of maximum read length). Also accepts percentages or an explicit list, same forms as `--spades2-kmer` |
| `--spades2-kmer` | `60%,75%` | rnaSPAdes k-mer(s) for the spadeshigh assembly — percentages of maximum read length, an explicit comma-separated list of odd sizes under 128, or `auto`. Percentages resolve per dataset (61,75 at 101bp reads; 89,113 at 150bp), clamped to rnaSPAdes' 127 ceiling |
| `--transabyss-kmer` | `32` | Trans-ABySS k-mer |
| `--transabyss-mpi` | `off` | `auto`/`on` run Trans-ABySS's single-threaded unitig assembly under MPI (`mpirun -np <its cores> ABYSS-P`). Preflight starts two ranks first: `on` refuses to run if that fails, `auto` falls back to threaded. `reports/<run>.transabyss.mode` records which ran |
| `--max-parallel` | `2` | `1` runs the short jobs that would otherwise run beside a long one (see [Parallel task management](#parallel-task-management) above) one after the other |
| `--no-cleanup` | off | Keep every file a run produces, uncompressed, including each assembler's working directory -- for debugging (see [What a finished run leaves behind](#what-a-finished-run-leaves-behind) below). Alias: `--keep-intermediates` |
| `--pytransrate-args` | none | Extra arguments passed verbatim to both pyTransRate runs, as one quoted string (see [Tuning pyTransRate](#tuning-pytransrate) below) |
| `--dir` | current directory | Working directory |
| `--version` | — | Print the installed ORP version and exit |
| `--help` | — | Print this same flag reference and exit |

### All flags (`chowder.py`)

There are no assembler flags — no k-mers, no `--strand`, no `--normalize-reads` — since a merge never runs an assembler. Everything below `--corrected-reads` means exactly what it does under `oyster.py`.

| Flag | Default | Description |
|------|---------|-------------|
| `--assemblies` | *(required)* | Two or more assembly fasta(.gz) files to merge |
| `--read1` | *(required)* | Path to R1 fastq(.gz) — the library the assemblies were built from |
| `--read2` | *(required)* | Path to R2 fastq(.gz) |
| `--labels` | from the filenames | One label per assembly, used as the contig-name prefix and in the report |
| `--assembly-order` | `shuffled` | `shuffled` orders the assemblies by `--seed`, independently of how you typed them; `given` uses your order, best first (see [Merging assemblies you already have](#merging-assemblies-you-already-have-chowderpy) above) |
| `--seed` | `23894` | Seed for `--assembly-order shuffled`; change it to measure what the order is worth on your data |
| `--corrected-reads` | off | Treat `--read1`/`--read2` as already trimmed and error-corrected, and skip both steps |
| `--mem` | `110` | Memory in GB |
| `--cpu` | `16` | CPU threads |
| `--busco-threads` | same as `--cpu` | BUSCO threads |
| `--runout` | `USER_RUN` | Run name prefix: letters, digits, `_` and `-` only |
| `--lineage` | `eukaryota_odb12.2` | BUSCO lineage |
| `--tpm-filt` | `0` | TPM filter threshold |
| `--min-cov` | `auto` | coverage filter threshold (NumReads x mean fragment length / contig length); a contig must also pass `--tpm-filt`. `auto` scales it with read pairs, `max(1, 2 x (pairs / 60M) ^ 0.85)`; a number is used as given -- see [The expression filter](#the-expression-filter---min-cov---tpm-filt) |
| `--max-parallel` | `2` | `1` runs the short jobs that would otherwise run beside a long one (see [Parallel task management](#parallel-task-management) above) one after the other |
| `--no-cleanup` | off | Keep every file the run produces, uncompressed -- for debugging. Alias: `--keep-intermediates` |
| `--pytransrate-args` | none | Extra arguments passed verbatim to both pyTransRate runs (see [Tuning pyTransRate](#tuning-pytransrate) above) |
| `--dir` | current directory | Working directory |
| `--version` | — | Print the installed ORP version and exit |
| `--help` | — | Print this same flag reference and exit |

### All flags (`assembly_eval.py`)

| Flag | Default | Description |
|------|---------|-------------|
| `--assembly` | *(required)* | Assembly fasta to score; must not be gzipped |
| `--read1` | *(required)* | Path to R1 fastq(.gz), used as given |
| `--read2` | *(required)* | Path to R2 fastq(.gz) |
| `--runout` | assembly filename up to the first `.` | Run name used in output file names |
| `--outdir` | `./<runout>_eval` | Output directory |
| `--cpu` | all cores | Total threads, split across the steps |
| `--mem` | `32` | Total memory in GB, split across the steps |
| `--lineage` | `eukaryota_odb12.2` | BUSCO lineage |
| `--skip` | none | Steps not to run: any of `pytransrate`, `busco`, `strandeval` |
| `--force` | off | Redo steps already marked done |
| `--no-cleanup` | off | Keep strandeval's BAM and bwa index |
| `--pytransrate-args` | none | Extra pyTransRate arguments as one quoted string (see [Tuning pyTransRate](#tuning-pytransrate)) |
| `--busco-config` | `software/config.ini` if present | BUSCO `config.ini` |
| `--busco-download-path` | `busco_dbs/` beside the script if present | BUSCO lineage dataset directory |
| `--examine-strand` | `scripts/examine_strand.pl` | Path to `examine_strand.pl` |
| `--no-conda` | off | Run tools from `PATH` instead of the conda envs; preflight checks they are all there |
| `--env` | `orp` | Conda env with pyTransRate and samtools |
| `--env-busco` | `orp_busco` | Conda env with BUSCO |
| `--env-trinity` | `orp_trinity` | Conda env with Trinity, bwa, seqtk and hist |
| `--help` | — | Print this same flag reference and exit |
