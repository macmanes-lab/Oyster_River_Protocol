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

### Parallel task management

See [docs/pipeline-schedule.html](docs/pipeline-schedule.html) for a full DAG of execution order and concurrency (download and open locally, or view via [htmlpreview](https://htmlpreview.github.io/?https://github.com/macmanes-lab/Oyster_River_Protocol/blob/master/docs/pipeline-schedule.html)), and [docs/pipeline-steps.md](docs/pipeline-steps.md) for what each step actually reads and writes.

Trinity itself runs in two stages, using its documented [multi-stage execution support](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity#running-trinity-in-multiple-sequential-stages), and each stage is paired with a different assembler rather than sharing one fixed lane for the whole run:

- **Stage A**: Trinity Phase 1 (Inchworm + Chrysalis, building the whole-transcriptome graph and partitioning reads per gene component) runs alongside rnaSPAdes55 and rnaSPAdes75, splitting `--cpu`/`--mem` 25/75 (`TRINITY_PHASE1_SHARE`) -- Phase 1's own share barely matters (Inchworm is capped at a fixed thread count regardless, and Chrysalis's clustering is brief), so SPAdes gets the bulk of the machine since it's fast but does scale with cores. Each SPAdes assembly is immediately followed by its own diamond search rather than waiting for the orthofuser/merge stage below, since that search only ever needed its own assembly.
- **Stage B**: once Stage A finishes, Trinity Phase 2 -- the actual per-gene-component assembly, thousands of small independent jobs and by far Trinity's dominant cost -- runs alongside Trans-ABySS, splitting `--cpu`/`--mem` 95/5 (`TRINITY_PHASE2_SHARE`). Trans-ABySS's own dominant cost (the initial FASTQ read + De Bruijn graph build) is single-threaded regardless of CPU count, so it gets just enough cores to keep its own threaded sub-stages moving while Phase 2 takes the rest; its memory share isn't cut along with its CPU share, since its memory footprint doesn't shrink the same way. Trans-ABySS's diamond search runs immediately after it finishes, same as the Stage A assemblers.

This split is fixed and not affected by `--max-parallel`.

By default (`--max-parallel 2`), oyster.py runs up to 2 jobs at once within the other stages of the pipeline that benefit from it, splitting `--cpu`/`--mem` across however many jobs are running concurrently:

- the orthofuser branch vs. the merge/orthotransrate branch
- transrate vs. strandeval

CPU-bound stages that don't benefit from splitting cores — diamond, orp_diamond, salmon, and BUSCO — always run sequentially at the full `--cpu` count regardless of this flag.

Set `--max-parallel 1` to disable concurrency for those stages and run them one at a time (useful when debugging, or on a machine where you'd rather not split cores). Raise it above 2 to run more jobs at once within a stage, at the cost of each job getting a smaller slice of `--cpu`/`--mem`.

### All flags

| Flag | Default | Description |
|------|---------|-------------|
| `--read1` | *(required)* | Path to R1 fastq(.gz) |
| `--read2` | *(required)* | Path to R2 fastq(.gz) |
| `--mem` | `110` | Memory in GB |
| `--cpu` | `16` | CPU threads |
| `--busco-threads` | same as `--cpu` | BUSCO threads |
| `--runout` | `USER_RUN` | Run name prefix |
| `--strand` | unstranded (`""`) | Strand-specificity: `RF`, `FR`, or unset |
| `--lineage` | `eukaryota_odb12.2` | BUSCO lineage |
| `--normalize-reads` | off | Let Trinity normalize reads (default is `--no_normalize_reads`) |
| `--tpm-filt` | `0` | TPM filter threshold |
| `--spades1-kmer` | `55` | rnaSPAdes k-mer for the spades55 assembly |
| `--spades2-kmer` | `75` | rnaSPAdes k-mer for the spades75 assembly |
| `--transabyss-kmer` | `32` | Trans-ABySS k-mer |
| `--max-parallel` | `2` | Max concurrent jobs per stage (see [Parallel task management](#parallel-task-management) above) |
| `--dir` | current directory | Working directory |
| `--version` | — | Print the installed ORP version and exit |
| `--help` | — | Print this same flag reference and exit |
