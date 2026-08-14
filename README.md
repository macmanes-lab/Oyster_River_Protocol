# Oyster River Protocol

[![Join the chat at https://gitter.im/macmanes-lab/Oyster_River_Protocol](https://badges.gitter.im/macmanes-lab/Oyster_River_Protocol.svg)](https://gitter.im/macmanes-lab/Oyster_River_Protocol?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Official Repository of the Oyster River Protocol for Transcriptome Assembly

Please see http://oyster-river-protocol.readthedocs.io/en/latest/ and https://hackmd.io/s/SJhOQvkVm# for instructions about how to run and install. 

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

By default (`--max-parallel 2`), oyster.py runs up to 2 jobs at once within each stage of the pipeline that benefits from it, splitting `--cpu`/`--mem` across however many jobs are running concurrently:

- the 4 assemblers (Trinity, rnaSPAdes55, rnaSPAdes75, Trans-ABySS)
- the orthofuser branch vs. the merge/orthotransrate branch
- transrate vs. strandeval

CPU-bound stages that don't benefit from splitting cores — diamond, orp_diamond, salmon, and BUSCO — always run sequentially at the full `--cpu` count regardless of this flag.

Set `--max-parallel 1` to disable concurrency entirely and run every step one at a time (useful when debugging, or on a machine where you'd rather not split cores). Raise it above 2 to run more jobs at once within a stage, at the cost of each job getting a smaller slice of `--cpu`/`--mem`.

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
