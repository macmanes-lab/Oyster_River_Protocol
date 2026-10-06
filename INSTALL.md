# Installing the Oyster River Protocol

These directions cover installing ORP 3.0+ (the `oyster.py` pipeline) on Linux. Two paths are covered: the one-command installer, and the same steps done by hand if you want more control or need to troubleshoot a partial install.

## Prerequisites

- Linux (the installer's Anaconda bootstrap is Linux x86_64 specific)
- `git`, `curl`, `bash`, and a system `python3` already available. **`oyster.py`/`chowder.py` run under that interpreter, not under a conda env**, so it sets the language floor: 3.6, which is what a cluster's `/usr/bin/python3` still commonly is. Everything they orchestrate runs inside the conda environments and is unaffected by it.
- Internet access on the install machine — it pulls down Anaconda, the UniProt/Swiss-Prot diamond database, and the BUSCO lineage database
- A working internet connection to GitHub as well: Trinity is built from a git clone (with submodules)
- Several GB of free disk space (Swiss-Prot + Anaconda + 6 conda environments + BUSCO database adds up)

## Option A: one-command install

```bash
git clone https://github.com/macmanes-lab/Oyster_River_Protocol.git
cd Oyster_River_Protocol
make
```

`make` chains through everything: creates the conda environments, builds Trinity from source, builds the diamond search database, downloads the BUSCO lineage database, and appends any needed PATH entries to `~/.profile`/`~/.bash_profile`. Each step checks whether it's already done and skips itself if so, so re-running `make` after a partial or failed install picks up where it left off rather than starting over.

When it finishes, run:

```bash
source ~/.profile
```

## Option B: manual step-by-step install

Useful if you want to see exactly what's happening, or if `make` fails partway through and you want to finish a specific step by hand. Run these from inside a cloned `Oyster_River_Protocol` directory, with `conda`/`mamba` already available on your system.

**1. Point mamba/conda at the right channels**

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda install mamba -n base -yc conda-forge
```

**2. Create the 4 isolated tool environments**

Only `orp_spades` pins `python=3.14` — that env has a confirmed bug (see [Known gotchas](#known-gotchas)) that requires forcing a modern Python. The others resolve whatever compatible Python their own recipe naturally wants; forcing 3.14 there bought no benefit (each env's Python is invisible outside that env) and broke `orp_transabyss`'s solve in practice.

```bash
mamba create -y -c bioconda -c conda-forge --override-channels --name orp_spades spades=4.3.0 python=3.14
mamba create -y -c bioconda -c conda-forge --override-channels --name orp_trinity compilers "cmake<4" make git kmer-jellyfish bowtie2 samtools salmon=1.10.3 "openjdk>=17" perl perl-db_file python numpy bwa=0.7.19 bashplotlib seqtk=1.5 libgomp zlib
mamba create -y -c bioconda -c conda-forge --override-channels --name orp_busco busco=6.1.0
mamba create -y -c bioconda -c conda-forge --override-channels --name orp_transabyss transabyss=2.0.1
```

`orp_trinity` holds Trinity's build tools and runtime dependencies but not Trinity itself, which is built from source in step 3. The bioconda `trinity=2.15.2` package (builds `_4` to `_6`, May 2025 on) ships a ParaFly compiled without OpenMP, so Trinity's phase 2 runs one component at a time whatever `--cpu` says: 34h27m instead of 2h03m on SRR1789336. Don't `mamba install trinity` into this env. `cmake` is pinned below 4 because Inchworm and Chrysalis declare a `cmake_minimum_required` that CMake 4 refuses.

**3. Build Trinity from source**

```bash
conda run --no-capture-output -n orp_trinity scripts/build_trinity.sh
```

This clones [trinityrnaseq/trinityrnaseq](https://github.com/trinityrnaseq/trinityrnaseq) at the commit pinned in the script (the `Trinity-v2.15.2` release, unmodified), builds it into `software/trinityrnaseq` with the env's compiler, and links `Trinity` into the env's `bin/`. It fails if ParaFly came out without OpenMP. You can check by hand, the count must be above 0:

```bash
nm -D software/trinityrnaseq/trinity-plugins/BIN/ParaFly | grep -c GOMP_
```

`oyster.py` repeats a timing version of that check at startup, so a serial ParaFly stops a run before it starts.

**4. Create the consolidated `orp` environment**

Everything else — rcorrector, trimmomatic, cd-hit, diamond, salmon (the pipeline's own, modern version), samtools, seqtk, sra-tools, blast, parallel, biopython, scipy, numpy, bashplotlib, pigz, plus pyTransRate and the snap-aligner it maps with — lives in one `orp` environment, defined in `orp_env.yml`:

```bash
mamba env create -f orp_env.yml
```

**5. Diamond's Swiss-Prot search database**

```bash
mkdir -p software/diamond
cd software/diamond
curl -LO ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gzip -d uniprot_sprot.fasta.gz
conda run -n orp diamond makedb --in uniprot_sprot.fasta -d swissprot
cd ../..
```

**6. BUSCO's eukaryota lineage database**

```bash
mkdir -p busco_dbs
conda run -n orp_busco busco --download eukaryota_odb12.2 --download_path busco_dbs
```

## Updating an existing `orp` environment

If `orp_env.yml` has changed (e.g. after pulling a newer version of this repo) and `mamba env create -f orp_env.yml` fails with `CondaValueError: prefix already exists`, update the environment in place instead of recreating it:

```bash
mamba env update -f orp_env.yml --prune
```

`--prune` removes anything currently in the env that's no longer listed in the yml, so it ends up matching the file exactly. Drop `--prune` if you've added extra packages by hand that you want to keep.

To start that environment fully fresh instead:

```bash
mamba env remove -n orp
mamba env create -f orp_env.yml
```

## Verifying the install

```bash
python3 oyster.py --version
```

Then run the pipeline once on real (or the bundled sample) data — `oyster.py`'s own preflight `check()` step verifies that every tool it runs is reachable in its env, and that the Swiss-Prot diamond database, the Swiss-Prot fasta and the `--lineage` BUSCO dataset are installed. If anything is missing it lists all of it and stops, rather than failing partway through a multi-hour run:

```bash
python3 oyster.py --read1 R1.fq.gz --read2 R2.fq.gz --runout myrun --cpu 24 --mem 110 --strand RF
```

`chowder.py` — the merge-only entry point, for assemblies you already have — shares that preflight, less the three assemblers it never runs. It still needs the `orp_trinity` env, for strandeval's `bwa`, `seqtk` and `hist`:

```bash
python3 chowder.py --assemblies a.fasta b.fasta --read1 R1.fq.gz --read2 R2.fq.gz --runout mymerge --cpu 24 --mem 110
```

## Uninstalling

```bash
make clean
```

Removes the conda install and downloaded software directories.

## Known gotchas

- **"Python version 3.6.8 is not supported"** from SPAdes: a [known SPAdes bug](https://github.com/ablab/spades/issues/1319) where it detects a stray system Python instead of its own environment's interpreter. Fixed by `orp_spades`'s explicit `python=3.14` pin above — if you hit this anyway, double check nothing earlier in your `PATH` (an HPC module system, for example) is shadowing the conda environment's own Python.
- **`nothing provides _python_rc needed by python-3.14.0rc1-...`**: means mamba is resolving to a Python 3.14 *release candidate* build because no stable 3.14 build exists yet for that particular package's dependency set. This is why only `orp_spades` pins `python=3.14` — pinning it on older/less-actively-maintained packages (like `transabyss=2.0.1`) can hit exactly this wall.
- **snap-aligner dies with SIGFPE (`Floating point exception`) partway through `score_pool` on a large merge.** This is [amplab/snap#171](https://github.com/amplab/snap/issues/171), an upstream integer divide-by-zero on snap's CIGAR-writing path, not a memory problem and not something ORP's flags can steer around: the crash needs a read with a secondary alignment that back-clips it to exactly half its length, and the output buffer to run out partway through writing it, at which point the flush-and-retry retains the back clipping and clips the original alignment to zero bases. The scale variable is therefore the size of the BAM rather than the size of the assembly — a merge producing ~160 GB of output refills that buffer more or less continuously, so a per-read-rare event becomes a certainty, while the same pipeline and flags against a smaller assembly are clean. The fix is two functional lines, sitting unreleased on snap's `dev` branch (`0e0997b`, `2.0.6.dev.2`), so bioconda's latest — the `snap-aligner=2.0.5` in `orp_env.yml` — still has the bug. Until upstream tags 2.0.6, the workaround is to build v2.0.5 with that commit's `SNAPLib/ReadWriter.cpp` and put the result ahead of the env's copy on `PATH`: pyTransRate resolves the aligner with a plain `which("snap-aligner")` and gates nothing on its version, so a patched binary is a drop-in needing no change to `orp_env.yml`, pyTransRate or ORP. Do not overwrite the binary inside the env — conda owns that file and a later `mamba env update` will silently revert it. Note that this changes the assembly, not just the report: the run gets past the crash, and the contig scores it then produces are what `pick_best_contigs.py` selects on.

## Updating Trinity

Change `TRINITY_COMMIT` in `scripts/build_trinity.sh` and re-run the step 3 command; it fetches, checks out the new commit and rebuilds.

To move an existing install off bioconda's Trinity, install everything the `orp_trinity` list in step 2 names, then remove the package and build. Do these in this order, and install the whole list rather than only the build tools: bioconda's `trinity` was what pulled in `kmer-jellyfish`, `bowtie2`, `samtools`, `openjdk` and `perl`, so removing it takes them away unless they are named.

```bash
mamba install -n orp_trinity -c bioconda -c conda-forge --override-channels compilers "cmake<4" make git kmer-jellyfish bowtie2 samtools salmon=1.10.3 "openjdk>=17" perl perl-db_file python numpy bwa=0.7.19 bashplotlib seqtk=1.5 libgomp zlib
mamba remove -n orp_trinity trinity
conda run --no-capture-output -n orp_trinity scripts/build_trinity.sh
```

Don't do this while a run is using the env: the build replaces `Trinity`, and a run in its Trinity stage needs its files.
