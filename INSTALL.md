# Installing the Oyster River Protocol

These directions cover installing ORP 3.0+ (the `oyster.py` pipeline) on Linux. Two paths are covered: the one-command installer, and the same steps done by hand if you want more control or need to troubleshoot a partial install.

## Prerequisites

- Linux (the installer's Anaconda bootstrap is Linux x86_64 specific)
- `git`, `curl`, `bash`, and a system `python3` already available
- Internet access on the install machine — it pulls down Anaconda, the UniProt/Swiss-Prot diamond database, and the BUSCO lineage database
- Several GB of free disk space (Swiss-Prot + Anaconda + 6 conda environments + BUSCO database adds up)

## Option A: one-command install

```bash
git clone https://github.com/macmanes-lab/Oyster_River_Protocol.git
cd Oyster_River_Protocol
make
```

`make` chains through everything: creates the conda environments (including OrthoFinder's), unpacks orp-transrate, builds the diamond search database, downloads the BUSCO lineage database, and appends any needed PATH entries to `~/.profile`/`~/.bash_profile`. Each step checks whether it's already done and skips itself if so, so re-running `make` after a partial or failed install picks up where it left off rather than starting over.

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

**2. Create the 5 isolated tool environments**

Only `orp_spades` pins `python=3.14` — that env has a confirmed bug (see [Known gotchas](#known-gotchas)) that requires forcing a modern Python. The others resolve whatever compatible Python their own recipe naturally wants; forcing 3.14 there bought no benefit (each env's Python is invisible outside that env) and broke `orp_transabyss`'s solve in practice.

`orp_orthofinder` is kept isolated rather than folded into the shared `orp` env because OrthoFinder 3.x requires `diamond <2.2`, which conflicts with `orp`'s own `diamond=2.2.5` (used for the Swiss-Prot search, a separate purpose from OrthoFinder's internal one).

```bash
mamba create -y -c bioconda -c conda-forge --override-channels --name orp_spades spades=4.3.0 python=3.14
mamba create -y -c bioconda -c conda-forge --override-channels --name orp_trinity trinity=2.15.2 bwa=0.7.19 bashplotlib seqtk=1.5 salmon=1.10.3
mamba create -y -c bioconda -c conda-forge --override-channels --name orp_busco busco=6.1.0
mamba create -y -c bioconda -c conda-forge --override-channels --name orp_transabyss transabyss=2.0.1
mamba create -y -c bioconda -c conda-forge --override-channels --name orp_orthofinder orthofinder=3.1.5
```

**3. Create the consolidated `orp` environment**

Everything else — rcorrector, trimmomatic, cd-hit, diamond, salmon (the pipeline's own, modern version), samtools, seqtk, mcl, sra-tools, blast, parallel, biopython, scipy, numpy, bashplotlib — lives in one `orp` environment, defined in `orp_env.yml`:

```bash
mamba env create -f orp_env.yml
```

**4. OrthoFinder**

Already created in step 2 above (`orp_orthofinder`) — nothing further needed here.

**5. orp-transrate** (already bundled in the repo — just unpack it)

```bash
cd software
tar -zxf orp-transrate.tar.gz
cd ..
```

**6. Diamond's Swiss-Prot search database**

```bash
mkdir -p software/diamond
cd software/diamond
curl -LO ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gzip -d uniprot_sprot.fasta.gz
conda run -n orp diamond makedb --in uniprot_sprot.fasta -d swissprot
cd ../..
```

**7. BUSCO's eukaryota lineage database**

```bash
mkdir -p busco_dbs
conda run -n orp_busco busco --download eukaryota_odb12.2 --download_path busco_dbs
```

**8. Add orp-transrate to your PATH**

```bash
echo "export PATH=\$PATH:$(pwd)/software/orp-transrate" >> ~/.profile
source ~/.profile
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

Then run the pipeline once on real (or the bundled sample) data — `oyster.py`'s own preflight `check()` step verifies every tool it needs is actually reachable and fails fast with a clear message naming whatever's missing, rather than failing partway through a multi-hour run:

```bash
python3 oyster.py --read1 R1.fq.gz --read2 R2.fq.gz --runout myrun --cpu 24 --mem 110 --strand RF
```

## Uninstalling

```bash
make clean
```

Removes the conda install and downloaded software directories.

## Known gotchas

- **"Python version 3.6.8 is not supported"** from SPAdes: a [known SPAdes bug](https://github.com/ablab/spades/issues/1319) where it detects a stray system Python instead of its own environment's interpreter. Fixed by `orp_spades`'s explicit `python=3.14` pin above — if you hit this anyway, double check nothing earlier in your `PATH` (an HPC module system, for example) is shadowing the conda environment's own Python.
- **`nothing provides _python_rc needed by python-3.14.0rc1-...`**: means mamba is resolving to a Python 3.14 *release candidate* build because no stable 3.14 build exists yet for that particular package's dependency set. This is why only `orp_spades` pins `python=3.14` — pinning it on older/less-actively-maintained packages (like `transabyss=2.0.1`) can hit exactly this wall.
