# Full ORP over the 108 pytransrate test assemblies

Array job over an existing manifest, then one table of metrics. Outputs are named
by SRR.

## Manifest

```
./make_manifest.sh
wc -l /mnt/home/macmaneslab/macmanes/compare/manifest.tsv
```

Walks `$COMPARE/pairs/` and writes `$COMPARE/manifest.tsv`, one headerless
tab-separated line per sample. Sample directories are those whose names start
with `tsa_` (the TSA-derived samples) or `orp_` (the ones assembled here);
anything else under `pairs/` is ignored rather than warned about. Override with
`SAMPLE_PREFIXES="tsa_ orp_ other_" ./make_manifest.sh`.

```
tsa_GADU <TAB> assembly.fasta.gz <TAB> SRR527277_1.fastq.gz <TAB> SRR527277_2.fastq.gz <TAB> SRR527277
```

Column 1 is the sample directory name, prefix included; the collector's `code`
column is the same thing with the prefix stripped, for joining against tables
that carry the bare identifier.

Column 2 is the pre-existing assembly. `oyster.py` has no input for it --
bringing your own assembly is `chowder.py`'s job -- so the full-ORP run reads and
ignores it, and a sample missing one is still listed, with that column empty.

Column 5 is the run name: what `--runout` gets and what the run directory is
called. It is settled here rather than parsed out of the R1 filename inside each
array task, because only this script sees all 108 samples at once -- a task
cannot tell that its accession belongs to two of them.

### What it does with the awkward cases

**No usable read pair** -- left out, since there is nothing to run, and the
message lists what the `reads/` directory actually contains, so an empty
directory and an unrecognised naming convention are distinguishable at a glance.

**Reads under 76 bp** -- left out. The length is the median of the first 100
records of each mate, and both mates have to clear the floor, so a 100/50 bp
library is excluded on its reverse read. The median rather than the maximum
because a submission trimmed before upload has a tail of short reads that should
not disqualify it, and one surviving full-length read should not qualify it
either. Raise or lower with `MIN_READ_LEN=100 ./make_manifest.sh`; sample more
records with `SAMPLE_READS=1000`.

The floor is there because the assemblers' k-mers make short reads a poor
bargain: rnaSPAdes runs at k=55 and k=75 in this pipeline, so a 50 bp library
contributes nothing to the second and little to the first.

**Two samples on one accession** -- three outcomes, decided by what the read
files really are rather than by their paths:

| what the reads are | outcome |
| --- | --- |
| resolve to the same files (symlink/hardlink) | folded to one run; the pairing recorded in `manifest.folded.tsv` |
| different paths, identical sizes | both run, suffixed with the sample directory (`SRR888888_orp_DUPY`), with a warning that they may be duplicates |
| genuinely different | both run, suffixed with the sample directory |

Folding only happens on proof that it is one set of files. Everything else runs
both ways, because spending the compute twice is recoverable and silently
dropping a sample is not.

## 0. Check the cluster checkout

The metrics only come out of pytransrate if the cluster's copy of the ORP is on a
branch that has the swap. `master` is back at 3.1.0 content:

```
cd $HOME/Oyster_River_Protocol && git rev-parse --abbrev-ref HEAD && cat version.txt
```

Expect `byo-assemblies` (or `pytransrate`) and `4.0.0`, not `master` / `3.1.0`.

## 1. Run

```
./submit_all_orp.sh /mnt/home/macmaneslab/macmanes/compare/manifest.tsv \
                    /mnt/home/macmaneslab/macmanes/compare/orp_runs
```

Both paths are given on the command line; a third argument sets the throttle,
which defaults to 6 concurrent tasks. `ORP=/path/to/oyster.py` overrides which
pipeline runs, defaulting to `$HOME/Oyster_River_Protocol/oyster.py`.

```
./submit_all_orp.sh manifest.tsv /scratch/orp_runs 12
```

Relative paths are resolved before being handed over, because a task resolves
neither against its own cwd. The runs directory is created if it does not exist.

Before submitting it checks the manifest is 5 tab-separated columns with unique
run names, that every R1/R2 in it exists, that `orp_array.sbatch` is the same
vintage as itself, and that the ORP checkout is on a branch with pytransrate.
Each of those otherwise fails one task at a time, at 24 cpus apiece.

It also creates `<runs>/logs` before calling sbatch. Slurm opens the `--output`
file when a task launches, *before* the job script runs, so submitting with a
missing log directory kills every task instantly, writing nothing anywhere. If a
run vanishes without a trace, check that directory first.

Six samples at a time, each 24 cpus and 120G: 144 cores and 720G in flight. The
allocation is stated once, in `orp_array.sbatch`'s `#SBATCH` directives; the job
reads `--cpu` and `--mem` back out of what slurm granted, so changing a directive
changes what ORP is told it has.

Each sample runs in `<runs>/<run name>/` with `--runout <run name>`, so the run
trees are fully isolated and every file inside carries that name. Logs are
`<runs>/logs/orp_<jobid>_<task>.log`, symlinked as both `<run>.log` and
`<tsa>.log`.

Reruns are cheap: finished samples exit immediately on the `qualreport.<run>.done`
guard, and ORP itself resumes a part-finished run from where it stopped. To retry
only the failures, submit the same command again.

### When a task fails

```
sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed,NodeList,Reason
```

A task that produced no log at all failed before the job script started: log
directory missing, or the node rejected the allocation. A task whose log stops
after the `=== task N of M ===` line failed in `module load` or manifest parsing.
Past that line the log names the sample, and ORP's own output follows.

## 2. Collect

```
cd $HOME/Oyster_River_Protocol/experiments/compare108
./collect_metrics.py -o metrics.csv
```

Paths default off `$COMPARE`, so there is usually nothing to pass. One row per
manifest line, in manifest order, identified four ways:

| column | example | what it is |
| --- | --- | --- |
| `run` | `SRR527277` | the run name: `--runout`, and the run directory |
| `srr` | `SRR527277` | the accession, from the R1 filename |
| `code` | `GADU` | the bare 4-letter code |
| `tsa` | `tsa_GADU` | the directory under `pairs/` |

`run` and `srr` differ only for the samples that shared an accession and were
disambiguated: there `run` is `SRR807358_tsa_BBBA` while `srr` stays `SRR807358`.

Then the metrics: every pytransrate column from
`reports/transrate_<run>/assemblies.csv` (including `score` and `optimal_score`),
BUSCO C/S/D/F/M/n from `reports/run_<run>.ORP/short*.txt`, unique SwissProt genes
and proper-pair rate. Columns come from the csv header by name, not by position.

Safe to run while the array is still going, and built for it: a run that has not
produced both pytransrate and BUSCO yet is left out of the csv entirely and
listed on stderr instead, so every row in the table is a row with numbers in it.

```
wrote 64 rows to metrics.csv (58 finished, 6 with metrics but still running)
left out 44 of 108 runs, no metrics yet:
  SRR500475: no busco dir; no transrate csv
  SRR544868: not started
```

`status` then separates a run that reached the end (`complete`) from one whose
metrics are in but whose last steps are not (`still running`) -- those write
`qualreport.<run>.done` after pytransrate and BUSCO. Pass `--include-partial` to
get a row for every manifest line regardless, cells empty.

### unique_genes_ORP and proper_pairs

These come from `assemblies/working/<run>.unique.ORP.txt` and `<run>.flagstat`,
and ORP's end-of-run cleanup deletes both -- so for exactly the runs that
finished, reading them directly returns nothing. The collector falls back to
`reports/qualreport.<run>`, which reportgen writes before the cleanup runs and
which survives it. A live run's own files stay authoritative where they exist.
