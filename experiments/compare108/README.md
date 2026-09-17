# Full ORP over the 108 pytransrate test assemblies

Array job over an existing manifest, then one table of metrics. Outputs are named
by SRR.

## Manifest

```
./make_manifest.sh
wc -l /mnt/home/macmaneslab/macmanes/compare/manifest.tsv
```

Walks `$COMPARE/pairs/tsa_XXXX/` and writes `$COMPARE/manifest.tsv`, one
headerless tab-separated line per sample:

```
tsa_GADU <TAB> assembly.fasta.gz <TAB> SRR527277_1.fastq.gz <TAB> SRR527277_2.fastq.gz <TAB> SRR527277
```

Column 2 is the TSA assembly. `oyster.py` has no input for it -- bringing your own
assembly is `chowder.py`'s job on this branch -- so the full-ORP run reads and
ignores it, and a sample missing one is still listed, with that column empty.

Column 5 is the run name: what `--runout` gets and what the run directory is
called. It is settled here rather than parsed out of the R1 filename inside each
array task, because only this script sees all 108 samples at once -- a task
cannot tell that its accession belongs to two of them.

### What it does with the awkward cases

**No usable read pair** -- left out, since there is nothing to run, and the
message lists what the `reads/` directory actually contains, so an empty
directory and an unrecognised naming convention are distinguishable at a glance.

**Two samples on one accession** -- three outcomes, decided by what the read
files really are rather than by their paths:

| what the reads are | outcome |
| --- | --- |
| resolve to the same files (symlink/hardlink) | folded to one run; the pairing recorded in `manifest.folded.tsv` |
| different paths, identical sizes | both run, under `SRR..._tsa_XXXX`, with a warning that they may be duplicates |
| genuinely different | both run, under `SRR..._tsa_XXXX` |

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
./submit.sh
```

That is the whole thing -- nothing needs to be exported first. `submit.sh` holds
the one definition of `RUNS`, `MANIFEST` and `ORP`, passes them into the job with
`--export`, sets the array range from the manifest and the throttle to 6, and
creates `$RUNS/logs` before submitting.

That last part matters: slurm opens the `--output` file when a task launches,
*before* the job script runs, so submitting with a missing log directory kills
every task instantly, writing nothing anywhere. If a run vanishes without a trace,
check that directory first.

Paths default off `COMPARE=/mnt/home/macmaneslab/macmanes/compare`: the manifest
is `$COMPARE/manifest.tsv` and the runs go in `$COMPARE/orp_runs/`, which
`submit.sh` creates. Override any of them on the command line:

```
COMPARE=/some/other/place THROTTLE=4 ./submit.sh
MANIFEST=/elsewhere/subset.tsv ./submit.sh      # e.g. to rerun a handful
```

Six samples at a time, each 24 cpus and 120G: 144 cores and 720G in flight.

Each sample runs in `$RUNS/<SRR>/` with `--runout <SRR>`, so the run trees are
fully isolated and every file inside carries the SRR. Logs are
`$RUNS/logs/orp_<jobid>_<task>.log`, symlinked as both `<SRR>.log` and `<tsa>.log`.

Reruns are cheap: finished samples exit immediately on the `qualreport.<SRR>.done`
guard, and ORP itself resumes a part-finished run from where it stopped. To retry
only the failures, resubmit -- `./submit.sh` again is safe.

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
./collect_metrics.py --runs $RUNS --manifest $RUNS/manifest.tsv -o metrics.csv
```

One row per manifest line, in manifest order: every pytransrate column from
`reports/transrate_<SRR>/assemblies.csv` (including `score` and `optimal_score`),
BUSCO C/S/D/F/M/n from `reports/run_<SRR>.ORP/short*.txt`, unique SwissProt genes
and proper-pair rate. Columns come from the csv header by name, not by position.

Samples still running or failed get blank cells and a reason in `status`, and are
listed on stderr, so it is safe to run this mid-array to see progress.
