# dev26 plan: phase 2 grows onto Trans-ABySS's cores

**Status 2026-10-12:** step 1 passed (Premise job 1321673), step 2 is
4.1.0-dev26, step 3's cluster runs are next. NOTES.md 2026-10-12.

Written 2026-10-11 at 4.1.0-dev25 (pipeline code = dev24). Background and
every number below: NOTES.md, entries from 2026-10-08 to 2026-10-11.

## Why

Under `--transabyss-mpi`, the assemblers run as two lanes on fixed core
counts: Trans-ABySS (TA) on a share of `--cpu` from the start, and Trinity
(phase 1 || SPAdes, then phase 2) on the rest. dev24 set TA's share from the
size of the corrected R1, fitted to the first datasets. Eight more runs broke
that fit:

- TA's work is predictable: **3.5-5.6 core-hours per GB** of uncompressed
  corrected R1 across 11 runs.
- Trinity's is not: **3.8-22.3 core-hours per GB**, no trend with size. It
  follows transcriptome complexity, which file size can't show.

So any up-front split misses in one direction or the other. Whichever lane
finishes first leaves its cores idle until the other is done. Measured
cost on the dev24 runs ("if balanced" = both lanes finishing together):

| run | R1 GB | slower lane | TOTAL | if balanced |
|---|---|---|---|---|
| SRR868662 | 25.54 | Trinity by 3h57m | 11h49m | 10h02m |
| SRR802670 | 53.42 | TA by 1h53m (0.45 cap) | 15h27m | 14h24m |
| DRR053698 | 33.23 | TA by 1h45m (0.45 cap) | 11h29m | 10h31m |
| SRR1139191 | 10.90 | Trinity by 2h47m | 9h33m | 8h43m |
| SRR1138707 | 6.24 | Trinity by 50m | 5h27m | 5h14m |
| SRR796576 | 3.73 | Trinity by 49m | 3h16m | 3h05m |
| DRR30368 | 21.83 | TA by 25m | 7h31m | 7h16m |
| DRR046632 | 18.44 | Trinity by 1m | 6h42m | 6h41m |

TOTAL wall time is the metric that matters (it is what users see).

## The change

1. **Phase 2 grows when TA finishes.** Phase 2's ParaFly runs with a core
   count fixed at launch, but ParaFly logs each finished command to
   `<trinity_out_dir>/recursive_trinity.cmds.completed` and skips those when
   re-run. So when the TA lane finishes while phase 2 is still running,
   `oyster.py` stops phase 2 and relaunches it on all of `--cpu`. Trinity
   then gets TA's cores instead of leaving them idle.
2. **A fixed MPI share of 0.5.** With (1), the share only has to keep TA
   from being the long pole: any cores TA doesn't need come back to Trinity.
   The largest balanced share measured is 0.511, so 0.5 covers every run so
   far. The R1-size fit (`TRANSABYSS_MPI_SHARE_BASE`/`_PER_GB`/`_MIN`/`_MAX`)
   goes. Threaded mode stays at 0.25 (TA is the long pole there whatever it
   gets, so growth never triggers).

Expected: each run lands near its "if balanced" TOTAL above, ~5h over these
eight. Caveats: that estimate assumes Trinity scales linearly with cores
(roughly true of phase 2); Trinity runs ~16-30% slower beside MPI TA
whatever the split; phase 1 gets fewer cores at first (lane 20 at --cpu 40,
phase 1 ~13), and it scales less than phase 2 (Inchworm capped at 10).

## Step 1: verify ParaFly resumes (before writing any pipeline code)

The plan rests on this. Test it on the cluster with a real phase 1:

1. Run Trinity phase 1 only, as oyster.py does (`_trinity_base_cmd` +
   `--no_distributed_trinity_exec`, `--normalize-reads`), on SRR1138704's
   corrected reads (its phase 2 is ~1h20m on 30 cores: long enough to
   interrupt, short enough to finish). Use a scratch `--output` dir.
2. Start phase 2 (same command, no stop flag, no `--full_cleanup`) at
   `--CPU 8`, in its own process group/session.
3. After ~15 minutes, record `wc -l recursive_trinity.cmds` and
   `wc -l recursive_trinity.cmds.completed`, then stop the whole group
   (SIGTERM, then SIGKILL after a grace period). Confirm no Trinity,
   ParaFly, Butterfly (java) or Chrysalis processes survive.
4. Relaunch phase 2 at `--CPU 40` (or the node's cores). Check:
   - it starts from the completed count, not zero (ParaFly prints how many
     commands it is skipping, and the `.completed` count only grows);
   - no command runs twice (its output dir isn't rebuilt from scratch);
   - it finishes and writes `<outdir>.Trinity.fasta`;
   - the run is not broken by components that were mid-assembly when
     killed (partial `read_partitions/.../*.out` dirs).
5. Compare with an uninterrupted phase 2 of the same phase 1 (copy the
   phase-1 dir before step 2): transcript count, and gene count after the
   ORP diamond step if convenient. Byte identity is not expected (Trinity
   depends on run dir and JVM hash, NOTES.md 2026-10-05), so judge by
   counts. Also note the restart overhead: how long from relaunch to
   ParaFly running jobs again.

If ParaFly does not resume cleanly, stop and report: the plan changes
(fallback ideas: oversubscribe phase 2 with all of `--cpu` from the start,
which costs TA when TA is the long pole; or Trinity's `--grid_exec` with a
local launcher).

## Step 2: implement in oyster.py (dev26)

- `run_trinity_phase2`: launch Trinity with `subprocess.Popen(...,
  start_new_session=True)` (through `conda run --no-capture-output -n
  orp_trinity`, like `conda_run`) so the whole tree can be signalled, and
  return a handle the lane can watch, instead of blocking in `self.run`.
  Keep retries/`--full_cleanup`/output handling as they are.
- `run_assemblers.trinity_lane`: start phase 2 on the lane's cores (or all
  of `--cpu` if TA is already done, as now). Poll the process and the
  existing `transabyss_done` event. When TA finishes and phase 2 is still
  running, read `recursive_trinity.cmds` and `.completed`; if enough is left
  (start with: at least 15% of commands not completed -- make it a named
  constant), stop the group (SIGTERM, grace period, SIGKILL) and relaunch
  phase 2 on `self.cpu`/`self.mem`. Once only.
- Log it: `=== Stage B: run_trinity_phase2 grows to 40 cpu (Trans-ABySS
  done; 12345/45678 components left) ===`, and record the restart in the
  timing log so balance.py and NOTES can see it (e.g. a
  `run_trinity_phase2_grow` entry, or extend `step()`'s timing).
- A failure in the relaunched phase 2 must fail the step exactly as today
  (marker left, run aborts after the other lane, `_lane_failed`).
- Interrupted oyster.py while phase 2 runs (SIGTERM to the job): the
  Trinity session must die with it -- it no longer shares oyster.py's
  process group. Handle that (signal handler or atexit that kills the
  session).
- `transabyss_share()`: MPI -> `TRANSABYSS_MPI_SHARE = 0.5`; keep logging the
  R1 size (balance.py and future analysis use it); drop the fit constants.
  Update the comments above `TRANSABYSS_SHARE`.
- `balance.py`: its balanced share assumes fixed lanes. With growth, report
  per run whether phase 2 grew and when; TOTAL against dev24's for the same
  dataset is the measure.
- Docs: README lane paragraph, `docs/pipeline-steps.md`,
  `docs/pipeline-schedule.html` (Stage B node and the concurrency list),
  `changelog.md`, NOTES.md, `version.txt` -> 4.1.0-dev26.

## Step 3: test

- Laptop: extend the stub harness used for dev19-dev24 (drive
  `run_assemblers` with stubbed steps; see NOTES 2026-10-08) to cover: TA
  done before Stage A ends (phase 2 starts on all cores, no restart); TA
  done mid-phase 2 with plenty left (one restart, on all cores); TA done
  near the end (no restart); relaunched phase 2 failing; threaded mode
  (never grows). Then `tests/release_check.py --tier quick`.
- Cluster: `tests/release_check.sbatch` (full tier, 22 cases), then dev26
  MPI runs of SRR868662 and SRR1139191 (Trinity-bound) and DRR053698 and
  SRR802670 (TA-bound), each in a fresh `--dir`. Compare TOTAL with dev24's
  and with "if balanced" above, and the qualreports with dev24's.
- Only then consider `--transabyss-mpi auto` as the default (quality
  already passed on three datasets; NOTES 2026-10-10).
