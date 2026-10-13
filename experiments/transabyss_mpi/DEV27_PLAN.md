# dev27 plan: make the phase 2 restart safe, then decide whether growth pays

Written 2026-10-13 from the dev26 cluster runs (NOTES.md 2026-10-13). Named
dev27 as the follow-up to dev26; `version.txt` is already 4.1.0-dev28 from
notes-only pushes, so the code ships as the next bump after that.

## What dev26 showed

Six dev26 runs finished, three failed (table and timings in NOTES.md
2026-10-13). In short:

1. **The grow restart breaks runs.** Of the six runs whose phase 2 grew, three
   failed (SRR1138704, SRR796576, SRR1138707). In each, one component that
   was in flight at the SIGTERM had Trinity's inchworm checkpoints
   (`.iworm.{18..24}.SR.ok`, `iworm_renamed.*.SR.ok`) written but not the
   `inchworm.DS.fa.SR.N` files the next step reads, and no
   `inchworm.DS.fa.finished`. On relaunch the component skips to that step,
   fails on `cat ... inchworm.DS.fa.SR.20: No such file or directory`, fails
   identically on every retry (~11 per run), and ParaFly exits 256 at the end
   of phase 2. All the broken files carry the stop's timestamp to the second.
   Step 1 of DEV26_PLAN.md interrupted 8 components and missed this window.
   An outage or a Slurm requeue that kills phase 2 can leave the same state,
   so this is not only a growth problem.
2. **Phase 2 barely scales past ~20-25 cores here.** Same dataset, dev24 ->
   dev26: SRR1139191 4h58m on 28 cores -> 4h32m on 40 (-9%); SRR866209 1h24m
   on 29 -> 1h14m mostly on 40 (-12%); SRR868662 5h24m on 22 -> 5h23m on 20
   then 40 (no change); DRR046632 2h38m on 24 -> 2h39m on 20 then 40. ParaFly
   is the repo's OpenMP build (6 GOMP symbols in
   `software/trinityrnaseq/trinity-plugins/BIN/ParaFly`), and the relaunch
   ran `-CPU 40`, so it is not the serial-bioconda problem. Unexplained.
3. **Phase 1 on 13 cores costs time on some datasets:** SRR1139191 +48m (19
   cores in dev24), DRR030368 +27m (15), SRR866209 +18m (19). Others were
   faster on 13 (SRR868662 -41m, DRR031870 -17m): node noise is as large as
   the effect.
4. **TOTAL is roughly a wash.** Over the six finished runs: +22m, -40m, +10m,
   -1h00m, +4m, +41m vs dev24, net -23m. The two largest moves line up with
   node138 (slow; DRR031870's dev24 and DRR030368's dev26 ran there). The ~5h
   saving DEV26_PLAN.md expected assumed phase 2 scales linearly (point 2).
5. **Quality is unchanged** in all six: BUSCO the same (two S/D splits moved
   by 2 of 125), unique genes within tens, pyTransRate within ~0.005.
6. **The 0.5 MPI share balanced the lanes well** where it was tested without
   growth: DRR031870 (TA 3h51m vs Trinity lane 3h53m), DRR030368 (5h30m vs
   5h23m).

## Step 1: make phase 2 safe to restart (code, small)

In `run_trinity_phase2` (`oyster.py`), before **every** launch of phase 2,
not only the grow relaunch:

- Read `recursive_trinity.cmds` and `recursive_trinity.cmds.completed`
  (`trinity_phase2_left()` already reads both). For each command not in
  `.completed`, take its `--output` dir (`.../cN.trinity.reads.fa.out`) and
  remove it if it exists. Those components restart from scratch; each is
  seconds of work.
- Log how many were removed, e.g. `[run_trinity_phase2] reset 3 interrupted
  components`.
- Doing it on every launch also makes a re-run of `oyster.py` in the same
  `--dir` after an outage or requeue safe. On a fresh run nothing matches.
- Match commands exactly as ParaFly writes them to `.completed`. Check one
  real `.completed` line against its `.cmds` line on Premise first, since
  whitespace or quoting differences would make every component look
  unfinished and wipe finished work.

Test:
- `tests/stub_assemblers.py`: add a component that writes its checkpoint
  files and is killed before its output, then a relaunch that must succeed.
  Keep the existing 6 cases passing. Then `tests/release_check.py --tier
  quick`.
- Cluster: the three failed dev26 runs still have their Trinity dirs
  (`~/SRR1138704_dev26`, `~/SRR796576_dev26`, `~/SRR1138707_dev26`). Re-run
  `oyster.py` in each `--dir` with the fix. It should reset the broken
  components, finish phase 2 and the run. That tests the fix on the real
  failure state, at no extra cost.

## Step 2: find out why phase 2 doesn't scale (measurement, no code)

Whether growth is worth keeping depends on this. On one dataset with a
known phase 1 (SRR1138704's, as in DEV26_PLAN.md step 1 /
`parafly_resume.sbatch`), run phase 2 alone at `--CPU` 10, 20, 30 and 40:

- with `--exclusive` (rule out other jobs sharing the node: the sbatch
  scripts don't ask for it and the partitions allow it);
- with `cpu_sampler.py` on the ParaFly tree, to see how many cores phase 2
  actually keeps busy, minute by minute;
- with the phase 2 dir on node-local disk (`/tmp` or `$TMPDIR`) for one run,
  against GPFS for the others: per-component work is many small files.

Read the result as: linear with `--exclusive` -> sharing was the cause, and
growth may pay; flat on GPFS but linear on local disk -> I/O is the cause,
and staging phase 2 locally is the bigger win; flat everywhere -> a long tail
of large components bounds phase 2, and growth can't help.

## Step 3: decide (after steps 1-2)

- Phase 2 scales -> keep growth (now safe) and the 0.5 share, and re-run the
  three failed datasets plus SRR802670 and DRR053698 if they need it.
- It doesn't -> drop growth (stop/relaunch adds risk for no gain). Keep the
  0.5 MPI share, which balanced the lanes on its own. Consider giving phase
  1 more of the Trinity lane's cores (SPAdes gets 7 of 20), since phase 1 is
  where dev26 lost time.
- Either way, keep step 1's cleanup on every phase 2 launch.
- Then `tests/release_check.sbatch` (full tier), and the `--transabyss-mpi
  auto` default decision.

## Housekeeping

- Premise's checkout lacks commit `0ba2e6c` (`--exclude=node117,node118` in
  `oyster_mpi.sbatch`). `git pull --ff-only` there before submitting from the
  repo copy, or pass `--exclude=node117,node118` on the `sbatch` line.
- balance.py's `DEV24` table keys DRR030368 as `DRR30368` (the dev24 run's
  name), so dev26's `DRR030368` gets no "vs dev24" line; same for
  watch_dev26.sh's qualreport pairing (`DRR30368_dev24`, `DRR31870_dev24`).
  Map the zero-padded names.
