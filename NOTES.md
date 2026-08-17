# Working notes (python_convert branch)

Cross-machine scratchpad so a session on either machine can pick up where
the other left off. Keep entries short; newest on top. Delete/trim once
stale.

## 2026-08-17 (3)

- User asked whether `oyster.py` lost any functionality vs. `oyster.mk` at
  tag `2.3.3`. Went target-by-target/method-by-method (both files read in
  full, plus `scripts/` diffed against `2.3.3`). Every pipeline rule from
  2.3.3 has a corresponding, functionally-equivalent method in `oyster.py`
  -- the many differences found are all changes already documented in
  `changelog.md` under 2.4.0/3.0.0 (env consolidation, BUSCO ODB v12.2,
  `build_list5.py`/`pick_best_contigs.py` replacing shell loops, per-step
  timing, concurrency, etc.), not accidental regressions. Variables in
  2.3.3 that `oyster.py` doesn't carry (`TRINITY_KMER`, `BUSCODB`, `START`,
  `LOWEXPFILE`, `RCORR`/`RCORRDIR`, `BUSCO`/`BUSCODIR`) were already dead/
  unused in 2.3.3 itself.
- Two genuine gaps found, both worth a decision rather than silently
  carrying forward:
  - **`clean:` target has no equivalent.** 2.3.3's `clean` (still present
    in the current `oyster.mk`, line ~452) deletes a `RUNOUT`'s
    intermediate files so a rerun starts fresh. `oyster.py` has no
    `--clean`/equivalent; its resumability (skip steps whose outputs are
    newer than their inputs) is a different mechanism aimed at resuming
    after failure, not deliberately discarding a run's outputs.
  - **`preprocess`/`update_merge` partial-pipeline targets have no
    equivalent.** 2.3.3 let you run just trimming+correction
    (`preprocess`) or resume from the merge stage onward
    (`update_merge`) as separate `make` invocations. `oyster.py`'s
    `main()` always runs the full pipeline top to bottom -- resumability
    covers the `update_merge` use case implicitly (already-done early
    steps get skipped), but there's no way to ask it to stop after
    preprocessing only.

## 2026-08-17 (2)

- User hit `rnaspades.py` failing with SPAdes's known "Python version 3.6.8
  is not supported!" bug (ablab/spades#1319) under the 80/20 lane split.
  Traced it rather than assuming the lane split caused it:
  - Concurrency between Trinity and the short-assembler lane (which is
    where SPAdes runs) predates the lane split -- assemblers have run
    concurrently since [2e7ee5b](oyster.py) (2026-08-13), a full day before
    the existing `python=3.14` fix for this exact bug landed in
    `orp_spades` ([b05da12](Makefile), 2026-08-14). So the lane-ratio
    change on 2026-08-16 isn't a new variable for this failure mode.
  - `conda run -n orp_spades python --version` -> 3.14.6 (pin is in place).
  - User's hypothesis that TransAByss's env was leaking into SPAdes (they
    run sequentially in the same short lane) doesn't hold up:
    `conda run -n orp_transabyss python --version` -> 2.7.15, which matches
    neither the reported bad version (3.6.8) nor `orp_spades`'s own
    (3.14.6). 3.6.8 is the cluster's system default Python, so SPAdes is
    finding *that* somehow, not confused by a neighboring conda env.
  - `conda run --no-capture-output -n orp_spades rnaspades.py --version`
    succeeds cleanly (reports 4.3.0, no complaint) -- but this likely just
    short-circuits before whatever internal check does the real
    python-version validation, so it doesn't prove the real-run path is
    clean. The user confirmed the real failure happens right at the start
    of an actual run, consistent with that being the first point the real
    validation path executes, not evidence of a mid-run env corruption.
  - **Not yet confirmed / still open.** Never got the actual traceback or
    SPAdes log (`spades.log`/`params.txt` in the run's `.spades_k*`
    working dir) from a real failed run -- that's the one thing that would
    show definitively which file/binary does the bad python lookup, and
    whether it's a fixable `PATH`-order issue or (per other reports on
    ablab/spades#1319) a python path baked into a compiled SPAdes
    component at bioconda build time, which wouldn't be fixable from this
    repo at all.

- Added general step-retry + failure-visibility handling in
  [oyster.py](oyster.py), prompted by the above (a step failing used to
  kill the whole run immediately, and even now, one lane failing stays
  silent until the other lane -- often Trinity, tens of hours out --
  finishes):
  - `run()` ([oyster.py:164](oyster.py#L164)) now retries a failed
    subprocess up to `STEP_RETRIES` (2) times with `STEP_RETRY_DELAY`
    (60s) between attempts before propagating, for transient cluster
    failures (node preemption, filesystem hiccups). Does nothing for
    deterministic failures like the SPAdes bug above -- those just fail
    the same way 3 times.
  - Added `retry_cleanup=` (path or paths to `rmtree` before each retry),
    wired into `run_spades()` and `run_transabyss()` -- both fail outright
    on a non-empty `-o`/`--outdir` from a prior attempt rather than
    resuming, so a bare retry would hit a different, unrelated error
    instead of actually re-attempting the assembly. Trinity is exempt from
    cleanup (and now from retries entirely, next bullet) since it resumes
    from its own checkpoints in-place.
  - `run_trinity()` now passes `retries=0` -- the blanket default would
    otherwise let a deterministically-failing ~37h Trinity step retry
    twice more (up to ~111h) before finally giving up.
  - `trinity_lane()`/`short_assembler_lane()` in `main()` now log
    immediately (`*** [lane] lane failed ... ***`) the moment either
    raises, since `ThreadPoolExecutor`'s context manager still blocks on
    `shutdown(wait=True)` for the other lane before the exception actually
    propagates and the run exits -- previously a short-lane failure gave
    no signal at all until Trinity separately finished.
  - Not yet real-run validated (no cluster access from this session).

## 2026-08-17

- New timing run from the user (weekend, home machine), SRR1789336,
  `--max-parallel 1 --normalize-reads`, TOTAL 43:39:49 vs the 38:28:11
  `--max-parallel 2` (no `--normalize-reads`) baseline in the entry below.
  Trinity/TransAByss/merge/transrate/busco all ~flat between the two runs;
  the entire regression is `run_orthofuser`: 5:15 -> 30:57 (~6x). Ruled out
  parallelism as the cause -- checked `run_parallel()`
  ([oyster.py:194](oyster.py#L194)): a solo job under `--max-parallel 1`
  gets the *full* `--cpu`/`--mem` (`workers=1`), not less, so orthofuser
  actually had *more* CPU this run (40 vs ~20 under mp=2) and was still 6x
  slower.
  User then mentioned they also upgraded OrthoFinder locally sometime
  between the two runs. Confirmed via `git show a7e5d62` (2026-08-14) that
  this repo's `orp_orthofinder` env was bumped 2.5.2 -> 3.1.5 in that
  window (changelog.md:7). The `-d -I 12 -f ... -og -t -a` CLI invocation
  in `run_orthofuser()` ([oyster.py:437](oyster.py#L437)) is unchanged
  across that bump, but a major-version algorithm change is a much more
  plausible cause of a 6x wall-time jump that *doesn't* respond to more
  CPU than either `--normalize-reads` or `--max-parallel` are. Leading
  theory now: the OrthoFinder 3.1.5 upgrade, not the other two flags,
  explains the orthofuser regression.
  **Not yet confirmed.** Three variables differ between these two runs
  (`--max-parallel`, `--normalize-reads`, OrthoFinder version) with no
  isolated comparison. To actually confirm, need a run holding OrthoFinder
  version and `--normalize-reads` fixed while only varying one thing at a
  time -- or at minimum, checking OrthoFinder's own log/timing output from
  the two runs (2.5.2 vs 3.1.5 report their internal stage timing
  differently and might show directly where the time went).

## 2026-08-16 (6)

- User suspects read normalization (Trinity's `--no_normalize_reads` /
  `--normalize-reads`) explains the 2.4.0-vs-3.0.0 Trinity diff in entry
  below, and is going to test it themselves. Checked the code first: both
  `oyster.mk` (`NORMALIZE_READS := FALSE` by default) and `oyster.py`
  (`--normalize-reads` is `action="store_true"`, default off) resolve to
  the same `--no_normalize_reads` flag by default, so this isn't a
  default-value mismatch between the two pipelines -- it would only apply
  if the two specific runs being compared used different actual flags/
  commands, which hasn't been confirmed either way. Waiting on the user's
  test result; don't assume this is or isn't the cause yet.

## 2026-08-16 (5)

- Investigating small quality-metric differences between ORP 2.4.0
  (oyster.mk) and 3.0.0 (oyster.py) on SRR1789336:
  ```
                          2.4.0                 3.0.0 (pre-lane-fix baseline)
  BUSCO                   C:88.8% M:3.2%        C:88.0% M:4.8%
  TRANSRATE               0.51558               0.51827
  TRANSRATE OPTIMAL       0.52685               0.53225
  UNIQUE GENES ORP        13810                 13772
  UNIQUE GENES TRINITY    13265                 13231
  UNIQUE GENES SPADES55   13493                 13493   (identical)
  UNIQUE GENES SPADES75   12809                 12809   (identical)
  UNIQUE GENES TRANSABYSS 12632                 12632   (identical)
  READS MAPPED PROPER     95.82%                96.06%
  ```
  SPAdes55/75 and Trans-ABySS are byte-identical between versions, which
  proves the trimmed/corrected reads feeding all four assemblers were
  identical and rules out an rcorrector/Trimmomatic difference. Only
  Trinity differs, and everything else that differs (ORP/BUSCO/transrate/
  mapping rate) is downstream of Trinity's output via the merge stage.
  Confirmed via `git log -p -- Makefile` that `orp_trinity`'s env
  (`trinity=2.15.2`, `bwa=0.7.19`, `seqtk=1.5`, `salmon=1.10.3`) hasn't
  changed across this whole window, and oyster.py's `run_trinity()` passes
  the same flags as oyster.mk's TRINITY rule -- except `--CPU`. oyster.mk
  ran Trinity sequentially at the full `--cpu`; the 3.0.0 run above predates
  tonight's two-lane fix (246804a) and gave Trinity only half `--cpu` via
  the old even `run_parallel()` split. Leading theory: Trinity's Chrysalis/
  Butterfly stages aren't perfectly reproducible across different thread
  counts (known Trinity behavior), so the `--CPU` difference alone plausibly
  explains a ~0.3% divergence with no actual bug involved.
  **To confirm**: compare tonight's in-progress validation run (Trinity now
  at ~80% `--cpu`, see entry below) against these two baselines -- if its
  Trinity/BUSCO/etc. numbers land closer to the 2.4.0 column than the 3.0.0
  column above, that's strong confirmation. Not yet checked.

## 2026-08-16 (4)

- Real-run validation of the two-lane assembler split (commit 246804a) is
  in progress, same dataset as before (SRR1789336). Time to beat: the old
  even-split TOTAL of 38:28:11. Not yet known whether it finished or what
  the new TOTAL/per-step timing looks like -- check with the user or look
  for an updated timing report / NOTES entry before assuming either way.

## 2026-08-16 (3)

- Updated README.md ("Parallel task management" section) and changelog.md
  (Unreleased) to describe the two-lane assembler split and the
  STEP_TIME_HINTS fix from the entry below -- both docs still described the
  old even `--max-parallel` split across all 4 assemblers.

## 2026-08-16 (2)

- Diagnosed why half the cores sit idle for ~36 hours during the assembler
  stage: `run_parallel()` split `--cpu`/`--mem` evenly across whichever jobs
  were pending *at group start*, and that split was fixed for a job's whole
  lifetime -- once TransAByss/SPAdes55/SPAdes75 finish (~62 min in), Trinity
  runs alone at half-CPU for the rest of its ~37h, since its thread count is
  fixed at launch (no way to hand it more CPU once it starts).
- Replaced the assembler `run_parallel()` call in [oyster.py](oyster.py)
  `main()` with two fixed resource "lanes" run concurrently via
  `ThreadPoolExecutor(max_workers=2)`:
  - **Trinity lane**: gets `TRINITY_LANE_SHARE` (0.8) of `--cpu`/`--mem` for
    the whole run.
  - **Short-assembler lane**: gets the remainder, runs TransAByss ->
    SPAdes55 -> SPAdes75 (slowest first), and fires each assembly's diamond
    search (`diamond_transabyss`/`diamond_spades75`/`diamond_spades55`)
    immediately after it finishes, rather than waiting for the full
    orthofuser/merge stage that follows (those diamond jobs only ever
    needed their own assembly fasta, not Trinity or the merge -- confirmed
    by tracing `diamond_jobs()` deps). Only `diamond_orthomerged` and
    `diamond_trinity` remain in the post-merge diamond block, since those
    two genuinely depend on the merge stage / Trinity.
  - `STEP_TIME_HINTS` no longer has assembler entries (dead now that
    assemblers don't go through `run_parallel()`); `--max-parallel` help
    text updated to say it no longer covers the assemblers.
- `TRINITY_LANE_SHARE = 0.8` is a flat constant, not a CLI flag -- decided
  against a flag for now since we don't have data yet on whether 80/20 is
  actually the right ratio.
- Not yet done: no real-run validation of this change (no cluster access
  from this session). Next real run should confirm: (a) `--max-parallel`
  still governs the orthofuser/merge and transrate/strandeval lanes
  correctly, (b) Trinity's wall time actually drops with ~80% instead of
  50% of the cores (there's no guarantee -- Trinity's Butterfly stage has
  diminishing returns past some thread count), and (c) memory headroom is
  fine with `TRINITY_LANE_SHARE=0.8` on real (larger) data, not just the
  sample dataset.
- Local-only so far, not yet pushed to `origin/python_convert` -- ask
  before pushing (see previous entry below).

## 2026-08-16

- Fixed `STEP_TIME_HINTS` in [oyster.py](oyster.py) (~line 43). These hints
  decide submission order inside the three `run_parallel()` groups
  (assemblers; `orthofuser_branch` vs `merge_branch`; `transrate` vs
  `strandeval`) -- the longer job should be submitted first so it isn't left
  running alone after its faster sibling finishes.
- Got a real timing report from a `--max-parallel 2` run (see table below).
  It showed two of the three groups were ordered backwards:
  - `merge_branch` (26:37) is actually much slower than `orthofuser_branch`
    (6:08) -- hints had it the other way.
  - `transrate` (15:33) is much slower than `strandeval` (1:34) -- hints
    had `strandeval` sorting first.
  - Assemblers group was roughly right except spades55 (9:24) vs spades75
    (7:48) were swapped.
  - Also dropped two dead dict entries (`run_orthofuser`, `orthotransrate`)
    that were never actually looked up -- `run_parallel()` only keys on the
    names passed to it (`orthofuser_branch`/`merge_branch`), not the
    sub-steps inside those branch functions.
- Not yet done: rerun the pipeline end-to-end with the corrected hints and
  confirm the new submission order actually shortens wall time vs. the
  38:28:11 TOTAL below (should mainly help the merge/orthofuser and
  transrate/strandeval pairs; assembler group changes are minor since
  trinity dominates regardless of order).

### Reference timing (this repo, `--max-parallel 2`)

```
run_trimmomatic  00:04:21
run_rcorrector   00:14:25
run_transabyss   04:46:42
run_spades75     00:07:48
run_spades55     00:09:24
run_trinity      37:11:49
run_filtershort  00:00:07
run_orthofuser   00:05:15
orthofuser_branch 00:06:08
orthotransrate   00:26:36
merge_branch     00:26:37
makeorthout      00:00:49
orthofusing      00:05:35
diamond_orthomerged 00:00:08
diamond_transabyss 00:00:07
diamond_spades75 00:00:06
diamond_spades55 00:00:06
diamond_trinity  00:00:08
make_list5       00:00:00
posthack         00:00:04
cdhit            00:00:34
orp_diamond      00:00:07
salmon_index     00:00:21
salmon           00:00:29
secondfilter     00:00:54
busco            00:05:34
strandeval       00:01:34
transrate        00:15:33
TOTAL            38:28:11
```
