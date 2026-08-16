# Working notes (python_convert branch)

Cross-machine scratchpad so a session on either machine can pick up where
the other left off. Keep entries short; newest on top. Delete/trim once
stale.

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
