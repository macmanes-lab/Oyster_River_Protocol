# Working notes (python_convert branch)

Cross-machine scratchpad so a session on either machine can pick up where
the other left off. Keep entries short; newest on top. Delete/trim once
stale.
	

## 2026-09-19

- **PATH cannot reach OrthoFinder's diamond. config.json can.** The dev4
  shim resolved correctly in the parent (`command -v diamond` -> the shim)
  and was still never called: `/proc/<pid>/environ` on a running search
  showed
  `PATH=<env>/bin:<env>/bin/src/orthofinder/bin:...:<run>/shim:...`
  -- OrthoFinder prepends its environment's bin and its own bundled bin at
  startup, so anything put in front from outside ends up behind them.
  Observed `-p 1` and no `--tmpdir` on all four searches with the shim
  sitting at position 9.
  - `-p 1` is a literal in the `search_cmd` template in
    `<env>/bin/src/orthofinder/run/config.json`. `orthofinder --help` has
    no `--config`, only `-S <txt>` to pick a program by name, so the
    install copy is the only lever.
  - ORP adds `diamond_orp_<threads>` -- the stock entry with `-p` set --
    and runs `-S diamond_orp_<threads>`. Stock entry untouched; original
    backed up once; thread count in the name so concurrent runs at
    different `--cpu` do not fight over one key; nothing run-specific in
    the file, since the cluster shares it.
  - **Three approaches, one lesson:** `env=` lost to conda's activation,
    an exported PATH lost to OrthoFinder's own prepending, and neither
    failure said anything in the log. Each was only caught by `ps` on the
    node. Whatever the mechanism, the step has to print what it actually
    resolved -- that is worth more than the mechanism being clever.
  - Memory, meanwhile, is still unmeasured at the peak: RSS held flat at
    ~1.5 GB per search (the loaded .dmnd) for the first three minutes of
    the 4-way run, ~6 GB against a 670 GB budget. Run 2's kills came at
    2h40m. Sample to the end before touching
    ORTHOFINDER_GB_PER_QUERY_GB.

- **`env=` cannot put anything ahead of a conda environment's own bin.**
  The dev2 diamond shim was written, PATH was set, the step ran -- and the
  shim was never called. `conda run -n X` activates X, and activation
  prepends `$CONDA_PREFIX/bin` to whatever PATH it inherited, so the shim
  ended up *behind* the real diamond. Measured on the node two hours in:
  `ps -o pid,nlwp,pcpu,rss -C diamond` gave `NLWP 2`, `%CPU 98` on all four
  processes, and `-p 1` with no `--tmpdir` on the command line -- four cores
  of forty, on a step that had been told to use ten threads each.
  - Fixed by exporting PATH inside the activated environment
    (`conda run -n X bash -c 'export PATH=shim:$PATH; exec orthofinder ...'`)
    rather than from the parent process.
  - Reproduced both directions before pushing: with the env bin prepended
    afterwards the shim loses, with the export inside it wins.
  - **The real lesson is that it was silent.** Nothing in the log
    distinguished a working shim from a bypassed one, which is why it ran
    for two hours before `ps` was the thing that caught it.
    `report_diamond_in_use()` now resolves `diamond` through the same PATH
    OrthoFinder will use and prints the answer before the searches start.
  - First real memory numbers, 5 minutes into the 4-way run: RSS 3.1, 3.1,
    2.7, 2.7 GB, ~11.6 GB total. Far under the 159 GB/search the model
    budgets -- but the run-2 kills came at 2h40m, not at 5 minutes, so this
    says nothing about the peak yet. Sample it to the end before re-fitting
    ORTHOFINDER_GB_PER_QUERY_GB; 96 GB/GB is still one dmesg line from a
    pre-masking run, and it is the number most likely to be wrong here.

- **The masking was necessary and not sufficient: the same run failed the
  same way, and the reason is that the concurrency cap could never bind.**
  Second attempt at 380C_0C5D_001Fv3_955, same node, `--cpu 40 --mem 670`,
  with `mask_search_input` doing its job (2,885,497 and 2,538,640 N->X for
  spades75/spades55, 0 and 0 for transabyss/trinity). Still 7 of 16
  searches lost, still every failure a SPAdes query: Species0 lost all four
  of its searches, Species1 three of four, Species2 and Species3 none.
  - **Four died on -9, three on exit 1 having written a 20-byte gzip** (an
    empty `Blast*.txt.gz`). All seven between 09:42 and 10:03, after the
    step started at 07:01 -- i.e. two and a half hours in, all at once --
    and the nine survivors finished between 13:05 and 16:08, once there was
    room. That is a node running out of memory collectively, not seven
    inputs each being individually poisonous.
  - **`ORTHOFINDER_GB_PER_SEARCH = 12` had no effect and never could.**
    `searches = min(cpu, mem // 12)` = `min(40, 55)` = 40; OrthoFinder runs
    `n_assemblies^2` = 16 diamonds and no more, so the cap sat above the
    real job count and 16 ran at once regardless. On a 670 GB node the
    constant has to exceed 670/16 = 42 GB before it changes anything at
    all. Every piece of advice the error messages gave -- raise it, lower
    `--cpu`, lower `--max-parallel` -- was advice about a knob that was not
    connected. (`--max-parallel` governs ORP's own step lanes and has never
    had anything to say about OrthoFinder's internal fan-out.)
  - **Sizing is now per GB of the largest search input**, at 96 GB per GB,
    which is the worst measured search (~145 GB peak RSS on a ~1.5 GB
    query) rather than the mean of the five. On this run that is 144 GB
    apiece, so 4 concurrent against 670 GB instead of 16.
  - **Lowering `-t` no longer costs cores.** OrthoFinder hands every
    diamond `-p 1` whatever `-t` says, so `-t 4` on a 40-core node would
    have run four diamonds on four cores -- which is why `-t` was never the
    knob anyone reached for. `write_diamond_shim()` puts a `diamond` ahead
    of the real one on PATH that rewrites `-p 1` to `cpu // searches`, so
    this run is 4 x 10 threads = the same 40 cores at a quarter of the peak
    memory. Only `blastp`/`blastx`/`blastn` are touched; `makedb` passes
    through.
  - **`--orthofinder-searches N`** pins concurrency when the model is wrong
    for a node. The thread count still follows from it.
  - Still unexplained: *why* a SPAdes query is expensive once its N runs
    are gone. Masking removed the poly-asparagine and the ranking did not
    change, so something else in those assemblies -- homopolymer runs read
    as poly-Ala/Gly/Thr is the obvious candidate -- is still making far
    more seed hits than Trinity's or TransAByss's contigs do. The fix above
    bounds the blast radius without knowing the answer; if it matters
    later, `diamond blastp --masking` is where to look.
  - **Recovering a part-done all-vs-all:** OrthoFinder will not recompute a
    search it can see an output file for, so the 20-byte and truncated
    `Blast*.txt.gz` must be deleted, not left in place, before a resume.

- **OrthoFinder searches DNA with `diamond blastp`, so SPAdes' N-gaps arrive
  as poly-asparagine and the searches that carry them OOM.** A chowder run
  (380C_0C5D_001Fv3_955, 4 assemblies, 5,354,958 contigs merged, `--cpu 40
  --mem 670` in a 720 GB cgroup) lost 7 of its 16 all-vs-all searches.
  Every one had a SPAdes assembly as its query; not one transabyss or
  trinity search failed.
  - **`-d` does not switch OrthoFinder to a nucleotide searcher.**
    `Log.txt` says `Search program: diamond` and the command it builds is
    `diamond blastp --ignore-warnings ...` -- the flag exists precisely so
    `makedb` will accept ACGT as protein. Every base is read as an amino
    acid, and **N is asparagine**. rnaSPAdes gap-fills scaffolds with N;
    Trinity and TransAByss emit none: 160,372 and 189,448 contigs with a
    >=10bp N run against 0 and 0. So each SPAdes assembly arrived with
    ~175K poly-asparagine tracts, and each one seeds against every other
    tract in the database.
  - **The kernel settled it.** `dmesg` on the batch host: five
    `constraint=CONSTRAINT_MEMCG, oom_memcg=/slurm/uid_46343/job_1314755,
    task=diamond` kills, at 92, 129, 133, 138 and 145 GB resident -- 638 GB
    between the five, against `ORTHOFINDER_GB_PER_SEARCH`'s 12 GB per
    process. Eleven times the model. The five kills match the five `-9`
    errors in the log exactly; the two `code 1` errors are the same
    pressure by another route.
  - **What the failure looks like from downstream is nothing at all.** All
    16 `Blast*.txt.gz` were valid gzip -- `printf "" | gzip -c` is 20 bytes
    and passes `gzip -t`, so an empty result is indistinguishable from an
    honest one to anything that does not look inside. OrthoFinder printed
    its `ERROR: external program returned code` lines, **exited 0**, and
    built orthogroups from the 9 searches that lived. Both SPAdes
    assemblies had lost their self-comparison *and* their comparison with
    each other, so they survived in the clustering only through hits found
    by the two assemblies that still worked. `run_orthofuser`'s existing
    guard checks that Orthogroups.txt exists, postdates the marker and is
    non-empty, and cannot see any of this.
  - **Four things were ruled out before the cause was found**, all worth
    not re-checking: the filesystem (345 TB free, no quota -- diamond's
    `--tmpdir` defaults to the output directory and OrthoFinder never
    passes `-t`); `MaxRSS 692.6 GiB` as evidence about diamond (it is a
    job-lifetime peak and `orthotransrate`'s 487 GB accumulator phase is an
    equally good candidate); OrthoFinder tearing down siblings after the
    first failure (the kills are spread 23:17-23:48 and the healthy
    searches ran on to 06:31 -- a teardown kills at once); and a first
    `dmesg` that showed nothing, which had been run on a login node rather
    than `BatchHost=node142`. Its one hit, a `global_oom` on `conda` at 494
    GB on Sep 14, is uid 46427 -- a different user, a different node.
  - **`write_search_inputs` is the fix**: OrthoFinder now reads N->X copies
    from `orthofuse/search/`. X is the unknown residue and diamond will not
    seed on it, so the tracts leave the search without a contig being
    shortened or renamed. **Deflines are copied byte for byte** -- a contig
    named `NODE_1_..._NNN` still has to answer to that name in
    Orthogroups.txt. Isolated Ns are translated too: an ambiguous base
    carries nothing to match on, and translating the lot beats deciding
    what counts as a run.
  - **The masked copies are deliberately not the merged ones.** `merge()`
    concatenates `short_fasta_paths()` into merged.fasta, which is what
    pytransrate scores and what `orthofusing` pulls final sequence out of
    via filter.py. Masking in place would edit the output assembly and
    invalidate a scoring run that takes twelve hours. Only the clustering
    sees an X.
  - **`check_orthofinder_searches` catches it next time, whatever the
    cause.** Two rules on the artifacts rather than on an exit status:
    every search produced at least one hit, and each `Blast{i}_i` found at
    least `BLAST_SELF_HIT_FLOOR` (0.5) of its own sequences, since every
    sequence aligns to itself. The floor is not 1.0 because diamond masks
    low-complexity before seeding and a fully masked sequence reports no
    self-hit -- ~13% of each SPAdes assembly here. 0.5 sits far below that
    legitimate shortfall and far above a failure: the worst surviving
    self-comparison in this run held about 0.2% of its input. Off-diagonal
    cells get no floor beyond one hit, because two assemblies can
    legitimately share very little; that misses `Blast1_3` and does not
    matter, since the diagonal already condemns both SPAdes species.
  - **`ORTHOFINDER_GB_PER_SEARCH` was left at 12** rather than fitted to
    the 145 GB measured. Fitting it would drop every run to 4 concurrent
    searches to insure against an input class that no longer reaches
    diamond. Its comment did have to go: it claimed the cap was
    "structurally incapable" of OOMing a node this size, which this run
    disproves, and it blamed an earlier OOM on snap's index on the other
    branch. Per-process cost is set by what is in the query and nothing
    there can see it.
  - **This probably also explains the 6x `run_orthofuser` regression**
    pinned on the OrthoFinder 2.5.2 -> 3.1.5 bump in the 2026-08-16 entry
    below. Same step, same poly-N seeding, not yet severe enough to cross
    the ceiling. Still not confirmed -- it would need a run of that sample
    with and without the masking.

## 2026-09-17

- **`clear_transrate_outdir` was deleting the two most expensive things in
  `-o` on every failed attempt: the BAM and `salmon/`.** It kept `logs/`
  and any directory with a `GenomeIndex` marker and deleted everything
  else, so a step that failed after mapping threw away 159 GB of snap
  output and a finished `quant.sf` -- the same two artifacts the 380C run
  proved were still on disk and still good when the retry fired. Three
  attempts at a failing step therefore paid for three mappings to reach the
  same failure.
  - **Keeping them unconditionally would have been worse than deleting
    them.** pytransrate 2.1.0 -- what ORP pinned when this was found --
    reuses both on existence alone: `if os.path.exists(self.bam): return
    self.bam`, and the same for `quant.sf`. So a half-written copy of
    either is picked straight back up. The BAM at least dies loudly on the
    missing BGZF EOF marker; a `quant.sf` cut at a line boundary does not
    fail at all, it scores the assembly off whichever contigs made it into
    the file, and those scores are what `pick_best_contigs.py` selects on.
    2.2.0 (now pinned) checks the BAM and **still** reuses `quant.sf`
    blind, so the `quant.sf` check stays load-bearing and the BAM check
    becomes a question of disk rather than of correctness -- keeping a
    truncated BAM is keeping hundreds of gigabytes nothing will read.
  - **`quant.sf` is only as good as the BAM it came from**, so `salmon/` is
    kept only beside a BAM that was kept, and then only when it holds one
    whole row per contig in the assembly. That is why the function now
    takes the assembly path -- it is the only thing that says how many rows
    a complete `quant.sf` has.
  - **Keep `*-read_count.txt` with the BAM.** Without it the reuse path
    falls to `_load_read_count`'s fallback, which counts lines in the fastq
    with a plain `open(path, "rb")` -- a full pass over the library, and
    nonsense if the reads are ever handed over gzipped.
  - **Every check fails towards deleting**, so being wrong about one costs
    a recomputation, never a wrong score. Worth keeping that property in
    mind before adding a fifth thing to spare.
  - **The pin moved to `v2.2.0` in the same change.** Checked first, not
    assumed: `score.py` and `output.py` are byte-identical across
    2.1.0..2.2.0, so the formulas and the CSV column order ORP reads
    positionally are untouched; the CLI is a strict superset, so nothing
    either call site passes was dropped; `requires-python` and the three
    runtime deps are unchanged, so only the tag in `orp_env.yml` moved.
    Scores still move, through snap's input rather than the scorer --
    contig padding 2000 -> 1000 on every assembly, and pipes no longer cut
    out of deflines, which is inert for oyster.py and a real fix for
    chowder on ENA/TSA input.
  - **What is still ahead of the pin.** The commits that fix pytransrate's
    own half of this -- `.align.done`, partial BAMs moved to `.partial`,
    BAM deletion deferred to the end of a successful run -- are on
    pytransrate master and **not pushed** (`master...origin/master [ahead
    2]` as of 2026-09-17), so there is nothing to pin to. Push and tag
    them, then this is worth revisiting. `clear_transrate_outdir` already
    keeps `.align.done` beside a kept BAM so the move needs no ORP change:
    drop it and that pytransrate would move a good BAM aside and map
    again.

## 2026-09-16

- **The snap SIGFPE is amplab/snap#171: an upstream bug, diagnosed and
  fixed fourteen months ago, never released.** The backtrace settled in one
  run what the log structurally could not. `AffineGapVectorized.cpp:351`,
  with `patternLen = 0`: `numVec = (patternLen + 7) / 8` is then 0 and
  `(patternLen - 1) % numVec` divides by zero. It is on the CIGAR path at
  *write* time, not in the aligner -- `SimpleReadWriter::writePairs` ->
  `BAMFormat::writePairs` -> `computeCigarOps` -> `SAMFormat::computeCigar`
  -> `computeGlobalScoreNormalized`. Issue #171 was opened 2025-06-06 by
  someone running transrate's Ruby flags against a Trinity assembly: same
  provenance as ours, reached independently.
  - **Bolosky's diagnosis (2025-07-21)**: the read must have multiple
    alignments, one secondary alignment must back-clip it to exactly half
    its length, and the writer must then run out of output buffer partway
    through writing it. The flush-and-retry retains the back clipping, which
    clips the *original* alignment to zero bases. His own note that
    "different lengths would result in incorrect alignments but no crash"
    is the part worth remembering: wherever this did not divide by zero it
    was quietly writing wrong alignments.
  - **The fix is `0e0997b` (2.0.6.dev.2) on `dev`, and it is two functional
    lines** -- a `setAdditionalBackClipping(0)` reset at the top of each
    alignment loop in `SimpleReadWriter::writePairs`. It has never been
    released: master's last code change is the v2.0.5 release of May 2025,
    there is no 2.0.6 tag, and so bioconda's latest is still 2.0.5. It is
    unreleased only because the reporter said they would test it and went
    quiet, which makes confirming it the cheapest lever anyone has on this.
  - **The fix appears to carry a typo.** Its rewritten loop indexes
    `reads[whichRead]` but leaves the subscript on
    `clippingForReadAdjustment` at `[0]`, where 2.0.5 used `[0]` for read 0
    and `[1]` for read 1 -- so read 1 takes read 0's front-clipping
    adjustment. Raised as a question on the issue, not asserted.
  - **`--max-parallel 1` worked; this was never the same failure.** The OOM
    is fixed -- the run mapped for 17m44s and passed 100 GB of BAM before it
    died. SIGFPE is not a memory symptom: x86-64 masks FP exceptions, so
    signal 8 is always integer divide-by-zero or `idiv` overflow, and an
    OOM would be SIGKILL.
  - **The scale variable is the size of the BAM, not the size of the
    assembly.** A full library against this merge produces ~160 GB of
    output, so the writer refills its buffer more or less continuously and a
    per-read-rare event becomes a certainty. Same pipeline, same flags, a
    smaller assembly: clean. It also explains why the crash reproduces every
    time at roughly the same elapsed point without being deterministic in
    *which* read it kills.
  - **MULTI_ALIGNMENT_SETTINGS' diagnosis was right and stays right; its
    remedy is wrong.** "The fault is in the multiple-alignment path" is
    exactly correct -- the bug requires a secondary alignment -- so last
    session's bisect was measuring the real mechanism. But "try lowering
    --max-alignments-per-pair and --max-seed-hits", which mapper.py's SIGFPE
    message also says, does not work: we are already at `-H 4000 -D 2 -om 2`
    against the Ruby's `-H 300000 -D 5 -om 5` and it still dies. Only
    dropping `-om`/`-omax` outright would help, and that is not reachable
    from the CLI. Both need correcting.
  - **Three hypotheses from the session were wrong**, and are recorded so
    they are not re-derived. (1) The >4 Gbase / `locationSize 5` story:
    CIGAR generation does not care about location size, and
    `doesGenomeIndexHave64BitLocations` gates nothing on this path. (2) "The
    alignment start lands in the inter-contig padding" is a genuine *second*
    route to `patternLen == 0` -- on the hangs-off-the-end branch the
    algebra reduces to `patternLen = (contig end) - (alignment start)` --
    but it is not the one that fires; ours reaches zero through `dataLength`
    itself. (3) Contig count is not the exposure driver.
  - **`--padding` is irrelevant here**, which closes that question a second
    time and by a simpler argument than 507ab63's measurement: nothing on
    the retained-back-clipping path touches padding at all.
  - **The retry loop costs ~35 minutes an attempt for nothing** -- 14.5 min
    of contig metrics, then ~17 min to the identical SIGFPE, three times
    over. pytransrate exits 1 for everything (`cli.py:573`), so `run()`
    cannot tell a deterministic crash from a transient failure. A distinct
    exit code for "snap died on a signal", with `retries=0` against it,
    pays for itself on the first failed merge.

- **Taking the patched snap costs this run nothing, and changes
  pytransrate's dependency story for everyone else.** pytransrate resolves
  the aligner with `which("snap-aligner")` (`mapper.py:445`) and there is no
  version gate anywhere -- `2.0.5` appears only in comments and docstrings,
  and oyster.py's preflight checks presence, not version. So a patched
  binary ahead of the env's copy on PATH is a drop-in: no `orp_env.yml`
  change, no pytransrate change, no chowder change. (Overwriting in place
  inside the env works but conda owns that file and a later install can
  revert it silently.)
  - **The whole delta from v2.0.5 to the fix is the fix.** Checked rather
    than assumed: `git diff v2.0.5 0e0997b` is ReadWriter.cpp, the version
    string, one help-text typo and one longer error message in
    `validateCigarString` -- 2.0.4.dev.1 and 2.0.5.dev.1 predate the v2.0.5
    tag and are already in it. ReadWriter.cpp is untouched by anything else
    since v2.0.5, so `git checkout 0e0997b -- SNAPLib/ReadWriter.cpp` onto
    v2.0.5 is exactly the fix and conflicts with nothing; building `dev` tip
    is equivalent. Everything that was not crashing comes out identical
    either way.
  - **The cost that cannot be packaged away**: pytransrate's requirement
    becomes "snap with #171 fixed", and no public release satisfies it.
    Anyone installing standalone takes 2.0.5 from bioconda and hits this on
    a large enough assembly. A lab conda package of 2.0.6.dev.2 would keep
    `conda env create` declarative and is the sane stopgap, but it is a
    binary we own until upstream tags 2.0.6. The real fix is the tag.
  - **`-G-` (disable affine gap) was considered and rejected.** It routes
    CIGAR generation to Landau-Vishkin, which has no such division, and
    needs no rebuild -- but it is not reachable from the CLI today, and it
    is assembly-changing: `contigs.csv` column 9 from orthotransrate is what
    `pick_best_contigs.py` uses to choose each orthogroup's representative,
    so different alignments mean a different `.ORP.fasta`.
    `--pytransrate-args` reaches both pytransrate call sites, so it cannot
    be confined to the report either. And with no pytransrate version stamp,
    adding it to an existing run directory would silently pair affine-gap
    picks with non-affine-gap scores -- the case the 2.1.0 entry below said
    to revisit when a release finally moves scores.
  - **The patched binary clears the reproduction.** v2.0.5 with
    `0e0997b`'s ReadWriter.cpp and the subscript restored to `[whichRead]`,
    built and installed over the env's copy as `2.0.5+snap171`, run against
    the standalone snap command that died at 17m44s on stock 2.0.5: past it.
    The full chowder run is going out on this binary, so this run's
    `.ORP.fasta` and its transrate score are the first produced with a
    patched aligner -- the `+snap171` banner in `merged/logs/snap.log` is
    the only trace of that, which is why the version string was changed.
  - **Outstanding**: post the confirmation on #171; correct mapper.py's
    SIGFPE message and MULTI_ALIGNMENT_SETTINGS' remedy paragraph; decide on
    packaging.

- **The chowder banner landed halfway down the 380C log, and the cause was
  stdout buffering, not print order.** 5251238 already calls `welcome()`
  before `check()`, and it does; what was wrong is that Python
  block-buffers stdout in 4-8 KB chunks whenever it is not a terminal --
  which on a cluster it never is -- while every tool we launch inherits the
  same descriptor and writes to it directly, unbuffered. So the whole
  parent-side narrative sat in our buffer while hours of OrthoFinder output
  streamed past it.
  - **The timestamps in the log prove it**: `run_filtershort` records
    `15:53:16`, OrthoFinder's first line is `15:54:51`, and yet every
    chowder line -- banner, `=== step -- start ===`, the `+ conda run`
    echoes, the retry warnings -- appears *after* OrthoFinder's. Not
    cosmetic: it makes a log read as though steps ran in an order they did
    not, and it puts a failure's explanation somewhere other than next to
    the failure.
  - `line_buffer_stdio()` in oyster.py, called first thing in both entry
    points' `main()`. Line buffering also guarantees we have flushed before
    a child we spawn writes anything, so the interleaving is correct rather
    than merely closer. **Not `sys.stdout.reconfigure()`** -- 3.7+, and the
    cluster is on 3.6.8; the fallback wraps `detach()` rather than `.buffer`
    so the discarded wrapper cannot close the descriptor out from under the
    replacement when it is collected. Both paths verified against a parent
    that prints and then spawns children, redirected to a file.

- **The 380C_0C5D_001F chowder run died OOM; four oyster.py bugs came out of
  the post-mortem, and only one of them is the OOM.** `sacct` on the job:
  `OUT_OF_MEMORY`, `MaxRSS 751,425,356K` = **716.6 GiB against a 720 GiB
  ReqMem**, 99.5% of the wall. Not close -- into it.
  - **The memory is snap's, not diamond's, and 494e0a8 tuned the wrong half
    of the node.** That commit's premise -- each diamond sizing its block
    against whatever memory looks free when it starts -- is simply not what
    diamond does. Its default block size is a fixed `-b2.0` and the manual's
    own rule is "roughly six times this number of memory (in GB)", so ~12 GB
    per process whatever the node. (`--more-sensitive` does not change it;
    only `--very-sensitive`/`--ultra-sensitive` do, to `-b0.4`.)
  - **And concurrency is capped by n_assemblies^2 before it is capped by
    anything in oyster.py.** Four assemblies = 16 searches total (the log
    says `16/16`), so `-t 20` could never have put more than 16 diamonds on
    that node: **16 x 12 = ~192 GB**, comfortably inside the orthofuser
    branch's 335 G budget. Which leaves **~525 GiB for the merge branch** --
    73% of the whole node, against the same 335 G budget, and
    `merge_branch` does not so much as pass `mem` to `orthotransrate`.
    `ORTHOFINDER_GB_PER_SEARCH` corrected 8 -> 12 so the constant is at
    least diamond's real number, but on a node over ~200 GB this cap is
    structurally incapable of being what OOMs a four-assembly run. It is
    not the lever and cannot be made into one.
  - **Diamond died first but is the victim, not the cause.** Killing a 12 GB
    diamond frees 12 GB, so the kernel has to keep going -- hence 27
    oom-kill events, picking off what it could actually reclaim while the
    process actually holding the memory ran on for another 8.4 hours.
    `ERROR: Blast1_1.txt is corrupted` is a red herring: that "offending
    line" is well-formed BLAST6, it is the last line read before the
    *truncated* file gave out.
  - From log ordering the peak looks like the snap **index build**
    (16:11->16:49), not the alignment -- every diamond ERROR block lands
    immediately after the "index built" line, and the 8.4 h alignment that
    followed drew no further kills. Inference from interleaving, not a
    measurement.
  - **Measured off the dead run's BAM header, and it closes the `--padding`
    question in the negative.** `5,354,958` contigs, `5,664,049,229` real
    bases, padding `5,354,959,000` = **1000 per contig, not the 2000
    2ac5d8d assumed**. Mean contig 1058 bp. (The same header reported
    `EOF marker is absent`, so the BAM was indeed truncated -- restarting
    rather than resuming was right on the facts.)
    - **`--location-size 5` is mandatory here, not a tuning choice.** Real
      sequence *alone* is 5.664 Gbp = **1.32x the 4.295 Gbp four-byte
      ceiling**, so even `--padding 0` cannot get the genome back under it.
      That was the main prize and it is not available.
    - **Padding is worth ~1% and nothing else.** It is written as `N` and
      seeds containing `N` are skipped, so it grows the 1 byte/base genome
      array and touches nothing else: `--padding 0` saves 5.35 GB out of a
      ~525 GiB branch. What does not shrink is the location table --
      ~5.55e9 seed positions x 5 B = ~27.7 GB -- because that is set by
      real sequence. **So drop `--padding`; 2ac5d8d's premise that it
      "lowers what snap counts as genome" as a memory lever does not
      survive the measurement, and its help text is corrected here.**
  - **Still unaccounted for: genome array (11 GB) + location table (28 GB)
    + hash overhead is ~50-70 GB for the finished index, against a merge
    branch measured at ~525 GiB.** Do not invent a story for the gap. The
    honest candidates are the index *build* (sorting 5.55e9 seed entries
    needs several times the finished index) and Slurm `MaxRSS` under cgroup
    v1 counting page cache -- this job wrote a 159 GB BAM plus a ~50 GB
    index. The next run settles it: `merged/logs/snap.log` now survives,
    and snap reports its own index and alignment memory there.
  - **The lever that is left is the size of the merge itself.** 5.35M
    contigs at 1058 bp mean is ~1.34M contigs / 1.42 Gbp *per input
    assembly* -- genome-scale for a transcriptome, and the reason snap is
    at the edge of what it can index. The knobs that would actually move it
    are `long.seq.py`'s threshold (currently 200 bp; short contigs dominate
    the count, which drives both padding total and per-contig overhead),
    pre-filtering inputs by expression, or merging fewer assemblies. All
    three change the science, so they are the user's call -- but that is
    where the order of magnitude is.
  - **Immediate mitigation regardless: `--max-parallel 1`.** Serialising the
    branches makes the peak the larger of them alone rather than their sum:
    orthofuser ~192 GB (528 GiB headroom), merge ~525 GiB (195 GiB
    headroom). Both fit on this node; together they do not. Costs the
    overlap -- roughly 19-20 h instead of ~10 h -- and note it is a fix for
    *this* assembly, not a general one: at ~525 GiB the merge branch alone
    is already at 73% of a 720 GiB node, so a somewhat larger merge OOMs
    with nothing to serialise against.
  - **chowder on a smaller dataset ran clean start to finish.** So none of
    this is a code-path problem; it is purely resource sizing at scale, and
    the scale variable is the merged assembly snap has to index.

- **`unlink(missing_ok=True)` in `clear_transrate_outdir` -- 3.8+, and the
  cluster launches oyster.py under the system python3, 3.6.8.** It fires
  only on the retry path, so it converts a retryable failure into a crash
  whose traceback buries the failure that caused it. That is what the
  `TypeError` at the bottom of the 380C log is.
  - **It also destroyed the evidence, and the run directory proves it.**
    `iterdir()` reached `logs/` (a directory with no `GenomeIndex` marker)
    and rmtree'd it, then hit the first *file* and raised -- so zero files
    were deleted. Confirmed against the real directory: `logs/snap.log`
    gone, snap index and the 159 GB BAM still there. snap.log was the only
    thing that would have said why snap took SIGFPE, and it is
    unrecoverable.
  - **837e23f claimed logs/ was already preserved; it was not.** The commit
    message says "It took merged/logs/snap.log with it too ... the evidence
    was destroyed by the retry that needed it", but the code only ever
    spared directories carrying a `GenomeIndex` marker. `logs/` is now
    spared by name.
  - **Next resume on the unfixed cluster code dies in seconds.** On the
    failed run the pre-step clear returned early through
    `if not outdir.is_dir()`. `merged/` now exists and holds files, so that
    same call -- the one *outside* the retry loop -- raises immediately.
    Deploying the fix is a precondition for resuming, not an improvement.

- **OrthoFinder exits 0 after printing its own fatal errors**, so
  `orthofuser.done` was written for a run that produced nothing.
  `needs_run` would then skip the 10-hour step on every later resume and
  hand `makeorthout` either nothing or a stale `Orthogroups.txt`, and the
  run would go on to build a final assembly off an orthogroup set that was
  never computed -- silently. `run_orthofuser` now drops a marker before
  launching and requires a non-empty `Orthogroups.txt` newer than it, so the
  sentinel comes from the artifact rather than from the exit status.
  - `find_orthogroups_txt` took `rglob`'s first match. OrthoFinder never
    reuses a results directory -- it makes a fresh `Results_<Mon><Day>`,
    then `_1`, `_2` -- so after a failed attempt there are several and the
    order is the filesystem's. Now newest by mtime.
  - **The existing `orthofuser.done` in the 380C directory is already
    poisoned.** The gate protects future runs, not that directory: delete
    it and `Results_Sep15` by hand before resuming.

## 2026-09-09

- **Branched `byo-assemblies` off `pytransrate` and added `chowder.py`:
  bring your own assemblies, merge them with the ORP.** Decision was
  shared-engine-plus-thin-entry-point over either a `--skip-assembly` flag
  on `oyster.py` or a standalone copy. The merge half is where every
  assembly-changing subtlety lives, so two copies of it would drift
  invisibly; but half of `oyster.py`'s flags (k-mers, `--strand`,
  `--normalize-reads`) are meaningless without assemblers, so one CLI would
  have been half-inert. `Chowder(Pipeline)` overriding `main()` gets both.
  - **The refactor that made it possible was the risky part, so it was
    measured, not eyeballed.** Two harnesses in the scratchpad: one dumps
    every assembly-derived path, order and report line; the other traces the
    entire step graph with the tools stubbed out. Both were captured from
    the pre-refactor `oyster.py` first and diffed after. Result: the 40-step
    graph, every declared input and output, `posthack`'s cat order,
    `build_list5.py`'s priority order and the `qualreport` line text are all
    unchanged. Worth keeping those harnesses in mind for the next
    structural change -- this repo has no test suite at all.
  - **There were three different orders of the same four assemblers**, and
    that was the trap. Concatenation order (sp55, sp75, ta, trinity) reaches
    cd-hit-est through `merged.fasta` and `posthack`, where it breaks length
    ties; diamond order (ta, sp75, sp55, trinity) is a real preference
    ranking because `build_list5.py` keeps the first hit per gene; report
    order (trinity, sp55, sp75, ta) is cosmetic. A tidy-minded "let's just
    sort them" here would have quietly changed assemblies. They are now
    named constants with the reason attached.
    - A *fourth* order existed in `diamond_uniq`'s dict literal (trinity,
      sp75, sp55, ta) and was inert -- four independent reads, four
      independent writes. Collapsed onto report order.
  - **Contig-name prefixing is not cosmetic and is the one thing that would
    have silently corrupted a merge.** Two Trinity assemblies of one library
    both start at `TRINITY_DN0_c0_g1_i1`; OrthoFinder, `contigs.csv` and
    `filter.py` all join on contig name. Ingest renames to
    `<label>_<original>` and refuses an input with duplicate names inside
    it.
  - **Assembly order is shuffled, not chosen -- with a fixed seed.** MM's
    call, to stop the order of the command line being a scientific choice
    nobody meant to make. Implemented as sort-by-label then permute with
    `random.Random(23894)` (strandeval's existing seed, rather than a second
    arbitrary constant), so the order is a function of the *set* of
    assemblies: all 24 typing orders of four assemblies give one merge
    order, and colliding labels are numbered by source path so that holds
    even for two files both called `trinity.fasta`.
    - **Not a plain `random.shuffle`, and the difference is the point.**
      Unseeded, the same command would give a different assembly on a
      different day and a resumed run could disagree with the run it was
      resuming -- trading a decision nobody made for one nobody can
      reproduce. Seeded, the order is arbitrary but stable and recorded
      (printed at startup, written to `<run>.ingest.done` with the seed).
    - **It does not make the pipeline order-independent, and the docs say
      so.** The picks still depend on the order; what changed is that the
      order no longer depends on typing. Genuine independence means breaking
      cd-hit-est's ties and the rescue ranking on merit instead of position
      -- the rescue in particular has a real preference to express, since
      `build_list5.py` picking the first hit per gene is a chance to prefer
      a better assembly that a shuffle throws away. Worth revisiting if a
      seed sweep on real data shows the choice is worth anything: `--seed`
      exists precisely to measure that, and `--assembly-order given` keeps
      the old behaviour for anyone who wants to rank them by hand.
  - `--corrected-reads` (renamed from `--reads-are-corrected`) symlinks the user's pair into `rcorr/` rather
    than copying tens of GB. `cleanup()` now skips symlinks entirely --
    without that it would have reported someone's own reads as "left
    uncompressed" and, worse, been one edit away from unlinking them.
  - **Not yet run against real data.** Same standing item as pytransrate
    itself: needs the cluster. Fold a chowder run into the same trip --
    two of the four assemblies from a finished ORP run are the obvious
    input, since the merge of those should land near that run's own
    `.ORP.fasta` and gives a sanity check with a known answer.
  - Noticed in passing and left alone: `oyster.py` has never checked for
    `bwa`, which `strandeval` needs from the `orp_trinity` env. It gets away
    with it because the Trinity check proves that env exists. chowder's
    preflight checks `bwa` directly since it drops the Trinity check.

- **Bumped the pytransrate pin to `v2.1.0`** (`orp_env.yml`). Verified before
  bumping rather than after: the tag is on the remote at 53fe488, `cli.py`'s
  only diff v2.0.0..v2.1.0 is passing `threads=args.threads` through to
  `read_metrics`, and `tests/test_output.py` still pins `score`/`optimal_score`
  at indices 36/37 of `assemblies.csv` and the contig score at column 9. So no
  ORP-side code change: both call sites already pass `-t <cpu>`, which is now
  what divides the scoring step too.
  - **Deliberately did *not* add a `stamp_tool_version("orp", "pytransrate", ...)`**
    alongside the salmon one. The salmon stamp exists because 2.7.0 *rejects*
    an index an older salmon wrote -- a hard failure a resumed run walks into.
    A pytransrate upgrade has no such edge: the artifacts stay readable, and
    the only thing a stamp would buy is forcing a rescore. Across 2.0.0 ->
    2.1.0 that rescore would cost hours of `orthotransrate` to reproduce the
    same numbers to the fifteenth decimal. Revisit at the next pytransrate
    release that actually moves scores -- that one wants the stamp, and wants
    it added *with* the bump so an in-flight run directory is invalidated.
  - `p_seq_true` is the one number that moves (<=3.3e-15, and only because it
    was made exact and thread-independent). `contigs.csv` rounds to six
    decimals, so the realistic blast radius is nil; noting it so it is not
    mistaken for drift if a rescored run differs in the last digit.
  - The stale-docstring item from the entry below is closed upstream: 2.1.0
    carries "Correct compare_orthogroup_picks' account of where groups come
    from", and its `--pick-best` now imports this repo's `best_in_group` when
    the target is ORP 4.0.0+, rather than mirroring the rule.

- **Did #2: the `.groups` round-trip is deleted.** `makelist`/`makegroups`
  are gone, `scripts/pick_best_contigs.py` takes `Orthogroups.txt` instead of
  a directory of `*.groups`. Went with option B (keep the script, change its
  input) over folding it into `Pipeline`: the only thing A saved was one
  `conda run` activation, a second or two, against a round-trip worth minutes
  -- and keeping the picker runnable on its own matters for a step that makes
  a scientific choice you may want to re-run by hand.
  - **Ordering was the whole risk and it is preserved.** The old glob sorted
    filenames, so lexicographic (`1, 10, 100, 2, ...`), not numeric; that
    order reaches `cd-hit-est` through `good.<run>.list` and
    `orthomerged.fasta`, where it breaks length ties. `good.<run>.list` is
    byte-identical old-vs-new on synthetic sets n=1..1111 (ties, zero and
    negative scores, contigs missing from `contigs.csv`, duplicate rows,
    blank lines), and the test asserts its own data distinguishes
    lexicographic from numeric order so it would actually catch a regression.
  - `compare_orthogroup_picks.py` in pytransrate was never at risk -- its
    `--orthogroups` mode already rebuilt groups from `Orthogroups.txt`.
    **Its docstring is now stale though**: it says "makeorthout deletes the
    *.groups files", when ORP no longer writes them at all, and its
    `oyster.py:580`/`oyster.py:601` line references have drifted. Worth a
    small commit in that repo; `--groups DIR` is now dead in practice.
  - `makelist`'s `<run>.list` went too: written, declared as `makegroups`'s
    input, never opened by anything.

- **Trinity's `--full_cleanup` was costing a resumed run both phases (~35h);
  fixed with a sentinel outside the directory it deletes.** `cmds.ok` was
  Phase 1's declared output *and* Phase 2's declared input, and Phase 2
  deletes it -- so a resume re-ran Phase 1, which rewrote `cmds.ok` newer than
  the finished `.Trinity.fasta`, which dragged Phase 2 along with it. Both now
  hang off `assemblies/<run>.trinity.phase1.done`. Verified against the real
  `needs_run()` across six states (fresh, phase-1-only, Stage-B-done, the two
  pre-sentinel migration cases, and stale corrected reads).
  - `already_complete()` does **not** cover this: the window is a run that
    finished Stage B and then died later, which is a walltime kill near the
    end of a long run -- exactly when the job gets resubmitted.
  - `seed_trinity_phase1_sentinel()` copies the *mtime* of whatever proves
    Phase 1 ran instead of stamping `now`. Stamping `now` would make the
    sentinel newer than `.Trinity.fasta` and re-trigger the same 34h re-run it
    exists to prevent.
  - `cleanup()` removes the sentinel on purpose. Keeping it would leave a run
    directory where Phase 1 looks done but the plain `.Trinity.fasta` is gone
    (only the `.gz` remains), so a forced re-run would skip Phase 1 and send
    Phase 2 in without its checkpoints.

- **Runs now clean up after themselves (`cleanup`, `reclaim_trimmed_reads`,
  `compress_async`, `already_complete`).** A finished run keeps `reports/`,
  `.ORP.fasta`, the four individual assemblies and the corrected read pair
  (the last two gzipped) and reclaims everything else. Written and unit-smoke-
  tested against a synthetic run directory; **not yet exercised on a real
  run** -- fold it into the same full run that first exercises pytransrate.

- **The compression is deliberately decoupled from the deletion**, and that's
  the whole design. A file stops being *written* long before it stops being
  *read*: `c1`/`c2` feed every assembler and every alignment step through to
  `strandeval`, and the four assemblies are read again at `run_filtershort`,
  `diamond_*` and `posthack`. So `compress_async()` builds the `.gz` beside
  the original as soon as the producing step returns, on a two-worker
  background pool, and `cleanup()` at the very end only unlinks. Net effect:
  the gzip cost lands in parallel with an assembler instead of being added to
  the end of the run. Peak disk goes up by roughly the `.gz` size (~25% of the
  reads) for the duration, which the trimmed-read reclaim below more than
  covers.
  - `pigz=2.8` added to `orp_env.yml`, with a `shutil.which` -> env -> plain
    `gzip` fallback chain. Thread count is capped at `min(4, cpu//8)` on
    purpose: this is background work sharing a machine an assembler already
    owns, not a stage of its own.

- **Two resumability traps, both handled, both worth remembering.** Deleting
  intermediates interacts badly with `needs_run()`, which decides everything
  from output presence/mtime:
  1. Deleting the `TRIM_*.fastq` makes `run_trimmomatic` look permanently out
     of date. Fixed with a `rcorr/<run>.trim.done` sentinel that `main()`
     swaps in as the step's declared output whenever the corrected pair is
     present and current. Old working dirs (no sentinel, TRIM files still
     there) take the original branch and behave exactly as before.
  2. Deleting most steps' outputs makes a *re-invocation* of a finished run
     reassemble from scratch instead of no-opping -- the failure mode a
     resubmitted cluster job would hit. `already_complete()` short-circuits
     `main()` on `reports/<run>.cleanup.done` + an `.ORP.fasta` newer than the
     raw reads. It has to return *before* `timing_init()`, which truncates
     `reports/<run>.timing.log` unconditionally and would otherwise wipe the
     finished run's timing report on the way past.

- **Open items for the first real run:** confirm the corrected-read gzip
  actually finishes inside Stage A (it should -- Stage A is hours, and pigz on
  4 threads does tens of GB in minutes -- but the fallback is a plain `gzip`
  on an env without pigz, which is the case worth watching); confirm nothing
  in `reports/` turns out to depend on something `cleanup` removes (`reportgen`
  runs before it and reads `assemblies/diamond/*.unique.*`, `assemblies/
  <run>.flagstat` and the transrate CSV, all of which is why `cleanup` is
  ordered last); and record the reclaimed total from
  `reports/<run>.cleanup.done` in `sampledata/benchmarks.md` alongside the
  timing, since "how much disk does a run leave behind" is now a number worth
  tracking.

## 2026-09-07

- **Opened branch `pytransrate` to swap the bundled Ruby orp-transrate for
  [pytransrate](https://github.com/macmanes-lab/pytransrate) v2.0.0** (local
  checkout at `~/transrate`, clean, tagged, pushed). **Swap is done in code**
  -- both call sites, the preflight check, the env, and the whole install
  surface. What follows is the survey it was done from; the open items are
  collected at the end. Nothing has been run yet.

- **The CLI is drop-in and the CSV column contract is preserved on purpose.**
  pytransrate keeps `-a/--assembly`, `-o/--output`, `-t/--threads`,
  `--left`, `--right`, and its `tests/test_output.py` pins the exact indices
  ORP reads: `score`/`optimal_score` at 36/37 of `assemblies.csv`
  ([oyster.py:949](oyster.py#L949)) and contig score at column 9 of
  `contigs.csv` ([scripts/pick_best_contigs.py](scripts/pick_best_contigs.py)).
  So neither `reportgen` nor the picker needs touching.
  - One layout difference, and it lands safely: with a single `-a`,
    pytransrate writes both CSVs directly into `-o`, where the Ruby put
    `contigs.csv` in a per-assembly subdirectory. The two consumers use
    `rglob` ([oyster.py:593](oyster.py#L593),
    [oyster.py:947](oyster.py#L947)) so either layout resolves. The one
    place that hardcodes a path is the step's declared output for
    resumability, `reports/transrate_<run>/assemblies.csv`
    ([oyster.py:1064](oyster.py#L1064)) -- that assumes `assemblies.csv`
    sits at the top of `-o`, which both implementations do. Worth
    re-checking on the first run rather than trusting it, since a miss here
    silently re-runs the step forever instead of erroring.
  - `--reference` is parsed but unimplemented; ORP never passes it.

- **Call sites to change (2):** [oyster.py:580](oyster.py#L580)
  `orthotransrate()` and [oyster.py:839](oyster.py#L839) `transrate()`. Each
  swaps `makedir/software/orp-transrate/transrate` for the `pytransrate`
  console script. The `rglob("*.bam")` unlink loops after both become
  no-ops -- pytransrate deletes the BAM on success unless `--keep-bam` --
  harmless to leave, cleaner to drop.

- **Install check to change (1):** [oyster.py:331](oyster.py#L331) tests
  `os.access` on the unpacked binary. Becomes a `which_in_env(<env>,
  "pytransrate")` call like every other tool above it.

- **Env: superseded -- see the 2.7.0 bump below. Original reasoning kept
  because the constraint it names is real.**
  pytransrate needs `snap-aligner=2.0.5` and `salmon=2.7.0`, but the `orp`
  env pins `salmon=2.5.1` ([orp_env.yml](orp_env.yml)) and oyster.py runs the
  real quantification against it ([oyster.py:737](oyster.py#L737),
  [oyster.py:746](oyster.py#L746)). Bumping salmon there would change
  `quant.sf` for reasons that have nothing to do with this swap. A separate
  env matches the existing per-tool pattern (`orp_spades`, `orp_trinity`,
  `orp_busco`, `orp_transabyss`, `orp_orthofinder`).
  - Proposed Makefile line, alongside the others in the `orp:` target:
    `mamba create -y -c bioconda -c conda-forge --override-channels --name
    orp_transrate python=3.11 numpy scipy pysam snap-aligner=2.0.5
    salmon=2.7.0 pip`, then
    `pip install git+https://github.com/macmanes-lab/pytransrate.git@v2.0.0`
    into it. Pin the tag -- scores move between versions.

- **Decision: bump the `orp` env to salmon 2.7.0 and run pytransrate out of
  that same env, rather than building a separate `orp_transrate`.** Checked
  every ORP salmon flag against the 2.x migration notes first; consequences
  below. `orp_env.yml` now carries `salmon=2.7.0`, plus `snap-aligner=2.0.5`
  and `pysam` for pytransrate. `orp_trinity`'s own `salmon=1.10.3` is a
  separate env and is untouched.

- **No ORP salmon invocation breaks on 2.x.** Both call sites were checked
  option by option:
  - `salmon index` ([oyster.py:732](oyster.py#L732)): `-t`, `-i`, `-k 31`,
    `--threads` all carried forward unchanged.
  - `salmon quant` ([oyster.py:741](oyster.py#L741)): `-p`, `-i`,
    `--seqBias`, `--gcBias`, `--libType A`, `-1/-2`, `-o` all unchanged.
  - `--validateMappings` is **ignored** in 2.x -- selective alignment is the
    default now -- so it was a no-op that read as though it still switched
    something on. Dropped.
  - `--no-version-check` is also a silent no-op (2.x never contacts the
    network). Kept: harmless, and still meaningful if the env ever resolves
    to a 1.x salmon.
  - `quant.sf` columns are unchanged (Name, Length, EffectiveLength, TPM,
    NumReads), so `filter_tpm`'s `cols[0]`/`cols[3]`
    ([oyster.py:762](oyster.py#L762)) is safe, as is pytransrate's own
    5-column assertion.

- **The one real operational hazard: 2.7.0 requires index format v2 and
  rejects every older index on load.** ORP's resumability is mtime-based
  ([oyster.py:221](oyster.py#L221)), and `salmon_index` declares
  `<run>.ortho.idx` as its output -- so **resuming a run whose index was
  built by 2.5.1 skips the rebuild and then fails in `salmon quant`**. It
  fails loudly rather than silently, but `run()` will burn its
  `STEP_RETRIES` attempts on it first. Anyone resuming an in-flight run
  across this upgrade must delete `quants/<run>.ortho.idx` by hand. Fresh
  runs are unaffected.
  - **Handled** (`stamp_tool_version`, [oyster.py:234](oyster.py#L234)).
    `quants/salmon.version` records the salmon version and is declared as an
    *input* to `salmon_index`, so an upgrade invalidates the index through
    the ordinary mtime path rather than a special case. The stamp is
    rewritten only when the version actually changes, or a rebuild would
    fire every run. `salmon_index` also `rmtree`s the index directory before
    rebuilding, since the rebuild is usually over one salmon has already
    refused. Checked against a scratch harness on all six cases: cold start,
    clean resume, the 2.5.1 -> 2.7.0 upgrade, the run after that rebuild, an
    unreadable version, and a regenerated intermediate fasta.

- **Quantification numbers move, and that is mostly a win.** 2.6.0 made
  deterministic quantification the default, so the same reads and assembly
  now give the same TPMs run to run -- ORP's salmon step stops being a
  source of run-to-run drift. 2.7.0 itself is byte-identical to 2.6.0, so
  all of the change lands in the 2.5.1 -> 2.6.0 step.
  - Downstream, TPM only reaches the assembly through `filter_tpm`, and
    **with the default `--tpm-filt 0` that path is inert**: `low` is written
    only when `tpm < 0`, which never happens, so LOWEXP stays empty and
    `secondfilter` copies the intermediate through unchanged
    ([oyster.py:765](oyster.py#L765)). Default runs get identical output.
    Only `--tpm-filt > 0` users see membership shift near the threshold.
  - 2.6.0 also stopped emitting decoys in `quant.sf`. ORP indexes its own
    assembly with no decoys, so no effect.

- **Still unverified, and only checkable on the cluster:** that
  `mamba env create -f orp_env.yml` actually solves with
  `salmon=2.7.0 + snap-aligner=2.0.5 + pysam` alongside the existing exact
  pins. No conda on the laptop.

- **Install surface to retire:** `transrate` var
  [Makefile:13](Makefile#L13), the `all` prerequisite
  [Makefile:20](Makefile#L20), the unpack target
  [Makefile:77](Makefile#L77), `postscript`
  [Makefile:84](Makefile#L84), `clean` [Makefile:100](Makefile#L100);
  `software/orp-transrate.tar.gz` itself; the PATH exports in
  [Dockerfile/Dockerfile:52](Dockerfile/Dockerfile#L52) and
  [Dockerfile/Dockerfile:54](Dockerfile/Dockerfile#L54); INSTALL.md steps 5
  and 8 plus the `make` summary at [INSTALL.md:20](INSTALL.md#L20); the two
  tool cells in [docs/pipeline-steps.md](docs/pipeline-steps.md). No PATH
  entry is needed at all now -- the console script lives in the env.

- **Scores will move, and that is expected, not a regression.** pytransrate's
  CHANGELOG documents it: the assembly score drops 0.008-0.070 across three
  assemblies of one library, driven by a genuine soft-clip coverage fix
  (SNAP 2.0 clips; the old `bam-read` advanced the reference cursor over
  clipped bases, shifting coverage rightward). Good-mapping *rate* agrees to
  +/-0.001, so the two implementations agree about the reads and disagree
  about the contigs.
  - **The part that reaches the assembly is ordering, not level.**
    14.5-19% of contig pairs order oppositely, so `makeorthout` will pick
    different representatives for some orthogroups and the final ORP.fasta
    will differ. Measure it rather than assume: pytransrate ships
    `scripts/compare_orthogroup_picks.py`, which runs ORP's own selection
    against two `contigs.csv` files and names the groups that change winner.
    Feed it the **orthotransrate** CSV over `merged.fasta`, not a run over
    the finished ORP.fasta.
  - Consequence for the record: the transrate/orthotransrate numbers in
    `sampledata/benchmarks.md` stop being comparable across this boundary.
    Needs a fresh baseline run, and a note in the benchmarks file marking
    where the evaluator changed.

- **Open, and all of it needs the cluster:**
  1. That `mamba env update -f orp_env.yml --prune` solves with
     `salmon=2.7.0`, `snap-aligner=2.0.5`, `pysam` and the pytransrate pip
     line against the existing exact pins. Nothing here is verifiable on the
     laptop -- no conda.
  2. A full run, which is the first real exercise of pytransrate under ORP.
     Watch two things specifically: that `reports/transrate_<run>/` really
     does get `assemblies.csv` at its top level, since
     [oyster.py:1096](oyster.py#L1096) hardcodes that path as the step's
     output and a miss there re-runs the step forever rather than erroring;
     and that `orthofuse/<run>/merged/contigs.csv` lands where `makeorthout`
     `rglob`s for it.
  3. `compare_orthogroup_picks.py` against the old and new orthotransrate
     `contigs.csv` for the same dataset, to put a number on how much of the
     assembly actually changes. Needs an old-run CSV kept aside before
     re-running anything.
  4. A fresh `sampledata/benchmarks.md` baseline -- the transrate numbers
     there are Ruby-era and no longer comparable.

- **Runtime is unknown and `STEP_TIME_HINTS["transrate"] = 16`
  ([oyster.py:48](oyster.py#L48)) is a Ruby-era measurement.** It only sets
  submission order against `strandeval` (2 min), so it is very unlikely to
  invert, but re-time it on the first full run.

## 2026-08-24

- **`--cpu 80` oversubscription run: net regression, 38:46:04 vs the
  `_955parallel` baseline's 37:06:19 (+1h39m45s, +4.5%).** Full writeup with
  per-step deltas in `sampledata/benchmarks.md` (2026-08-24 entry); the
  investigation doc has been updated to match.
  - Critical-path decomposition (from start timestamps, reconciles to the
    second): preprocessing **-4m28s**, Stage A **-13m41s**, Stage B
    **+1h42m56s**, post-Trinity tail **+14m58s**.
  - **Section 3.3 of the investigation doc is falsified.** Phase 2 throughput
    *fell* from 35.67 to 33.98 jobs/min when ParaFly went from 38 to 76 slots
    on 40 physical cores. Phase 2 is **CPU-bound**, not I/O-latency-bound, so
    the extra concurrency bought nothing but context-switching. Don't retry
    oversubscription anywhere.
  - **Section 3.4 (working dir off GPFS) downgraded** by the same evidence: a
    filesystem-starved Phase 2 would have *gained* from more in-flight jobs.
    Inference, not direct measurement -- a `%iowait` reading during any future
    Phase 2 would close it for free.
  - **Interacts with the 2026-08-22 closure entry below**, which listed three
    still-generalizable candidates: `--normalize_max_read_cov 50`,
    oversubscribing ParaFly's `-CPU`, and `--min_kmer_cov 2`. This run kills
    the second of those. Not reopening the `--grid_exec` decision -- it was
    ruled out on generalizability, which this run doesn't speak to. It does
    mean that with `--grid_exec` off the table, nothing order-of-magnitude
    remains on the list at all.
  - **Section 2d confirmed.** Phase 1 dropped 1:22:49 -> 1:09:07 (-16.5%) on
    double the slots with `--inchworm_cpu` still pinned at 10, so the gain is
    coming from the `-t $CPU` Chrysalis stages, as predicted. Nothing about
    inchworm changed.
  - The tax landed on every step already saturating 40 cores: orthotransrate
    +40%, orthofusing +81%, busco +45%, transrate +26%, strandeval +42%.
    Trans-ABySS lost 32m48s but is off the critical path (done 19:02 vs Phase 2
    running to 01:29 next day), so it cost nothing.
  - **Next move, and it fits the 2026-08-22 generalizability bar better than
    anything else left: `--cpu 40` with `TRINITY_PHASE1_SHARE = 0.5`**
    ([oyster.py:61](oyster.py#L61)). One constant in our own code -- no Trinity
    flag, no cluster-specific setup, no change to assembly output, so no fresh
    BUSCO/TransRate pass needed. This run showed 52m38s of idle SPAdes headroom
    in Stage A and SPAdes flat between 30 and 60 slots, so the cores are there
    to move. Worth ~13 min. Not yet implemented; would ride along with any
    future run rather than needing one.
  - **Loose end:** Phase 1 ran at 20 threads instead of 10, and inchworm output
    is thread-count-dependent (doc 2c), so this run's Phase 2 may not have had
    the baseline's 73,737 jobs. `wc -l recursive_trinity.cmds` on that run's
    output dir settles whether part of the +5.0% is extra work rather than
    worse throughput. Doesn't rescue oversubscription either way.
  - **Also missing:** the run's actual command line wasn't captured. The
    benchmarks entry assumes it was `_955parallel`'s with `--cpu 40` swapped
    for `--cpu 80` and flags the assumption -- worth confirming from the run
    log, since if `--mem` moved too the Phase 2 comparison isn't clean.
## 2026-08-22

- Closed the Trinity Phase 1/2 speedup investigation
  ([docs/trinity-speedup-investigation.md](docs/trinity-speedup-investigation.md),
  entry below). Decision: not pursuing `--grid_exec`, the one option with
  order-of-magnitude upside. It only works given HpcGridRunner plus a
  scheduler-specific config (SLURM/SGE/PBS/LSF) tuned to a particular
  cluster's queues and fair-share policy -- a real win on Premise, but not
  something that generalizes into ORP for end users, who'd each need their
  own cluster-specific setup rather than a flag that just works. The
  remaining candidates (`--normalize_max_read_cov 50`, oversubscribing
  ParaFly's `-CPU`, `--min_kmer_cov 2`) are plain Trinity flags with no
  environment-specific setup and would generalize fine, but none is
  order-of-magnitude and two change the assembly output, requiring a fresh
  BUSCO/TransRate validation pass. Calling this optimization effort done
  for now rather than chasing diminishing returns.
- Tagged this as ORP 3.1.0 and pushed to GitHub.

## 2026-08-21

- Researched ways to speed up Trinity Phase 1 and Phase 2. Full writeup with
  source citations and line refs in
  [docs/trinity-speedup-investigation.md](docs/trinity-speedup-investigation.md)
  -- **research only, nothing implemented, nothing benchmarked** (no cluster
  access from this session). Headlines:
  - Phase 2 is **92.9%** of the `_955parallel` run (34:27:02 of 37:06:19);
    Phase 1 is 3.7%. Deleting Phase 1 entirely buys under 4%, so Stage A
    tuning (the 25/75 -> 75/25 question) is bounded at ~45min minus whatever
    SPAdes gives back. Not the lever.
  - **Raising `--inchworm_cpu` past 10 will not help.** Inchworm's kmer
    parsing is hard-capped at 6 threads in the binary
    (`Inchworm/src/IRKE_run.cpp:29`, only ever clamped *down*), and the part
    that does use all our threads (contig building, `IRKE.cpp:466`) is
    exactly the shared-hash-contention code Brian Haas cites when explaining
    the cap (trinityrnaseq issue #648). We're already at 10, past Trinity's
    documented default of 6.
  - Same code path **explains the open reproducibility question in the
    2026-08-16 (5) entry below**: `PARALLEL_IWORM` (on by default) walks an
    unsorted kmer list while threads zap kmers from a shared hash, so
    inchworm output is genuinely thread-count-dependent. The theory in that
    entry was right; this is the mechanism. It also means `--inchworm_cpu` is
    not a free knob -- changing it moves our quality metrics.
  - **Correction to the 2026-08-17 entry below**: it calls the 38:28:11 run
    the "no `--normalize-reads`" baseline, but `sampledata/benchmarks.md`
    shows that run's command line *did* pass the flag -- as did all four
    SRR1789336 runs. Combined with Trinity 2.15's `normalize_max_read_cov`
    default of **200x** (`Trinity:214`), which discards almost nothing at
    typical depth, this means **normalization has never actually been tested
    as a speed lever here**. A real test at 50x is the best cheap experiment
    available; oyster.py has no flag for it yet.
  - Biggest single lever if Premise policy allows it: **`--grid_exec`**.
    Phase 2 is 73,737 independent 1-core/1GB jobs currently squeezed through
    38 slots on one node. Adding the flag touches only `run_trinity_phase2`
    and cannot perturb the per-component command list (already written and
    `.ok`-checkpointed by Phase 1).
  - Also worth testing, both cheap: oversubscribing ParaFly past physical
    core count (jobs are 1GB each against 670GB), and getting the Trinity
    working dir off GPFS onto node-local disk (73,737 job dirs created and
    torn down = metadata storm). These two conflict with `--grid_exec`;
    `%iowait` during Phase 2 tells us which regime we're in.
  - Suggested first moves, none of which need a rerun: read the Phase 1
    per-stage split straight off the `.ok` file mtimes in a finished run, and
    read `%iowait` + the ParaFly rate during any live Phase 2. See section 4
    of the doc.

## 2026-08-20 (1)

- First full-scale (non-sample-CI) validation run of the 95/5 Stage
  A/Stage B restructure from entry, 2026-08-19 (2). Timing so far:
  - Stage A started 23:20 (`run_trinity_phase1` + the SPAdes lane
    together).
  - SPAdes lane (`run_spades55`/`75`, chained with their diamond
    searches) finished 23:40 -- **20min**.
  - `run_trinity_phase1` finished 00:40 -- **1h20m**, later than the
    SPAdes lane. Stage A therefore converges Phase-1-bound here, the
    opposite of `samplerun3`'s sample-CI result (SPAdes-bound there) --
    confirms the full-size prediction flagged in the `samplerun3`
    writeup (`sampledata/benchmarks.md`, 2026-08-19).
  - Stage B (`run_transabyss` + `run_trinity_phase2`) started 00:40, at
    the exact moment Stage A converged -- pairing confirmed working at
    real scale, same as `samplerun3`.
  - `run_transabyss` finished after **5h9m** (00:40 -> 05:49). Longer
    than the 4.5h Trans-ABySS took under the old 50/50 short-lane split
    (entry, 2026-08-19 (1), ~20 cores there) despite now running on only
    ~5% of `--cpu` -- consistent with that entry's finding that
    Trans-ABySS's dominant read-in stage is core-count-invariant, so the
    ~40min increase is coming from its smaller threaded sub-stages, not
    the bottleneck stage. Supports the 95/5 split's core assumption:
    squeezing Trans-ABySS's share doesn't blow up its wall time.
  - `run_trinity_phase2`/Butterfly running at **38.3 jobs/min** on 95%
    `--cpu`, vs. the 36.2/min measured at 100% `--cpu` (40 cores) in
    entry 2026-08-19 (2) -- slightly *faster* despite fewer cores
    available, most likely run-to-run noise in per-component job-size
    distribution rather than a real effect, but worth noting since it
    means the 95% share isn't visibly costing throughput here. Assuming
    the same 73,737-process total from that entry (same dataset), this
    revises the `T100` extrapolation down slightly: 73737/38.3 ≈ 32h6m,
    vs. the earlier 33h57m estimate -- projected Phase 2 finish ≈ 08:46
    on 2026-08-21 (from its 00:40 start).
  - **Direct progress comparison against the 50/50 baseline**, both
    read from ParaFly's `succeeded(N)` counter mid-run: the 50/50 run
    (`TIME2_SRR1789336_norm_py_5050parallel`) shows `succeeded(55682)
    75.5143% completed` at the 32.5h Phase-2 mark; the 95/5 run shows
    `succeeded(33244) 45.0919% completed` at the 17h mark. Both back out
    to the same 73,737-process total (55682/0.755143 ≈ 73737;
    33244/0.450919 ≈ 73737, matching entry 2026-08-19 (2)'s figure) --
    confirms same dataset/job count, so the two `%completed` figures are
    directly comparable.
    - Average rate so far: 50/50 = 55682/32.5h ≈ **28.6/min**; 95/5 =
      33244/17h ≈ **32.6/min**. Both are well below their respective
      early-run instantaneous readings (36.2/min for 50/50 in entry
      2026-08-19 (2), 38.3/min for 95/5 above) -- rate decelerates over
      the run in *both* designs, so that's a property of Butterfly's
      per-component job-size distribution (slower components running
      later), not an artifact of the 95/5 pairing.
    - Linear extrapolation from each run's own average-so-far:
      50/50 Phase 2 total ≈ 32.5h/0.755143 ≈ **43.0h**; 95/5 Phase 2
      total ≈ 17h/0.450919 ≈ **37.7h**. If these hold, the 95/5 design's
      Phase 2 alone finishes ~5.3h (~12%) faster than the 50/50 run's
      Phase 2 alone -- and that's *before* accounting for the 50/50
      design's extra ~5h20m of Phase 2 sitting idle while the short lane
      ran sequentially first (entry, 2026-08-19 (1)), which the 95/5
      design avoids entirely by starting Phase 2 concurrently with
      Trans-ABySS at 00:40. Still extrapolation, not a finished number --
      confirm against actual finish times once both complete.
  - **Correction, now that the 50/50 run has actually finished** (entry,
    2026-08-19 (1), results in `sampledata/benchmarks.md` 2026-08-21):
    real `run_trinity_phase2` was **35:57:48**, well under the 43.0h
    linear extrapolation above (off by ~7h, ~19% -- the checkpoint's
    28.6/min average-so-far undershot because the remaining 24.5% of
    jobs after the 32.5h mark actually ran much faster than the run's
    own average, ~87/min). So the "95/5 finishes Phase 2 ~12% faster"
    claim above doesn't hold up: 35.96h (50/50, actual) vs. 37.7h (95/5,
    still just a same-flawed-method extrapolation) -- if the 95/5 rate
    also picks up in its back half the way the 50/50 run's did, its real
    finish could easily undercut 37.7h too, but there's no way to know
    from a mid-run linear extrapolation alone. Treat both runs'
    mid-run % complete as directional only; wait for the 95/5 run's
    actual Phase 2 finish time before drawing any speed conclusion.
    What *does* still hold from the finished 50/50 run: it lost to the
    no-split baseline (38:28:11) by 4h00m25s overall, so the old 50/50
    design was a net regression on this dataset regardless of how the
    95/5 comparison lands -- see entry 2026-08-19 (1) for the full
    writeup.
  - **Run finished** (results, `sampledata/benchmarks.md` 2026-08-21):
    `run_trinity_phase2` took **34:27:02** -- close to the 50/50 run's
    actual Phase 2 (35:57:48), only 1h30m46s shorter despite running on
    95% instead of 100% `--cpu` for its whole duration, confirming
    Trans-ABySS's slack really was large enough to absorb the 5% cut
    without meaningfully costing Phase 2. TOTAL wallclock **37:06:19**
    -- the fastest of all four SRR1789336 designs tested, beating the
    no-split `_parallel` baseline (38:28:11) by 1h21m52s (~3.6%) and the
    50/50 split (42:28:36) by 5h22m17s (~12.7%). The win traces to Stage
    A converging in just 1h22m49s vs. the 50/50 design's ~5h sequential
    short lane, so Phase 2 starts ~3.7h earlier here -- not from Phase 2
    itself running meaningfully faster. **This closes out the
    validation**: the 95/5 Stage A/Stage B restructure is confirmed both
    correct (scheduling, all sample-CI and full-scale runs) and a real
    wall-clock win on the one full-scale dataset tested so far, unlike
    the 50/50 split it replaced. Still outstanding: quality-metric
    comparison (BUSCO/TransRate/gene counts) for this run wasn't
    reported yet, and only SRR1789336 has been tested at full scale --
    the original motivating dataset (entry, 2026-08-18 (1), where the
    short lane's idle extrapolated to ~110h) hasn't been re-run under
    the new design.

## 2026-08-19 (2)

- **Implemented** the cross-pairing restructure from entry (1) below,
  based on real numbers from the live `TIME2_SRR1789336_norm_py_5050parallel`
  run: Phase 2 processes at 36.2/min * 40 cores, 73,737 processes total ->
  `T100 = 73737/36.2 ~= 33.95h`. User's proposal (pair Phase 1 with
  SPAdes55/75 instead of the whole short lane, pair Phase 2 with
  Trans-ABySS instead of making it wait) sidesteps the SPAdes-starvation
  risk that made a blanket short-lane CPU reservation risky -- SPAdes gets
  a generous share where it doesn't matter (finishes in minutes either
  way), and only Trans-ABySS (structurally CPU-insensitive for its
  dominant cost, see entry (1)) gets squeezed. Pushed the ratio to 95/5
  (past the user's suggested 90/10, which only broke even) since
  Trans-ABySS's real slack is large enough to be safe.
  - `oyster.py`: replaced `TRINITY_LANE_SHARE=0.5` with
    `TRINITY_PHASE1_SHARE=0.25` (Stage A: `run_trinity_phase1` vs.
    `run_spades55`/`run_spades75`, sequential, each -> diamond) and
    `TRINITY_PHASE2_SHARE=0.95` (Stage B: `run_trinity_phase2` vs.
    `run_transabyss` -> `diamond_transabyss`). Trans-ABySS's mem share is
    *not* cut by `TRINITY_PHASE2_SHARE` -- kept at Stage A's `spades_mem`
    level instead, since its memory footprint doesn't shrink with fewer
    cores the way SPAdes/Phase 1's do. `run_filtershort` (needs all 4
    assemblies) now runs after both stages, unchanged position otherwise.
  - Added wall-clock timestamps throughout: `Pipeline._ts()` helper,
    `step()`/`run_parallel()`'s `run_one()` now print start+done
    timestamps to stdout and write them into the timing log
    (`name\telapsed\tstart_ts`, third column); `timing_report()` parses
    and displays them. Requested so live `squeue`-style monitoring of the
    new Stage A/Stage B overlap can be correlated against real wall time
    without cross-referencing separate elapsed-seconds math by hand.
  - Updated `docs/pipeline-schedule.html` (the DAG), `docs/pipeline-steps.md`,
    `README.md`'s Parallel task management section, and `changelog.md`'s
    pending 3.1.0 entry to match. **Not committed yet.**
  - **Scheduling confirmed correct** on a real run: `samplerun3`
    (`--cpu 20 --mem 100`, sample CI dataset, see
    `sampledata/benchmarks.md` 2026-08-19) shows `run_spades55` and
    `run_trinity_phase1` starting together, then `run_transabyss` and
    `run_trinity_phase2` starting together at the exact instant Stage A's
    `spades_lane` converges (15:38:25 for both), and `run_filtershort`
    waiting on Stage B's slower lane rather than Phase 2 alone. The
    Phase-2-waits-for-Trans-ABySS behavior from the old design is
    confirmed gone.
  - **Confirmed again** on `samplerun4` (`--cpu 40 --mem 600`, same
    sample dataset, different CPU/mem config, see
    `sampledata/benchmarks.md` 2026-08-19): same pairing behavior, and
    this time `run_filtershort`'s gate flipped to the *other* Stage B
    lane -- `run_trinity_phase2` (56s) outlasted `diamond_transabyss`
    this run, vs. `diamond_transabyss` being the slower one in
    `samplerun3`. `run_filtershort` correctly waited for whichever was
    slower both times, confirming the gate isn't hardcoded to one
    specific step.
  - **Still not validated at real scale.** This was the tiny sample
    dataset (seconds per step, `--cpu 20`), so it doesn't stress-test:
    (a) whether the actual node has enough free RAM for Trinity Phase 2's
    real, enforced `--max_memory` cap (95% of `--mem`) running
    concurrently with Trans-ABySS's own memory footprint, which
    oyster.py doesn't actually bound at all -- `run_transabyss()`'s `mem`
    parameter is accepted but unused, since Trans-ABySS's CLI has no
    memory flag to pass it to (only `--threads` is real);
    `transabyss_mem` in `main()` is therefore inert bookkeeping, not an
    enforced reservation, so the genuine open question is node headroom,
    not a percentage-sum bug in the code. (b) Trans-ABySS's threaded
    sub-stages don't become a real bottleneck at only ~2 cores (40 * 0.05
    = 2 at `--cpu 40`), (c) whether the total run actually beats,
    matches, or loses to the old 38:28:11 baseline on the SRR1789336
    dataset (see the caveat in entry (1) about this being the "wrong"
    dataset to validate the split's original motivation on). Needs a
    full-size run to confirm any of these.

## 2026-08-19 (1)

- Real-run validation (in progress) of the 50/50 Phase1/short-lane split
  from entry (1) on 2026-08-18: `TIME2_SRR1789336_norm_py_5050parallel`,
  `--cpu 40 --mem 670`. Trinity Phase 1 finished in **90 minutes**,
  confirming the "largely CPU-insensitive" claim -- Inchworm's hard
  `--inchworm_cpu 10` cap and Chrysalis's brief clustering mean Phase 1
  doesn't need its full 20-CPU share.
- **Trans-ABySS is not meaningfully CPU-scalable, and the reason is
  structural, not just algorithmic.** Traced `abyss-pe`'s binary dispatch
  (`bin/abyss-pe` in [bcgsc/abyss](https://github.com/bcgsc/abyss)): the
  initial De Bruijn graph build (`%-1.fa`, which includes reading the
  FASTQ files) only threads via `-j` if Bloom-filter mode (`-b`/`B`) is
  requested, or runs distributed via MPI if `np` is set (`mpirun -np
  $(np) ABYSS-P`). Neither `oyster.py` nor `transabyss`'s
  `dbg_assembly()`/`contig_assembly()` ever sets either, so both graph
  builds (stage 1 and, again, at the start of stage 3) fall to the plain
  `ABYSS` binary -- confirmed zero `omp`/`pthread`/`std::thread` anywhere
  in that binary or the `Assembly`/`Common`/`DataLayer` code it links
  against (including the FASTA/FASTQ reader itself). Also confirmed the
  middle `unitig_assembly` graph-simplification stage
  (`unbraid()`/`walk()`, pure single-process igraph) is single-threaded
  by design regardless of `--useblat` (which oyster.py doesn't set
  anyway).
  - Real numbers from this run: Trans-ABySS took **4.5h total**, of which
    the user measured **~3.2h as the single-threaded FASTQ read-in**
    stage specifically. Confirms the prediction: 20 CPU (this run) barely
    beat the old 8-CPU SRR1789336 benchmark's 4h46m (entry, 2026-08-17)
    -- most of the wall time is core-count-invariant.
  - Implication: tuning `TRINITY_LANE_SHARE` further (e.g. 25/75) has
    limited upside for the short lane specifically, since Trans-ABySS
    dominates it and can't use the extra cores for its dominant stage.
- SPAdes55/75 just started (sequential after Trans-ABySS in
  `short_assembler_lane()`); Trinity Phase 2 still gated on the short
  lane finishing entirely (`main()`'s `as_completed([...])` join before
  the `run_trinity_phase2` step) -- idle since Phase 1 finished at the
  90-minute mark, so ~3h+ of idle Trinity-lane capacity so far.
  - Considered decoupling `run_trinity_phase2` from the short lane (it
    only actually depends on Phase 1's `recursive_trinity.cmds.ok`
    checkpoint, not on Trans-ABySS/SPAdes output) so Phase 2 could start
    immediately and overlap with Trans-ABySS's mostly-idle-of-CPU tail.
    **Not free**: Trinity's `--CPU` is fixed at process launch for its
    entire run (same constraint behind the original Phase1/Phase2
    split), and Phase 2 is the dominant cost by a wide margin (~30h
    extrapolated from the old 37h/32-CPU benchmark at 40 CPU) -- so
    starting early at a reduced share trades a one-time ~3.5-4h idle
    window against a permanent CPU discount applied across Phase 2's
    entire runtime. Breakeven is roughly a 89%+ CPU share for Phase 2
    during the overlap; below that, early-start is a net loss. Decided
    to wait for this run's actual Phase 2 wall-clock time at 100% CPU
    (the missing data point) before committing to that restructure --
    too easy to make a multi-day job slower on a bad guess.
  - Update: short lane (Trans-ABySS + SPAdes55/75 + their diamond
    searches) finished at the **5h20m** mark; `run_trinity_phase2`
    started then. Confirms idle Trinity-lane window empirically: 90min
    (Phase 1 done) to 320min (Phase 2 start) = **3h50m idle**, matching
    the ~3.5-4h estimate used in the early-start breakeven math above.
    With that idle figure fixed, breakeven CPU share for an early-start
    Phase 2 is `T100 / (T100 + 3.83h)` -- e.g. ~88.7% if `T100` lands
    near the ~30h extrapolation, still needing the real number below.
- **Early-start restructure: decided against it.** User reported Phase 2
  processing 36.2 processes/min at 40 cores, 73,737 processes total ->
  `T100 = 73737/36.2 ≈ 2037min ≈ 33.95h`, higher than the ~30h
  extrapolation used above. That pushes the early-start breakeven share
  to `33.95/(33.95+3.83) ≈ 89.9%` -- realistic upside is ~2h off a ~38h
  run (~5%), only achievable by reserving as little as ~2 cores for the
  short lane during the overlap, which is tight enough that
  underestimating the short lane's real minimum (Trans-ABySS's threaded
  steps, SPAdes) could erase the gain. Given Trinity's CPU is locked in
  at launch for the full 34h+ duration, not worth the risk for a ~5%
  upside -- current wait-then-100%-CPU design stands.
  - Caveat worth revisiting: the original Phase1/Phase2 split (entry,
    2026-08-18) was motivated by a *different*, larger dataset where the
    short lane's idle extrapolated to ~110h against Trinity. On
    SRR1789336, Trinity already dominates so heavily (~34h vs.
    Trans-ABySS's 4.5h) that the old fixed 80/20 design's idle cost was
    small to begin with -- this run may land close to the old
    38:28:11 baseline (entry, 2026-08-17) rather than clearly beating
    it. Compare final totals once this run finishes to check whether the
    split actually helped on *this* dataset specifically.
- **Run finished** (results, `sampledata/benchmarks.md` 2026-08-21):
  `run_trinity_phase2` took **35:57:48**, close to the 33.95h estimate
  above (~5.6% under). TOTAL wallclock **42:28:36** -- confirms the
  worry two paragraphs up: this **lost** to the no-split `_parallel`
  baseline (38:28:11 from entry, 2026-08-17) by 4h00m25s, and only beat
  the `_NOparallel` baseline (43:39:49) by 1h11m13s. The 50/50 split, as
  designed, made SRR1789336 slower than not splitting at all -- the
  short lane's sequential-before-Phase-2 cost outweighed whatever
  concurrency benefit it bought elsewhere. Direct motivation for the
  Stage A/Stage B pairing restructure below (entry, 2026-08-19 (2))
  that lets Phase 2 start immediately instead of waiting on the short
  lane.

## 2026-08-18 (2)

- Real cluster run hit `strandeval()` failing with `Can't exec "samtools":
  No such file or directory` from `SAM_reader.pm` (part of
  `examine_strand.pl`, called via Trinity's own PerlLib). Root cause: that
  one call ran bare `perl` (`self.run(["perl", ...])`) instead of going
  through `conda_run()`/`conda run -n <env>` like literally every other
  subprocess in [oyster.py](oyster.py) -- confirmed by grepping every
  `self.run(` call in the file; this was the only one not wrapped (or
  wrapped only internally, like the `bash -c` bwa|samtools pipe a few
  lines above it). It inherited the bare process `PATH`, which per entry
  (5) on 2026-08-17 no longer has any conda env active (the `conda
  activate orp` line was deliberately dropped from the SLURM script, on
  the assumption every subprocess call already wrapped itself -- this one
  slipped through that audit).
  - Fix: `self.conda_run("orp_trinity", "perl", "-I", perllib, ...)`,
    matching the `bwa index` call right above it in the same function.
    Picked `orp_trinity` over `orp` (which also has samtools) because the
    `-I perllib` path is Trinity's own PerlLib, pulled from `orp_trinity`'s
    install (`trinity_perllib_dir()`) -- running it under a different env's
    perl risks a version/XS mismatch for no reason, and Trinity's own
    conda package almost certainly already depends on samtools internally
    (Trinity uses it in its own pipeline).
  - Not yet confirmed on a real run that `orp_trinity` actually has
    samtools available -- inferred from Trinity needing it internally, not
    observed. If this is wrong, the fallback is `orp` instead (confirmed
    to have both perl and samtools, since `run_rcorrector.pl` already runs
    successfully under `orp`).

## 2026-08-18 (1)

- Real-run evidence that the 80/20 lane split (entry (2) below, 2026-08-16)
  still leaves a lot of the machine idle even though it fixed the original
  "Trinity runs alone at half `--cpu`" bug: on the user's current dataset,
  the short-assembler lane (20% of `--cpu`) finished in ~15h, and at that
  point Trinity's Butterfly/Phase-2 stage was only ~12% complete. Linear
  extrapolation from that (`15h / 0.12 ≈ 125h` total) puts the idle window
  at ~110h, not the ~20h a naive read of the old ~37h `SRR1789336`-based
  benchmark (entry (4)/(3), 2026-08-16) would suggest -- that benchmark was
  from a much smaller dataset and never applied here.
- Root cause: Trinity's own architecture is two phases -- Phase 1
  (Inchworm + Chrysalis: builds the whole-transcriptome graph, partitions
  reads per gene component) then Phase 2 (thousands of small, independent,
  single-threaded per-component assembly jobs dispatched via ParaFly).
  Phase 2 is the dominant cost by a wide margin and scales close to
  linearly with core count, but the old `TRINITY_LANE_SHARE=0.8` capped it
  at 80% of `--cpu` for the *entire* run, including the long stretch after
  the short lane was done and its 20% was sitting unused.
- Fix, using Trinity's documented [multi-stage execution
  support](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity#running-trinity-in-multiple-sequential-stages)
  (confirmed against the actual `Trinity-v2.15.2` source, not just the
  wiki): split `run_trinity()` in [oyster.py](oyster.py) into
  `run_trinity_phase1()` (`--no_distributed_trinity_exec`, stops right
  after Phase 1) and `run_trinity_phase2()` (same command, no flag --
  Trinity resumes from Phase 1's on-disk checkpoints straight into Phase 2,
  per the docs). `main()`'s lane split now only covers Phase 1, at a new
  `TRINITY_LANE_SHARE=0.5` (was 0.8) -- Phase 1 is largely insensitive to
  its CPU share (Inchworm is hard-capped at `--inchworm_cpu` regardless,
  Chrysalis's clustering is brief next to Phase 2), so a 50/50 split gets
  the short lane through its own work faster without meaningfully slowing
  Trinity's prep. Once both lanes join, Phase 2 runs alone at the full
  `--cpu`/`--mem` budget instead of staying capped at a fixed share.
  Updated [changelog.md](changelog.md) (new `ORP Version 3.1.0 <- 3.0.0
  (pending release)` section, later renamed from an initial `Unreleased`
  heading at the user's request -- the `3.0.0` tag is already published at
  an earlier commit, so the tagged section itself wasn't touched) and
  [README.md](README.md)'s "Parallel task management" section to match.
- Not yet done: no real-run validation (no cluster access from this
  session). Next real run should confirm: (a) `recursive_trinity.cmds.ok`
  is actually the right checkpoint file to gate the phase1->phase2
  transition on (confirmed from reading the Trinity source, not observed
  in an actual output dir), (b) Phase 1 really is short/CPU-insensitive
  enough that 50/50 doesn't meaningfully slow it down, (c) Phase 2 at
  100% `--cpu` doesn't blow past `--mem` now that it's not sharing with the
  short lane.

## 2026-08-17 (5)

- **SPAdes python-3.6.8 bug (see (2) below) -- root cause found, fix
  confirmed.** Not a lane-split or `oyster.py` issue at all: the user's
  SLURM submission script did
  `source ~/.bashrc; module purge; module load anaconda/colsa; conda
  activate orp` before launching `oyster.py`. `conda activate orp`, run
  *after* `module load anaconda/colsa` has already altered `PATH`, was
  the thing breaking python resolution -- confirmed by reproducing it
  standalone (`conda activate orp_spades` then running `rnaspades.py`
  directly fails; `conda run -n orp_spades rnaspades.py` from an
  unactivated shell doesn't).
  - Audited every subprocess call in `oyster.py` (all ~20, including the
    ones not going through the `conda_run()` helper -- lines 283, 548,
    726, 749, 792, 817, 838, plus the `bwa | samtools` pipe in
    `strandeval()`): every single one already explicitly wraps its
    command in `conda run -n <env> ...` itself. Nothing in `oyster.py`
    ever assumes an env is pre-activated, so dropping `conda activate
    orp` from the job script costs nothing functionally.
  - Fix: remove the `conda activate orp` line from the SLURM script;
    keep `module purge`/`module load anaconda/colsa` as-is. Confirmed
    working: with the module loaded but no `conda activate`, `conda run
    -n orp_spades python --version` correctly reports 3.14.6 (the pin).
  - Residual risk noted but not yet hit: `conda_run()` invokes bare
    `"conda"` via `subprocess.run()` (no shell), so it does a raw `PATH`
    lookup for the `conda` *executable*, not the bash *function*
    `.bashrc`'s `conda.sh` defines -- `module load anaconda/colsa` could
    in principle still shadow ORP's own private conda independent of the
    `activate` line. Not observed in practice here (`orp_spades`
    resolved correctly), but if this bug resurfaces on a differently
    configured node, check `which conda` and `conda run -n orp_spades
    python --version` first.
  - This closes the "still open" status from entry (2) below -- no code
    change needed in this repo, purely a job-script fix on the user's
    end.

## 2026-08-17 (4)

- Pre-3.0.0-release cleanup, per user request:
  - Folded changelog's "Unreleased" section (two-lane split, `STEP_TIME_HINTS`
    fix, OrthoFinder 3.1.5, `orp_orthofinder` preflight check, dead
    `orthofuser.py` Makefile cruft removal) into the `ORP Version 3.0.0 <-
    2.4.0` section -- `version.txt` already said 3.0.0, so these were the
    actual 3.0.0 changes, not a future release.
  - Removed `scripts/for_loop.sh`, `scripts/numbers.sh`, `scripts/tpm.sh` --
    unreferenced anywhere (pipeline, docs, changelog) and predate the
    current pipeline entirely (`orp.mk`, a `shannon` assembler, `SAMP=`,
    `eukaryota_odb9` -- none of which exist anymore).
  - Removed `oyster.mk` -- fully superseded by `oyster.py` per changelog,
    not referenced by README/INSTALL.
  - Did NOT touch: the still-open SPAdes python-3.6.8 investigation (see
    entry below), or the `clean`/`preprocess`/`update_merge` functionality
    gap noted in the oyster.py-vs-2.3.3 comparison two entries up -- both
    are still open decisions, not yet acted on.

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
  `--max-parallel 2` baseline in the entry below. **Correction
  (2026-08-21):** this entry originally described that baseline as "no
  `--normalize-reads`". That is wrong -- `sampledata/benchmarks.md` shows
  both runs passed the flag, so `--normalize-reads` is not a variable
  between them at all. See the 2026-08-21 entry above.
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
