# Trinity phase 2: why it took 34 hours, and engineering fixes

> The changes now live in a fork: https://github.com/macmanes-lab/trinityrnaseq_2026
> (`master`, with Chrysalis and Butterfly on branches `chrysalis-2026` and
> `butterfly-2026`; write-up in its `PERFORMANCE_2026.md`). A built, verified
> checkout is at `~/trinity_eng/fork_check` on Premise. This directory is the
> working record; the patches here predate two build fixes made in the fork.

Work done 2026-10-05 on Premise against a fresh clone of
trinityrnaseq/trinityrnaseq (master 653e43a = v2.15.2 + 1 docs commit), using
the existing SRR1789336 phase-1 output in
`~/assemblies/TIME2_SRR1789336_norm_py_5050parallel.trinity/` (73,737
`recursive_trinity.cmds`, the read partitions, and that run's per-component
outputs as the reference). Working tree on Premise: `~/trinity_eng/`.

## Summary

| phase 2, SRR1789336, all 73,737 components, one 40-core node | wall |
| --- | --- |
| production (`orp_trinity`, bioconda build `_6`), from benchmarks | 34h27m |
| same conda Trinity, same launch, ParaFly built with OpenMP | **2h03m** |
| Trinity + Chrysalis + Butterfly patches below | **58m** |

1. **Phase 2 ran one component at a time.** Bioconda's `trinity=2.15.2` builds
   `_4`, `_5`, `_6` (May 2025 on; `orp_trinity` was created 2026-08-14 with
   `_6`) ship a ParaFly compiled without OpenMP, so `ParaFly -CPU 40` runs its
   commands serially. Builds `_0`-`_3` are fine. Every SRR1789336 benchmark
   in NOTES.md since 2026-08-14 measured this.
2. **Trinity's per-component work can be cut ~2x** with patches that leave
   the assembly unchanged.
3. **Trinity's output is not reproducible as shipped**: the SuperTranscripts
   polish depends on the run directory's path, and Butterfly on the JVM's
   identity hash sequence (JDK, locale, flags). Fixed; needed anyway to
   validate the speed patches byte-for-byte.

## Method

`harness/replay.py` copies real phase-2 commands (a seeded sample, or all of
them) into a sandbox, runs them with a given Trinity under ParaFly, and
compares every component's `.out.Trinity.fasta` with a reference. Sample:
1,500 of 73,737 components (seed 1; 2.0% of jobs, 2.3% of read bytes). Each
run: one exclusive 40-core macmanes node, GPFS, 38 slots (sample) or 40
(full). `patches/profiling-only/` + `harness/TrinProf.pm` add an opt-in log
(`TRINITY_PROF_LOG`) of wall and child CPU for every command Trinity runs;
`harness/javashim`, `pyshim` time each java/python call;
`harness/analyze.py` aggregates.

## 1. Serial ParaFly

- A replay launched exactly like production (`conda run --no-capture-output
  -n orp_trinity ParaFly -c cmds -CPU 40 -v -shuffle`) left the node 97% idle:
  one `perl`, one `java`, ParaFly `Threads: 1`. First job finished at 186 s.
- `nm -D .../trinity-plugins/BIN/ParaFly | grep -c GOMP_`: build `_6` 0,
  `_5` 0, `_4` 0, `_3` 6, `_2` 6, `_1` 6, `_0` 6, source build 6.
- Cause (`harness/pf_fix_test.sh` reproduces it):
  - bioconda's `makefile.patch` replaces ParaFly's `CFLAGS="-fopenmp"
    CXXFLAGS="-fopenmp"` with the recipe's flags. Until the 2025-05-29 rebuild
    (bioconda-recipes #56483) those flags included `${LDFLAGS}`, which holds
    `-fopenmp`, so it worked by accident; since then they don't.
  - ParaFly's `configure.ac` would catch this, but `AC_OPENMP` compiles its
    probe with `LDFLAGS`, so `-fopenmp` there reads as "none needed" while the
    object compiles never see `LDFLAGS`: the `#pragma omp parallel for` is
    dropped and libgomp is still linked.
  - Separately, `AC_SUBST([AM_CXXFLAGS], [-m64 $OPENMP_CXXFLAGS])` is unquoted,
    so when `AC_OPENMP` does find `-fopenmp` configure executes it as a
    command (`configure: line 3103: -fopenmp: command not found`) and drops
    it. That is presumably why `trinity-plugins/Makefile` hard-codes
    `CXXFLAGS="-fopenmp"` for ParaFly.
- This explains the earlier phase-2 observations in
  `docs/trinity-speedup-investigation.md` and NOTES.md: ~36 jobs/min is about
  1.7 s per component run serially; `--cpu 80` could not help; 95% vs 100% of
  `--cpu` made no difference. Phase 2 was not CPU-bound on 40 cores.
- Fixed ParaFly, same conda Trinity, all 73,737 components: 7,394 s (2h03m).
  65,617 of 66,440 outputs identical to production; the 823 that differ are
  section 3's path dependence (the replay runs in a different directory).

## 2. Where a component's time goes

Baseline (det-fix, shipped Butterfly), 1,500-component sample: median
component 3.4 s wall; median read file 1.7 KB (about 7 reads); 99% of jobs
under 15 s. Almost all of it is fixed per-run cost.

| stage | CPU share | cause |
| --- | --- | --- |
| Butterfly | 42% | JVM C2 JIT threads (~1.7 CPU-s per call on a 40-core node, `cpu/wall` 1.8); debug strings built and thrown away |
| GraphFromFasta | 10% | allocates, zeroes, sorts and frees a dense 4^12-slot k-mer table (~450 MB) per component |
| bowtie2-build, bowtie2 + samtools | 12% | python/perl wrappers; run even with fewer than 2 contigs to scaffold (74% of components) |
| SuperTranscripts polish | 7% | `import numpy` (~0.1 s) for one `argmax` and one `zeros` |
| shell-outs, tool probes, Perl startup | ~12% | `sh -c command -v` x5, `samtools --version`, `jellyfish --version`, `sort --help`, `touch`/`mkdir -p` via bash, every component |

Butterfly, JFR on mid-size graphs: 60-75% of samples in
`SeqVertex.getWeightAvg()` via `toString()`, called while building
`debugMes(...)` arguments that `debugMes` then discards (`DEBUG` is a
compile-time `true`; the level test is inside the call; 775 call sites).

The five largest components in the whole run (4.9-6.7 MB of reads) take
20-27 s each, so there is no long tail.

## 3. Reproducibility

- **SuperTranscripts polish** (`Compact_graph_pruner.pop_small_bubbles`) keeps
  `bubble_node_list.pop()` from a list built off a `set` of `TNode`s, which
  hash by `id()`: the base kept at a SNP-sized bubble depends on memory
  layout. The same `Trinity.tmp.fasta` gave four different outputs depending
  only on the input/output path strings. The code has `TODO: should sort by
  expression`.
- **Butterfly** keys HashMaps/HashSets (jung keys its edge maps by
  `SimpleEdge`) with eight classes that have no `hashCode`, so order follows
  the JVM identity hash. For ~0.15% of components (2 of 1,346 in the sample)
  the output flips between two alternatives with the identity hash sequence:
  `-XX:hashCode` mode, the locale (it changes how many identity hashes the
  JVM uses at startup; Trinity sets `LC_ALL=C` in `COMMON.pm`), and the
  debug-message patch all move it.
- **Butterfly source**: Trinity's `Butterfly/Butterfly` submodule points at
  the 2019 `master` ("init"); the shipped `Butterfly.jar` (classes dated
  2022-03-10) matches `origin/devel` (d3cb730). Building from the submodule as
  pinned gives a different Butterfly. Patches here are on `devel`.
- `scaffold_iworm_contigs.pl` breaks count ties in Perl's randomized hash
  order (no downstream effect seen; fixed anyway).
- Not an issue: locale for `sort`. `COMMON.pm` sets `LC_ALL=C` for the run; an
  early apparent locale effect was the path dependence above.

## Patches (`patches/`)

| repo | patch | effect | validation |
| --- | --- | --- | --- |
| trinityrnaseq | 0001 det-fix | bubble representative = node in most isoforms, ties by loc id; scaffold tie-break | same output at 8 paths; two replays at different paths: 1346/1346 |
| trinityrnaseq | 0002 Butterfly C1 JIT, `--bfly_full_jit` restores | Butterfly CPU -58%; 5 largest components: wall within +/-15%, CPU 2.3-3.3x lower | identical |
| trinityrnaseq | 0003 numpy-free polish | -0.1 s per component | identical |
| trinityrnaseq | 0004 fewer processes per component | PATH lookup in Perl; no version probes in `--trinity_complete` runs; one sort probe per tree; Perl touch/mkdir/checkpoints; `bowtie2-build-s`/`bowtie2-align-s` directly; no scaffolding with <2 contigs | identical |
| trinityrnaseq | 0005 ParaFly configure | OpenMP whatever CFLAGS/CXXFLAGS/LDFLAGS (configure regenerated with automake 1.16, hence the size; the change is 8 lines of configure.ac) | `pf_fix_test.sh`: plain / bioconda-like / Trinity-Makefile flags give 6/6/6 GOMP symbols (before 0/0/6) |
| Chrysalis | 0001 sparse k-mer table | GraphFromFasta 0.20 -> 0.02 s, 464 -> 6 MB on a small component; dense above bound/16 positions (phase 1) | dense/sparse/auto identical; `CHRYSALIS_KMER_TABLE=dense\|sparse` forces either |
| Butterfly (devel) | 0001 reproducible hash codes | creation counter through MurmurHash3 fmix32 | c1122, c3226: one output across 2 locales x 4 hash modes x guard on/off; sample vs shipped jar: 1344/1346 identical, +0.09% transcripts |
| Butterfly (devel) | 0002 build debug messages only when printed | Butterfly 2.4 -> 0.6 s (c3986), 2.3 -> 0.8 (c1128), 0.55 -> 0.31 (c77) | identical given 0001 |

Butterfly jars are built by `harness/build_bfly.sh` (javac `--release 8`,
classes overlaid on the shipped fat jar); `devel` built unmodified
reproduces the shipped jar.

Sample, 1,500 components, 38 slots:

| run | ParaFly wall | job CPU | median job | Butterfly CPU |
| --- | --- | --- | --- | --- |
| base5: det-fix + Butterfly 0001 | 157 s | 4,036 s | 3.4 s | 1,659 s |
| perf5: + every other patch | 75 s | 1,527 s | 1.4 s | 449 s |

perf5 vs base5: 1346/1346 identical (154 components produce no transcripts in
either).

Full phase 2, 73,737 components, 40 slots:

| run | wall | outputs |
| --- | --- | --- |
| base5 | 7,545 s (2h06m) | reference |
| perf5 | 3,499 s (58m) | 66,440/66,440 identical to base5; none failed |

## Upstream actions (for the user)

1. **bioconda-recipes `recipes/trinity`** (affects every bioconda Trinity
   install since May 2025): in `makefile.patch` give ParaFly
   `CFLAGS="${CFLAGS} -fopenmp" CXXFLAGS="${CXXFLAGS} -fopenmp"`, bump the
   build number, and add a test that fails on a serial ParaFly
   (`nm -D $PREFIX/bin/trinity-plugins/BIN/ParaFly | grep -q GOMP_`, or two
   `sleep 1` commands at `-CPU 2` finishing in under ~1.8 s).
2. **trinityrnaseq**: 0001-0005; independent of each other except that 0002's
   flag lives in the Trinity script. Also: point the Butterfly submodule at
   the commit the jar is built from.
3. **trinityrnaseq/Chrysalis** 0001 and **trinityrnaseq/Butterfly** (devel)
   0001-0002, then a jar rebuild.
4. Seen in passing: `Analysis/SuperTranscripts/pylib` has non-raw regex
   strings that Python 3.14 warns about and a later Python will reject.

## ORP-side options (not done; need a decision)

- Pin `trinity=2.15.2=*_3` in `Makefile`/`Dockerfile` (a `mamba create
  --dry-run` with ORP's other `orp_trinity` packages and `python=3.14` solves,
  openjdk 23.0.2).
- Preflight in `oyster.py`: fail if Trinity's ParaFly runs two `sleep 1`
  commands at `-CPU 2` in more than ~1.8 s.
- Revisit the lane design (`TRINITY_PHASE1_SHARE`/`PHASE2_SHARE`, pairing
  phase 2 with Trans-ABySS): it was tuned against a 34-hour serial phase 2.
  At 1-2 h, Trans-ABySS (4.5 h, mostly single-threaded) is the long pole.
