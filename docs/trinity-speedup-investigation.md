# Trinity Phase 1 / Phase 2 speedup investigation

**Status:** written 2026-08-21 as research only. **Updated 2026-08-24** with
results from the first run that actually tested any of it -- the `--cpu 80`
oversubscription run (`sampledata/benchmarks.md`, 2026-08-24). Section 3.3 is
now **falsified** and 3.4 is **downgraded** on the strength of that run; 2d
gained a confirmation. Everything else in here remains untested source-reading.

**This investigation was closed on 2026-08-22 with ORP 3.1.0** (NOTES.md), which
ruled out 3.1 (`--grid_exec`) on the grounds that it needs cluster-specific
setup and so would not generalize to ORP's end users. The 2026-08-24 updates
were folded in afterwards because the run was already in flight; they record a
result against a closed investigation rather than reopening it.

Original framing, still true of the untested sections: written from a session
with no cluster access, every number either read out of our own finished runs
(`sampledata/benchmarks.md`) or read out of the Trinity source.

Trinity version this applies to: **2.15.2**, the version pinned for the
`orp_trinity` env in [Makefile](../Makefile) (line 50). The `Trinity` driver
script is byte-identical between the `Trinity-v2.15.2` tag and current `master`
(4278 lines both), so line numbers cited below are valid for either.
`Inchworm/` is a git submodule -- its sources live at
<https://github.com/trinityrnaseq/Inchworm>, not in the main trinityrnaseq repo.

---

## 1. The time budget we're actually optimizing

From `TIME2_SRR1789336_norm_py_955parallel`, our current best design
(`--cpu 40 --mem 670`, see [sampledata/benchmarks.md](../sampledata/benchmarks.md)):

| stage | wall time | % of TOTAL |
| --- | --- | --- |
| `run_trinity_phase1` | 1:22:49 | 3.7% |
| `run_trinity_phase2` | 34:27:02 | **92.9%** |
| everything else combined | ~1:16 | 3.4% |
| **TOTAL** | **37:06:19** | |

Consequences, before any mechanism discussion:

- Deleting Phase 1 **entirely** buys 3.7%. Any Phase 1 tuning is bounded by
  1h22m, and realistically by a fraction of that.
- Making Phase 2 2x faster buys ~46% of total wall clock.

This is the single most important framing in this document. The 25/75 -> 75/25
Stage A question is real but small; Phase 2 is where the run lives.

**Second data point, 2026-08-24** (`--cpu 80` oversubscription run, same
dataset and node):

| stage | cpu 80 | cpu 40 | delta |
| --- | --- | --- | --- |
| preprocessing (trim + rcorrector) | 0:18:52 | 0:23:20 | **-4m28s** |
| Stage A (Phase 1-bound) | 1:09:08 | 1:22:49 | **-13m41s** |
| Stage B (Phase 2-bound) | 36:09:58 | 34:27:02 | **+1h42m56s** |
| tail (filtershort -> transrate) | 1:08:06 | 0:53:08 | **+14m58s** |
| **TOTAL** | **38:46:04** | **37:06:19** | **+1h39m45s (+4.5%)** |

The decomposition reconciles to the second. It reinforces the framing above
rather than changing it: the two segments that improved are worth ~18 min
combined and are structurally capped, while the single segment that regressed
is the one holding 93% of the run.

---

## 2. Phase 1: raising `--inchworm_cpu` past 10 will not help

We currently pass `--inchworm_cpu 10` ([oyster.py:436](../oyster.py#L436)), and
Phase 1 gets `TRINITY_PHASE1_SHARE = 0.25` of `--cpu`
([oyster.py:61](../oyster.py#L61)), i.e. 10 cores at `--cpu 40`. So inchworm is
currently getting every core Stage A gives Phase 1.

Three independent reasons not to raise it:

### 2a. The binary hard-caps kmer parsing at 6 and only ever clamps *down*

`Inchworm/src/IRKE_run.cpp:29`:

```cpp
int IRKE_COMMON::MAX_THREADS_READ_PARSING = 6;  // going higher leads to decreased performance in Inchworm due to thread collisions.
```

`IRKE_run.cpp:292-293` only lowers it (`if (MAX_THREADS_READ_PARSING >
NUM_THREADS) MAX_THREADS_READ_PARSING = NUM_THREADS;`), never raises it.
`IRKE.cpp:89` then pins the kmer-loading loop to that value regardless of what
`--inchworm_cpu` said. That stage is capped at 6 threads no matter what we pass.

### 2b. Contig building *does* use all our threads -- and that's the contended part

`Inchworm/src/IRKE.cpp:466-470`:

```cpp
if (PARALLEL_IWORM) {
    omp_set_num_threads(IRKE_COMMON::NUM_THREADS);
} else {
    omp_set_num_threads(1); // turn off multithreading for the contig building.
}
```

`$PARALLEL_IWORM_FLAG = 1` by default (`Trinity:204`), so contig building really
is running on our 10 threads. But the loop is
`#pragma omp parallel for schedule(dynamic,1000)` over the kmer list, with
threads zeroing out kmers in a single shared hash guarded by
`#pragma omp critical (HashMap)` (`KmerCounter.cpp:406`, `:482`).

That is exactly the contention Brian Haas describes when asked why inchworm is
capped, in <https://github.com/trinityrnaseq/trinityrnaseq/issues/648>:

> The performance of inchworm doesn't scale well beyond 6 threads, so it's
> capped at this level. It uses shared memory and there are parts where the
> threads compete for resources.

The [Running Trinity wiki](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity)
repeats it: "Inchworm alone will be internally capped at 6 threads, since
performance will not improve for this step beyond that setting." Trinity's own
default is `my $inchworm_cpu = 6` (`Trinity:119`), clamped only downward against
`--CPU` (`Trinity:1232-1233`). **We are already at 10, past the documented cap.**

### 2c. It is not a free knob -- it changes the assembly

Because `PARALLEL_IWORM` iterates an *unsorted* kmer list (`IRKE.cpp:451`
comments this explicitly) and threads zap kmers out of the shared hash as they
consume them, which kmers win as seeds depends on thread interleaving. So
inchworm contig output is thread-count-dependent.

**This is the mechanism behind the open reproducibility question in
[NOTES.md](../NOTES.md) 2026-08-16 (5)** (2.4.0 vs 3.0.0 giving slightly
different Trinity output, theory being "Trinity isn't perfectly reproducible
across thread counts"). That theory is correct and this is the reason. It also
means changing `--inchworm_cpu` would perturb our BUSCO / unique-gene numbers
for no speed gain.

### 2d. What a 75/25 flip *would* actually affect

Not inchworm -- the Chrysalis stages, which key off `--CPU` rather than
`--inchworm_cpu`:

| stage | Trinity line | threading |
| --- | --- | --- |
| `jellyfish count` | 2612 | `-t $CPU` |
| `GraphFromFasta` | 2180 | `-t $CPU` |
| `ReadsToTranscripts` | 2240 | `-t $CPU` |
| `sort` of read->component map | 2257 | **single-threaded**, `-S $max_memory` |
| `inchworm` | 2686-2687 | `--num_threads $inchworm_cpu` (see above) |

So flipping Stage A's split would speed up three of five sub-stages, with an
absolute ceiling of ~45 min, minus whatever SPAdes gives back by dropping from
30 cores to 10. Probably worth well under half an hour. Matches the intuition
that this isn't the lever.

> **Confirmed 2026-08-24.** The `--cpu 80` run gave Phase 1 twenty slots instead
> of ten with `--inchworm_cpu` still pinned at 10, and Phase 1 dropped
> **1:22:49 -> 1:09:07 (-13m42s, -16.5%)**. Since inchworm's own thread count
> did not change, that gain has to be coming from the three `-t $CPU` Chrysalis
> stages in the table above, exactly as predicted. Stage A is Phase-1-bound, so
> it shrank by the full amount.
>
> **But do not buy it with a global `--cpu` bump** -- that run lost far more
> elsewhere (see 3.3). Buy it with `TRINITY_PHASE1_SHARE`. The slack is
> demonstrably there: in that same run the SPAdes lane finished at 12:26:53 and
> Phase 1 ran until 13:19:38, leaving **52m38s** of idle SPAdes headroom, and
> SPAdes was flat between 30 and 60 slots (8:30 vs 8:36). Concretely:
> **`--cpu 40` with `TRINITY_PHASE1_SHARE = 0.5`**
> ([oyster.py:61](../oyster.py#L61)) should capture most of the 13m42s with no
> oversubscription penalty anywhere else in the run. This is now the
> recommended Stage A change, and it supersedes the 75/25 framing above.

**Cheaper and lower-risk than the CPU flip: decouple mem share from cpu share in
Stage A** ([oyster.py:1075-1078](../oyster.py#L1075)). Today both are
`* TRINITY_PHASE1_SHARE`. That single-threaded `sort -T . -S $max_memory` runs
in the Trinity output dir on GPFS -- more RAM keeps it from spilling temp files
onto the parallel filesystem, more cores do nothing for it. SPAdes55/75 do not
need 500 GB. Giving Phase 1 a large mem share while keeping its CPU share small
costs nothing and may take a bite out of the sort.

### 2e. Free diagnostic -- get the real Phase 1 breakdown before touching anything

Trinity leaves timestamped checkpoint files. In the trinity output dir of any
completed run:

```bash
ls -la --time-style=full-iso .jellyfish_count.* .jellyfish_dump.* .jellyfish_histo.* .iworm.25..ok chrysalis/*.ok partitioned_reads.files.list.ok recursive_trinity.cmds.ok
```

Successive mtimes give the per-stage split of the 1h22m for free, with no rerun.
Do this before spending a run on any Phase 1 change -- it tells us whether the
time is in jellyfish, inchworm, GraphFromFasta, ReadsToTranscripts, or the sort,
and only some of those respond to cores at all.

(`Trinity.timing` exists but is only written under `--monitoring`, which shells
out to collectl and is Linux-only -- `Trinity:171`, `:347`, `:1071`. The mtime
trick needs no rerun and no extra flag.)

---

## 3. Phase 2: where the 34.5 hours actually are

Phase 2 is:

```
ParaFly -c recursive_trinity.cmds -CPU $CPU -v -shuffle
```

(`Trinity:3652`), over **73,737** independent mini-Trinity jobs, each generated
with `--CPU $grid_node_CPU --max_memory $grid_node_max_memory`
(`Trinity:3669`), whose defaults are **1 core and 1G** (`Trinity:219-220`).

So: 73,737 single-core, 1 GB jobs, currently squeezed through 38 slots on one
node, against a machine with 670 GB of RAM. Options below are ranked by expected
payoff.

### 3.1 `--grid_exec` -- the only order-of-magnitude lever

> **Ruled out by decision, 2026-08-22** (NOTES.md, ORP 3.1.0): not being
> pursued. It needs HpcGridRunner plus a scheduler-specific config tuned to a
> particular cluster, so it would not generalize to ORP's end users, each of
> whom would need their own setup rather than a flag that just works. That is a
> scope judgement about ORP, not a performance one, and the 2026-08-24 run does
> not bear on it.
>
> Recorded here because the mechanics still hold and the section is the
> reference if that scope judgement is ever revisited: with 3.3 falsified and
> 3.4 downgraded, this is the only identified option that changes the order of
> magnitude. A CPU-bound Phase 2 on a single node cannot be fixed on that node;
> it needs more nodes.

`run_partitioned_cmds` branches (`Trinity:3644-3652`): if `--grid_exec` is set it
runs `$grid_exec_toolname $cmds_file` instead of ParaFly. Trinity's recommended
tool is HpcGridRunner, with reported support for LSF/SGE/SLURM/PBS.

This is mechanically **clean** given how oyster.py already splits the phases:

- `recursive_trinity.cmds` is written and `.ok`-checkpointed during Phase 1
  (`Trinity:2511-2513`), so adding a flag at Phase 2 time cannot change the
  per-component command list -- it's already on disk and won't be regenerated.
- `--grid_exec` is stripped out of the per-component commands by the arg filter
  at `Trinity:3706` (the alternation includes `grid`, which substring-matches
  `--grid_exec`, and the filter consumes its value too). So there's no risk of
  each mini-job recursively trying to grid-dispatch.
- Change is confined to `run_trinity_phase2` ([oyster.py:454](../oyster.py#L454)).

Scale precedent: GACRC document a Trinity Phase 2 of **2,020,460** mini-assemblies
finishing in ~19h via HpcGridRunner at `cmds_per_node`/`max_nodes` = 100
(<https://wiki.gacrc.uga.edu/wiki/Trinity-HpcGridRunner>). We have 27x fewer jobs
taking 1.8x longer. At 160 slots instead of 38 this is roughly 34h -> ~9h.

Open questions: whether we can `sbatch` from inside an allocation on Premise, and
what the fair-share cost of a 73k-job campaign is.

### 3.2 Our normalization result is a no-op -- and NOTES.md has it wrong

`Trinity:214`:

```perl
my $normalize_max_read_cov = 200; # better for polymorphic transcriptomes
```

Trinity 2.15's in-silico normalization defaults to a **200x** coverage cap, which
for typical RNA-seq depth discards almost nothing. That is very likely why
Trinity came out flat in the `--normalize-reads` comparison.

Worse, the record is wrong: [NOTES.md](../NOTES.md) 2026-08-17 describes the
38:28:11 run as the "no `--normalize-reads`" baseline, but
[sampledata/benchmarks.md](../sampledata/benchmarks.md) line 9 shows its command
line *did* include `--normalize-reads`. So did all four SRR1789336 runs
(`_parallel`, `_NOparallel`, `_5050parallel`, `_955parallel`).

**Conclusion: normalization has never actually been tested as a speed lever
here.** Every run normalized, and every run normalized at a cap high enough to be
nearly inert. Phase 2 cost scales with reads-per-component, so a genuine test at
`--normalize_max_read_cov 50` (the classic Trinity setting) is the highest-value
cheap experiment available.

oyster.py currently gives no way to set it -- `--normalize-reads` is a bare
boolean ([oyster.py:1260](../oyster.py#L1260)) and `_trinity_base_cmd`
([oyster.py:422](../oyster.py#L422)) never passes `--normalize_max_read_cov`.
Would need a new flag.

Caveat: this changes the assembly, so it needs the usual BUSCO / TransRate /
unique-genes comparison, not just a stopwatch.

### 3.3 Oversubscribe ParaFly -- **FALSIFIED 2026-08-24, do not retry**

ParaFly's `-CPU` is just a concurrency count -- it runs each command via
`system()` in its own slot (`ParaFly.cpp:61`, `:175`). Each of our jobs is
1 core / 1 GB against 670 GB of RAM, so memory is nowhere near binding.

If Phase 2 is I/O-latency-bound rather than CPU-bound (likely, see 3.4), running
`phase2_cpu` at 1.5-2x physical cores is close to free throughput. Today
`phase2_cpu = round(self.cpu * 0.95)` ([oyster.py:1080](../oyster.py#L1080)) can
never exceed `--cpu`.

**This is testable against a live run in ~20 minutes** by comparing ParaFly's
`succeeded(N) ...% completed` rate at two different `-CPU` values -- we already
read that counter for progress monitoring (NOTES.md 2026-08-20 (1)).

> **Result: wrong. Tested 2026-08-24 at `--cpu 80` (76 ParaFly slots on 40
> physical cores).**
>
> | | 76 slots | 38 slots |
> | --- | --- | --- |
> | `run_trinity_phase2` | 36:09:58 | 34:27:02 |
> | throughput | **33.98 jobs/min** | **35.67 jobs/min** |
>
> Doubling ParaFly's concurrency made Phase 2 **1h42m56s slower (+5.0%)**. The
> premise of this section -- "if Phase 2 is I/O-latency-bound" -- is false.
> Phase 2 is **CPU-bound**, and the extra slots bought nothing but
> context-switching. The memory reasoning above was correct and irrelevant:
> memory was never the binding constraint.
>
> The run also showed the tax lands everywhere else already core-saturated:
> orthotransrate +40%, orthofusing +81%, busco +45%, transrate +26%, strandeval
> +42%, for +14m58s on the post-Trinity tail. Net for the whole run
> **+1h39m45s (+4.5%)**. Full decomposition in `sampledata/benchmarks.md`.
>
> One caveat not closed: Phase 1 ran at 20 threads instead of 10, and per 2c
> inchworm output is thread-count-dependent, so this run's component count may
> not be the baseline's 73,737. `wc -l recursive_trinity.cmds` would settle
> whether part of the +5% is extra work rather than worse throughput. Even if
> it is, nothing here argues for oversubscription.

### 3.4 Get the Trinity working dir off GPFS -- **downgraded 2026-08-24**

Our assemblies live under `/mnt/gpfs01/home/macmaneslab/...`. Phase 2 creates
73,737 job directories, each doing dozens of small file creates/writes/deletes,
with `--full_cleanup` tearing them down as it goes. That's a metadata storm, and
it's the standard reason Phase 2 crawls on shared HPC storage -- Purdue RCAC
explicitly recommend node-local disk for Trinity for this reason
(<https://www.rcac.purdue.edu/news/555>).

Node-local NVMe (or `/dev/shm`, we have 670 GB) for `--output`, then copy the
resulting fasta back. Phase 1 and Phase 2 are separate oyster.py steps but run in
the same process on the same node, so the dir persists between them.

**Conflicts with 3.1** -- grid dispatch needs shared storage. Pick one.

Diagnose which regime we're in by checking `%iowait` (`sar`, `top`) during
Phase 2. High iowait -> 3.4 and 3.3. Pegged CPU -> 3.1 is the only real answer.

> **Downgraded 2026-08-24.** The 3.3 result answers that diagnostic from the
> other direction: if Phase 2 were starved waiting on GPFS metadata, doubling
> the number of in-flight jobs would have soaked up the idle time and gone
> *faster*. It went slower. So Phase 2 is CPU-bound and this section's premise
> -- that the filesystem is the binding constraint -- is probably wrong too.
>
> Not fully closed: this is an inference from the 3.3 result, not a direct
> measurement, and a metadata-server bottleneck could in principle degrade
> under added concurrency rather than absorb it. A direct `%iowait` reading
> during any future Phase 2 settles it for free. But it is no longer worth a
> run of its own, and 3.1 does not depend on the answer.

### 3.5 `--min_kmer_cov 2`

Default is 1 (`Trinity:117`), applied as `jellyfish dump -L $min_kmer_cov`
(`Trinity:2625`). Raising to 2 drops singleton kmers, which means fewer spurious
inchworm contigs, which means fewer and smaller components, which means fewer
Phase 2 jobs. Helps both phases.

The usual objection is loss of sensitivity for lowly-expressed transcripts. That
objection is weaker for us specifically because ORP already runs rcorrector
before Trinity, so error kmers are largely corrected rather than merely rare.

Like 3.2, this changes the assembly and needs the quality battery to sign off.

### 3.6 Considered and deprioritized

- **`--max_reads_per_graph`** (default 200000, `Trinity:144`): caps the biggest
  components. Our observed rate profile doesn't show a long slow tail -- the
  final 24.5% of the 50/50 run ran at ~87/min vs a 28.6/min run average
  (NOTES.md 2026-08-20 (1)) -- so the big jobs aren't stranded at the end. Not a
  priority.
- **Butterfly JVM tuning** (`-Xmx10G -Xms1G -Xss1G`, `Trinity:2373`, defaults at
  `Trinity:102-103`): 38 concurrent JVMs each committing 1 GB init heap is
  trivial against 670 GB. Lowering `--bflyHeapSpaceInit` would also shrink
  `-Xss`, which Butterfly needs large because it recurses. Low reward, nonzero
  risk. Skip.
- **Dynamic core handoff when Trans-ABySS finishes**: Trans-ABySS completes 5h09m
  into Stage B, leaving 5% of cores idle for the remaining ~29h. Same class of
  problem as the one fixed in NOTES.md 2026-08-16 (2). But measured Phase 2
  throughput at 95% vs 100% `--cpu` was within noise (38.3 vs 36.2 jobs/min), so
  the recoverable time is small. Would fall out for free from 3.1 or from an
  oyster-controlled dispatcher anyway.

---

## 4. Suggested order of work

**Revised 2026-08-24** after the `--cpu 80` run.

Context: this investigation was formally closed on 2026-08-22 with ORP 3.1.0,
having ruled out `--grid_exec` on generalizability grounds (3.1). What follows
is what remains actionable within that decision, not an argument to reopen it.

Free, no rerun required:

1. `wc -l recursive_trinity.cmds` on the `--cpu 80` run. Confirms whether its
   Phase 2 had the same 73,737 jobs as the baseline, which is the one loose end
   in the 3.3 falsification.
2. Pull the Phase 1 `.ok` mtimes from a finished run (section 2e). Now more
   useful than before: we know Stage A responds to cores, and this says which
   sub-stage is actually paying.

Then, one run each, in this order:

3. **`--cpu 40` with `TRINITY_PHASE1_SHARE = 0.5`** (2d). One constant, no new
   flags, ~13 min expected, and it rides along with any other run -- it does not
   need a dedicated one.
4. **`--normalize_max_read_cov 50`** (3.2). One flag, real mechanism, needs a
   new oyster.py option. Still the highest value per run spent, and completely
   untouched by the `--cpu 80` result.
5. `--min_kmer_cov 2` (3.5), if 4 disappoints.

Ruled out, do not spend runs on: oversubscribing `--cpu` anywhere (3.3, tested
2026-08-24 and lost 1h40m), node-local scratch (3.4, downgraded by inference
from the same run), `--inchworm_cpu` changes (2a-2c), and the 75/25 Stage A CPU
flip as originally framed -- superseded by item 3. Out of scope by decision:
`--grid_exec` (3.1).

Worth stating plainly: with 3.1 out of scope and 3.3 dead, nothing
order-of-magnitude remains identified. Items 3 and 4 are worth roughly 13
minutes and an unknown-but-probably-hours amount respectively, against a 37-hour
run. Section 5 is the honest place to look next.

---

## 5. The thing worth saying out loud

Trinity contributes 13,231 unique genes to ORP's 13,772 on SRR1789336, and
SPAdes55 alone contributes 13,493 -- for 93% of the pipeline's wall clock. The
largest available speedup is not inside Trinity. That's a scientific call about
what ORP is for, not a performance one, but it should be stated explicitly
somewhere in this repo rather than left implicit in the timing table.

---

## Sources

- Trinity driver script, v2.15.2 tag: <https://github.com/trinityrnaseq/trinityrnaseq/blob/Trinity-v2.15.2/Trinity>
- Inchworm submodule sources: <https://github.com/trinityrnaseq/Inchworm> (`src/IRKE.cpp`, `src/IRKE_run.cpp`, `src/KmerCounter.cpp`)
- Inchworm thread cap, dev response: <https://github.com/trinityrnaseq/trinityrnaseq/issues/648>
- Running Trinity wiki: <https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity>
- Trinity + HpcGridRunner at scale: <https://wiki.gacrc.uga.edu/wiki/Trinity-HpcGridRunner>
- Trinity on shared vs local filesystems: <https://www.rcac.purdue.edu/news/555>
- Trinity RNA-Seq Assembler Performance Optimization (XSEDE12): <https://wwwpub.zih.tu-dresden.de/~mlieber/publications/2012_XSEDE12_Trinity_post_print.pdf>
