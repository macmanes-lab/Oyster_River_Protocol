#!/usr/bin/env python3
"""Python port of oyster.mk - the Oyster River Protocol assembly pipeline.

Usage:
    oyster.py --read1 R1.fq.gz --read2 R2.fq.gz --mem 110 --cpu 24 \\
              --runout myrun --strand RF

Runs the full pipeline end to end. Steps whose output files already exist
and are newer than their inputs are skipped, so a failed/interrupted run
can be re-invoked to resume where it left off.
"""

import argparse
import csv
import gzip
import io
import json
import math
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import threading
import time
import concurrent.futures
from functools import partial
from pathlib import Path
from typing import NamedTuple

HERE = Path(__file__).resolve().parent
RED = "\033[31m"
RESET = "\033[0m"

class Assembly(NamedTuple):
    """One input assembly, and the four names the pipeline knows it by.

    The names are irregular because they are the ones oyster.mk used and the
    ones every existing run directory on disk already carries: the file is
    `<runout>.trinity.Trinity.fasta` but its diamond output is
    `<runout>.trinity.diamond.txt`, while SPAdes' unique-gene count lands in
    `<runout>.unique.sp75.txt` and not `...spades75...`. Renaming any of
    them would silently invalidate every resumable run directory that
    exists, so they are carried as data instead of being derived.

    The SPAdes assemblies are a deliberate exception. They changed which k
    values they run at, so a directory holding the old `spades55.fasta` no
    longer describes what this code would produce from the same reads.
    Keeping the old name would let step() find that stale file and skip the
    assembly, reporting an auto-k run while serving k=55 output -- the same
    silent wrongness the paragraph above is guarding against, arriving from
    the other direction. Renaming costs a re-assembly; not renaming costs a
    wrong answer that looks right.
    """

    fasta_name: str      # assemblies/<runout>.<fasta_name>
    diamond_label: str   # assemblies/diamond/<runout>.<diamond_label>.diamond.txt
    unique_label: str    # assemblies/diamond/<runout>.unique.<unique_label>.txt
    report_label: str    # reportgen's "UNIQUE GENES <report_label>" line


SPADES_AUTO = Assembly("spadesauto.fasta", "spadesauto", "spauto", "SPADESAUTO")
SPADES75 = Assembly("spades75.fasta", "spades75", "sp75", "SPADES75")
TRANSABYSS = Assembly("transabyss.fasta", "transabyss", "transabyss", "TRANSABYSS")
TRINITY = Assembly("trinity.Trinity.fasta", "trinity", "trinity", "TRINITY")

# oyster.py's own four assemblers, in the three different orders oyster.mk
# used them in. The three are not interchangeable and two of them reach the
# assembly, so they are spelled out rather than sorted:
#
#   ASSEMBLY_ORDER    concatenation order. Sets contig order in
#                     orthofuse/merged.fasta and in posthack's `cat` of the
#                     assemblies, which flows through to cd-hit-est -- where
#                     input order breaks length ties and so decides which
#                     representative survives into .ORP.fasta.
#   DIAMOND_PRIORITY  search order. build_list5.py keeps the *first* diamond
#                     hit per gene in the order it is given the files, so
#                     this is a preference ranking between assemblies for
#                     the contigs the orthogroup pass missed.
#   REPORT_ORDER      the order the UNIQUE GENES lines appear in
#                     reports/qualreport.<run>. Cosmetic, but people diff
#                     those reports across runs.
ASSEMBLY_ORDER = (SPADES_AUTO, SPADES75, TRANSABYSS, TRINITY)
DIAMOND_PRIORITY = (TRANSABYSS, SPADES75, SPADES_AUTO, TRINITY)
REPORT_ORDER = (TRINITY, SPADES_AUTO, SPADES75, TRANSABYSS)

# Preflight, in the order it prints. Everything here is shelled out to at
# some point in a full run, and finding it missing hours in -- at
# orthotransrate, or at the assembler that was going to run overnight -- is
# the thing this list exists to prevent. snap-aligner is on it because
# pytransrate maps with it.
SPADES_TOOL = ("orp_spades", "rnaspades.py", "SPADES")
TRINITY_TOOL = ("orp_trinity", "Trinity", "TRINITY")
TRANSABYSS_TOOL = ("orp_transabyss", "transabyss", "TRANSABYSS")
# The three an entry point that doesn't assemble has no use for.
ASSEMBLER_TOOLS = (SPADES_TOOL, TRINITY_TOOL, TRANSABYSS_TOOL)
#: The pytransrate this pipeline needs, checked at preflight rather than
#: assumed from orp_env.yml. The pin in that file describes the environment
#: as built; it says nothing about the environment as it actually is, and the
#: two diverge the moment anyone installs by hand or reuses an older env. The
#: gap is expensive in exactly one direction: 2.2.0 caps the read-metrics
#: workers against a memory budget and 2.2.1 stops a failed run deleting the
#: BAM, so an env still holding 2.1.0 runs for nine hours, is OOM-killed,
#: throws away the BAM, and does it again on the retry -- with a log whose
#: only sign of the problem is a version number in the banner. Checked in
#: seconds instead.
PYTRANSRATE_MIN_VERSION = "2.2.1"

#: --max-memory's own spellings in pytransrate, both of which it accepts. A
#: user who set one in --pytransrate-args means it, so pytransrate_memory_args
#: stands aside rather than passing the flag twice.
PYTRANSRATE_MEMORY_FLAGS = ("--max-memory", "--mem")

CHECK_TOOLS = (
    ("orp", "salmon", "SALMON"),
    ("orp", "pytransrate", "PYTRANSRATE"),
    ("orp", "seqtk", "SEQTK"),
    ("orp_busco", "busco", "BUSCO"),
    ("orp", "mcl", "MCL"),
    SPADES_TOOL,
    TRINITY_TOOL,
    ("orp", "trimmomatic", "TRIMMOMATIC"),
    TRANSABYSS_TOOL,
    ("orp", "run_rcorrector.pl", "RCORRECTOR"),
    ("orp_orthofinder", "orthofinder", "ORTHOFINDER"),
    ("orp", "snap-aligner", "SNAP-ALIGNER"),
)

# Reference profile (minutes, from a representative run at --max-parallel 2)
# used only to decide submission order within the two remaining
# run_parallel() concurrent groups (orthofuser_branch vs. merge_branch;
# transrate vs. strandeval) -- run the historically slow step first so it
# isn't left waiting behind a quick one. The assemblers no longer go
# through run_parallel (see TRINITY_PHASE1_SHARE/TRINITY_PHASE2_SHARE and
# the assembly-lane pairings in
# main()), so they have no entries here. Each dataset is normally assembled
# only once, so this is a fixed relative ranking rather than something
# learned per-run. Names with no entry sort after every hinted step, in the
# order they were given.
STEP_TIME_HINTS = {
    "merge_branch": 27,
    "orthofuser_branch": 6,
    "transrate": 16,
    "strandeval": 2,
}

# Trinity's own CPU/mem is fixed at launch for however long that stage
# runs, so which assembler it's paired against matters more than a single
# fixed lane ratio -- see main() for the two-stage pairing.
#
# Stage A: Phase 1 (Inchworm + Chrysalis prep, see run_trinity_phase1) pairs
# with SPAdes55/SPAdes75. Phase 1 is largely insensitive to its CPU share
# above its own --inchworm_cpu=10 cap (Chrysalis's clustering is brief next
# to Phase 2), while SPAdes is fast but does scale with cores -- so SPAdes
# gets the bulk of the machine and Phase 1 gets just enough to clear its cap.
TRINITY_PHASE1_SHARE = 0.25
#
# Stage B: Phase 2 (thousands of independent per-gene-component ParaFly
# jobs, see run_trinity_phase2 -- by far Trinity's dominant cost, and scales
# close to linearly with cores) pairs with Trans-ABySS instead of running
# alone after every assembler finishes. Trans-ABySS's own dominant cost
# (the initial FASTQ read + De Bruijn graph build) runs single-threaded no
# matter how many cores it's given -- traced to abyss-pe falling through to
# the plain, unthreaded `ABYSS` binary whenever neither Bloom-filter mode
# nor MPI is requested, which oyster.py never does (see NOTES.md
# 2026-08-19). Real-run evidence (TIME2_SRR1789336_norm_py_5050parallel,
# 2026-08-19): at 40 cores, Trans-ABySS took 4.5h (3.2h of it
# single-threaded and CPU-invariant) while Phase 2 alone took ~34h -- so
# Phase 2 gets the large majority of the machine, leaving Trans-ABySS just
# enough cores to keep its own threaded sub-stages moving. Mem is NOT split
# by this same ratio (see main()): Trans-ABySS's memory footprint doesn't
# shrink with its CPU share.
TRINITY_PHASE2_SHARE = 0.95
#
# Not yet validated on a real run -- next run should confirm Trans-ABySS
# doesn't OOM at its reduced mem share and doesn't become the long pole of
# Stage B (it shouldn't: even a large slowdown from fewer cores is still
# small next to Phase 2's ~34h).

# A step that fails on a cluster is often transient (node preemption,
# filesystem hiccup, scheduler blip) rather than a real bug, so retry before
# giving up -- previously any single failed subprocess call killed the whole
# run immediately, discarding hours of unrelated concurrent work in the
# other lane. This does not help steps that fail deterministically (e.g. a
# missing dependency), which will just fail the same way on every attempt.
STEP_RETRIES = 2
STEP_RETRY_DELAY = 60

# OrthoFinder's -t is the count of *concurrent diamond processes* it launches
# for its all-vs-all -- n_assemblies^2 of them, each given `-p 1` -- so set
# from cores alone, -t 20 on a four-assembly run puts 16 `--more-sensitive`
# diamonds on the node at once. Cap it by memory as well as by cores, at
# roughly one concurrent search per this many GB.
#
# The figure is diamond's own: its default block size is a fixed -b2.0, and
# "the program can be expected to use roughly six times this number of memory
# (in GB)" -- so ~12 GB per process, whatever the node. (--more-sensitive does
# not change it; only --very-sensitive and --ultra-sensitive do, to -b0.4.)
# An earlier version of this comment had each diamond sizing its block against
# whatever memory looked free when it started. That is not what diamond does,
# and it is worth being precise about: a fixed per-process cost is one this
# cap can actually model.
#
# Concurrency is capped by the searches that exist, n_assemblies^2, before it
# is capped by anything here: a four-assembly run cannot put more than 16
# diamonds on a node whatever -t says.
#
# This comment used to conclude from that that search memory "cannot exceed
# ~16 * 12 = 192 GB", and so that on anything bigger than a ~200 GB node the
# cap was structurally incapable of being what OOMs the job. The arithmetic
# was right and the premise was wrong. A four-assembly chowder run on a 720 GB
# cgroup was killed five times over by the memory cgroup's OOM killer, every
# victim a diamond, at 92, 129, 133, 138 and 145 GB resident -- eleven times
# the figure below, and 638 GB between the five of them. Per-process cost is
# not fixed. It is set by what is in the query, and nothing here can see that.
#
# What made those five expensive was poly-asparagine: OrthoFinder searches DNA
# with `diamond blastp`, SPAdes gap-fills with N, and N is asparagine. That
# particular cause is dealt with upstream now, in write_search_inputs, which
# is the reason this constant was left at 12 rather than raised to fit the
# measurements above -- fitting it to them would cost every run wall time to
# insure against an input class that no longer reaches diamond. The number is
# still a guess about a cost it cannot measure, so treat it as a floor and not
# a guarantee: check_orthofinder_searches is what actually catches an
# all-vs-all that lost searches, whatever the reason.
#
# 2026-09-19: the masking landed and the run failed again, identically. Same
# four assemblies, same node, 670 GB, `--cpu 40`: seven of sixteen searches
# lost, every one of them with a SPAdes assembly as its query. Species0 lost
# all four of its searches and Species1 three of four; Species2 and Species3
# lost none. Four died on returncode -9 and three on returncode 1 having
# written an empty file, all seven between 09:42 and 10:03 after two and a
# half hours of running, and the nine survivors then finished between 13:05
# and 16:08 -- the shape of a node that ran out of memory all at once and of
# survivors that only had room once the kernel had made some.
#
# So masking the N runs was necessary and was not sufficient. What it did not
# do is bound anything: it removed one known source of seed hits from one
# input class, and the cost of a search is still set by how many seed hits
# the query actually makes, which nothing here can see in advance.
#
# The bug this constant had all along is that it could never bind. `searches`
# is min(cpu, mem // this), OrthoFinder runs n_assemblies^2 diamonds and no
# more, and 16 is below min(40, 670//12 = 55) -- so the cap was computed,
# logged, and then had no effect on a four-assembly run, which is every
# chowder run and every ORP run. At 670 GB it takes a figure above 670/16 =
# 42 GB before the cap changes a single thing. 12 was not a cautious estimate
# that turned out low; it was an estimate that was never consulted.
#
# Measured, not inferred. Three searches run alone on the node, 40 threads,
# `/usr/bin/time -v`, each a self-comparison (the worst case: every long
# contig aligns against itself at full length):
#
#     Species0 spades55    1.61 GB  1,072,398 seqs  >10kb=19,387  143.0 GiB
#     Species3 trinity     1.61 GB  1,902,240 seqs  >10kb= 2,139  101.9 GiB
#     Species2 transabyss  1.18 GB  1,325,909 seqs  >10kb=   471   51.2 GiB
#
# Three things fall out, in order of how much they matter.
#
# **Cost is quadratic in query size.** Species2 to Species3 is 1.36x the
# bytes for 1.99x the memory -- an exponent of 2.22, and 2.19 after
# correcting for their different tails. That is not a curve fit looking for
# a shape: diamond's work goes as query x database, so a self-comparison is
# q^2, and the measured exponent agrees with the mechanism. It is the only
# figure here with support from something other than three points.
#
# **Sequence count is not it.** Species3 carries 77% more sequences than
# Species0 at the same file size and costs 29% less. Count was the obvious
# candidate and it is dead.
#
# **The long-contig tail is a real modifier, and only a modifier.** Same
# bytes, 9.1x the contigs over 10 kb, 29% more memory. It is not in the
# formula below: `40 * GB^2 + 0.0025 * (contigs > 10kb)` fits all three
# within 5-11%, but the tail needs a pass over the assemblies to count, and
# two parameters on three points is how the last two guesses here went
# wrong. It goes in when there is a fourth and fifth measurement, not
# before.
#
# So: GiB per GB-of-query squared, sized off the largest search input.
# Against the three measurements it predicts +5%, +48% and +58% -- always
# conservative, tightest on the expensive one, which is the right way round.
#
# What this replaced was a linear 96 GB per GB, fitted at ~1.5 GB and
# accurate only there: +8% on Species0, +121% on Species2. Worse, it fell
# the dangerous way as inputs grew. At 3 GB assemblies linear reads 288 GiB
# and would have planned two concurrent searches that each want 522 -- the
# same OOM this whole constant exists to prevent, arrived at by trusting a
# straight line outside the range it was fitted in.
ORTHOFINDER_GB_PER_SEARCH_GB2 = 58

# Floor under the above, for inputs small enough that the linear term says
# less than one diamond's fixed cost. diamond's default block size is -b2.0
# and it "can be expected to use roughly six times this number of memory (in
# GB)", so ~12 GB before a single seed hit is stored.
ORTHOFINDER_GB_PER_SEARCH_FLOOR = 12

# OrthoFinder's own ceiling on -a: its documented default is "16 or t/8
# (whichever lower)". Worth keeping, because the rest of that default is not
# usable here -- see orthofinder_analysis_threads.
ORTHOFINDER_MAX_ANALYSIS = 16

# Lowering OrthoFinder's -t lowers the core count with it: OrthoFinder hands
# every diamond `-p 1` whatever -t says, so -t 4 on a 40-core node runs four
# diamonds on four cores and leaves thirty-six idle. That is why -t was never
# lowered -- the only way to make this run was to make it slow.
#
# It is not the only way. Concurrency and core count are separate knobs in
# diamond and only OrthoFinder ties them together, so ensure_diamond_program()
# unties them: a `diamond_orp_<threads>` entry in OrthoFinder's config.json,
# a copy of its own diamond entry with `-p` set to `cpu // searches`, chosen
# with `-S`. Four diamonds at ten threads each is the same forty cores as
# sixteen at one, at a quarter of the peak memory.
#
# This is not a nicety. The measured search is 14.8 core-hours, so at `-p 1`
# it is 14.8 hours of wall time on one core; four waves of that is 59 hours
# before the cheap searches are counted. With the threads it is the same
# forty cores throughout and the step is hours, not days. Untying -p from -t
# is what makes a memory-safe concurrency affordable at all.
#
# PATH was tried first and cannot work. A shim ahead of the real diamond is
# overtaken by OrthoFinder itself, which prepends its environment's bin and
# its own bundled bin at startup: measured on the node, the shim sat at
# position 9 of a PATH whose first entry was the real diamond, and the
# searches ran `-p 1` with the shim present and unused. There is no
# `--config` flag either, so the install copy of config.json is the only
# place this can be said from.

# Floor for a self-comparison in check_orthofinder_searches: the fraction of
# an assembly's sequences that must show up as queries in its own Blast{i}_i.
# Every sequence aligns to itself, so the honest value is ~1.0 less whatever
# diamond's low-complexity masking removes from seeding entirely. See that
# method for why 0.5 is both far above a dead search and far below a real one.
BLAST_SELF_HIT_FLOOR = 0.5


def line_buffer_stdio():
    """Make our own output appear where it happened in a redirected log.

    Python block-buffers stdout in 4-8 KB chunks when it is not a terminal,
    which on a cluster it never is. Every tool we launch, though, inherits
    the same file descriptor and writes to it directly, unbuffered. So the
    pipeline's own narrative -- the banner, the `=== step -- start ===`
    lines, the `+ <command>` echoes, the retry warnings -- sits in our
    buffer while hours of OrthoFinder and pytransrate output stream past it,
    and only lands when the buffer happens to fill.

    It is not a cosmetic problem. On the 380C_0C5D_001F run the banner and
    the whole of the first stage appeared *after* OrthoFinder's 15:54:51
    output even though run_filtershort's own timestamp says 15:53:16, which
    makes a log read as though steps ran in an order they did not, and puts
    a failure's explanation somewhere other than next to the failure.

    Line buffering also guarantees we have flushed before a child we spawn
    writes anything, so the interleaving is right and not merely closer.

    Not sys.stdout.reconfigure(): that is 3.7+, and the cluster launches
    this under the system python3, which is 3.6.8. detach() rather than
    wrapping .buffer directly, so replacing sys.stdout cannot leave the old
    wrapper to close the descriptor out from under the new one when it is
    collected.
    """
    for name in ("stdout", "stderr"):
        stream = getattr(sys, name, None)
        if stream is None:
            continue
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(line_buffering=True)
            continue
        if not hasattr(stream, "detach") or not hasattr(stream, "encoding"):
            continue  # already replaced by something that isn't a text stream
        encoding, errors = stream.encoding, stream.errors
        try:
            detached = stream.detach()
        except (AttributeError, ValueError):
            continue
        setattr(sys, name, io.TextIOWrapper(
            detached, encoding=encoding, errors=errors, line_buffering=True))


def awk_first_field(src: Path, dst: Path) -> None:
    with open(src) as inf, open(dst, "w") as outf:
        for line in inf:
            fields = line.split()
            outf.write((fields[0] if fields else "") + "\n")


def extract_gene_ids(diamond_txt: Path) -> set:
    ids = set()
    with open(diamond_txt) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            parts = cols[1].split("|")
            if len(parts) < 3:
                continue
            ids.add(parts[2].split("_")[0])
    return ids


def parse_unique_count(diamond_txt: Path) -> int:
    return len(extract_gene_ids(diamond_txt))


def write_sorted(path: Path, ids) -> None:
    with open(path, "w") as f:
        for i in sorted(ids):
            f.write(i + "\n")


def human_size(nbytes) -> str:
    size = float(nbytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024


def path_size(path: Path) -> int:
    """Bytes on disk under `path`, whether it's one file or a whole tree."""
    path = Path(path)
    if not path.is_dir():
        try:
            return path.lstat().st_size
        except OSError:
            return 0
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += (Path(root) / name).lstat().st_size
            except OSError:
                pass
    return total


def is_gzip(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


#: The 28-byte empty BGZF block that closes every complete BAM (SAM spec
#: 4.1). A BAM not ending in it was still being written when whatever was
#: writing it stopped.
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003" + "00" * 9)

#: Columns in a salmon quant.sf: Name, Length, EffectiveLength, TPM,
#: NumReads. pytransrate rejects any other count as a version mismatch.
QUANT_COLUMNS = 5


def version_below(version: str, minimum: str) -> bool:
    """Is ``version`` older than ``minimum``, comparing release numbers only?

    Numeric components, left to right, shorter padded with zeros, so 2.10.0
    beats 2.9.0 rather than losing to it as a string compare would have it.

    A pre-release suffix is deliberately ignored: 2.2.2.dev1 compares equal to
    2.2.2 rather than below it, because those builds are how a fix is tested
    on the cluster before it is tagged, and a check that rejected them would
    make the version gate an obstacle to the work it exists to support. The
    exact string, suffix and all, is what gets printed -- the gate catches an
    install that is plainly old, and the printed version answers everything
    finer-grained than that.
    """
    def parts(text):
        numbers = re.match(r"\d+(?:\.\d+)*", text)
        return [int(n) for n in numbers.group(0).split(".")] if numbers else []

    have, want = parts(version), parts(minimum)
    if not have or not want:
        return False
    width = max(len(have), len(want))
    have += [0] * (width - len(have))
    want += [0] * (width - len(want))
    return have < want


def bam_is_complete(path: Path) -> bool:
    """Whether `path` is a BAM that was written all the way to the end.

    The same check pytransrate 2.2.0 makes before reusing one, made here
    because the decision to *keep* a BAM is ORP's: a truncated one is
    hundreds of gigabytes that no attempt can use, and 2.1.0 -- which ORP
    pinned until 4.0.0 -- would reuse it without asking.
    """
    try:
        if path.stat().st_size < len(BGZF_EOF):
            return False
        with open(path, "rb") as fh:
            fh.seek(-len(BGZF_EOF), os.SEEK_END)
            return fh.read(len(BGZF_EOF)) == BGZF_EOF
    except OSError:
        return False


def count_sequences(path: Path) -> int:
    """Records in a fasta, counted in binary chunks.

    Chunks rather than lines because this runs on merged.fasta, which is
    four assemblies concatenated -- millions of records and gigabytes of
    sequence -- and the answer is only wanted to compare against a row
    count. The one-byte `tail` carries a chunk boundary that falls between
    the newline and the ">"; seeding it with a newline is what makes the
    header on the very first line count.
    """
    opener = gzip.open if is_gzip(path) else open
    total = 0
    tail = b"\n"
    with opener(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            total += (tail + chunk).count(b"\n>")
            tail = chunk[-1:]
    return total


def quant_sf_is_complete(path: Path, expected: int) -> bool:
    """Whether a salmon quant.sf holds one whole row per contig.

    salmon writes quant.sf in a single pass at the end of its run, so a run
    killed during that pass leaves a short file -- and pytransrate reuses
    quant.sf on existence alone, without counting rows, at every version to
    date including the 2.2.0 that does check its BAM. That would not fail
    loudly: it would quietly score the assembly off whichever contigs made
    it into the file before the kill, and those scores are what
    pick_best_contigs.py then selects on.

    Both halves are needed. The row count catches a file cut at a line
    boundary; the width of the last row catches the commoner case of a file
    cut in the middle of one. A short count also catches a quant.sf left
    over from a different assembly.
    """
    try:
        rows = 0
        last = ""
        with open(path) as fh:
            next(fh, None)  # header, skipped as load_expression skips it
            for line in fh:
                if line.strip():
                    rows += 1
                    last = line
    except OSError:
        return False
    return (
        rows == expected
        and len(last.rstrip("\n").split("\t")) == QUANT_COLUMNS
    )


def read_length_stats(path: Path, n_records: int = 1000):
    """(mean, max) read length over the first n_records reads."""
    opener = gzip.open if is_gzip(path) else open
    lengths = []
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i >= n_records * 4:
                break
            if i % 4 == 1:
                lengths.append(len(line.strip()))
    if not lengths:
        return 0, 0
    return sum(lengths) // len(lengths), max(lengths)


def average_read_length(path: Path, n_records: int = 100) -> int:
    return read_length_stats(path, n_records)[0]


MAX_SPADES_K = 127  # rnaSPAdes requires every k to be odd and strictly below 128


def parse_kmer_spec(value):
    """Parse a --spadesN-kmer value into None ("auto") or a list of k-mer sizes.

    None makes run_spades() omit -k entirely, so rnaSPAdes picks its own two
    k values from the observed read length -- approximately 1/3 and 1/2 of the
    maximum. That is the documented default and the configuration the
    rnaSPAdes authors recommend: they warn that smaller k-mer sizes typically
    produce chimeric transcripts, and ORP forcing a single k per run has
    always deviated from it. Two single-k runs are not one two-k run.

    rnaSPAdes requires every k to be odd, below 128, and in ascending order,
    so reject violations here rather than after the reads have already been
    trimmed and corrected.
    """
    if value.strip().lower() == "auto":
        return None
    try:
        kmers = [int(x) for x in value.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{value!r}: expected 'auto' or a comma-separated list of integers"
        )
    for k in kmers:
        if k % 2 == 0 or not 0 < k <= MAX_SPADES_K:
            raise argparse.ArgumentTypeError(
                f"{value!r}: rnaSPAdes k-mer sizes must be odd and less than 128 (got {k})"
            )
    if kmers != sorted(kmers) or len(set(kmers)) != len(kmers):
        raise argparse.ArgumentTypeError(
            f"{value!r}: rnaSPAdes k-mer sizes must be distinct and in ascending order"
        )
    return kmers


def hostname_suffix() -> str:
    parts = socket.gethostname().split(".")
    return ".".join(parts[2:5])


class Pipeline:
    # What a run calls itself in its quality report. chowder.py overrides it:
    # a merge-only run writes a qualreport in exactly the same layout as a
    # full one, so without this the file is the one place a finished run
    # can't be told apart from an ORP that ran its own assemblers.
    RUN_DESCRIPTION = "the ORP"

    def __init__(self, args):
        self.dir = Path(args.dir).resolve() if args.dir else Path.cwd()
        self.makedir = HERE
        self.cpu = args.cpu
        self.busco_threads = args.busco_threads or self.cpu
        self.mem = args.mem
        # Assembler-only settings. An entry point that doesn't assemble
        # (chowder.py) has no flags for these and never reads them back.
        self.spades1_kmer = getattr(args, "spades1_kmer", None)
        self.spades2_kmer = getattr(args, "spades2_kmer", [75])
        self.transabyss_kmer = getattr(args, "transabyss_kmer", 32)
        self.read1 = Path(args.read1)
        self.read2 = Path(args.read2)
        self.runout = args.runout
        self.lineage = args.lineage
        self.strand = getattr(args, "strand", "")
        self.normalize_reads = getattr(args, "normalize_reads", False)
        self.tpm_filt = args.tpm_filt
        self.max_parallel = max(1, args.max_parallel)
        self.keep_intermediates = args.keep_intermediates
        # Appended to both pytransrate invocations. shlex so a value can be
        # quoted, and so the flags arrive as separate argv entries rather
        # than one string pytransrate would reject.
        self.pytransrate_args = shlex.split(getattr(args, "pytransrate_args", "") or "")
        self.orthofinder_searches = getattr(args, "orthofinder_searches", None) or 0
        self.orthofinder_analysis = getattr(args, "orthofinder_analysis", None) or 0

        # Everything from run_filtershort onwards works on "the assemblies"
        # rather than on four named assemblers, so a caller that brings its
        # own (chowder.py) supplies them here and shares the whole merge
        # half. It gives one order and means it for all three uses; only
        # oyster.py's own four carry the historical split between them (see
        # ASSEMBLY_ORDER / DIAMOND_PRIORITY / REPORT_ORDER).
        supplied = getattr(args, "assemblies", None)
        if supplied:
            self.assemblies = tuple(supplied)
            self.diamond_priority = self.assemblies
            self.report_order = self.assemblies
        else:
            self.assemblies = ASSEMBLY_ORDER
            self.diamond_priority = DIAMOND_PRIORITY
            self.report_order = REPORT_ORDER

        self.version = (self.makedir / "version.txt").read_text().strip()
        self.busco_config = self.makedir / "software" / "config.ini"
        os.environ["BUSCO_CONFIG_FILE"] = str(self.busco_config)
        self.diamond_db = self.makedir / "software" / "diamond" / "swissprot"

        self.rcorr_dir = self.dir / "rcorr"
        self.assemblies_dir = self.dir / "assemblies"
        self.assemblies_working = self.assemblies_dir / "working"
        self.diamond_dir = self.assemblies_dir / "diamond"
        self.reports_dir = self.dir / "reports"
        self.orthofuse_dir = self.dir / "orthofuse" / self.runout
        self.orthofuse_working = self.orthofuse_dir / "working"
        # OrthoFinder gets its own copies of the filtered assemblies rather
        # than the originals -- see write_search_inputs.
        self.orthofuse_search = self.orthofuse_dir / "search"
        self.quants_dir = self.dir / "quants"

        self.timing_log = self.reports_dir / f"{self.runout}.timing.log"
        # Quoted, so the line a run prints and logs is the line you can paste
        # back to repeat it: a read path with a space in it is otherwise
        # recorded as two arguments.
        self.run_cmd = " ".join(
            shlex.quote(a) for a in [Path(sys.argv[0]).name] + sys.argv[1:]
        )
        self.steps = []
        self._timing_lock = threading.Lock()
        # Background gzip of the files a finished run keeps -- see
        # compress_async(). Two workers is enough: the six files queued over
        # a run are queued in pairs, and this is meant to run *beside* an
        # assembler, not to compete with one.
        self._compress_pool = concurrent.futures.ThreadPoolExecutor(
            max_workers=2, thread_name_prefix="compress")
        self._compressions = {}
        self._compress_lock = threading.Lock()
        self._compress_cmd = None

    # -- process helpers -------------------------------------------------

    def run(self, cmd, cwd=None, retries=STEP_RETRIES, retry_delay=STEP_RETRY_DELAY,
            retry_cleanup=None, **kwargs):
        """retry_cleanup: what to clear before each retry, for tools (SPAdes,
        TransAByss) that refuse to reuse a non-empty output dir rather than
        resuming, so a bare retry would just fail differently instead of
        actually re-attempting the work. Not needed for tools like Trinity
        that resume from their own checkpoints in-place.

        A path or iterable of paths is rmtree'd. A callable is called
        instead, for a step where "clear the output directory" is too blunt
        and something in it has to survive the retry -- see
        clear_transrate_outdir.
        """
        printable = " ".join(str(c) for c in cmd)
        for attempt in range(retries + 1):
            print("+", printable)
            try:
                subprocess.run(cmd, check=True, cwd=str(cwd or self.dir), **kwargs)
                return
            except subprocess.CalledProcessError as e:
                if attempt == retries:
                    raise
                print(
                    f"*** step failed (exit {e.returncode}), retrying in {retry_delay}s "
                    f"[attempt {attempt + 2}/{retries + 1}] ***"
                )
                if callable(retry_cleanup):
                    retry_cleanup()
                elif retry_cleanup is not None:
                    paths = [retry_cleanup] if isinstance(retry_cleanup, (str, Path)) else retry_cleanup
                    for path in paths:
                        shutil.rmtree(path, ignore_errors=True)
                time.sleep(retry_delay)

    def conda_run(self, conda_env, *cmd, **kwargs):
        """`conda run -n conda_env -- cmd`, with kwargs going on to subprocess.

        The first parameter is `conda_env` and not `env` because `env` is
        subprocess's own name for the environment block, and a caller that
        wants to set one would otherwise be handing this method two values
        for the same parameter -- a TypeError raised at the call, hours into
        a run. (The caller that wanted one, run_orthofuser, no longer does:
        setting PATH from out here could not beat OrthoFinder's own
        prepending. The name stays right regardless.)
        """
        self.run(["conda", "run", "--no-capture-output", "-n", conda_env,
                  *[str(c) for c in cmd]], **kwargs)

    # -- background compression / reclaim -----------------------------------

    def _rel(self, path):
        try:
            return str(Path(path).relative_to(self.dir))
        except ValueError:
            return str(path)

    def _resolve_compressor(self):
        """argv prefix that writes a gzip stream of its file argument to stdout.

        pigz wherever one is available: these jobs are hidden behind a stage
        that already owns most of the machine, but a corrected read pair is
        tens of GB and single-threaded gzip can still be running long after
        the stage that was meant to hide it has ended. The thread count is
        deliberately small for the same reason -- this is background work,
        not a stage of its own. Falls back to plain gzip so an `orp` env
        built before pigz was added to orp_env.yml keeps working.
        """
        if self._compress_cmd is None:
            threads = str(max(1, min(4, self.cpu // 8)))
            if shutil.which("pigz"):
                self._compress_cmd = ["pigz", "-c", "-p", threads]
            elif self.which_in_env("orp", "pigz"):
                self._compress_cmd = ["conda", "run", "--no-capture-output", "-n", "orp",
                                      "pigz", "-c", "-p", threads]
            else:
                self._compress_cmd = ["gzip", "-c"]
            print(f"[compress] using: {' '.join(self._compress_cmd)}")
        return self._compress_cmd

    def compress_async(self, path):
        """Queue `path` for gzipping in the background, original left in place.

        Called when a file stops being *written*, which is much earlier than
        when it stops being *read*: the corrected reads feed every assembler
        and every alignment step through to `strandeval`, and the four
        assemblies are read again at `run_filtershort`, `diamond_*` and
        `posthack`. So the .gz is built alongside the original while the
        pipeline runs, and cleanup() at the end only has to unlink -- which
        is what puts the compression cost in parallel with an assembler
        instead of on the end of the run, where it would be pure added wall
        time.

        A no-op if the .gz is already at least as new as its source, so a
        resumed run doesn't recompress work the previous one finished.
        """
        path = Path(path)
        if self.keep_intermediates or not path.exists():
            return
        gz = path.with_suffix(path.suffix + ".gz")
        if gz.exists() and gz.stat().st_mtime >= path.stat().st_mtime:
            return
        with self._compress_lock:
            if path in self._compressions:
                return
            self._resolve_compressor()
            self._compressions[path] = self._compress_pool.submit(self._compress, path, gz)

    def _compress(self, src, gz):
        # _compress_cmd was resolved by compress_async() under _compress_lock
        # before this job was submitted, so workers only ever read it.
        part = gz.with_name(gz.name + ".part")
        start = time.time()
        try:
            with open(part, "wb") as out:
                subprocess.run(self._compress_cmd + [str(src)], check=True,
                               stdout=out, cwd=str(self.dir))
            part.replace(gz)
        except BaseException:
            # Never leave a truncated .gz behind that a later run would
            # mistake for a finished one on mtime alone.
            if part.exists():
                part.unlink()
            raise
        print(f"[compress] {self._rel(gz)} written in {int(time.time() - start)}s "
              f"({human_size(path_size(src))} -> {human_size(path_size(gz))})", flush=True)

    def compression_done(self, path) -> bool:
        """Block until `path`'s .gz is finished; True only if it really is.

        The precondition for deleting the uncompressed original. A failed
        compression is reported and returns False rather than raising --
        losing a background gzip is a reason to keep the plain file, not to
        fail a run whose actual work is already done.
        """
        path = Path(path)
        future = self._compressions.get(path)
        if future is not None:
            try:
                future.result()
            except Exception as e:
                print(f"*** compressing {self._rel(path)} failed ({e}); "
                      "keeping the uncompressed file ***")
                return False
        gz = path.with_suffix(path.suffix + ".gz")
        return gz.exists() and gz.stat().st_size > 0

    def finish_compression(self):
        """Join the background pool. Called on every exit path, including failure."""
        pending = [p for p, f in self._compressions.items() if not f.done()]
        if pending:
            print("\n=== waiting on background compression: "
                  + ", ".join(self._rel(p) for p in pending) + " ===", flush=True)
        self._compress_pool.shutdown(wait=True)

    def reclaim_trimmed_reads(self):
        """Delete the trimmed-but-uncorrected reads, once rcorrector has read them.

        Trimmomatic writes four files (both paired mates and both unpaired
        ones) and nothing past run_rcorrector ever opens any of them again --
        every assembler and every alignment step reads the corrected pair.
        They're the same order of magnitude as the raw input, so this is both
        the largest reclaim in the run and the earliest one available, which
        is why it doesn't wait for cleanup() at the end.

        The sentinel is what keeps this from costing a resumed run a re-trim:
        main() asks for the TRIM files back as outputs only while the
        corrected pair is missing or stale.
        """
        if self.keep_intermediates:
            return
        freed = 0
        for suffix in ("1P", "2P", "1U", "2U"):
            p = self.rcorr_dir / f"{self.runout}.TRIM_{suffix}.fastq"
            if p.exists():
                freed += path_size(p)
                p.unlink()
        (self.rcorr_dir / f"{self.runout}.trim.done").touch()
        if freed:
            print(f"[cleanup] reclaimed {human_size(freed)} of trimmed reads "
                  f"(rcorr/{self.runout}.TRIM_*.fastq); the corrected pair is what "
                  "everything downstream reads")

    def run_inputs(self):
        """The files a finished run is checked for staleness against.

        The raw read pair for oyster.py; chowder.py adds the assemblies it
        was handed, since replacing one of those is as much a new run as
        replacing the reads.
        """
        return [self.read1, self.read2]

    def already_complete(self) -> bool:
        """True if a previous run finished *and* cleanup() already ran on it.

        needs_run() decides each step from its outputs, and cleanup deletes
        most of those -- so without this guard, re-invoking oyster.py on a
        finished, cleaned run directory (a resubmitted cluster job, say)
        would find nearly every stage 'out of date' and quietly reassemble
        from scratch over a completed run, where today it no-ops. The marker
        plus a .ORP.fasta still newer than the raw reads is the whole
        condition; anything else (new reads, a deleted assembly) falls
        through to the normal resume path.
        """
        marker = self.reports_dir / f"{self.runout}.cleanup.done"
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        if not marker.exists() or self.needs_run([orp_fasta], self.run_inputs()):
            return False
        print(f"\n=== {self.runout} already finished, and its intermediate files "
              "have been reclaimed ===")
        print(f"    assembly:  {orp_fasta}")
        print(f"    reports:   {self.reports_dir}")
        print(f"    reclaimed: {marker}")
        print("\n    Nothing to do. Assemble these reads again under a different "
              "--runout/--dir,\n    or delete the marker above to force a full "
              "re-run in place.\n")
        return True

    def cleanup(self):
        """Reclaim everything a finished run doesn't need any more.

        What survives: reports/, the final .ORP.fasta, the four individual
        assemblies and the corrected read pair -- the last two as the .gz
        compress_async() has been building in the background since each was
        written, so this step only unlinks. A file whose compression didn't
        finish is kept uncompressed instead of being deleted.

        Everything removed here is reproducible from what's kept: the
        orthofuse tree (OrthoFinder's all-vs-all output plus pytransrate's
        scoring of the pooled fasta -- normally the largest directory in the
        run), the diamond hits and the list1-list7 set algebra built from
        them, the salmon index and quantification, and the chain of working
        assemblies between orthofusing and .ORP.fasta. Every number any of
        it contributed is already in reports/qualreport.<run>.
        """
        if self.keep_intermediates:
            print("[cleanup] --keep-intermediates given; leaving intermediates in place")
            return

        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        freed = 0
        kept = [f"{self._rel(orp_fasta)}  (the assembly)",
                f"{self._rel(self.reports_dir)}/  (all reports)"]
        removed = []

        for src in [self.cor1(), self.cor2()] + self.assembly_fasta_paths():
            gz = src.with_suffix(src.suffix + ".gz")
            if src.is_symlink():
                # Not ours to reclaim: chowder.py points the corrected pair
                # straight at the user's own reads under
                # --reads-are-corrected, and compressing or unlinking those
                # is not what "reclaim this run's intermediates" means.
                kept.append(f"{self._rel(src)}  (symlink to a file this run did not create)")
                continue
            if src.exists() and self.compression_done(src):
                freed += path_size(src)
                src.unlink()
            elif src.exists():
                kept.append(f"{self._rel(src)}  (left uncompressed -- gzip did not finish)")
                continue
            if gz.exists():
                kept.append(f"{self._rel(gz)}  ({human_size(path_size(gz))})")

        for path in (
            self.dir / "orthofuse",
            self.assemblies_working,
            self.diamond_dir,
            self.quants_dir,
            # Trinity's --full_cleanup normally removes this itself; a run
            # that was interrupted and resumed can still leave it behind.
            self.trinity_out_dir(),
            self.assemblies_dir / f"{self.runout}.orthomerged.fasta",
            self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta",
            self.assemblies_dir / f"{self.runout}.ORP.diamond.txt",
            self.assemblies_dir / f"{self.runout}.flagstat",
            self.assemblies_dir / f"{self.runout}.filter.done",
            self.trinity_phase1_done(),
        ):
            if not path.exists():
                continue
            size, is_dir = path_size(path), path.is_dir()
            if is_dir:
                shutil.rmtree(path, ignore_errors=True)
            else:
                path.unlink()
            freed += size
            removed.append(f"{self._rel(path)}{'/' if is_dir else ''}  ({human_size(size)})")

        lines = [f"Command: {self.run_cmd}", "",
                 f"Reclaimed {human_size(freed)} of intermediate files "
                 f"at {self._ts()}.", "", "kept:"]
        lines += [f"  {k}" for k in kept]
        lines += ["", "removed:"]
        lines += [f"  {r}" for r in removed] or ["  (nothing left to remove)"]
        text = "\n".join(lines) + "\n"
        (self.reports_dir / f"{self.runout}.cleanup.done").write_text(text)
        print("\n" + text)

    # -- resumability ------------------------------------------------------

    def needs_run(self, outputs, inputs) -> bool:
        if not outputs:
            return True
        outputs = [Path(o) for o in outputs]
        if not all(o.exists() for o in outputs):
            return True
        existing_inputs = [Path(i) for i in inputs if i and Path(i).exists()]
        if not existing_inputs:
            return False
        newest_input = max(i.stat().st_mtime for i in existing_inputs)
        oldest_output = min(o.stat().st_mtime for o in outputs)
        return newest_input > oldest_output

    def tool_version(self, env, binary):
        """``<binary> --version`` in ``env``, as a bare version string or None.

        Tools print anything from "2.2.1" to "pytransrate 2.2.1" to a banner,
        so the first thing on the first line that looks like a version is what
        comes back. None when the tool cannot be run or says nothing usable,
        which every caller here treats as "cannot tell" rather than as bad
        news -- see check_pytransrate_version.
        """
        try:
            result = subprocess.run(
                ["conda", "run", "-n", env, binary, "--version"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
            )
        except OSError:
            return None
        text = (result.stdout.strip() or result.stderr.strip()).splitlines()
        if not text:
            return None
        match = re.search(r"\d+(?:\.\d+)+(?:\.?(?:dev|a|b|rc)\d*)?", text[0])
        return match.group(0) if match else None

    def stamp_tool_version(self, env, binary, stamp):
        """Record a tool's version beside its artifacts; return the stamp path.

        needs_run only ever compares mtimes, so an artifact that is still
        newer than its inputs looks up to date even when the tool that has to
        read it can no longer do so. salmon 2.7.0 is exactly that case: it
        rejects any index built by an earlier salmon, so on a resumed run a
        stale <run>.ortho.idx would be kept, salmon_index skipped, and
        salmon quant left to fail against an index it cannot read.

        Declaring this stamp as an input turns a version change into an
        ordinary out-of-date input, which is the machinery every other step
        already uses. The file is rewritten only when the version actually
        changes -- rewriting unconditionally would force a rebuild on every
        run.
        """
        try:
            result = subprocess.run(
                ["conda", "run", "-n", env, binary, "--version"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
            )
        except FileNotFoundError:
            return stamp
        version = result.stdout.strip() or result.stderr.strip()
        if not version:
            # Can't tell. Leave any existing stamp alone rather than writing a
            # placeholder that would itself look like a version change later.
            return stamp
        if not stamp.exists() or stamp.read_text() != version:
            stamp.parent.mkdir(parents=True, exist_ok=True)
            stamp.write_text(version)
        return stamp

    def _ts(self, t=None):
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t if t is not None else time.time()))

    def step(self, name, outputs, inputs, func, timed=True):
        if not self.needs_run(outputs, inputs):
            print(f"[{name}] up to date, skipping")
            return
        start = time.time()
        print(f"\n=== {name} -- start {self._ts(start)} ===")
        func()
        elapsed = int(time.time() - start)
        print(f"=== {name} -- done {self._ts()} ({elapsed}s) ===")
        if timed:
            self._record_timing(name, elapsed, start)

    def _record_timing(self, name, elapsed, start):
        with self._timing_lock:
            self.steps.append((name, elapsed))
            with open(self.timing_log, "a") as f:
                f.write(f"{name}\t{elapsed}\t{self._ts(start)}\n")

    def run_parallel(self, jobs, max_workers=2):
        """Run independent (name, outputs, inputs, func) jobs concurrently.

        Each pending job's func is called as func(cpu=<split>, mem=<split>),
        where self.cpu/self.mem are divided across however many jobs actually
        run at once (capped at max_workers). Already up-to-date jobs are
        skipped and don't count toward the split.
        """
        pending = []
        for name, outputs, inputs, func in jobs:
            if self.needs_run(outputs, inputs):
                pending.append((name, func))
            else:
                print(f"[{name}] up to date, skipping")
        if not pending:
            return

        # Longest-processing-time-first: submit the slowest known jobs first so
        # they start immediately instead of waiting behind quick ones, using
        # a fixed reference profile (STEP_TIME_HINTS) rather than this run's
        # own history -- each dataset is normally only ever assembled once,
        # so there's no prior run of *this* data to learn from. Jobs with no
        # hint keep their given relative order, after every hinted job.
        pending.sort(key=lambda item: STEP_TIME_HINTS.get(item[0], -1), reverse=True)

        workers = min(max_workers, len(pending), max(1, self.cpu))
        job_cpu = max(1, self.cpu // workers)
        job_mem = max(1, self.mem // workers)
        if workers > 1:
            print(f"\n=== running {len(pending)} step(s), {workers} at a time (cpu={job_cpu}, mem={job_mem}G each) ===")

        def run_one(name, func):
            start = time.time()
            print(f"\n=== {name} (cpu={job_cpu}, mem={job_mem}G) -- start {self._ts(start)} ===")
            func(cpu=job_cpu, mem=job_mem)
            elapsed = int(time.time() - start)
            print(f"=== {name} -- done {self._ts()} ({elapsed}s) ===")
            self._record_timing(name, elapsed, start)

        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            futures = [ex.submit(run_one, name, func) for name, func in pending]
            for future in concurrent.futures.as_completed(futures):
                future.result()

    # -- setup / preflight -------------------------------------------------

    def setup(self):
        for d in (
            self.assemblies_dir, self.rcorr_dir, self.reports_dir,
            self.orthofuse_dir, self.quants_dir, self.diamond_dir, self.assemblies_working,
        ):
            d.mkdir(parents=True, exist_ok=True)

    def timing_init(self):
        self.reports_dir.mkdir(parents=True, exist_ok=True)
        self.timing_log.write_text(f"Command: {self.run_cmd}\n\n")

    def which_in_env(self, env, binary):
        try:
            result = subprocess.run(
                ["conda", "run", "-n", env, "which", binary],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
            )
        except FileNotFoundError:
            return None
        path = result.stdout.strip()
        if result.returncode == 0 and path and os.path.basename(path) == binary:
            return path
        return None

    def required_tools(self):
        """(env, binary, label) for every tool this entry point shells out to.

        Order is the order preflight prints them in. An entry point that
        doesn't assemble overrides this rather than demanding assemblers it
        will never run (see chowder.py).
        """
        return CHECK_TOOLS

    def check(self):
        """Verify every tool is present, saying nothing when they all are.

        A dozen "installed" lines is noise on every successful run; the only
        news preflight has is a tool that is missing.
        """
        for env, binary, label in self.required_tools():
            if not self.which_in_env(env, binary):
                sys.exit(f"*** {label} is not installed, must fix ***")
        self.check_pytransrate_version()
        self.log_provenance()

    def log_provenance(self):
        """Print what a later post-mortem needs and cannot recover.

        `sacct` is the only place a peak memory figure for a finished job
        exists, and it is keyed on a job ID that the log never carried --
        so the run that raised the memory question could not be asked about
        memory afterwards. The scheduler puts the ID in the environment;
        writing it down costs a line and is the difference between
        measuring a failure and arguing about it.

        The diamond version goes here for the same reason it matters: the
        one that runs OrthoFinder's all-vs-all comes in as an unpinned
        dependency of the orthofinder package, not from orp_env.yml, and
        several of the memory fixes in diamond's own ChangeLog land in
        specific versions (the hash join stage in 2.1.11, very long
        queries in 2.0.1). Which one is installed is not knowable from
        this repository.
        """
        host = socket.gethostname()
        job = os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOBID")
        array = os.environ.get("SLURM_ARRAY_JOB_ID")
        task = os.environ.get("SLURM_ARRAY_TASK_ID")
        print(f"[provenance] host {host}, pid {os.getpid()}")
        if job:
            label = f"{array}_{task}" if array and task else job
            print(f"[provenance] slurm job {label}")
            print(f"[provenance] after this run, peak memory is: "
                  f"sacct -j {job} --units=G "
                  f"--format=JobID,JobName%20,State,ExitCode,ReqMem,MaxRSS,MaxDiskWrite,Elapsed")
        else:
            print("[provenance] no SLURM_JOB_ID in the environment -- if this is a "
                  "batch job, peak memory will not be recoverable afterwards")
        diamond = self.tool_version("orp_orthofinder", "diamond")
        print(f"[provenance] diamond in orp_orthofinder: {diamond or 'unknown'} "
              "(unpinned -- it arrives as an orthofinder dependency)")

    def check_pytransrate_version(self):
        """Refuse to start on a pytransrate older than the pipeline needs.

        Present-and-runnable is the wrong question for this one tool: the
        version that matters is the difference between a 16-hour failure and
        a run that finishes, and nothing else in the pipeline notices which
        one is installed. See PYTRANSRATE_MIN_VERSION.

        The version is printed either way, because the second half of the
        problem is a fix that was installed and did not work being
        indistinguishable, in the log, from a fix that was never installed.
        The log now carries the answer at the top, before the hours.

        A version that cannot be read is not fatal. This check exists to
        catch a known-old install, not to become a new way for the run to
        refuse to start.
        """
        version = self.tool_version("orp", "pytransrate")
        if version is None:
            print("[preflight] could not read the pytransrate version; "
                  f"carrying on (this pipeline needs >= {PYTRANSRATE_MIN_VERSION})")
            return
        print(f"[preflight] pytransrate {version}")
        if version_below(version, PYTRANSRATE_MIN_VERSION):
            sys.exit(
                f"\n*** pytransrate {version} is installed and this pipeline "
                f"needs at least {PYTRANSRATE_MIN_VERSION}. ***\n\n"
                "    Older versions size the read-metrics step against the\n"
                "    whole machine rather than the memory budget, and delete\n"
                "    the BAM when a run fails -- so the failure costs a full\n"
                "    remap on every retry. Update the orp environment:\n\n"
                "      conda run -n orp pip install --upgrade --force-reinstall \\\n"
                "        --no-deps \\\n"
                "        'pytransrate @ git+https://github.com/macmanes-lab/"
                f"pytransrate.git@v{PYTRANSRATE_MIN_VERSION}'\n"
            )

    def welcome(self):
        print(RED)
        print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("                    _.-~~~~~-._")
        print("                 .-~   o   o   ~-.")
        print("                (   .-'~~~~~'-.   )        OYSTER RIVER PROTOCOL")
        print(f"                 '-.___________.-'         version {self.version}")
        print("                     '-.___.-'")
        print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + RESET + "\n")

    def readcheck(self):
        if not self.read1.exists():
            sys.exit("\n\n\n\n ERROR: YOUR READ1 FILE DOES NOT EXIST AT THE LOCATION YOU SPECIFIED\n\n\n\n ")
        if not self.read2.exists():
            sys.exit("\n\n\n\n ERROR: YOUR READ2 FILE DOES NOT EXIST AT THE LOCATION YOU SPECIFIED\n\n\n\n ")
        len1 = average_read_length(self.read1)
        len2 = average_read_length(self.read2)
        # Only k values we picked need checking against read length. An "auto"
        # assembly derives its k from the reads themselves, so it cannot be
        # too large by construction.
        explicit = [k for spec in (self.spades1_kmer, self.spades2_kmer) if spec for k in spec]
        if not explicit:
            return
        max_k = max(explicit)
        if not (len1 > max_k and len2 > max_k):
            sys.exit(
                f"\n\n\n\n IT LOOKS LIKE YOUR READS ARE NOT AT LEAST {max_k} BP LONG,\n "
                'PLEASE EDIT YOUR COMMAND USING THE "--spades1-kmer"/"--spades2-kmer" FLAGS,\n'
                " SETTING EACH ASSEMBLY KMER LENGTH TO AN ODD NUMBER LESS THAN YOUR READ LENGTH,\n"
                ' OR TO "auto" TO LET rnaSPAdes PICK FROM THE READS \n\n\n\n'
            )

    # -- trimming / correction ----------------------------------------------

    def trim1(self):
        return self.rcorr_dir / f"{self.runout}.TRIM_1P.fastq"

    def trim2(self):
        return self.rcorr_dir / f"{self.runout}.TRIM_2P.fastq"

    def cor1(self):
        return self.rcorr_dir / f"{self.runout}.TRIM_1P.cor.fq"

    def cor2(self):
        return self.rcorr_dir / f"{self.runout}.TRIM_2P.cor.fq"

    def run_trimmomatic(self):
        baseout = self.rcorr_dir / f"{self.runout}.TRIM.fastq"
        clip = self.makedir / "barcodes" / "barcodes.fa"
        pe_args = [
            "PE", "-threads", str(self.cpu), "-baseout", str(baseout),
            str(self.read1), str(self.read2),
            "LEADING:3", "TRAILING:3",
            f"ILLUMINACLIP:{clip}:2:30:10:8:TRUE", "MINLEN:25",
        ]
        if hostname_suffix() == "bridges.psc.edu":
            trimmomatic_home = os.environ.get("TRIMMOMATIC_HOME", "")
            jar = f"{trimmomatic_home}/trimmomatic-0.36.jar"
            self.run(["java", f"-Xmx{self.mem}G", "-jar", jar, *pe_args])
        else:
            env = dict(os.environ, _JAVA_OPTIONS=f"-Xmx{self.mem}G")
            self.run(
                ["conda", "run", "--no-capture-output", "-n", "orp", "trimmomatic", *pe_args],
                env=env,
            )

    def run_rcorrector(self):
        self.conda_run(
            "orp", "run_rcorrector.pl", "-t", self.cpu, "-k", "31",
            "-1", self.trim1(), "-2", self.trim2(), "-od", self.rcorr_dir,
        )

    # -- assemblers ----------------------------------------------------------

    def trinity_out_dir(self):
        return self.assemblies_dir / f"{self.runout}.trinity"

    def trinity_phase1_done(self):
        """Phase 1's completion sentinel, deliberately outside trinity_out_dir().

        Phase 2 passes --full_cleanup, which deletes the whole
        <run>.trinity/ working directory -- including
        recursive_trinity.cmds.ok, the file that used to serve as both
        phase 1's declared output and phase 2's declared input. So a
        finished Trinity erased its own evidence that it had run: on the
        next invocation needs_run() found phase 1's output missing and
        re-ran it (~90min), which rewrote cmds.ok *newer* than
        <run>.trinity.Trinity.fasta, which in turn made phase 2 look out of
        date and re-run against an assembly that was already complete
        (~34h). The window is a resumed run whose Stage B had succeeded and
        which then failed or was killed later -- i.e. a walltime kill near
        the end of a long run, which is exactly when a job gets
        resubmitted.
        """
        return self.assemblies_dir / f"{self.runout}.trinity.phase1.done"

    def seed_trinity_phase1_sentinel(self):
        """Back-fill the sentinel for run directories created before it existed.

        Without this, the very first resume after upgrading would hit the
        bug the sentinel exists to prevent, once. Seeds from whatever
        already proves phase 1 ran -- cmds.ok if Trinity's working dir is
        still there, otherwise the finished assembly -- and copies that
        file's mtime rather than stamping 'now', so the ordering needs_run()
        compares stays exactly what it was.
        """
        sentinel = self.trinity_phase1_done()
        if sentinel.exists():
            return
        for evidence in (self.trinity_out_dir() / "recursive_trinity.cmds.ok",
                         self.assemblies_dir / f"{self.runout}.trinity.Trinity.fasta"):
            if evidence.exists():
                sentinel.touch()
                st = evidence.stat()
                os.utime(sentinel, (st.st_atime, st.st_mtime))
                print(f"[resume] seeded {self._rel(sentinel)} from "
                      f"{self._rel(evidence)} (run directory predates it)")
                return

    def _trinity_base_cmd(self, cpu, mem):
        cmd = ["Trinity"]
        if self.strand == "RF":
            cmd += ["--SS_lib_type", "RF"]
        elif self.strand == "FR":
            cmd += ["--SS_lib_type", "FR"]
        cmd += ["--no_version_check", "--bypass_java_version_check"]
        if not self.normalize_reads:
            cmd.append("--no_normalize_reads")
        cmd += [
            "--seqType", "fq",
            "--output", str(self.trinity_out_dir()),
            "--max_memory", f"{mem}G",
            "--left", str(self.cor1()), "--right", str(self.cor2()),
            "--CPU", str(cpu), "--inchworm_cpu", "10",
        ]
        return cmd

    def run_trinity_phase1(self, cpu=None, mem=None):
        # Inchworm + Chrysalis prep only, via Trinity's documented staged-
        # execution flag (https://github.com/trinityrnaseq/trinityrnaseq/wiki/
        # Running-Trinity#running-trinity-in-multiple-sequential-stages):
        # builds the whole-transcriptome contig graph, partitions reads per
        # gene component, writes the Phase-2 command list, then stops --
        # doesn't touch the actual per-component assembly (run_trinity_phase2
        # below). Runs alongside SPAdes55/SPAdes75 (Stage A) at
        # TRINITY_PHASE1_SHARE of --cpu/--mem.
        cpu = self.cpu if cpu is None else cpu
        mem = self.mem if mem is None else mem
        cmd = self._trinity_base_cmd(cpu, mem) + ["--no_distributed_trinity_exec"]
        self.conda_run("orp_trinity", *cmd, retries=0)
        self.trinity_phase1_done().touch()

    def run_trinity_phase2(self, cpu=None, mem=None):
        # Same command, no stop flag: per the docs above, Trinity resumes
        # from Phase 1's on-disk checkpoints straight into Phase 2 -- the
        # thousands of small, independent, single-threaded per-component
        # assembly jobs dispatched via ParaFly, and by far Trinity's
        # dominant cost. Runs alongside Trans-ABySS (Stage B) at
        # TRINITY_PHASE2_SHARE of --cpu/--mem instead of waiting for it to
        # finish -- see TRINITY_PHASE2_SHARE above for why that pairing is
        # safe.
        cpu = self.cpu if cpu is None else cpu
        mem = self.mem if mem is None else mem
        out = self.assemblies_dir / f"{self.runout}.trinity.Trinity.fasta"
        cmd = self._trinity_base_cmd(cpu, mem) + ["--full_cleanup"]
        # No retries: this step's wall time dwarfs every other (hours to
        # days), so blindly retrying a deterministic failure could multiply
        # the wall time before finally giving up. It also resumes from its
        # own checkpoints in-place (same as Phase 1 above), so a manual
        # re-run of oyster.py after a transient failure loses little anyway.
        self.conda_run("orp_trinity", *cmd, retries=0)
        tmp = out.with_suffix(".fa")
        awk_first_field(out, tmp)
        tmp.replace(out)
        for f in self.assemblies_dir.glob("*gene_trans_map"):
            f.unlink()

    def run_spades(self, kmers, outname, workdir_suffix, cpu=None, mem=None):
        """kmers: a list of k-mer sizes, or None to leave -k off the command
        line so rnaSPAdes selects its own (see parse_kmer_spec)."""
        cpu = self.cpu if cpu is None else cpu
        mem = self.mem if mem is None else mem
        out = self.assemblies_dir / f"{self.runout}.{outname}.fasta"
        workdir = self.assemblies_dir / f"{self.runout}.spades_k{workdir_suffix}"
        cmd = ["rnaspades.py"]
        if self.strand == "RF":
            cmd.append("--ss-rf")
        elif self.strand == "FR":
            cmd.append("--ss-fr")
        cmd += [
            "--only-assembler", "-o", str(workdir),
            "--threads", str(cpu), "--memory", str(mem),
        ]
        if kmers:
            cmd += ["-k", ",".join(str(k) for k in kmers)]
        cmd += ["-1", str(self.cor1()), "-2", str(self.cor2())]
        self.conda_run("orp_spades", *cmd, retry_cleanup=workdir)
        shutil.move(str(workdir / "transcripts.fasta"), str(out))
        shutil.rmtree(workdir, ignore_errors=True)

    def run_spadesauto(self, cpu=None, mem=None):
        self.run_spades(self.spades1_kmer, "spadesauto", "auto", cpu=cpu, mem=mem)

    def run_spades75(self, cpu=None, mem=None):
        self.run_spades(self.spades2_kmer, "spades75", "75", cpu=cpu, mem=mem)

    def run_transabyss(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        out = self.assemblies_dir / f"{self.runout}.transabyss.fasta"
        workdir = self.assemblies_dir / f"{self.runout}.transabyss"
        cmd = ["transabyss"]
        if self.strand in ("RF", "FR"):
            cmd.append("--SS")
        cmd += [
            "--threads", str(cpu), "--outdir", str(workdir),
            "--kmer", str(self.transabyss_kmer), "--length", "250",
            "--name", f"{self.runout}.transabyss.fasta",
            "--pe", str(self.cor1()), str(self.cor2()),
        ]
        self.conda_run("orp_transabyss", *cmd, retry_cleanup=workdir)
        final = workdir / f"{self.runout}.transabyss.fasta-final.fa"
        awk_first_field(final, out)
        shutil.rmtree(workdir, ignore_errors=True)

    # -- orthofuse merge -------------------------------------------------------

    def assembly_fasta(self, assembly):
        return self.assemblies_dir / f"{self.runout}.{assembly.fasta_name}"

    def assembly_fasta_paths(self):
        return [self.assembly_fasta(a) for a in self.assemblies]

    def diamond_txt(self, assembly):
        return self.diamond_dir / f"{self.runout}.{assembly.diamond_label}.diamond.txt"

    def unique_txt(self, assembly):
        return self.diamond_dir / f"{self.runout}.unique.{assembly.unique_label}.txt"

    def short_fasta_paths(self):
        return [self.orthofuse_working / f"{self.assembly_fasta(a).name}.short.fasta"
                for a in self.assemblies]

    def run_filtershort(self):
        self.orthofuse_working.mkdir(parents=True, exist_ok=True)
        for a in self.diamond_priority:
            fasta = self.assembly_fasta(a)
            outp = self.orthofuse_working / f"{fasta.name}.short.fasta"
            self.conda_run("orp", "python", self.makedir / "scripts" / "long.seq.py", fasta, outp, "200")

    def search_fasta_paths(self):
        return [self.orthofuse_search / p.name for p in self.short_fasta_paths()]

    def write_search_inputs(self):
        """Copy the filtered assemblies for OrthoFinder, with N runs neutered.

        OrthoFinder's `-d` does not switch it to a nucleotide searcher: it
        runs `diamond blastp` over the DNA, which is why the command it
        builds carries `--ignore-warnings` -- that is what lets `makedb`
        accept ACGT as protein. Every base is therefore read as an amino
        acid, and `N` is asparagine.

        rnaSPAdes gap-fills scaffolds with runs of N; Trinity and TransAByss
        emit none. In the run this was written for that was 160,372 and
        189,448 contigs carrying a >=10bp N run, against zero and zero -- so
        each SPAdes assembly arrived with ~175K poly-asparagine tracts, and
        a poly-asparagine tract seeds against every other one in the
        database. diamond's memory went with it: the kernel killed five
        searches at 92-145 GB resident apiece, against the ~12 GB per
        process that the sizing model then assumed. Every one of the
        seven failed searches had a SPAdes assembly as its query; not one
        transabyss or trinity search failed.

        N becomes X, the unknown residue, which diamond will not seed on.
        That removes the tracts from the search without shortening a contig
        or renaming one, so the orthogroups still refer to the same IDs.
        Isolated Ns go with them: an ambiguous base carries no information
        to match on either, and translating the lot is both simpler and
        safer than deciding what counts as a run.

        These are written to their own directory and are not the assemblies
        that go on to be merged. merge() concatenates short_fasta_paths()
        into merged.fasta, which is what pytransrate scores and what
        orthofusing pulls the final sequence out of -- masking in place
        would edit the output assembly and invalidate a scoring run that
        takes twelve hours. Only the clustering sees an X.
        """
        self.orthofuse_search.mkdir(parents=True, exist_ok=True)
        table = bytes.maketrans(b"Nn", b"XX")
        for src, dst in zip(self.short_fasta_paths(), self.search_fasta_paths()):
            masked = 0
            tmp = dst.with_suffix(dst.suffix + ".partial")
            with open(src, "rb") as inf, open(tmp, "wb") as outf:
                for line in inf:
                    # Deflines are copied byte for byte: a contig whose name
                    # contains an N still has to answer to that name in
                    # Orthogroups.txt.
                    if not line.startswith(b">"):
                        masked += line.count(b"N") + line.count(b"n")
                        line = line.translate(table)
                    outf.write(line)
            tmp.replace(dst)
            print(f"    {dst.name}: {masked} N -> X")

    def orthofinder_search_plan(self, cpu, mem):
        """(concurrent diamonds, threads each, GB apiece) for the all-vs-all.

        Sized off the largest search input, because the searches run
        together and it is the biggest of them that decides when the node
        runs out: a plan that fits the mean fits nothing on the run where
        one assembly is twice its neighbours. See
        ORTHOFINDER_GB_PER_SEARCH_GB2 for where the figure comes from
        and why the cap it feeds had no effect before this.

        `--orthofinder-searches` overrides the memory term and nothing else.
        The thread count still follows from it, so pinning concurrency on a
        node whose memory this model has wrong does not also mean giving up
        the cores.
        """
        inputs = self.search_fasta_paths()
        jobs = max(1, len(inputs) ** 2)
        biggest = max((p.stat().st_size for p in inputs if p.is_file()), default=0)
        per_search = max(
            ORTHOFINDER_GB_PER_SEARCH_FLOOR,
            math.ceil(ORTHOFINDER_GB_PER_SEARCH_GB2 * (biggest / 1e9) ** 2),
        )
        if self.orthofinder_searches:
            searches = max(1, min(cpu, jobs, self.orthofinder_searches))
        else:
            searches = max(1, min(cpu, jobs, mem // per_search))
        return searches, max(1, cpu // searches), per_search

    def orthofinder_config(self):
        """OrthoFinder's config.json, the only place its diamond line lives.

        `-p 1` is written into the `search_cmd` template in that file.
        Nothing outside the process can move it: PATH was tried and
        OrthoFinder re-prepends its own two bin directories at startup, so
        a shim put in front of them lands behind them by the time the
        searches run. Measured, on the node, at position 9 of a PATH whose
        first entry was the environment's real diamond.

        There is no `--config` flag, so the install copy is the one that
        counts. Located from the `orthofinder` on PATH in its own env
        rather than by hard-coding a layout: it sits at
        `<env>/bin/src/orthofinder/run/config.json` in this install, and
        the globs below cover the other shapes a pip or conda install
        leaves behind.
        """
        exe = self.which_in_env("orp_orthofinder", "orthofinder")
        if exe is None:
            return None
        bindir = Path(exe).resolve().parent
        candidates = [bindir / "src" / "orthofinder" / "run" / "config.json"]
        candidates += sorted(bindir.glob("src/*/run/config.json"))
        candidates += sorted(bindir.parent.glob("lib/python*/site-packages/orthofinder/run/config.json"))
        for c in candidates:
            if c.is_file():
                return c
        return None

    def ensure_diamond_program(self, threads):
        """Add a `diamond_orp_<threads>` search program, and return its name.

        A copy of OrthoFinder's own `diamond` entry with `-p` set, added
        beside it rather than over it: the stock entry keeps working for
        anything else using this environment, and `-S diamond` still means
        exactly what it meant before.

        The name carries the thread count because this file is shared by
        every run on the cluster that uses this env. Two runs at different
        `--cpu` want different `-p`, and one entry per thread count lets
        them coexist instead of overwriting each other; a second run at the
        same thread count finds its entry already there and writes nothing.
        Nothing run-specific goes in -- no `--tmpdir`, no paths -- because a
        shared file must not carry one run's directories into another's.

        The write is temp-and-rename, then read back: two jobs adding
        different entries at the same moment is a lost update, and the
        read-back is what notices. One retry, because the loser of a race
        is not likely to lose twice.

        Returns None if the file cannot be read or written, which is a slow
        all-vs-all and not a wrong one -- run_orthofuser falls back to the
        stock program and says what that costs.
        """
        config = self.orthofinder_config()
        if config is None:
            return None
        name = f"diamond_orp_{threads}"
        for attempt in range(2):
            try:
                with open(config) as f:
                    data = json.load(f)
            except (OSError, ValueError) as e:
                print(f"    cannot read {config}: {e}")
                return None
            stock = data.get("diamond")
            if not isinstance(stock, dict) or "search_cmd" not in stock:
                print(f"    no usable 'diamond' entry in {config}")
                return None
            if name in data:
                return name
            entry = dict(stock)
            toks = entry["search_cmd"].split()
            if "-p" in toks:
                toks[toks.index("-p") + 1] = str(threads)
            else:
                toks += ["-p", str(threads)]
            entry["search_cmd"] = " ".join(toks)
            data[name] = entry
            backup = config.with_suffix(".json.orp-backup")
            first_backup = not backup.exists()
            try:
                if first_backup:
                    shutil.copy2(str(config), str(backup))
                tmp = config.with_suffix(f".json.orp-{os.getpid()}")
                with open(tmp, "w") as f:
                    json.dump(data, f, indent=4)
                    f.write("\n")
                os.replace(str(tmp), str(config))
            except OSError as e:
                print(f"    cannot write {config}: {e}")
                return None
            try:
                with open(config) as f:
                    if name in json.load(f):
                        print(f"    added search program '{name}' to {config}"
                              + (f" (original saved as {backup.name})" if first_backup else ""))
                        return name
            except (OSError, ValueError):
                pass
            print(f"    '{name}' did not survive the write -- another run writing "
                  f"the same file? retrying ({attempt + 1}/2)")
        return None

    def orthofinder_analysis_threads(self, cpu):
        """`-a`: OrthoFinder's workers for the algorithm phase after the searches.

        Upstream's default is "16 or t/8 (whichever lower)", deriving -a from
        -t. That is reasonable where -t means "cores you have" and wrong here,
        because this pipeline lowers -t to fit diamond's memory -- a
        constraint the algorithm phase does not share, since by the time it
        runs the searches have exited and their 500-odd GB with them. Deriving
        -a from the throttled -t is how a 40-core node ended up running that
        phase on one worker; it stalled and took the run with it.

        So: upstream's shape and ceiling, but computed from `cpu`, and capped
        at the number of species. That last cap is free and exact -- "Initial
        processing of each species" has one task per species, so a fifth
        worker on a four-assembly run has nothing to do.

        **There is deliberately no memory term.** Sizing this by memory needs
        a per-worker figure, and there is not one: upstream documents no RAM
        guidance for -a, and this phase has never been measured here. Its
        whole input is the Blast files -- 523 MB gzipped, ~4 GiB of text, on
        the run this was written for, against 143 GiB for a single search --
        so it is very unlikely to be what runs a node out of memory. That is
        an expectation, not a measurement, and it is the reason
        `--orthofinder-analysis` exists. Fit a memory term when there is a
        number to fit it to, the way ORTHOFINDER_GB_PER_SEARCH_GB2 was fitted
        and then validated; not before.

        Note the asymmetry with the searches, which is why this errs high
        where that errs low: too few searches costs wall time, too many loses
        the step to the OOM killer. Too few analysis workers is what trips
        OrthoFinder's 200s stall watchdog; too many, on present evidence,
        costs nothing.
        """
        if self.orthofinder_analysis:
            return max(1, self.orthofinder_analysis)
        species = max(1, len(self.search_fasta_paths()))
        return max(1, min(cpu // 8, ORTHOFINDER_MAX_ANALYSIS, species))

    def run_orthofuser(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        mem = self.mem if mem is None else mem
        # -t is the count of concurrent diamonds and so a memory knob before
        # it is a core count -- see ORTHOFINDER_GB_PER_SEARCH_GB2. -a, the
        # analysis threads, is RAM-hungry in its own right and OrthoFinder's
        # own default is t/8; it used to be handed the whole core count here,
        # which under -og buys nothing at all, since that run stops at
        # orthogroups and never reaches the MSA/tree work -a exists to
        # parallelise.
        searches, threads, per_search = self.orthofinder_search_plan(cpu, mem)
        analysis = self.orthofinder_analysis_threads(cpu)
        jobs = max(1, len(self.search_fasta_paths()) ** 2)
        why = ("--orthofinder-searches" if self.orthofinder_searches
               else f"{mem}G / {per_search}G per search")
        print(f"    all-vs-all: {jobs} searches, {searches} at a time ({why}), "
              f"{threads} thread(s) each; -a {analysis} for the algorithm phase")
        if not self.orthofinder_searches and per_search > mem:
            # mem // per_search is 0 here and max(1, ...) floors it to one
            # search, so the run proceeds having already computed that the
            # search does not fit. Memory and assembly size are both known
            # before anything starts; this is answerable now rather than as
            # an OOM four hours in. Not fatal: the estimate is a conservative
            # proxy (it predicted 159 GiB where 143.0 was measured), so a
            # budget just under it may still be survivable, and refusing
            # would turn a warning into a new way for the run not to start.
            biggest = max((q.stat().st_size for q in self.search_fasta_paths()
                           if q.is_file()), default=0)
            print(f"    *** one search is estimated at {per_search}G and the budget is "
                  f"{mem}G, so even a single search may not fit. The estimate is sized "
                  f"off the largest input ({biggest / 1e9:.1f} GB) and is deliberately "
                  "conservative, so this may still run -- but if diamond is OOM-killed "
                  "(returncode -9), the fix is more memory or smaller assemblies, not "
                  "--orthofinder-searches, which is already at its floor of 1. ***")
        program = self.ensure_diamond_program(threads)
        if program is None:
            program = "diamond"
            print(f"    *** could not set diamond's thread count: OrthoFinder's own `-p 1` "
                  f"stands, so this step gets {searches} of {cpu} cores. Correct, but slow. "
                  "Raise --orthofinder-searches to trade memory for cores. ***")
        # OrthoFinder reports its own fatal errors and then exits 0. A run
        # whose diamonds were OOM-killed leaves truncated Blast*.txt behind,
        # prints "ERROR: Blast1_1.txt is corrupted" and "ERROR: An error
        # occurred", and still returns success -- so conda_run is happy, the
        # sentinel gets written, and needs_run skips this step on every
        # later resume. makeorthout is then handed either nothing or a stale
        # Orthogroups.txt from an earlier attempt, and the run goes on to
        # build a final assembly off an orthogroup set that was never
        # computed. Take the sentinel from the artifact instead of from the
        # exit status: drop a marker first, and require orthogroups newer
        # than it, so a stale result from a previous attempt cannot pass.
        marker = self.orthofuse_dir / "orthofuser.attempt"
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.touch()
        self.conda_run(
            "orp_orthofinder", "orthofinder",
            "-d", "-I", "12", "-f", self.orthofuse_search,
            "-og", "-t", searches, "-a", analysis, "-S", program,
        )
        groups = self.newest_orthogroups_txt()
        if groups is None or groups.stat().st_mtime < marker.stat().st_mtime:
            sys.exit(self.no_orthogroups_message(cpu, mem))
        if groups.stat().st_size == 0:
            sys.exit(f"orthofinder produced an empty {groups}")
        # The checks above catch an attempt that produced no orthogroups at
        # all. They cannot see an attempt that produced them from a partial
        # all-vs-all, which is the more dangerous failure because everything
        # downstream accepts it -- see check_orthofinder_searches.
        workdir = self.orthofinder_working_dir(groups)
        if workdir is None:
            sys.exit(f"no WorkingDirectory above {groups} -- cannot check orthofinder's all-vs-all")
        self.check_orthofinder_searches(workdir)
        (self.orthofuse_dir / "orthofuser.done").touch()

    def merge(self):
        out = self.orthofuse_dir / "merged.fasta"
        with open(out, "wb") as outf:
            for p in self.short_fasta_paths():
                with open(p, "rb") as inf:
                    shutil.copyfileobj(inf, outf)

    def newest_orthogroups_txt(self):
        """The most recently written Orthogroups.txt, or None.

        Newest rather than rglob's first: OrthoFinder never reuses a results
        directory, it makes a new Results_<Mon><Day> (then _1, _2, ...) per
        invocation, so a working directory that has seen a failed attempt
        holds several. rglob's order is the filesystem's, which on a resume
        after a failure is as likely to hand back the attempt that died as
        the one that succeeded -- and makeorthout would pick contigs from it
        without complaint.
        """
        matches = list(self.orthofuse_search.rglob("Orthogroups.txt"))
        if not matches:
            return None
        return max(matches, key=lambda p: p.stat().st_mtime)

    def find_orthogroups_txt(self):
        match = self.newest_orthogroups_txt()
        if match is None:
            sys.exit("Orthogroups.txt not found under orthofuse working directory")
        return match

    # -- all-vs-all validation ---------------------------------------------

    @staticmethod
    def _count_fasta_records(path):
        """Deflines in a FASTA, counted without decoding the sequence."""
        n = 0
        tail = b"\n"
        with open(path, "rb") as fh:
            while True:
                chunk = fh.read(8 << 20)
                if not chunk:
                    return n
                n += (tail[-1:] + chunk).count(b"\n>")
                tail = chunk

    @staticmethod
    def _count_lines(path, stop_at):
        """Lines in a plain or gzipped file, giving up once stop_at is reached.

        The healthy self-comparisons in a four-assembly merge run to tens of
        megabytes compressed, and nothing here needs their exact size -- only
        whether they cleared the floor. Short-circuiting keeps this check a
        few seconds on a directory that took ten hours to produce.
        """
        opener = gzip.open if path.suffix == ".gz" else open
        n = 0
        with opener(path, "rb") as fh:
            while n < stop_at:
                chunk = fh.read(8 << 20)
                if not chunk:
                    break
                n += chunk.count(b"\n")
        return n

    def orthofinder_working_dir(self, groups):
        """The WorkingDirectory holding the searches behind `groups`.

        Walked up from Orthogroups.txt rather than globbed for the newest
        Results_<Mon><Day>: a working directory that has seen a failed
        attempt holds several of those trees, and the searches that need
        checking are the ones belonging to *this* Orthogroups.txt.
        """
        for parent in groups.parents:
            candidate = parent / "WorkingDirectory"
            if candidate.is_dir():
                return candidate
            if parent == self.orthofuse_search:
                break
        return None

    def newest_working_dir(self):
        """The newest Results_*/WorkingDirectory, with or without orthogroups.

        check_orthofinder_searches reaches its working directory by walking
        up from an Orthogroups.txt. That is the right anchor when there is
        one and no help at all when there is not -- which is exactly the
        case where the searches most need looking at, because an attempt
        that produced no orthogroups may still have produced sixteen
        perfectly good Blast files that cost four and a half hours.
        """
        matches = [d for d in self.orthofuse_search.rglob("WorkingDirectory") if d.is_dir()]
        if not matches:
            return None
        return max(matches, key=lambda d: d.stat().st_mtime)

    def audit_orthofinder_searches(self, workdir):
        """(species, failures) for the all-vs-all in `workdir`, judging nothing.

        Split out of check_orthofinder_searches so the same audit can be run
        where a failure must not be fatal: when there is no Orthogroups.txt
        at all, the question "are the searches good?" decides whether the
        right advice is to delete this directory or to guard it with your
        life. Returns (None, []) when there are no Species*.fa to judge.
        """
        species = {}
        for fa in workdir.glob("Species*.fa"):
            m = re.fullmatch(r"Species(\d+)\.fa", fa.name)
            if m:
                species[int(m.group(1))] = fa
        if not species:
            return None, []
        counts = {i: self._count_fasta_records(fa) for i, fa in species.items()}
        failures = []
        for i in sorted(species):
            for j in sorted(species):
                stem = workdir / f"Blast{i}_{j}.txt"
                blast = next((p for p in (stem.with_suffix(".txt.gz"), stem) if p.is_file()), None)
                if blast is None:
                    failures.append(f"  Blast{i}_{j}: no output file")
                    continue
                floor = int(counts[i] * BLAST_SELF_HIT_FLOOR) if i == j else 1
                got = self._count_lines(blast, floor)
                if got < floor:
                    why = (f"{got} hits, expected at least {floor} "
                           f"({BLAST_SELF_HIT_FLOOR:.0%} of {counts[i]} sequences in {species[i].name})"
                           if i == j else "no hits at all")
                    failures.append(f"  Blast{i}_{j} ({blast.stat().st_size} bytes): {why}")
        return species, failures

    def no_orthogroups_message(self, cpu, mem):
        """What to say when an attempt produced no Orthogroups.txt.

        This used to say one thing: diamond was OOM-killed, lower the
        concurrency, delete the Results_* directory and start again. That is
        right when the searches died and catastrophic when they did not --
        the run this was rewritten for completed all sixteen searches in 4h36m
        and then failed in OrthoFinder's own algorithm phase ("Initial
        processing of each species", stalled at 3/4). Deleting on that advice
        throws away four and a half hours of good alignments to fix something
        that was never wrong.

        So audit the searches first and let them decide which advice this is.
        """
        planned = self.orthofinder_search_plan(cpu, mem)[0]
        head = ("orthofinder exited 0 but produced no Orthogroups.txt for this "
                "attempt -- read its ERROR lines above.\n")
        workdir = self.newest_working_dir()
        if workdir is None:
            return head + "\nNo WorkingDirectory to inspect, so the searches cannot be judged."
        species, failures = self.audit_orthofinder_searches(workdir)
        if species is None:
            return head + f"\nNo Species*.fa in {workdir}, so the searches cannot be judged."
        if failures:
            return (
                head
                + f"\n{len(failures)} of {len(species) ** 2} searches in {workdir} are bad:\n"
                + "\n".join(failures)
                + "\n\nThat is the usual cause: diamond OOM-killed mid-search (returncode -9) "
                  "leaves truncated Blast*.txt that OrthoFinder reports as corrupted. Run again "
                  f"with a lower --orthofinder-searches (this attempt used {planned}), and delete\n"
                  f"  {workdir.parent}\n"
                  "before resuming -- OrthoFinder will not recompute a search it can see an "
                  "output file for, so a truncated one left in place is a truncated one reused."
            )
        return (
            head
            + f"\nAll {len(species) ** 2} searches in\n  {workdir}\nare complete and non-empty, "
              "so the all-vs-all is not what failed and diamond is not what to fix.\n\n"
              "*** Do NOT delete that directory. *** It holds the finished alignments, which are "
              "the expensive part of this step; OrthoFinder will reuse them rather than recompute "
              "them. Lowering --orthofinder-searches would cost hours and change nothing.\n\n"
              "The failure is in OrthoFinder's own algorithm phase, after the searches. Look for "
              "its 'Initial processing of each species' or 'Stalled for' lines above, and check "
              f"{workdir.parent}/Log.txt."
        )

    def check_orthofinder_searches(self, workdir):
        """Fail on an all-vs-all whose searches died without saying so.

        OrthoFinder runs n_assemblies^2 diamonds and reports a child that
        failed by printing `ERROR: external program returned code` and then
        exiting 0 anyway. The run that prompted this lost six of sixteen
        searches -- two killed before they wrote a byte, two that wrote a
        zero-hit file, two that stopped a fraction of the way in -- and
        every one of the sixteen was still a valid gzip stream, so nothing
        downstream had any reason to object. `printf '' | gzip -c` is 20
        bytes and passes `gzip -t`: an empty result is indistinguishable
        from an honest one by any check that does not look inside.

        OrthoFinder then built orthogroups from what survived and wrote a
        perfectly well-formed Orthogroups.txt. Both SPAdes assemblies had
        lost their self-comparison and their comparison with each other, so
        they were present in the clustering only through hits found by the
        two assemblies that still worked -- a merge silently missing half
        its inputs, which makeorthout would have picked contigs from
        without complaint.

        Two rules, both on the artifacts rather than on an exit status:

        **Every search produced at least one hit.** No pair of assemblies
        from the same library has nothing whatsoever in common, so an empty
        Blast{i}_{j} is a dead search, whatever returncode it reported.

        **Each self-comparison found most of its own sequences.** Blast{i}_i
        aligns Species{i} against itself, so every sequence hits itself and
        the file should carry at least one line per record in the input.
        The floor is well under 1.0 because that is not quite guaranteed:
        diamond masks low-complexity before seeding, and a sequence masked
        end to end cannot seed, so it reports no self-hit. SPAdes N-gaps
        are exactly that case -- N is asparagine to `diamond blastp`, and
        ~13% of the contigs in each SPAdes assembly here carry a run of
        them. A floor of 0.5 sits far below that legitimate shortfall and
        far above a failure: the worst surviving self-comparison in the run
        this was written for held about 0.2% of its input.
        """
        species, failures = self.audit_orthofinder_searches(workdir)
        if species is None:
            sys.exit(f"orthofinder left no Species*.fa in {workdir} -- cannot check its all-vs-all")
        if failures:
            names = "\n".join(f"  {i}: {fa.name}" for i, fa in sorted(species.items()))
            sys.exit(
                "orthofinder exited 0 but its all-vs-all is incomplete -- "
                f"{len(failures)} of {len(species) ** 2} searches failed:\n"
                + "\n".join(failures)
                + "\n\nspecies:\n" + names
                + "\n\nOrthogroups.txt built on this is missing real relationships, "
                "so the merge would silently drop whatever those searches would have "
                "found. Read the diamond errors in this run's log -- OrthoFinder does "
                "not capture its children's stderr, so their own message is there and "
                "not in OrthoFinder's output. Fix what they report (memory "
                f"pressure from running too many of the {len(species) ** 2} diamonds "
                "at once is the usual cause, and --orthofinder-searches is the knob "
                f"for it), then delete {workdir.parent} before resuming."
            )

    @staticmethod
    def clear_transrate_outdir(outdir, assembly):
        """Clear pytransrate's -o of what a retry must not reuse, and only that.

        A retry has to start from a directory pytransrate can work in:
        it will not overwrite an existing assemblies.csv, and it reuses the
        BAM and the quant.sf it finds in -o on their existence alone, so a
        step killed part-way leaves half-written copies of both behind and
        every later attempt picks the same ones back up.

        The question to ask of each thing in -o is therefore not whether it
        is there but whether whatever wrote it finished -- because in this
        directory the expensive artifacts and the half-written ones are the
        same files. rmtree'ing the lot answers that question by throwing
        away the run, which is the wrong answer twice over: three attempts
        at a failing step meant three identical index builds, three
        identical mappings and three identical waits to reach the same
        failure. What survives does so on evidence that it is whole, and
        everything below is that evidence.

        **The snap index.** Building it is the longest single piece of work
        in the run -- the better part of an hour on a multi-million-contig
        merge, and two builds rather than one whenever the -locationSize
        sweep steps up. pytransrate keys its own reuse on the GenomeIndex
        marker snap writes when a build completes, and a build that died
        half way leaves its directory behind without one, so that is the
        marker checked here too: trusting a partial index yields a corrupt
        one.

        **A complete BAM.** Mapping is the other multi-hour step, and on a
        large merged assembly the BAM runs to hundreds of gigabytes. A BAM
        ending in the empty BGZF block of SAM spec 4.1 was closed by a
        writer that got to the end; one that does not was still being
        written when snap was OOM-killed, hit its wall clock, or died on
        amplab/snap#171. pytransrate 2.2.0 makes that same check before
        reusing one and remaps over a BAM that fails it, but the decision
        to keep the file is this function's, and keeping a truncated one
        would be keeping hundreds of gigabytes nothing can use. It is not
        kept as evidence either -- logs/snap.log is the evidence. A BAM
        that is kept keeps its .align.done marker beside it, which is what
        pytransrate reads after 2.2.0 to decide the same question: drop the
        marker and it would move a perfectly good BAM aside and map again.
        The <index>.index.lock file is kept for the same reason -- it sits
        beside the index rather than inside it precisely so an rmtree of a
        partial index cannot pull it out from under its own holder, and it
        belongs to the index either way.

        **The read count** that goes with it, `*-read_count.txt`. It is
        keyed on the read filenames and depends only on the reads, so it
        cannot go stale while those names hold. Keeping it matters more
        than its size suggests: it is what pytransrate reads when it reuses
        a BAM, and without it that path falls back to counting lines in the
        fastq itself.

        **salmon/**, but only beside a BAM that was kept, and only when
        quant.sf is whole. quant.sf is a completed quantification *of that
        BAM*, so keeping it when the BAM it was computed from has gone
        would score the assembly off numbers belonging to a file that no
        longer exists; and see quant_sf_is_complete for why existence is
        not enough on its own -- no pytransrate to date checks it.

        **logs/**, which holds snap.log, the file pytransrate points at
        when snap dies without explaining itself, so deleting it is
        deleting the evidence the retry exists to gather. pytransrate
        rewrites it per attempt, so what survives the last retry is the
        last attempt's output, which is the one worth reading.

        Everything else goes: assemblies.csv, contigs.csv and the score
        optimisation csv are the outputs being recomputed, and anything a
        future pytransrate leaves behind that this does not recognise is
        cleared rather than assumed safe.
        """
        outdir = Path(outdir)
        assembly = Path(assembly)
        if not outdir.is_dir():
            return

        keep = {outdir / "logs"}
        keep.update(outdir.glob("*.index.lock"))
        kept_bam = False
        for bam in outdir.glob("*.bam"):
            if bam_is_complete(bam):
                kept_bam = True
                keep.add(bam)
                keep.add(bam.with_name(bam.name + ".align.done"))

        if kept_bam:
            keep.update(outdir.glob("*read_count.txt"))
            quant_sf = outdir / "salmon" / "quant.sf"
            # Ordered so the assembly is only walked when there is a
            # quant.sf whose fate depends on the answer.
            if (quant_sf.is_file() and assembly.is_file()
                    and quant_sf_is_complete(quant_sf, count_sequences(assembly))):
                keep.add(outdir / "salmon")

        survived = []
        for path in outdir.iterdir():
            if path in keep:
                survived.append(path)
                continue
            if path.is_dir():
                if (path / "GenomeIndex").is_file():
                    survived.append(path)
                    continue
                shutil.rmtree(path, ignore_errors=True)
            else:
                # Not unlink(missing_ok=True): that keyword is 3.8+, and
                # oyster.py is launched by whatever system python3 the
                # cluster has -- 3.6.8 on ours. A TypeError raised here
                # fires only on the retry path, i.e. only once a step has
                # already failed, so it converts a retryable failure into a
                # crash whose traceback hides the failure that caused it.
                try:
                    path.unlink()
                except OSError:
                    pass

        # Said out loud because the decision is worth hours either way, and
        # because the log is the only place anyone can check it after the
        # fact: a retry that silently remapped and a retry that reused a
        # good BAM look identical from outside until the wall time comes in.
        print("[transrate] {}: kept from the last attempt: {}".format(
            outdir.name,
            ", ".join(p.name for p in sorted(survived)) or "nothing",
        ))

    def pytransrate_memory_args(self, mem):
        """``--max-memory`` for a pytransrate call, or nothing.

        pytransrate sizes the read-metrics step's shared accumulators by the
        assembly and multiplies them by ``--threads``: on a 5.8 Gbp merge
        that is 23 GB per worker, so ``-t 40`` asks for 928 GB. It caps the
        workers at what fits, but only against a budget it can find, and on
        an unconstrained login or interactive node the only figure available
        is the whole machine's free memory. A run given ``--mem 670`` on a
        1.5 TB node therefore passed its own cap and was OOM-killed after
        nine hours of mapping and quantifying -- twice, because the retry had
        no more reason to cap than the first attempt did.

        So --mem is forwarded. It is the number the user chose and the number
        every step above already splits between concurrent jobs; leaving it
        at this one boundary meant the most memory-hungry step in the run was
        the only one that never heard it. Requires pytransrate >= 2.2.0,
        which is what PYTRANSRATE_MIN_VERSION enforces at preflight.
        """
        if mem is None:
            return []
        if any(arg.split("=")[0] in PYTRANSRATE_MEMORY_FLAGS
               for arg in self.pytransrate_args):
            return []
        return ["--max-memory", f"{mem}G"]

    def orthotransrate(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        outdir = self.orthofuse_dir / "merged"
        merged = self.orthofuse_dir / "merged.fasta"
        # needs_run() re-runs this step whenever the corrected reads are
        # newer than merged/assemblies.csv -- not only when it is absent --
        # so a resumed run would abort on that csv unless it is cleared
        # first. retry_cleanup repeats the clear before each retry, because
        # the one below happens once, outside run()'s retry loop. See
        # clear_transrate_outdir for what survives it and why.
        self.clear_transrate_outdir(outdir, merged)
        self.conda_run(
            "orp", "pytransrate",
            "-o", outdir, "-t", cpu, "-a", merged,
            "--left", self.cor1(), "--right", self.cor2(),
            *self.pytransrate_memory_args(mem),
            *self.pytransrate_args,
            retry_cleanup=partial(self.clear_transrate_outdir, outdir, merged),
        )

    def makeorthout(self):
        """Pick the best-scoring contig per orthogroup.

        The picker used to consume a directory of one <i>.groups file per
        orthogroup, written by a makelist/makegroups pair here and unlinked
        again on the way out -- of order 1e5 small files created, globbed
        back in and deleted, purely to hand data between two Python
        processes. It reads Orthogroups.txt directly now; the group ordering
        that used to come out of sorting those filenames is reproduced
        inside the script, deliberately, because it reaches cd-hit-est and
        so the final assembly (see the note at the top of
        scripts/pick_best_contigs.py).
        """
        print("Picking the best contig per orthogroup")
        contigs_csv = next(self.orthofuse_dir.rglob("contigs.csv"), None)
        if contigs_csv is None:
            sys.exit("contigs.csv not found under orthofuse directory")
        good_list = self.orthofuse_dir / f"good.{self.runout}.list"
        self.conda_run(
            "orp", "python", self.makedir / "scripts" / "pick_best_contigs.py",
            contigs_csv, self.find_orthogroups_txt(), good_list,
        )

    def orthofusing(self):
        good_list = self.orthofuse_dir / f"good.{self.runout}.list"
        out = self.assemblies_dir / f"{self.runout}.orthomerged.fasta"
        with open(out, "w") as outf:
            subprocess.run(
                ["conda", "run", "--no-capture-output", "-n", "orp", "python",
                 str(self.makedir / "scripts" / "filter.py"), str(self.orthofuse_dir / "merged.fasta"), str(good_list)],
                check=True, stdout=outf, cwd=self.dir,
            )

    # -- diamond passes ----------------------------------------------------

    def diamond_jobs(self):
        return [
            (self.assemblies_dir / f"{self.runout}.orthomerged.fasta",
             self.diamond_dir / f"{self.runout}.orthomerged.diamond.txt"),
        ] + [(self.assembly_fasta(a), self.diamond_txt(a)) for a in self.diamond_priority]

    def run_diamond_one(self, query, out, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        self.conda_run(
            "orp", "diamond", "blastx", "--quiet", "-p", cpu,
            "-e", "1e-8", "--top", "0.1", "-q", query, "-d", self.diamond_db, "-o", out,
        )

    def diamond_uniq(self):
        for a in self.report_order:
            count = parse_unique_count(self.diamond_txt(a))
            self.unique_txt(a).write_text(f"{count}\n")

    def make_list1(self):
        ids = extract_gene_ids(self.diamond_dir / f"{self.runout}.orthomerged.diamond.txt")
        write_sorted(self.diamond_dir / f"{self.runout}.list1", ids)

    def make_list2(self):
        ids = set()
        for a in self.diamond_priority:
            ids |= extract_gene_ids(self.diamond_txt(a))
        write_sorted(self.diamond_dir / f"{self.runout}.list2", ids)

    def make_list3(self):
        list1 = set(line.rstrip("\n") for line in open(self.diamond_dir / f"{self.runout}.list1"))
        list2_path = self.diamond_dir / f"{self.runout}.list2"
        out = self.diamond_dir / f"{self.runout}.list3"
        with open(list2_path) as f, open(out, "w") as o:
            for line in f:
                if line.rstrip("\n") not in list1:
                    o.write(line)

    def make_list5(self):
        # build_list5.py keeps the first hit per gene in the order it is
        # given the files, so diamond_priority is a real preference ranking
        # here, not just an iteration order.
        self.conda_run(
            "orp", "python", self.makedir / "scripts" / "build_list5.py",
            self.diamond_dir / f"{self.runout}.list3", self.diamond_dir / f"{self.runout}.list5",
            *[self.diamond_txt(a) for a in self.diamond_priority],
        )

    def make_list6(self):
        fasta = self.assemblies_dir / f"{self.runout}.orthomerged.fasta"
        out = self.diamond_dir / f"{self.runout}.list6"
        with open(fasta) as f, open(out, "w") as o:
            for line in f:
                if line.startswith(">"):
                    o.write(line[1:])

    def make_list7(self):
        list6 = set(line.rstrip("\n") for line in open(self.diamond_dir / f"{self.runout}.list6"))
        list5_path = self.diamond_dir / f"{self.runout}.list5"
        out = self.diamond_dir / f"{self.runout}.list7"
        with open(list5_path) as f, open(out, "w") as o:
            for line in f:
                if line.rstrip("\n") not in list6:
                    o.write(line)

    def posthack(self):
        # Concatenation order, so ASSEMBLY_ORDER and not diamond_priority:
        # this reaches cd-hit-est, where input order breaks length ties.
        fastas = " ".join(str(p) for p in self.assembly_fasta_paths())
        list7 = self.diamond_dir / f"{self.runout}.list7"
        newbies = self.diamond_dir / f"{self.runout}.newbies.fasta"
        orthomerged = self.assemblies_dir / f"{self.runout}.orthomerged.fasta"
        working_out = self.assemblies_working / f"{self.runout}.orthomerged.fasta"
        filter_py = self.makedir / "scripts" / "filter.py"
        script = f"python {filter_py} <(cat {fastas}) {list7} >> {newbies}"
        self.run(["conda", "run", "--no-capture-output", "-n", "orp", "bash", "-c", script])
        with open(working_out, "wb") as outf:
            for p in (newbies, orthomerged):
                with open(p, "rb") as inf:
                    shutil.copyfileobj(inf, outf)

    # -- dedup / quantify -----------------------------------------------------

    def cdhit(self):
        src = self.assemblies_working / f"{self.runout}.orthomerged.fasta"
        out = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        self.conda_run(
            "orp", "cd-hit-est", "-M", self.mem * 1000, "-T", self.cpu,
            "-c", ".98", "-i", src, "-o", out,
        )

    def orp_diamond(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        src = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        out = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        self.conda_run(
            "orp", "diamond", "blastx", "--quiet", "-p", cpu,
            "-e", "1e-8", "--top", "0.1", "-q", src, "-d", self.diamond_db, "-o", out,
        )
        out.touch()

    def orp_uniq(self):
        diamond_txt = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        count = parse_unique_count(diamond_txt)
        (self.assemblies_working / f"{self.runout}.unique.ORP.txt").write_text(f"{count}\n")
        (self.assemblies_working / f"{self.runout}.unique.ORP.done").touch()

    def salmon_index(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        src = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        idx = self.quants_dir / f"{self.runout}.ortho.idx"
        # A rebuild here is usually a rebuild *over* an index salmon has
        # already refused to load, so clear it rather than writing into the
        # old directory alongside whatever format it was in.
        shutil.rmtree(idx, ignore_errors=True)
        self.conda_run(
            "orp", "salmon", "index", "--no-version-check", "-t", src,
            "-i", idx, "-k", "31", "--threads", cpu,
        )

    def salmon(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        idx = self.quants_dir / f"{self.runout}.ortho.idx"
        outdir = self.quants_dir / f"salmon_orthomerged_{self.runout}"
        self.conda_run(
            "orp", "salmon", "quant", "--no-version-check",
            "-p", cpu, "-i", idx, "--seqBias", "--gcBias", "--libType", "A",
            "-1", self.cor1(), "-2", self.cor2(), "-o", outdir,
        )

    def filter_tpm(self):
        quant = self.quants_dir / f"salmon_orthomerged_{self.runout}" / "quant.sf"
        high = self.assemblies_working / f"{self.runout}.HIGHEXP.txt"
        low = self.assemblies_working / f"{self.runout}.LOWEXP.txt"
        with open(quant) as f, open(high, "w") as hf, open(low, "w") as lf:
            next(f, None)
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 4:
                    continue
                try:
                    tpm = float(cols[3])
                except ValueError:
                    continue
                if tpm > self.tpm_filt:
                    hf.write(cols[0] + "\n")
                elif tpm < self.tpm_filt:
                    lf.write(cols[0] + "\n")
        (self.assemblies_dir / f"{self.runout}.filter.done").touch()
        print("\n\n\n\n PART: TPM_FILT MAKE LOW AND HIGH\n\n\n\n")

    def secondfilter(self):
        intermediate = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        low = self.assemblies_working / f"{self.runout}.LOWEXP.txt"
        high = self.assemblies_working / f"{self.runout}.HIGHEXP.txt"
        diamond_txt = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        filter_py = self.makedir / "scripts" / "filter.py"

        if low.exists() and low.stat().st_size > 0:
            before = self.assemblies_working / f"{self.runout}.ORP_BEFORE_TPM_FILT.fasta"
            shutil.copy(intermediate, before)

            highexp_fasta = self.assemblies_working / f"{self.runout}.ORP.HIGHEXP.fasta"
            with open(highexp_fasta, "w") as outf:
                subprocess.run(
                    ["conda", "run", "--no-capture-output", "-n", "orp", "python",
                     str(filter_py), str(intermediate), str(high)],
                    check=True, stdout=outf, cwd=self.dir,
                )

            low_ids = set(line.strip() for line in open(low) if line.strip())
            blasted = self.assemblies_working / f"{self.runout}.blasted"
            donotremove = self.assemblies_working / f"{self.runout}.donotremove.list"
            do_not_remove_ids = set()
            with open(diamond_txt) as f, open(blasted, "a") as bf:
                for line in f:
                    cols = line.rstrip("\n").split("\t")
                    if cols and cols[0] in low_ids:
                        bf.write(line)
                        do_not_remove_ids.add(cols[0])
            with open(donotremove, "a") as df:
                for i in sorted(do_not_remove_ids):
                    print(i)
                    df.write(i + "\n")

            saveme = self.assemblies_working / f"{self.runout}.saveme.fasta"
            with open(saveme, "w") as outf:
                subprocess.run(
                    ["conda", "run", "--no-capture-output", "-n", "orp", "python",
                     str(filter_py), str(intermediate), str(donotremove)],
                    check=True, stdout=outf, cwd=self.dir,
                )

            with open(orp_fasta, "wb") as outf:
                for p in (saveme, highexp_fasta):
                    with open(p, "rb") as inf:
                        shutil.copyfileobj(inf, outf)
            print("\n\n\n\n PART: FILTER LOW STUFF \n\n\n\n")
        else:
            shutil.copy(intermediate, orp_fasta)
            print("\n\n\n\n PART: THERE IS NO LOW STUFF \n\n\n\n")

    # -- QC / reporting ----------------------------------------------------

    def busco(self, cpu=None, mem=None):
        cpu = self.busco_threads if cpu is None else cpu
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        self.conda_run(
            "orp_busco", "busco", "--offline", "--lineage", self.lineage,
            "--download_path", self.makedir / "busco_dbs",
            "-i", orp_fasta, "-m", "transcriptome", "--cpu", cpu,
            "-o", f"run_{self.runout}.ORP", "--config", self.busco_config,
        )
        for p in self.dir.glob(f"run_{self.runout}*"):
            shutil.move(str(p), str(self.reports_dir / p.name))
        (self.reports_dir / f"{self.runout}.busco.done").touch()

    def transrate(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        outdir = self.reports_dir / f"transrate_{self.runout}"
        # See orthotransrate() and clear_transrate_outdir.
        self.clear_transrate_outdir(outdir, orp_fasta)
        self.conda_run(
            "orp", "pytransrate",
            "-o", outdir, "-a", orp_fasta,
            "--left", self.cor1(), "--right", self.cor2(), "-t", cpu,
            *self.pytransrate_memory_args(mem),
            *self.pytransrate_args,
            retry_cleanup=partial(self.clear_transrate_outdir, outdir, orp_fasta),
        )

    def trinity_perllib_dir(self):
        result = subprocess.run(
            ["conda", "run", "-n", "orp_trinity", "which", "Trinity"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True,
        )
        trinity_path = Path(result.stdout.strip()).resolve()
        return trinity_path.parent / "PerlLib"

    def strandeval(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        r1, r2 = self.cor1(), self.cor2()
        self.conda_run("orp_trinity", "bwa", "index", "-p", self.runout, orp_fasta)

        pipeline_script = (
            f'conda run --no-capture-output -n orp_trinity bash -c '
            f'"bwa mem -t {cpu} {self.runout} '
            f'<(seqtk sample -s 23894 {r1} 400000) <(seqtk sample -s 23894 {r2} 400000)" '
            f"| conda run --no-capture-output -n orp samtools view -@10 -Sb - "
            f"| conda run --no-capture-output -n orp samtools sort -T {self.runout} -O bam -@10 "
            f"-o {self.runout}.sorted.bam -"
        )
        self.run(["bash", "-o", "pipefail", "-c", pipeline_script])

        flagstat_out = self.assemblies_dir / f"{self.runout}.flagstat"
        with open(flagstat_out, "w") as outf:
            subprocess.run(
                ["conda", "run", "--no-capture-output", "-n", "orp", "samtools", "flagstat",
                 f"{self.runout}.sorted.bam"],
                check=True, stdout=outf, cwd=self.dir,
            )

        perllib = self.trinity_perllib_dir()
        self.conda_run(
            "orp_trinity", "perl", "-I", str(perllib),
            str(self.makedir / "scripts" / "examine_strand.pl"),
            f"{self.runout}.sorted.bam", self.runout,
        )

        dat_file = self.dir / f"{self.runout}.dat"
        hist_input = self.dir / f"{self.runout}.hist_input.txt"
        with open(dat_file) as f, open(hist_input, "w") as out:
            next(f, None)
            for line in f:
                cols = line.rstrip("\n").split()
                if len(cols) >= 5:
                    out.write(cols[4] + "\n")

        hist_result = subprocess.run(
            ["conda", "run", "--no-capture-output", "-n", "orp_trinity", "bash", "-c",
             f"hist -p '#' -c red {hist_input}"],
            check=True, cwd=self.dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
        )

        if hist_input.exists():
            hist_input.unlink()
        sorted_bam = self.dir / f"{self.runout}.sorted.bam"
        if sorted_bam.exists():
            sorted_bam.unlink()
        for ext in ("bwt", "pac", "ann", "amb", "sa", "dat"):
            p = self.dir / f"{self.runout}.{ext}"
            if p.exists():
                p.unlink()

        (self.reports_dir / f"{self.runout}.strandeval.done").touch()
        histogram_text = hist_result.stdout.rstrip("\n")
        summary = (
            "\n*****  STRAND EXAMINATION HISTOGRAM ***** \n"
            f"{histogram_text}\n"
            "\n*****  See the following link for interpretation ***** \n"
            "*****  https://oyster-river-protocol.readthedocs.io/en/latest/strandexamine.html ***** \n"
        )
        print(summary)
        (self.reports_dir / f"{self.runout}.strandeval_summary.txt").write_text(summary)

    def reportgen(self):
        runout = self.runout
        lines = []

        def emit(label, value):
            text = f"{label}      {value}"
            print(text)
            lines.append(text)

        header = f"*****  QUALITY REPORT FOR: {runout} using {self.RUN_DESCRIPTION} version {self.version} ****"
        print(f"\n\n{header}")
        lines.append(header)
        orp_fasta = self.assemblies_dir / f"{runout}.ORP.fasta"
        print(f"\n*****  THE ASSEMBLY CAN BE FOUND HERE: {orp_fasta} **** \n")

        busco_line = ""
        busco_dir = self.reports_dir / f"run_{runout}.ORP"
        for p in busco_dir.rglob("short*.txt"):
            for line in open(p):
                if re.match(r"^\s*C:[0-9]", line):
                    busco_line = line.strip()
        emit("*****  BUSCO SCORE ~~~~~~~~~~~~~~~~~~~~~~>", busco_line)

        csv_path = next((self.reports_dir / f"transrate_{runout}").rglob("assemblies.csv"), None)
        rows = list(csv.reader(open(csv_path))) if csv_path else []
        transrate_score = rows[1][36] if len(rows) > 1 and len(rows[1]) > 36 else ""
        transrate_optimal = rows[1][37] if len(rows) > 1 and len(rows[1]) > 37 else ""
        emit("*****  TRANSRATE SCORE ~~~~~~~~~~~~~~~~~~>     ", transrate_score)
        emit("*****  TRANSRATE OPTIMAL SCORE ~~~~~~~~~~>     ", transrate_optimal)

        def read_count(path):
            return path.read_text().strip() if path.exists() else ""

        def unique_genes(label):
            # The tildes pad each label out to the same column so the counts
            # line up under one another; with the four built-in assemblers
            # this reproduces the hand-written arrows exactly.
            return f"*****  UNIQUE GENES {label} " + "~" * max(1, 20 - len(label)) + ">     "

        emit(unique_genes("ORP"), read_count(self.assemblies_working / f"{runout}.unique.ORP.txt"))
        for a in self.report_order:
            emit(unique_genes(a.report_label), read_count(self.unique_txt(a)))

        proper_pairs = ""
        flagstat = self.assemblies_dir / f"{runout}.flagstat"
        if flagstat.exists():
            for line in open(flagstat):
                if "properly paired" in line:
                    fields = line.split()
                    if len(fields) > 5:
                        proper_pairs = fields[5].lstrip("(")
        emit("*****  READS MAPPED AS PROPER PAIRS ~~~~~>     ", proper_pairs)
        print(" \n")

        strandeval_summary_path = self.reports_dir / f"{runout}.strandeval_summary.txt"
        if strandeval_summary_path.exists():
            strandeval_summary = strandeval_summary_path.read_text()
            print(strandeval_summary)
            lines.append(strandeval_summary.rstrip("\n"))

        qualreport_path = self.reports_dir / f"qualreport.{runout}"
        qualreport_path.write_text("\n".join(lines) + "\n")
        (self.reports_dir / f"qualreport.{runout}.done").touch()

    def timing_report(self, wall_clock_seconds):
        with open(self.timing_log) as f:
            all_lines = f.readlines()
        header = all_lines[0].rstrip("\n") if all_lines else ""
        body = [l for l in all_lines if "\t" in l]

        out_lines = [header, "", f"*****  STEP TIMING for {self.runout} ***** ", ""]
        for line in body:
            parts = line.rstrip("\n").split("\t")
            name, secs = parts[0], int(parts[1])
            started = parts[2] if len(parts) > 2 else ""
            h, m, s = secs // 3600, (secs % 3600) // 60, secs % 60
            suffix = f"  (started {started})" if started else ""
            out_lines.append(f"{name:<16} {h:02d}:{m:02d}:{s:02d}{suffix}")
        # TOTAL is real wall-clock elapsed time, not a sum of the lines above:
        # concurrent jobs overlap (their durations shouldn't add), and a
        # branch's own entry already includes the sub-steps it calls, so
        # summing every line double-counts and overstates the true runtime.
        h, m, s = wall_clock_seconds // 3600, (wall_clock_seconds % 3600) // 60, wall_clock_seconds % 60
        out_lines.append(f"{'TOTAL':<16} {h:02d}:{m:02d}:{s:02d}")

        text = "\n".join(out_lines) + "\n"
        print(text)
        self.timing_log.write_text(text)
        print(f"\nFull timing log saved to: {self.timing_log}\n")

    # -- orchestration -------------------------------------------------------

    def main(self):
        pipeline_start = time.time()
        self.setup()
        if self.already_complete():
            return
        self.timing_init()
        self.check()
        self.welcome()
        self.readcheck()
        self.prepare_reads()
        self.run_assemblers()
        self.merge_and_report(pipeline_start)

    def prepare_reads(self):
        """Trim and error-correct the raw pair, and start compressing it.

        Split out of main() because every entry point needs it: the merge
        half scores, quantifies and strand-checks against the corrected
        pair, so a run that brings its own assemblies still comes through
        here.
        """
        t1, t2 = self.trim1(), self.trim2()
        trim_done = self.rcorr_dir / f"{self.runout}.trim.done"
        c1, c2 = self.cor1(), self.cor2()

        # Ask for the trimmed pair back as trimmomatic's outputs only while
        # the corrected pair it feeds is missing or stale: reclaim_trimmed_
        # reads() deletes those files as soon as rcorrector is done with them,
        # so on a resumed run they're legitimately gone and their sentinel,
        # not the files, is what records that trimming happened. A working
        # directory from before the sentinel existed has no sentinel and its
        # TRIM files still present, so it takes the first branch and skips
        # the same way it always did.
        trim_outputs = [t1, t2]
        if trim_done.exists() and not self.needs_run([c1, c2], [self.read1, self.read2]):
            trim_outputs = [trim_done]
        self.step("run_trimmomatic", trim_outputs, [self.read1, self.read2], self.run_trimmomatic)
        self.step("run_rcorrector", [c1, c2], [t1, t2], self.run_rcorrector)
        self.reclaim_trimmed_reads()
        # Every stage from here through strandeval reads c1/c2, so they can't
        # be replaced by their .gz until cleanup() -- but the compression
        # itself starts now, behind the assemblers, rather than being paid
        # for serially once the run is otherwise over.
        self.compress_async(c1)
        self.compress_async(c2)

    def run_assemblers(self):
        """Build the four assemblies this pipeline is named for.

        oyster.py's own half of the work, and the only half that is specific
        to a particular set of assemblers -- everything downstream of
        run_filtershort treats them as an unordered set of inputs.
        """
        c1, c2 = self.cor1(), self.cor2()
        trinity_fa = self.assembly_fasta(TRINITY)
        phase1_done = self.trinity_phase1_done()
        sp75 = self.assembly_fasta(SPADES75)
        spauto = self.assembly_fasta(SPADES_AUTO)
        ta = self.assembly_fasta(TRANSABYSS)
        diamond_ta = self.diamond_txt(TRANSABYSS)
        diamond_sp75 = self.diamond_txt(SPADES75)
        diamond_spauto = self.diamond_txt(SPADES_AUTO)

        # Two sequential stage-pairings rather than one lane split across all
        # four assemblers for the whole run -- see TRINITY_PHASE1_SHARE and
        # TRINITY_PHASE2_SHARE above for why each pairing is chosen the way
        # it is.
        phase1_cpu = max(1, round(self.cpu * TRINITY_PHASE1_SHARE))
        phase1_mem = max(1, round(self.mem * TRINITY_PHASE1_SHARE))
        spades_cpu = max(1, self.cpu - phase1_cpu)
        spades_mem = max(1, self.mem - phase1_mem)

        self.seed_trinity_phase1_sentinel()

        phase2_cpu = max(1, round(self.cpu * TRINITY_PHASE2_SHARE))
        phase2_mem = max(1, round(self.mem * TRINITY_PHASE2_SHARE))
        transabyss_cpu = max(1, self.cpu - phase2_cpu)
        # Not split by TRINITY_PHASE2_SHARE: Trans-ABySS's memory footprint
        # doesn't shrink along with its CPU share the way SPAdes/Phase 1's
        # do, so it keeps the same generous share Stage A used rather than
        # being squeezed down with its CPU. Revisit if this run OOMs or if
        # it turns out to be more mem than Trans-ABySS actually needs.
        transabyss_mem = spades_mem

        def _lane_failed(lane_name, e):
            # ThreadPoolExecutor.__exit__ (below) calls shutdown(wait=True) even
            # while an exception unwinds, so the *other* lane's still-running
            # future keeps blocking the run's exit -- Trinity in particular can
            # still be tens of hours out. Log the failure the moment it happens
            # rather than leaving it silent until the other lane finally finishes.
            print(
                f"\n*** [{lane_name}] lane failed ({e}) -- other lane keeps running; "
                "this run will still abort once it finishes ***",
                flush=True,
            )

        def trinity_phase1_lane():
            try:
                self.step(
                    "run_trinity_phase1", [phase1_done], [c1],
                    partial(self.run_trinity_phase1, cpu=phase1_cpu, mem=phase1_mem),
                )
            except Exception as e:
                _lane_failed("run_trinity_phase1", e)
                raise

        def spades_lane():
            # spadesauto (the lower-k assembly, whose smaller k dominates its
            # cost) first, since it has been the slower of the two --
            # historically the fixed k=55 run against the fixed k=75 one.
            # diamond_{spadesauto, spades75} depend only on their own assembly
            # (not on Trinity or the orthofuser merge below), so each fires as
            # soon as its assembly is done instead of waiting for the merge.
            try:
                for step_name, outputs, inputs, assemble, diamond_name, diamond_out in (
                    ("run_spadesauto", [spauto], [c1, c2], self.run_spadesauto, "spadesauto", diamond_spauto),
                    ("run_spades75", [sp75], [c1, c2], self.run_spades75, "spades75", diamond_sp75),
                ):
                    self.step(step_name, outputs, inputs, partial(assemble, cpu=spades_cpu, mem=spades_mem))
                    query = outputs[0]
                    self.compress_async(query)
                    self.step(
                        f"diamond_{diamond_name}", [diamond_out], [query],
                        partial(self.run_diamond_one, query, diamond_out, cpu=spades_cpu),
                    )
            except Exception as e:
                _lane_failed("spades", e)
                raise

        print(f"\n=== Stage A: run_trinity_phase1 ({phase1_cpu} cpu) || spadesauto/spades75 ({spades_cpu} cpu) -- start {self._ts()} ===")
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
            for f in concurrent.futures.as_completed([ex.submit(trinity_phase1_lane), ex.submit(spades_lane)]):
                f.result()
        print(f"=== Stage A done -- {self._ts()} ===")

        def trinity_phase2_lane():
            try:
                self.step(
                    "run_trinity_phase2", [trinity_fa], [phase1_done],
                    partial(self.run_trinity_phase2, cpu=phase2_cpu, mem=phase2_mem),
                )
                self.compress_async(trinity_fa)
            except Exception as e:
                _lane_failed("run_trinity_phase2", e)
                raise

        def transabyss_lane():
            try:
                self.step(
                    "run_transabyss", [ta], [c1, c2],
                    partial(self.run_transabyss, cpu=transabyss_cpu, mem=transabyss_mem),
                )
                self.compress_async(ta)
                self.step(
                    "diamond_transabyss", [diamond_ta], [ta],
                    partial(self.run_diamond_one, ta, diamond_ta, cpu=transabyss_cpu),
                )
            except Exception as e:
                _lane_failed("transabyss", e)
                raise

        print(f"\n=== Stage B: run_trinity_phase2 ({phase2_cpu} cpu) || transabyss ({transabyss_cpu} cpu) -- start {self._ts()} ===")
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
            for f in concurrent.futures.as_completed([ex.submit(trinity_phase2_lane), ex.submit(transabyss_lane)]):
                f.result()
        print(f"=== Stage B done -- {self._ts()} ===")

    def merge_and_report(self, pipeline_start):
        """Fuse the assemblies into one, then score and report on the result.

        Everything here is generic over `self.assemblies`: it is the half
        chowder.py reuses wholesale for assemblies it did not build.
        """
        c1, c2 = self.cor1(), self.cor2()
        assembly_fastas = self.assembly_fasta_paths()
        short_fastas = self.short_fasta_paths()
        orthofuser_done = self.orthofuse_dir / "orthofuser.done"
        merged_fasta = self.orthofuse_dir / "merged.fasta"
        merged_csv = self.orthofuse_dir / "merged" / "assemblies.csv"
        good_list = self.orthofuse_dir / f"good.{self.runout}.list"
        orthomerged_fasta = self.assemblies_dir / f"{self.runout}.orthomerged.fasta"
        diamond_outs = [o for _, o in self.diamond_jobs()]
        diamond_orthomerged = self.diamond_dir / f"{self.runout}.orthomerged.diamond.txt"
        uniq_outs = [self.unique_txt(a) for a in self.report_order]
        list1 = self.diamond_dir / f"{self.runout}.list1"
        list2 = self.diamond_dir / f"{self.runout}.list2"
        list3 = self.diamond_dir / f"{self.runout}.list3"
        list5 = self.diamond_dir / f"{self.runout}.list5"
        list6 = self.diamond_dir / f"{self.runout}.list6"
        list7 = self.diamond_dir / f"{self.runout}.list7"
        newbies = self.diamond_dir / f"{self.runout}.newbies.fasta"
        working_orthomerged = self.assemblies_working / f"{self.runout}.orthomerged.fasta"
        orp_intermediate = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        orp_diamond_txt = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        unique_orp_done = self.assemblies_working / f"{self.runout}.unique.ORP.done"
        ortho_idx = self.quants_dir / f"{self.runout}.ortho.idx"
        quant_sf = self.quants_dir / f"salmon_orthomerged_{self.runout}" / "quant.sf"
        filter_done = self.assemblies_dir / f"{self.runout}.filter.done"
        low_txt = self.assemblies_working / f"{self.runout}.LOWEXP.txt"
        high_txt = self.assemblies_working / f"{self.runout}.HIGHEXP.txt"
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        busco_done = self.reports_dir / f"{self.runout}.busco.done"
        transrate_csv = self.reports_dir / f"transrate_{self.runout}" / "assemblies.csv"
        strandeval_done = self.reports_dir / f"{self.runout}.strandeval.done"
        qualreport_done = self.reports_dir / f"qualreport.{self.runout}.done"
        cleanup_done = self.reports_dir / f"{self.runout}.cleanup.done"

        self.step("run_filtershort", short_fastas, assembly_fastas, self.run_filtershort)

        def orthofuser_branch(cpu=None, mem=None):
            self.step("mask_search_input", self.search_fasta_paths(), short_fastas,
                      self.write_search_inputs)
            # mem reaches run_orthofuser because OrthoFinder's search
            # concurrency is now capped against it; left unforwarded it would
            # cap against the whole machine while holding half of it.
            self.step("run_orthofuser", [orthofuser_done], self.search_fasta_paths(),
                      partial(self.run_orthofuser, cpu=cpu, mem=mem))

        def merge_branch(cpu=None, mem=None):
            self.step("merge", [merged_fasta], short_fastas, self.merge)
            self.step("orthotransrate", [merged_csv], [merged_fasta, c1, c2],
                      partial(self.orthotransrate, cpu=cpu, mem=mem))

        # run_orthofuser and merge->orthotransrate are independent chains that
        # both only need short_fastas; they join at makeorthout below.
        self.run_parallel(
            [
                ("orthofuser_branch", [orthofuser_done], short_fastas, orthofuser_branch),
                ("merge_branch", [merged_csv], short_fastas, merge_branch),
            ],
            max_workers=self.max_parallel,
        )
        self.step("makeorthout", [good_list], [orthofuser_done, merged_csv], self.makeorthout)
        self.step("orthofusing", [orthomerged_fasta], [good_list, merged_fasta], self.orthofusing)

        # Every assembly needs a diamond pass, and under oyster.py most of
        # them already had one: the assembler lanes fire each assembly's
        # diamond the moment that assembler returns, rather than leaving all
        # four to queue up here. Those steps are up to date by now and skip;
        # what is genuinely left is orthomerged, which depends on the merge
        # stage just above, and Trinity, whose lane only just finished. A run
        # that brought its own assemblies had no lanes, so all of them run
        # here -- which is why this is a loop over the set and not the two
        # named steps it used to be.
        print("\n\n\n\n Starting diamond \n\n\n\n")
        self.step(
            "diamond_orthomerged", [diamond_orthomerged], [orthomerged_fasta],
            partial(self.run_diamond_one, orthomerged_fasta, diamond_orthomerged),
        )
        for a in self.diamond_priority:
            fasta, out = self.assembly_fasta(a), self.diamond_txt(a)
            self.step(
                f"diamond_{a.diamond_label}", [out], [fasta],
                partial(self.run_diamond_one, fasta, out),
            )
        self.step("diamond_uniq", uniq_outs, diamond_outs, self.diamond_uniq)
        self.step("make_list1", [list1], [diamond_orthomerged], self.make_list1)
        self.step("make_list2", [list2], [self.diamond_txt(a) for a in self.assemblies], self.make_list2)
        self.step("make_list3", [list3], [list1, list2], self.make_list3)
        self.step("make_list5", [list5], [list3] + [self.diamond_txt(a) for a in self.diamond_priority], self.make_list5)
        self.step("make_list6", [list6], [orthomerged_fasta], self.make_list6)
        self.step("make_list7", [list7], [list6, list5], self.make_list7)
        self.step("posthack", [newbies, working_orthomerged], [list7], self.posthack)
        self.step("cdhit", [orp_intermediate], [working_orthomerged], self.cdhit)

        # orp_diamond is the same CPU-bound diamond blastx as above, paired
        # here with salmon_branch which is tiny (~2s); halving orp_diamond's
        # CPU to overlap with it costs more than the overlap saves, so both
        # run sequentially at full CPU instead.
        self.step("orp_diamond", [orp_diamond_txt], [orp_intermediate], self.orp_diamond)
        self.step("orp_uniq", [unique_orp_done], [orp_diamond_txt], self.orp_uniq)
        salmon_stamp = self.stamp_tool_version("orp", "salmon", self.quants_dir / "salmon.version")
        self.step("salmon_index", [ortho_idx], [orp_intermediate, salmon_stamp], self.salmon_index)
        self.step("salmon", [quant_sf], [ortho_idx, c1, c2], self.salmon)
        self.step("filter", [filter_done], [orp_intermediate, quant_sf, orp_diamond_txt], self.filter_tpm)
        self.step(
            "secondfilter", [orp_fasta],
            [filter_done, low_txt, high_txt, orp_intermediate, quant_sf, orp_diamond_txt],
            self.secondfilter,
        )
        # BUSCO dwarfs transrate/strandeval (minutes vs. seconds); splitting its
        # CPUs to overlap with them would slow it down far more than the
        # overlap could ever save, so it gets the full --cpu/--busco-threads
        # budget to itself. transrate and strandeval are independent of each
        # other and cheap, so they still run as a small concurrent pair.
        self.step("busco", [busco_done], [orp_fasta], self.busco)
        self.run_parallel(
            [
                ("transrate", [transrate_csv], [orp_fasta, c1, c2], self.transrate),
                ("strandeval", [strandeval_done], [orp_fasta], self.strandeval),
            ],
            max_workers=self.max_parallel,
        )
        self.step("reportgen", [qualreport_done], [unique_orp_done, orp_fasta], self.reportgen)
        # Last, because it deletes inputs several of the steps above declare.
        self.step("cleanup", [cleanup_done], [qualreport_done], self.cleanup)

        self.timing_report(int(time.time() - pipeline_start))


def parse_args():
    p = argparse.ArgumentParser(description="Python port of oyster.mk - the Oyster River Protocol pipeline.")
    version = (HERE / "version.txt").read_text().strip()
    p.add_argument("--version", action="version", version=f"Oyster River Protocol {version}")
    p.add_argument("--read1", required=True, help="path to R1 fastq(.gz)")
    p.add_argument("--read2", required=True, help="path to R2 fastq(.gz)")
    p.add_argument("--mem", type=int, default=110, help="memory in GB (default: 110)")
    p.add_argument("--cpu", type=int, default=16, help="CPU threads (default: 16)")
    p.add_argument("--busco-threads", type=int, default=None, help="BUSCO threads (default: same as --cpu)")
    p.add_argument("--runout", default="USER_RUN", help="run name prefix (default: USER_RUN)")
    p.add_argument("--strand", choices=["RF", "FR", ""], default="", help="strand-specificity (default: unstranded)")
    p.add_argument("--lineage", default="eukaryota_odb12.2", help="BUSCO lineage (default: eukaryota_odb12.2)")
    p.add_argument("--normalize-reads", action="store_true", help="let Trinity normalize reads (default: off, i.e. --no_normalize_reads)")
    p.add_argument("--tpm-filt", type=float, default=0, help="TPM filter threshold (default: 0)")
    p.add_argument("--spades1-kmer", type=parse_kmer_spec, default="auto",
                   help="rnaSPAdes k-mer(s) for the spadesauto assembly: 'auto' to let "
                        "rnaSPAdes pick its documented default pair from read length, or a "
                        "comma-separated list of odd sizes under 128 (default: auto)")
    p.add_argument("--spades2-kmer", type=parse_kmer_spec, default="75",
                   help="rnaSPAdes k-mer(s) for the spades75 assembly: 'auto' or a "
                        "comma-separated list of odd sizes under 128 (default: 75)")
    p.add_argument("--transabyss-kmer", type=int, default=32, help="Trans-ABySS k-mer (default: 32)")
    p.add_argument(
        "--orthofinder-searches", type=int, default=None, metavar="N",
        help="how many of OrthoFinder's n_assemblies^2 diamond searches may run "
             "at once. Default: as many as --mem allows, sized off the largest "
             "search input. --cpu is split across them, so lowering this costs "
             "memory rather than cores. Lower it if diamonds are OOM-killed "
             "(returncode -9 in the log)",
    )
    p.add_argument(
        "--orthofinder-analysis", type=int, default=None, metavar="N",
        help="OrthoFinder's -a, the workers for its algorithm phase after the "
             "searches. Default: min(--cpu/8, 16, number of assemblies). Not "
             "sized by memory -- that phase reads only the search output and has "
             "not been measured; raise it if it stalls, lower it if it runs a "
             "node out of memory",
    )
    p.add_argument(
        "--max-parallel", type=int, default=2,
        help="max concurrent jobs within each independent stage that benefits from "
             "it (orthofuser vs. merge/orthotransrate; transrate vs. strandeval), "
             "splitting --cpu/--mem across however many run at once; the 4 "
             "assemblers instead run as two sequential stage-pairings (see "
             "TRINITY_PHASE1_SHARE/TRINITY_PHASE2_SHARE), unaffected by this flag; "
             "other CPU-bound stages "
             "that don't benefit (diamond, orp_diamond, salmon, busco) always run "
             "sequentially at full --cpu; 1 disables concurrency for the remaining "
             "stages entirely (default: 2)",
    )
    p.add_argument(
        "--keep-intermediates", action="store_true",
        help="keep every file a run produces: skips both the end-of-run cleanup "
             "(orthofuse/, quants/, diamond/, the working assemblies) and the "
             "reclaim of the trimmed reads, and leaves the four assemblies and "
             "the corrected reads uncompressed. For debugging a run (default: off)",
    )
    p.add_argument(
        "--pytransrate-args", default="",
        help="extra arguments passed verbatim to both pytransrate runs, as one "
             "quoted string, e.g. --pytransrate-args '--location-size 5'. For "
             "the snap index tuning a large merge needs: --location-size skips "
             "the sweep when you already know four byte locations will not hold "
             "the genome. Note --padding is not the memory lever it looks like: "
             "snap writes it as N and skips seeds containing N, so it grows the "
             "1 byte/base genome array and nothing else -- on a 5.4M-contig, "
             "5.7 Gbp merge, dropping it entirely saved 1%% of the branch, while "
             "real sequence alone still exceeded the four-byte ceiling. Run "
             "`pytransrate --help` for the full set (default: none)",
    )
    p.add_argument("--dir", default=None, help="working directory (default: current directory)")
    return p.parse_args()


def main():
    line_buffer_stdio()
    args = parse_args()
    pipeline = Pipeline(args)
    try:
        pipeline.main()
    except subprocess.CalledProcessError as e:
        sys.exit(f"\n*** step failed: {' '.join(str(c) for c in e.cmd)} (exit {e.returncode}) ***")
    except FileNotFoundError as e:
        sys.exit(f"\n*** required command not found: {e.filename} ***")
    finally:
        # The pool's threads are not daemons, so a run that dies mid-stage
        # would otherwise sit at interpreter exit with no explanation of what
        # it's waiting for. The .gz files themselves are still worth
        # finishing: nothing has been deleted on this path.
        pipeline.finish_compression()


if __name__ == "__main__":
    main()
