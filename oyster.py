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
    `<runout>.unique.sphigh.txt` and not `...spadeshigh...`. Renaming any of
    them would silently invalidate every resumable run directory that
    exists, so they are carried as data instead of being derived.

    The SPAdes assemblies are a deliberate exception. They changed which k
    values they run at, so a directory holding the old `spades55.fasta` or
    `spades75.fasta` no longer describes what this code would produce from
    the same reads. Keeping the old names would let step() find those stale
    files and skip the assemblies, reporting an auto-k run while serving k=55
    output -- the same silent wrongness the paragraph above is guarding
    against, arriving from the other direction. Renaming costs a
    re-assembly; not renaming costs a wrong answer that looks right.
    """

    fasta_name: str      # assemblies/<runout>.<fasta_name>
    diamond_label: str   # assemblies/diamond/<runout>.<diamond_label>.diamond.txt
    unique_label: str    # assemblies/diamond/<runout>.unique.<unique_label>.txt
    report_label: str    # reportgen's "UNIQUE GENES <report_label>" line


SPADES_AUTO = Assembly("spadesauto.fasta", "spadesauto", "spauto", "SPADESAUTO")
SPADES_HIGH = Assembly("spadeshigh.fasta", "spadeshigh", "sphigh", "SPADESHIGH")
TRANSABYSS = Assembly("transabyss.fasta", "transabyss", "transabyss", "TRANSABYSS")
TRINITY = Assembly("trinity.Trinity.fasta", "trinity", "trinity", "TRINITY")

# oyster.py's own four assemblers, in the three different orders oyster.mk
# used them in. The three are not interchangeable and two of them reach the
# assembly, so they are spelled out rather than sorted:
#
#   ASSEMBLY_ORDER    concatenation order. Sets contig order in
#                     shuck/pool.fasta and in posthack's `cat` of the
#                     assemblies, which flows through to cd-hit-est -- where
#                     input order breaks length ties and so decides which
#                     representative survives into .ORP.fasta.
#   DIAMOND_PRIORITY  search order. build_list5.py keeps the *first* diamond
#                     hit per gene in the order it is given the files, so
#                     this is a preference ranking between assemblies for
#                     the genes the two-track selection missed.
#   REPORT_ORDER      the order the UNIQUE GENES lines appear in
#                     reports/qualreport.<run>. Cosmetic, but people diff
#                     those reports across runs.
ASSEMBLY_ORDER = (SPADES_AUTO, SPADES_HIGH, TRANSABYSS, TRINITY)
DIAMOND_PRIORITY = (TRANSABYSS, SPADES_HIGH, SPADES_AUTO, TRINITY)
REPORT_ORDER = (TRINITY, SPADES_AUTO, SPADES_HIGH, TRANSABYSS)

# Preflight. Everything here is shelled out to at some point in a full run,
# and finding it missing hours in -- at score_pool, or at the assembler that
# was going to run overnight -- is the thing this list exists to prevent.
# snap-aligner is on it because pyTransRate maps with it.
SPADES_TOOL = ("orp_spades", "rnaspades.py", "SPADES")
TRINITY_TOOL = ("orp_trinity", "Trinity", "TRINITY")
TRANSABYSS_TOOL = ("orp_transabyss", "transabyss", "TRANSABYSS")
# The three an entry point that doesn't assemble has no use for.
ASSEMBLER_TOOLS = (SPADES_TOOL, TRINITY_TOOL, TRANSABYSS_TOOL)
TRIMMOMATIC_TOOL = ("orp", "trimmomatic", "TRIMMOMATIC")
RCORRECTOR_TOOL = ("orp", "run_rcorrector.pl", "RCORRECTOR")
# The two a run handed already-corrected reads has no use for.
READ_PREP_TOOLS = (TRIMMOMATIC_TOOL, RCORRECTOR_TOOL)
#: The pyTransRate this pipeline needs, checked at preflight rather than
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
#: What orp_env.yml installs, and so what the upgrade hint above installs: an
#: env brought up to date by hand should match one built fresh. Can sit
#: above the minimum -- 2.2.2 only adds logging -- and moves with that file.
PYTRANSRATE_PINNED_VERSION = "2.2.2"

#: --max-memory's own spellings in pyTransRate, both of which it accepts. A
#: user who set one in --pytransrate-args means it, so pytransrate_memory_args
#: stands aside rather than passing the flag twice.
PYTRANSRATE_MEMORY_FLAGS = ("--max-memory", "--mem")

CHECK_TOOLS = (
    ("orp", "salmon", "SALMON"),
    ("orp", "pytransrate", "PYTRANSRATE"),
    ("orp", "snap-aligner", "SNAP-ALIGNER"),
    ("orp", "diamond", "DIAMOND"),
    ("orp", "cd-hit-est", "CD-HIT-EST"),
    # twotrack_select.py's distinct-copy test.
    ("orp", "blastn", "BLASTN"),
    ("orp", "makeblastdb", "MAKEBLASTDB"),
    ("orp", "samtools", "SAMTOOLS"),
    ("orp_busco", "busco", "BUSCO"),
    # strandeval maps a read sample with these, out of the Trinity env.
    ("orp_trinity", "bwa", "BWA"),
    ("orp_trinity", "seqtk", "SEQTK"),
    ("orp_trinity", "hist", "HIST (bashplotlib)"),
    SPADES_TOOL,
    TRINITY_TOOL,
    TRIMMOMATIC_TOOL,
    TRANSABYSS_TOOL,
    RCORRECTOR_TOOL,
)

# A job run *beside* another (see run_beside) gets this many threads on top
# of the main job's, never more than SIDE_JOB_MAX_CPU, plus at most
# SIDE_JOB_MAX_MEM GB taken out of the main job's memory budget.
SIDE_JOB_CPU_SHARE = 0.25
SIDE_JOB_MAX_CPU = 8
SIDE_JOB_MAX_MEM = 16

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

def line_buffer_stdio():
    """Make our own output appear where it happened in a redirected log.

    Python block-buffers stdout in 4-8 KB chunks when it is not a terminal,
    which on a cluster it never is. Every tool we launch, though, inherits
    the same file descriptor and writes to it directly, unbuffered. So the
    pipeline's own narrative -- the banner, the `=== step -- start ===`
    lines, the `+ <command>` echoes, the retry warnings -- sits in our
    buffer while hours of OrthoFinder and pyTransRate output stream past it,
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
#: NumReads. pyTransRate rejects any other count as a version mismatch.
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

    Chunks rather than lines because this runs on pool.fasta, which is
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
    killed during that pass leaves a short file -- and pyTransRate reuses
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
MIN_SPADES_K = 15


def parse_kmer_spec(value):
    """Parse a --spadesN-kmer value into one of three forms:

    - None                -- "auto": run_spades() omits -k entirely so
      rnaSPAdes picks its own two k values from the observed read length
      (approximately 1/3 and 1/2 of the maximum). This is the documented
      default and the configuration the rnaSPAdes authors recommend: they
      warn that smaller k-mer sizes typically produce chimeric transcripts,
      and ORP forcing a single k per run has always deviated from it. Two
      single-k runs are not one two-k run.
    - a list of ints      -- absolute k values, used as given.
    - a list of floats    -- fractions of maximum read length ("60%,75%"),
      resolved per-dataset by resolve_kmers() once the reads have been read.

    Fractions exist because an absolute k means completely different things at
    different read lengths: k=75 is half the read length at 150bp but 74% of
    it at 101bp, which is why ORP's historical fixed k values behaved so
    differently across datasets.

    Absolute and fractional entries cannot be mixed in one spec -- the two
    can't be ordered against each other until read length is known, and a spec
    that silently reorders itself per dataset would be worse than an error.
    rnaSPAdes' own constraints (odd, below 128, ascending, distinct) are
    enforced here for absolute values and in resolve_kmers() for fractions,
    rather than after the reads have already been trimmed and corrected.
    """
    text = value.strip()
    if text.lower() == "auto":
        return None
    parts = [x.strip() for x in text.split(",")]
    if not parts or any(not x for x in parts):
        raise argparse.ArgumentTypeError(f"{value!r}: empty k-mer entry")

    if all(x.endswith("%") for x in parts):
        try:
            fracs = [float(x[:-1]) / 100.0 for x in parts]
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"{value!r}: expected percentages like '60%,75%'"
            )
        for f in fracs:
            if not 0.0 < f < 1.0:
                raise argparse.ArgumentTypeError(
                    f"{value!r}: k-mer percentages must be above 0% and below 100% "
                    f"(got {f * 100:g}%)"
                )
        if fracs != sorted(fracs) or len(set(fracs)) != len(fracs):
            raise argparse.ArgumentTypeError(
                f"{value!r}: k-mer percentages must be distinct and in ascending order"
            )
        return fracs

    if any(x.endswith("%") for x in parts):
        raise argparse.ArgumentTypeError(
            f"{value!r}: cannot mix absolute k values and percentages in one spec"
        )

    try:
        kmers = [int(x) for x in parts]
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{value!r}: expected 'auto', a comma-separated list of integers, "
            "or a comma-separated list of percentages like '60%,75%'"
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


def resolve_kmers(spec, max_read_len, label=""):
    """Turn a parsed spec into the concrete k list to hand rnaSPAdes.

    None and absolute lists pass straight through. Fractions are multiplied by
    max_read_len and rounded *down* to the nearest odd number -- down rather
    than nearest, so a fraction can never round up into the read-length wall,
    where a read yields too few k-mers to be worth anything (at k=L a read
    yields one k-mer, and rnaSPAdes errors outright).

    Results are clamped to rnaSPAdes' 127 ceiling, which long reads hit
    easily: 75% of a 250bp read is 187. Clamping can collapse two fractions
    onto the same k, so duplicates are dropped and the collapse is reported
    rather than passed to rnaSPAdes as an invalid repeated -k.
    """
    if spec is None or not spec or isinstance(spec[0], int):
        return spec

    kmers, notes = [], []
    for f in spec:
        k = int(f * max_read_len + 0.5)
        if k % 2 == 0:
            k -= 1
        if k > MAX_SPADES_K:
            notes.append(f"{f * 100:g}% of {max_read_len} = {k}, clamped to {MAX_SPADES_K}")
            k = MAX_SPADES_K
        if k < MIN_SPADES_K:
            sys.exit(
                f"\n\n\n\n ERROR: {label}{f * 100:g}% OF A {max_read_len} BP READ IS K={k},\n"
                " WHICH IS TOO SMALL TO ASSEMBLE WITH. USE A LARGER PERCENTAGE OR AN\n"
                " EXPLICIT ODD K VALUE. \n\n\n\n"
            )
        if k not in kmers:
            kmers.append(k)
        else:
            notes.append(f"{f * 100:g}% of {max_read_len} also resolves to {k}, dropped as a duplicate")
    for n in notes:
        print(f"[kmer] {label}{n}")
    return kmers


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
        self.spades2_kmer = getattr(args, "spades2_kmer", [0.60, 0.75])
        self.transabyss_kmer = getattr(args, "transabyss_kmer", 32)
        self.read1 = Path(args.read1)
        self.read2 = Path(args.read2)
        self.runout = args.runout
        self.lineage = args.lineage
        self.strand = getattr(args, "strand", "")
        self.normalize_reads = getattr(args, "normalize_reads", False)
        self.tpm_filt = args.tpm_filt
        self.max_parallel = max(1, args.max_parallel)
        self.no_cleanup = args.no_cleanup
        # oyster.py spells it --trimmed-corrected-reads, chowder.py
        # --corrected-reads; both mean trimmomatic and rcorrector are done.
        self.corrected_reads = getattr(args, "corrected_reads", False)
        # Appended to both pyTransRate invocations. shlex so a value can be
        # quoted, and so the flags arrive as separate argv entries rather
        # than one string pyTransRate would reject.
        self.pytransrate_args = shlex.split(getattr(args, "pytransrate_args", "") or "")

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
        # The fasta the diamond database was built from: twotrack_select.py
        # reads protein lengths from it.
        self.sprot_fasta = self.makedir / "software" / "diamond" / "uniprot_sprot.fasta"

        self.rcorr_dir = self.dir / "rcorr"
        self.assemblies_dir = self.dir / "assemblies"
        self.assemblies_working = self.assemblies_dir / "working"
        self.diamond_dir = self.assemblies_dir / "diamond"
        self.reports_dir = self.dir / "reports"
        self.shuck_dir = self.dir / "shuck" / self.runout
        self.shuck_working = self.shuck_dir / "working"
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
        clear_pytransrate_outdir.
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
        a run.
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
        if self.no_cleanup or not path.exists():
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
        if self.no_cleanup:
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
        shuck tree (the pooled fasta and pyTransRate's scoring of it --
        normally the largest directory in the run), the diamond hits and the list1-list7 set algebra built from
        them, the salmon index and quantification, and the chain of working
        assemblies between shuck and .ORP.fasta. Every number any of
        it contributed is already in reports/qualreport.<run>.
        """
        if self.no_cleanup:
            # No cleanup.done either, so a later run of the same command
            # without --no-cleanup finds this step pending and does it then.
            print("[cleanup] --no-cleanup given; leaving intermediates in place")
            return

        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        freed = 0
        kept = [f"{self._rel(orp_fasta)}  (the assembly)",
                f"{self._rel(self.reports_dir)}/  (all reports)"]
        removed = []

        for src in [self.cor1(), self.cor2()] + self.assembly_fasta_paths():
            gz = src.with_suffix(src.suffix + ".gz")
            if src.is_symlink():
                # Not ours to reclaim: under --corrected-reads (chowder.py)
                # or --trimmed-corrected-reads (oyster.py) the corrected pair
                # points straight at the user's own reads, and compressing or
                # unlinking those is not what "reclaim this run's
                # intermediates" means.
                kept.append(f"{self._rel(src)}  (symlink to a file this run did not create)")
                continue
            if self.corrected_reads and src in (self.cor1(), self.cor2()):
                # A plain copy of the user's gzipped reads, made only so the
                # assemblers see content matching the name (see
                # use_corrected_reads). The .gz it came from is still where
                # the user left it, so there is nothing here worth keeping.
                if src.exists():
                    size = path_size(src)
                    src.unlink()
                    freed += size
                    removed.append(f"{self._rel(src)}  ({human_size(size)}, "
                                   "uncompressed copy of your reads)")
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
            self.dir / "shuck",
            # The same tree under its pre-4.1.0-dev8 name, left behind when a
            # run directory from before the rename was resumed under it.
            self.dir / "orthofuse",
            self.assemblies_working,
            self.diamond_dir,
            self.quants_dir,
            # Trinity's --full_cleanup normally removes this itself; a run
            # that was interrupted and resumed can still leave it behind.
            self.trinity_out_dir(),
            self.assemblies_dir / f"{self.runout}.shucked.fasta",
            self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta",
            # cd-hit-est's cluster report, written beside its -o; nothing
            # reads it.
            self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta.clstr",
            self.assemblies_dir / f"{self.runout}.ORP.diamond.txt",
            self.assemblies_dir / f"{self.runout}.flagstat",
            self.assemblies_dir / f"{self.runout}.filter.done",
            self.trinity_phase1_done(),
            # The rest are normally deleted by the step that made them, and
            # are only here after a --no-cleanup run (or an interrupted one).
            self.assemblies_dir / f"{self.runout}.transabyss",
            # Directories only: a chowder assembly may be labelled spades_k*.
            *(d for d in self.assemblies_dir.glob(f"{self.runout}.spades_k*") if d.is_dir()),
            *self.assemblies_dir.glob(f"{self.runout}.*gene_trans_map"),
            *self.strandeval_scratch(),
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
        marker = self.reports_dir / f"{self.runout}.cleanup.done"
        marker.write_text(text)
        print(f"\nReclaimed {human_size(freed)} of intermediate files; "
              f"what was kept and removed is in {marker}")

    # -- resumability ------------------------------------------------------

    def step_marker(self, name):
        """A file that exists while step `name` runs, and is left behind if it
        never finishes.

        needs_run() trusts any output that exists and is newer than its
        inputs, which is the wrong question for a step that died partway:
        trimmomatic, rcorrector and diamond all write their output as they
        go, so a run OOM-killed, killed at walltime or out of retries leaves
        a truncated file with a fresh mtime, and every later resume skips the
        step and builds on it. For rcorrector that is the whole assembly run
        off part of the reads, with the trimmed reads deleted behind it.

        Exit status can't be relied on to clear such a file up -- SIGKILL
        and a walltime kill never reach Python -- so the evidence is written
        before the work instead: a step whose marker is still there did not
        finish, whatever its outputs look like.
        """
        return self.reports_dir / f".{self.runout}.{name}.running"

    def is_pending(self, name, outputs, inputs) -> bool:
        """needs_run(), plus a re-run of any step an earlier attempt left unfinished."""
        marker = self.step_marker(name)
        if marker.exists():
            print(f"[{name}] an attempt started {self._ts(marker.stat().st_mtime)} "
                  "and never finished; re-running it rather than trusting what it "
                  "left behind")
            return True
        return self.needs_run(outputs, inputs)

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
        stale <run>.shucked.idx would be kept, salmon_index skipped, and
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
        if not self.is_pending(name, outputs, inputs):
            print(f"[{name}] up to date, skipping")
            return
        marker = self.step_marker(name)
        marker.touch()
        start = time.time()
        print(f"\n=== {name} -- start {self._ts(start)} ===")
        func()
        # Only on success: a step that raised keeps its marker, so the next
        # run redoes it (see step_marker).
        marker.unlink()
        elapsed = int(time.time() - start)
        print(f"=== {name} -- done {self._ts()} ({elapsed}s) ===")
        if timed:
            self._record_timing(name, elapsed, start)

    def _record_timing(self, name, elapsed, start):
        with self._timing_lock:
            self.steps.append((name, elapsed))
            with open(self.timing_log, "a") as f:
                f.write(f"{name}\t{elapsed}\t{self._ts(start)}\n")

    def side_job_budget(self, max_cpu=SIDE_JOB_MAX_CPU):
        """(cpu, mem) for a job run beside another; see run_beside."""
        cpu = max(1, min(max_cpu, int(self.cpu * SIDE_JOB_CPU_SHARE)))
        mem = max(1, min(SIDE_JOB_MAX_MEM, self.mem // 4))
        return cpu, mem

    def _run_job(self, name, func, cpu, mem):
        """One run_beside job: step()'s marker, banner and timing, with a budget."""
        marker = self.step_marker(name)
        marker.touch()
        start = time.time()
        print(f"\n=== {name} (cpu={cpu}, mem={mem}G) -- start {self._ts(start)} ===")
        func(cpu=cpu, mem=mem)
        marker.unlink()
        elapsed = int(time.time() - start)
        print(f"=== {name} -- done {self._ts()} ({elapsed}s) ===")
        self._record_timing(name, elapsed, start)

    def run_beside(self, main, side, side_cpu, side_mem):
        """Run two independent jobs at once, `side` on a few threads beside `main`.

        Each job is (name, outputs, inputs, func), and func is called as
        func(cpu=, mem=). `main` keeps every one of --cpu's cores and gives
        up `side_mem` GB of --mem; `side` runs on `side_cpu` threads on top.

        This replaced an even split of the cores, which was the wrong shape
        for every pair it was used on: one job of the pair is short
        (strandeval, minutes; a diamond pass) and the other long (the final
        pyTransRate scoring; score_pool, hours), so once the short one finished its half of the
        machine sat idle while the long one carried on at half speed.
        Oversubscribing by a few threads instead costs the main job a little
        while both run and nothing afterwards, and the side job mostly fills
        cores pyTransRate leaves idle in its serial phases (snap's index
        build, salmon). Memory is never oversubscribed: side_mem comes out
        of main's budget, which is what pyTransRate sizes its workers to.

        With --max-parallel 1, or when only one of them is pending, they run
        one after the other with the whole machine each.
        """
        pending = []
        for name, outputs, inputs, func in (main, side):
            if self.is_pending(name, outputs, inputs):
                pending.append((name, func))
            else:
                print(f"[{name}] up to date, skipping")
        if len(pending) < 2 or self.max_parallel < 2:
            for name, func in pending:
                self._run_job(name, func, self.cpu, self.mem)
            return
        side_cpu = max(1, min(side_cpu, self.cpu))
        side_mem = max(1, min(side_mem, self.mem - 1))
        main_mem = max(1, self.mem - side_mem)
        (main_name, main_func), (side_name, side_func) = pending
        print(f"\n=== {main_name} (cpu={self.cpu}, mem={main_mem}G) with {side_name} "
              f"beside it (cpu={side_cpu}, mem={side_mem}G) ===")
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
            futures = [ex.submit(self._run_job, main_name, main_func, self.cpu, main_mem),
                       ex.submit(self._run_job, side_name, side_func, side_cpu, side_mem)]
            for future in concurrent.futures.as_completed(futures):
                future.result()

    # -- setup / preflight -------------------------------------------------

    def setup(self):
        for d in (
            self.assemblies_dir, self.rcorr_dir, self.reports_dir,
            self.shuck_dir, self.quants_dir, self.diamond_dir, self.assemblies_working,
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

    def which_all_in_env(self, env, binaries):
        """The subset of `binaries` on PATH in `env`, from one `conda run`.

        One call per environment rather than per tool: each `conda run`
        costs a second or more of conda's own startup, and preflight checks
        a dozen and a half tools.
        """
        try:
            result = subprocess.run(
                ["conda", "run", "-n", env, "which", *binaries],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
            )
        except FileNotFoundError:
            return set()
        found = {os.path.basename(line.strip()) for line in result.stdout.splitlines()}
        return found & set(binaries)

    def required_tools(self):
        """(env, binary, label) for every tool this entry point shells out to.

        Order is the order preflight prints them in. An entry point that
        doesn't assemble overrides this rather than demanding assemblers it
        will never run (see chowder.py), and a run handed corrected reads
        drops trimmomatic and rcorrector the same way.
        """
        if self.corrected_reads:
            return tuple(t for t in CHECK_TOOLS if t not in READ_PREP_TOOLS)
        return CHECK_TOOLS

    def check(self):
        """Verify every tool is present, saying nothing when they all are.

        A dozen "installed" lines is noise on every successful run; the only
        news preflight has is a tool that is missing.
        """
        by_env = {}
        for env, binary, label in self.required_tools():
            by_env.setdefault(env, []).append((binary, label))
        missing = []
        for env, tools in by_env.items():
            found = self.which_all_in_env(env, [b for b, _ in tools])
            missing += [f"{label} (in the {env} env)" for b, label in tools if b not in found]
        if missing:
            sys.exit("*** not installed, must fix: " + ", ".join(missing) + " ***")
        self.check_databases()
        self.check_pytransrate_version()
        self.log_provenance()

    def busco_lineage_present(self) -> bool:
        """Whether --lineage names a BUSCO dataset this install has.

        Lenient on purpose: --lineage may be a path, and the directory a
        download leaves can carry a version suffix the flag does not (or the
        reverse), so any directory under busco_dbs/ whose name starts with
        the lineage less a trailing `.N` counts -- the same test the
        Makefile's busco_data target uses to decide it is installed.
        """
        if Path(self.lineage).is_dir():
            return True
        base = re.sub(r"\.\d+$", "", Path(self.lineage).name)
        root = self.makedir / "busco_dbs"
        for depth in ("*", "*/*", "*/*/*"):
            if any(d.is_dir() and d.name.startswith(base) for d in root.glob(depth)):
                return True
        return False

    def check_databases(self):
        """Refuse to start without the databases a run reads.

        Each of these used to be found missing late. The diamond database at
        the first diamond pass, which under oyster.py is in Stage A, hours
        in; the BUSCO lineage at the very end. The swissprot fasta was never
        found missing at all: twotrack_select.py quietly ranked contigs by
        aligned length instead of protein coverage, and so built a different
        assembly, with one line in the log to say so.
        """
        problems = []
        dmnd = self.diamond_db.with_suffix(".dmnd")
        if not dmnd.is_file():
            problems.append(f"{dmnd}\n    the swissprot diamond database; `make diamond_data` builds it")
        if not self.sprot_fasta.is_file():
            problems.append(f"{self.sprot_fasta}\n    the fasta that database was built from; "
                            "the two-track merge reads protein lengths from it. "
                            "`make diamond_data` downloads it")
        if not self.busco_lineage_present():
            problems.append(f"BUSCO lineage {self.lineage!r} under {self.makedir / 'busco_dbs'}\n"
                            f"    `conda run -n orp_busco busco --download {self.lineage} "
                            f"--download_path {self.makedir / 'busco_dbs'}`, or pass "
                            "--lineage a dataset you have")
        if problems:
            sys.exit("*** missing, must fix:\n  " + "\n  ".join(problems) + "\n***")

    def log_provenance(self):
        """Print what a later post-mortem needs and cannot recover.

        `sacct` is the only place a peak memory figure for a finished job
        exists, and it is keyed on a job ID that the log never carried --
        so the run that raised the memory question could not be asked about
        memory afterwards. The scheduler puts the ID in the environment;
        writing it down costs a line and is the difference between
        measuring a failure and arguing about it.
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

    def check_pytransrate_version(self):
        """Refuse to start on a pyTransRate older than the pipeline needs.

        Present-and-runnable is the wrong question for this one tool: the
        version that matters is the difference between a 16-hour failure and
        a run that finishes, and nothing else in the pipeline notices which
        one is installed. See PYTRANSRATE_MIN_VERSION.

        The version is printed either way, because the second half of the
        problem is a fix that was installed and did not work being
        indistinguishable, in the log, from a fix that was never installed.
        The log now carries the answer at the top, before the hours.

        Between the minimum and PYTRANSRATE_PINNED_VERSION it only warns:
        such a version works, and the pin moves for reasons as small as
        logging.

        A version that cannot be read is not fatal. This check exists to
        catch a known-old install, not to become a new way for the run to
        refuse to start.
        """
        version = self.tool_version("orp", "pytransrate")
        if version is None:
            print("[preflight] could not read the pyTransRate version; "
                  f"carrying on (this pipeline needs >= {PYTRANSRATE_MIN_VERSION})")
            return
        print(f"[preflight] pyTransRate {version}")
        if version_below(version, PYTRANSRATE_MIN_VERSION):
            sys.exit(
                f"\n*** pyTransRate {version} is installed and this pipeline "
                f"needs at least {PYTRANSRATE_MIN_VERSION}. ***\n\n"
                "    Older versions size the read-metrics step against the\n"
                "    whole machine rather than the memory budget, and delete\n"
                "    the BAM when a run fails -- so the failure costs a full\n"
                "    remap on every retry. Update the orp environment:\n\n"
                "      conda run -n orp pip install --upgrade --force-reinstall \\\n"
                "        --no-deps \\\n"
                "        'pytransrate @ git+https://github.com/macmanes-lab/"
                f"pytransrate.git@v{PYTRANSRATE_PINNED_VERSION}'\n"
            )
        if version_below(version, PYTRANSRATE_PINNED_VERSION):
            # Between the minimum and the pin: works, but not what a fresh
            # env would hold. Said, not enforced -- see PYTRANSRATE_PINNED_VERSION.
            print(f"[preflight] WARNING: pyTransRate {version} is older than the "
                  f"{PYTRANSRATE_PINNED_VERSION} this release pins. The run will "
                  "carry on and its scores are unaffected; to match a fresh env:\n"
                  "      conda run -n orp pip install --upgrade --force-reinstall "
                  "--no-deps \\\n"
                  "        'pytransrate @ git+https://github.com/macmanes-lab/"
                  f"pytransrate.git@v{PYTRANSRATE_PINNED_VERSION}'")

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
        avg1, max1 = read_length_stats(self.read1)
        avg2, max2 = read_length_stats(self.read2)
        max_read_len = max(max1, max2)

        def is_fractional(spec):
            return bool(spec) and isinstance(spec[0], float)

        frac1, frac2 = is_fractional(self.spades1_kmer), is_fractional(self.spades2_kmer)

        # Resolve percentage specs here, the one place read length is already
        # in hand, so the concrete k values are fixed and printed before any
        # of the assembly work starts rather than discovered mid-run.
        self.spades1_kmer = resolve_kmers(self.spades1_kmer, max_read_len, "spadesauto: ")
        self.spades2_kmer = resolve_kmers(self.spades2_kmer, max_read_len, "spadeshigh: ")
        mean_shown = str(avg1) if avg1 == avg2 else f"{min(avg1, avg2)}-{max(avg1, avg2)}"
        for name, spec in (("spadesauto", self.spades1_kmer), ("spadeshigh", self.spades2_kmer)):
            shown = "auto (rnaSPAdes picks from read length)" if spec is None else ",".join(str(k) for k in spec)
            print(f"[kmer] {name}: k = {shown}   (reads: mean {mean_shown} bp, max {max_read_len} bp)")

        # A k at or above read length is fatal -- rnaSPAdes extracts no k-mers
        # at all (ablab/spades#237). Split by where the k came from: an
        # absolute k that doesn't fit the data is a mistake in the command
        # line and should stop the run, while a resolved percentage was
        # derived from these very reads and is only worth warning about. The
        # warning still matters, because percentages resolve against *max*
        # read length: on variable-length input a k below the max can still
        # sit above the mean, where most reads contribute nothing.
        absolute = [k for spec, frac in ((self.spades1_kmer, frac1), (self.spades2_kmer, frac2))
                    if spec and not frac for k in spec]
        if absolute and not (avg1 > max(absolute) and avg2 > max(absolute)):
            sys.exit(
                f"\n\n\n\n IT LOOKS LIKE YOUR READS ARE NOT AT LEAST {max(absolute)} BP LONG,\n "
                'PLEASE EDIT YOUR COMMAND USING THE "--spades1-kmer"/"--spades2-kmer" FLAGS,\n'
                " SETTING EACH ASSEMBLY KMER LENGTH TO AN ODD NUMBER LESS THAN YOUR READ LENGTH,\n"
                ' OR TO "auto" TO LET rnaSPAdes PICK FROM THE READS \n\n\n\n'
            )

        resolved = [k for spec, frac in ((self.spades1_kmer, frac1), (self.spades2_kmer, frac2))
                    if spec and frac for k in spec]
        if resolved and not (avg1 > max(resolved) and avg2 > max(resolved)):
            print(
                f"\n{RED} WARNING: k={max(resolved)} was derived from a maximum read length of "
                f"{max_read_len} bp, but the mean read length is only {min(avg1, avg2)} bp.\n"
                " Reads shorter than k contribute no k-mers at all, so this assembly may "
                f"recover very little.{RESET}\n"
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
        cmd = self._trinity_base_cmd(cpu, mem)
        # Trinity writes <outdir>.Trinity.fasta either way; --full_cleanup
        # only decides whether <outdir>/ survives it.
        if not self.no_cleanup:
            cmd.append("--full_cleanup")
        # No retries: this step's wall time dwarfs every other (hours to
        # days), so blindly retrying a deterministic failure could multiply
        # the wall time before finally giving up. It also resumes from its
        # own checkpoints in-place (same as Phase 1 above), so a manual
        # re-run of oyster.py after a transient failure loses little anyway.
        self.conda_run("orp_trinity", *cmd, retries=0)
        tmp = out.with_suffix(".fa")
        awk_first_field(out, tmp)
        tmp.replace(out)
        if not self.no_cleanup:
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
        # Cleared first as well as between retries: rnaSPAdes won't start in a
        # non-empty output directory, and the one an interrupted attempt left
        # would otherwise fail every retry of the resume too.
        shutil.rmtree(workdir, ignore_errors=True)
        self.conda_run("orp_spades", *cmd, retry_cleanup=workdir)
        shutil.move(str(workdir / "transcripts.fasta"), str(out))
        if not self.no_cleanup:
            shutil.rmtree(workdir, ignore_errors=True)

    def run_spadesauto(self, cpu=None, mem=None):
        self.run_spades(self.spades1_kmer, "spadesauto", "auto", cpu=cpu, mem=mem)

    def run_spadeshigh(self, cpu=None, mem=None):
        self.run_spades(self.spades2_kmer, "spadeshigh", "high", cpu=cpu, mem=mem)

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
        shutil.rmtree(workdir, ignore_errors=True)  # see run_spades
        self.conda_run("orp_transabyss", *cmd, retry_cleanup=workdir)
        final = workdir / f"{self.runout}.transabyss.fasta-final.fa"
        awk_first_field(final, out)
        if not self.no_cleanup:
            shutil.rmtree(workdir, ignore_errors=True)

    # -- shuck: pool the assemblies, score the pool, pick from it --------------

    def assembly_fasta(self, assembly):
        return self.assemblies_dir / f"{self.runout}.{assembly.fasta_name}"

    def assembly_fasta_paths(self):
        return [self.assembly_fasta(a) for a in self.assemblies]

    def diamond_txt(self, assembly):
        return self.diamond_dir / f"{self.runout}.{assembly.diamond_label}.diamond.txt"

    def unique_txt(self, assembly):
        return self.diamond_dir / f"{self.runout}.unique.{assembly.unique_label}.txt"

    def short_fasta_paths(self):
        return [self.shuck_working / f"{self.assembly_fasta(a).name}.short.fasta"
                for a in self.assemblies]

    def run_filtershort(self):
        """Drop contigs of 200 bp or less from every assembly, all at once.

        One single-threaded Biopython pass per assembly, independent of each
        other, with nothing else running at this point -- so they go in
        parallel rather than one after another, each paying its own
        `conda run` startup in the same few seconds.
        """
        self.shuck_working.mkdir(parents=True, exist_ok=True)

        def filter_one(a):
            fasta = self.assembly_fasta(a)
            outp = self.shuck_working / f"{fasta.name}.short.fasta"
            self.conda_run("orp", "python", self.makedir / "scripts" / "long.seq.py", fasta, outp, "200")

        workers = max(1, min(len(self.assemblies), self.cpu))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            for future in [ex.submit(filter_one, a) for a in self.diamond_priority]:
                future.result()

    def build_pool(self):
        out = self.shuck_dir / "pool.fasta"
        with open(out, "wb") as outf:
            for p in self.short_fasta_paths():
                with open(p, "rb") as inf:
                    shutil.copyfileobj(inf, outf)

    @staticmethod
    def clear_pytransrate_outdir(outdir, assembly):
        """Clear pyTransRate's -o of what a retry must not reuse, and only that.

        A retry has to start from a directory pyTransRate can work in:
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
        sweep steps up. pyTransRate keys its own reuse on the GenomeIndex
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
        pyTransRate reads after 2.2.0 to decide the same question: drop the
        marker and it would move a perfectly good BAM aside and map again.
        The <index>.index.lock file is kept for the same reason -- it sits
        beside the index rather than inside it precisely so an rmtree of a
        partial index cannot pull it out from under its own holder, and it
        belongs to the index either way.

        **The read count** that goes with it, `*-read_count.txt`. It is
        keyed on the read filenames and depends only on the reads, so it
        cannot go stale while those names hold. Keeping it matters more
        than its size suggests: it is what pyTransRate reads when it reuses
        a BAM, and without it that path falls back to counting lines in the
        fastq itself.

        **salmon/**, but only beside a BAM that was kept, and only when
        quant.sf is whole. quant.sf is a completed quantification *of that
        BAM*, so keeping it when the BAM it was computed from has gone
        would score the assembly off numbers belonging to a file that no
        longer exists; and see quant_sf_is_complete for why existence is
        not enough on its own -- no pyTransRate to date checks it.

        **logs/**, which holds snap.log, the file pyTransRate points at
        when snap dies without explaining itself, so deleting it is
        deleting the evidence the retry exists to gather. pyTransRate
        rewrites it per attempt, so what survives the last retry is the
        last attempt's output, which is the one worth reading.

        Everything else goes: assemblies.csv, contigs.csv and the score
        optimisation csv are the outputs being recomputed, and anything a
        future pyTransRate leaves behind that this does not recognise is
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
        print("[pytransrate] {}: kept from the last attempt: {}".format(
            outdir.name,
            ", ".join(p.name for p in sorted(survived)) or "nothing",
        ))

    def pytransrate_memory_args(self, mem):
        """``--max-memory`` for a pyTransRate call, or nothing.

        pyTransRate sizes the read-metrics step's shared accumulators by the
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
        the only one that never heard it. Requires pyTransRate >= 2.2.0,
        which is what PYTRANSRATE_MIN_VERSION enforces at preflight.
        """
        if mem is None:
            return []
        if any(arg.split("=")[0] in PYTRANSRATE_MEMORY_FLAGS
               for arg in self.pytransrate_args):
            return []
        return ["--max-memory", f"{mem}G"]

    def score_pool(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        # The two-track path calls pool_branch() bare, so mem arrives as
        # None there -- which pytransrate_memory_args reads as "no budget".
        mem = self.mem if mem is None else mem
        outdir = self.shuck_dir / "pool"
        pool = self.shuck_dir / "pool.fasta"
        # needs_run() re-runs this step whenever the corrected reads are
        # newer than pool/assemblies.csv -- not only when it is absent --
        # so a resumed run would abort on that csv unless it is cleared
        # first. retry_cleanup repeats the clear before each retry, because
        # the one below happens once, outside run()'s retry loop. See
        # clear_pytransrate_outdir for what survives it and why.
        self.clear_pytransrate_outdir(outdir, pool)
        self.conda_run(
            "orp", "pytransrate",
            "-o", outdir, "-t", cpu, "-a", pool,
            "--left", self.cor1(), "--right", self.cor2(),
            *self.pytransrate_memory_args(mem),
            *self.pytransrate_args,
            retry_cleanup=partial(self.clear_pytransrate_outdir, outdir, pool),
        )

    def twotrack_select(self):
        """Choose the contigs to carry forward by swissprot gene, not orthogroup.

        See scripts/twotrack_select.py: one representative per gene (the
        longest of those with near-best protein coverage), distinct expressed
        copies of kept genes, and cd-hit-est over the contigs without a hit.
        Writes good.<run>.list, so shuck onwards is unchanged, and
        twotrack.<run>.tsv beside it saying what happened to every contig.
        """
        print("Selecting contigs by swissprot gene (two-track)")
        good_list = self.shuck_dir / f"good.{self.runout}.list"
        self.conda_run(
            "orp", "python", self.makedir / "scripts" / "twotrack_select.py",
            "--pool", self.shuck_dir / "pool.fasta",
            "--contigs-csv", self.shuck_dir / "pool" / "contigs.csv",
            "--diamond", *[self.diamond_txt(a) for a in self.diamond_priority],
            "--sprot", self.sprot_fasta,
            "--table", self.shuck_dir / f"twotrack.{self.runout}.tsv",
            "--threads", self.cpu, "--mem-mb", self.mem * 1000, "--out", good_list,
        )

    def shuck(self):
        good_list = self.shuck_dir / f"good.{self.runout}.list"
        out = self.assemblies_dir / f"{self.runout}.shucked.fasta"
        with open(out, "w") as outf:
            subprocess.run(
                ["conda", "run", "--no-capture-output", "-n", "orp", "python",
                 str(self.makedir / "scripts" / "filter.py"), str(self.shuck_dir / "pool.fasta"), str(good_list)],
                check=True, stdout=outf, cwd=self.dir,
            )

    # -- diamond passes ----------------------------------------------------

    def diamond_jobs(self):
        return [
            (self.assemblies_dir / f"{self.runout}.shucked.fasta",
             self.diamond_dir / f"{self.runout}.shucked.diamond.txt"),
        ] + [(self.assembly_fasta(a), self.diamond_txt(a)) for a in self.diamond_priority]

    def run_diamond_one(self, query, out, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        self.conda_run(
            "orp", "diamond", "blastx", "--quiet", "-p", cpu,
            "-e", "1e-8", "--top", "0.1", "-q", query, "-d", self.diamond_db, "-o", out,
        )

    def hits_for(self, fasta, out):
        """Write to `out` the swissprot hits, already in hand, of every contig in `fasta`.

        shucked.fasta and ORP.intermediate.fasta contain nothing but contigs
        lifted whole, under their own names, out of the assemblies: shuck
        filters the pool, posthack filters the assemblies, and cd-hit-est
        only chooses among them. Each assembly has already been through
        diamond blastx against the same database with the same settings,
        and diamond scores each query on its own, so a contig's hits are the
        same whichever file it is searched from. These two used to be blastx
        runs of their own; now they are a lookup.

        Grouped by assembly rather than in `fasta`'s order. Nothing reads
        them in order: make_list1 and orp_uniq take the set of genes, and
        secondfilter the set of contigs with a hit.
        """
        with open(fasta) as f:
            ids = {(line[1:].split() or [""])[0] for line in f if line.startswith(">")}
        part = out.with_name(out.name + ".part")
        with open(part, "w") as o:
            for a in self.diamond_priority:
                with open(self.diamond_txt(a)) as f:
                    for line in f:
                        if line.split("\t", 1)[0] in ids:
                            o.write(line)
        part.replace(out)

    def diamond_uniq(self):
        for a in self.report_order:
            count = parse_unique_count(self.diamond_txt(a))
            self.unique_txt(a).write_text(f"{count}\n")

    def make_list1(self):
        ids = extract_gene_ids(self.diamond_dir / f"{self.runout}.shucked.diamond.txt")
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
        fasta = self.assemblies_dir / f"{self.runout}.shucked.fasta"
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
        fastas = " ".join(shlex.quote(str(p)) for p in self.assembly_fasta_paths())
        list7 = self.diamond_dir / f"{self.runout}.list7"
        newbies = self.diamond_dir / f"{self.runout}.newbies.fasta"
        shucked = self.assemblies_dir / f"{self.runout}.shucked.fasta"
        working_out = self.assemblies_working / f"{self.runout}.shucked.fasta"
        filter_py = self.makedir / "scripts" / "filter.py"
        # `>`, not `>>`: a retry or a resumed run would otherwise add a second
        # copy of every rescued contig on top of the first.
        q = shlex.quote
        script = f"python {q(str(filter_py))} <(cat {fastas}) {q(str(list7))} > {q(str(newbies))}"
        self.run(["conda", "run", "--no-capture-output", "-n", "orp", "bash", "-c", script])
        with open(working_out, "wb") as outf:
            for p in (newbies, shucked):
                with open(p, "rb") as inf:
                    shutil.copyfileobj(inf, outf)

    # -- dedup / quantify -----------------------------------------------------

    def cdhit(self):
        src = self.assemblies_working / f"{self.runout}.shucked.fasta"
        out = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        self.conda_run(
            "orp", "cd-hit-est", "-M", self.mem * 1000, "-T", self.cpu,
            "-c", ".98", "-i", src, "-o", out,
        )

    def orp_diamond(self):
        self.hits_for(self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta",
                      self.assemblies_dir / f"{self.runout}.ORP.diamond.txt")

    def orp_uniq(self):
        diamond_txt = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        count = parse_unique_count(diamond_txt)
        (self.assemblies_working / f"{self.runout}.unique.ORP.txt").write_text(f"{count}\n")
        (self.assemblies_working / f"{self.runout}.unique.ORP.done").touch()

    def salmon_index(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        src = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        idx = self.quants_dir / f"{self.runout}.shucked.idx"
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
        idx = self.quants_dir / f"{self.runout}.shucked.idx"
        outdir = self.quants_dir / f"salmon_shucked_{self.runout}"
        self.conda_run(
            "orp", "salmon", "quant", "--no-version-check",
            "-p", cpu, "-i", idx, "--seqBias", "--gcBias", "--libType", "A",
            "-1", self.cor1(), "-2", self.cor2(), "-o", outdir,
        )

    def filter_tpm(self):
        quant = self.quants_dir / f"salmon_shucked_{self.runout}" / "quant.sf"
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
                # At-threshold counts as high: with `>` and `<` a contig whose
                # TPM equalled --tpm-filt landed in neither list and so was
                # dropped from .ORP.fasta whenever LOWEXP was non-empty.
                if tpm >= self.tpm_filt:
                    hf.write(cols[0] + "\n")
                else:
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
            # "w", not "a": IDs left from an earlier attempt would be kept
            # again, and one that is now in HIGHEXP would be written twice.
            with open(diamond_txt) as f, open(blasted, "w") as bf:
                for line in f:
                    cols = line.rstrip("\n").split("\t")
                    if cols and cols[0] in low_ids:
                        bf.write(line)
                        do_not_remove_ids.add(cols[0])
            with open(donotremove, "w") as df:
                for i in sorted(do_not_remove_ids):
                    df.write(i + "\n")
            print(f"[secondfilter] {len(do_not_remove_ids)} low-expression contigs kept "
                  "for their swissprot hit")

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
        name = f"run_{self.runout}.ORP"
        work, final = self.dir / name, self.reports_dir / name
        # BUSCO refuses to start over an existing -o, which an interrupted
        # attempt leaves behind. And shutil.move into an existing directory
        # moves *inside* it: a re-run used to land at
        # reports/<name>/<name>, where reportgen could read either summary,
        # and the run after that failed outright. Clear both ends.
        shutil.rmtree(work, ignore_errors=True)
        self.conda_run(
            "orp_busco", "busco", "--offline", "--lineage", self.lineage,
            "--download_path", self.makedir / "busco_dbs",
            "-i", orp_fasta, "-m", "transcriptome", "--cpu", cpu,
            "-o", name, "--config", self.busco_config,
            retry_cleanup=work,
        )
        shutil.rmtree(final, ignore_errors=True)
        shutil.move(str(work), str(final))
        (self.reports_dir / f"{self.runout}.busco.done").touch()

    def adopt_pre_rename_reports(self):
        """Carry reports/transrate_<run>/ across its rename to pytransrate_<run>/.

        Before 4.1.0-dev9 the final scoring step was called `transrate` and
        wrote there. Renaming the directory in place, rather than letting the
        step find its new output missing, keeps a resumed run from scoring
        the assembly all over again.
        """
        old = self.reports_dir / f"transrate_{self.runout}"
        new = self.reports_dir / f"pytransrate_{self.runout}"
        if old.is_dir() and not new.exists():
            old.rename(new)
            print(f"[resume] renamed {self._rel(old)} to {self._rel(new)}")

    def pytransrate(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        mem = self.mem if mem is None else mem
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        outdir = self.reports_dir / f"pytransrate_{self.runout}"
        # See score_pool() and clear_pytransrate_outdir.
        self.clear_pytransrate_outdir(outdir, orp_fasta)
        self.conda_run(
            "orp", "pytransrate",
            "-o", outdir, "-a", orp_fasta,
            "--left", self.cor1(), "--right", self.cor2(), "-t", cpu,
            *self.pytransrate_memory_args(mem),
            *self.pytransrate_args,
            retry_cleanup=partial(self.clear_pytransrate_outdir, outdir, orp_fasta),
        )

    def trinity_perllib_dir(self):
        result = subprocess.run(
            ["conda", "run", "-n", "orp_trinity", "which", "Trinity"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True,
        )
        trinity_path = Path(result.stdout.strip()).resolve()
        return trinity_path.parent / "PerlLib"

    def strandeval_scratch(self):
        """strandeval's working files: the sampled BAM, its bwa index, and
        the column hist was fed. Removed as soon as strandeval finishes,
        or by cleanup() after a --no-cleanup run."""
        paths = [self.dir / f"{self.runout}.hist_input.txt",
                 self.dir / f"{self.runout}.sorted.bam"]
        paths += [self.dir / f"{self.runout}.{ext}"
                  for ext in ("bwt", "pac", "ann", "amb", "sa", "dat")]
        return paths

    def strandeval(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        r1, r2 = self.cor1(), self.cor2()
        self.conda_run("orp_trinity", "bwa", "index", "-p", self.runout, orp_fasta)

        pipeline_script = (
            f'conda run --no-capture-output -n orp_trinity bash -c '
            f'"bwa mem -t {cpu} {self.runout} '
            f'<(seqtk sample -s 23894 {r1} 400000) <(seqtk sample -s 23894 {r2} 400000)" '
            f"| conda run --no-capture-output -n orp samtools view -@{cpu} -Sb - "
            f"| conda run --no-capture-output -n orp samtools sort -T {self.runout} -O bam -@{cpu} "
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

        if not self.no_cleanup:
            for p in self.strandeval_scratch():
                if p.exists():
                    p.unlink()

        (self.reports_dir / f"{self.runout}.strandeval.done").touch()
        histogram_text = hist_result.stdout.rstrip("\n")
        summary = (
            "\n*****  STRAND EXAMINATION HISTOGRAM ***** \n"
            f"{histogram_text}\n"
            "\n*****  See the following link for interpretation ***** \n"
            "*****  https://github.com/macmanes-lab/Oyster_River_Protocol/blob/master/docs/strandexamine.md ***** \n"
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

        csv_path = next((self.reports_dir / f"pytransrate_{runout}").rglob("assemblies.csv"), None)
        rows = list(csv.reader(open(csv_path))) if csv_path else []
        pytransrate_score = rows[1][36] if len(rows) > 1 and len(rows[1]) > 36 else ""
        pytransrate_optimal = rows[1][37] if len(rows) > 1 and len(rows[1]) > 37 else ""
        emit("*****  PYTRANSRATE SCORE ~~~~~~~~~~~~~~~~>     ", pytransrate_score)
        emit("*****  PYTRANSRATE OPTIMAL SCORE ~~~~~~~~>     ", pytransrate_optimal)

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
        self.timing_log.write_text(text)
        print(f"\nStep timing saved to: {self.timing_log}\n")

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
        if self.corrected_reads:
            return self.use_corrected_reads()
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

    def use_corrected_reads(self):
        """Stand the user's already-corrected pair in for rcorrector's output.

        A plain fastq is symlinked into place. A gzipped one is decompressed
        instead, because the corrected pair is named .cor.fq and the
        assemblers believe the name: Trinity and SPAdes decide whether to
        gunzip from the extension, so a symlink to gzip data under that name
        would be read as garbage. Neither is queued for compression -- the
        user already has the reads, and cleanup() keeps the symlinks and
        deletes the copies.
        """
        self.rcorr_dir.mkdir(parents=True, exist_ok=True)
        for src, dst in ((self.read1, self.cor1()), (self.read2, self.cor2())):
            if dst.exists() and not self.needs_run([dst], [src]):
                continue
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            if is_gzip(src):
                self.decompress(src, dst)
            else:
                dst.symlink_to(src.resolve())
        print(f"[reads] corrected reads given: using {self.read1} / {self.read2} "
              "as they are; trimmomatic and rcorrector skipped")

    def decompress(self, src, dst):
        """gunzip `src` to `dst`, via a .part so a killed run leaves no stub."""
        part = dst.with_name(dst.name + ".part")
        start = time.time()
        try:
            with open(part, "wb") as out:
                self.run(self._resolve_compressor() + ["-d", str(src.resolve())], stdout=out,
                         retries=0)
            part.replace(dst)
        except BaseException:
            if part.exists():
                part.unlink()
            raise
        print(f"[reads] {self._rel(dst)} decompressed in {int(time.time() - start)}s "
              f"({human_size(path_size(src))} -> {human_size(path_size(dst))})")

    def run_assemblers(self):
        """Build the four assemblies this pipeline is named for.

        oyster.py's own half of the work, and the only half that is specific
        to a particular set of assemblers -- everything downstream of
        run_filtershort treats them as an unordered set of inputs.
        """
        c1, c2 = self.cor1(), self.cor2()
        trinity_fa = self.assembly_fasta(TRINITY)
        phase1_done = self.trinity_phase1_done()
        sphigh = self.assembly_fasta(SPADES_HIGH)
        spauto = self.assembly_fasta(SPADES_AUTO)
        ta = self.assembly_fasta(TRANSABYSS)
        diamond_ta = self.diamond_txt(TRANSABYSS)
        diamond_sphigh = self.diamond_txt(SPADES_HIGH)
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
                    "run_trinity_phase1", [phase1_done], [c1, c2],
                    partial(self.run_trinity_phase1, cpu=phase1_cpu, mem=phase1_mem),
                )
            except Exception as e:
                _lane_failed("run_trinity_phase1", e)
                raise

        def spades_lane():
            # spadesauto (the lower-k assembly, whose smaller k dominates its
            # cost) first, since it has been the slower of the two --
            # historically the fixed k=55 run against the fixed k=75 one.
            # diamond_{spadesauto, spadeshigh} depend only on their own assembly
            # (not on Trinity or the shuck stage below), so each fires as
            # soon as its assembly is done instead of waiting for the merge.
            try:
                for step_name, outputs, inputs, assemble, diamond_name, diamond_out in (
                    ("run_spadesauto", [spauto], [c1, c2], self.run_spadesauto, "spadesauto", diamond_spauto),
                    ("run_spadeshigh", [sphigh], [c1, c2], self.run_spadeshigh, "spadeshigh", diamond_sphigh),
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

        print(f"\n=== Stage A: run_trinity_phase1 ({phase1_cpu} cpu) || spadesauto/spadeshigh ({spades_cpu} cpu) -- start {self._ts()} ===")
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
        pool_fasta = self.shuck_dir / "pool.fasta"
        pool_csv = self.shuck_dir / "pool" / "assemblies.csv"
        good_list = self.shuck_dir / f"good.{self.runout}.list"
        shucked_fasta = self.assemblies_dir / f"{self.runout}.shucked.fasta"
        diamond_outs = [o for _, o in self.diamond_jobs()]
        diamond_shucked = self.diamond_dir / f"{self.runout}.shucked.diamond.txt"
        uniq_outs = [self.unique_txt(a) for a in self.report_order]
        list1 = self.diamond_dir / f"{self.runout}.list1"
        list2 = self.diamond_dir / f"{self.runout}.list2"
        list3 = self.diamond_dir / f"{self.runout}.list3"
        list5 = self.diamond_dir / f"{self.runout}.list5"
        list6 = self.diamond_dir / f"{self.runout}.list6"
        list7 = self.diamond_dir / f"{self.runout}.list7"
        newbies = self.diamond_dir / f"{self.runout}.newbies.fasta"
        working_shucked = self.assemblies_working / f"{self.runout}.shucked.fasta"
        orp_intermediate = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        orp_diamond_txt = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        unique_orp_done = self.assemblies_working / f"{self.runout}.unique.ORP.done"
        shucked_idx = self.quants_dir / f"{self.runout}.shucked.idx"
        quant_sf = self.quants_dir / f"salmon_shucked_{self.runout}" / "quant.sf"
        filter_done = self.assemblies_dir / f"{self.runout}.filter.done"
        low_txt = self.assemblies_working / f"{self.runout}.LOWEXP.txt"
        high_txt = self.assemblies_working / f"{self.runout}.HIGHEXP.txt"
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        busco_done = self.reports_dir / f"{self.runout}.busco.done"
        pytransrate_csv = self.reports_dir / f"pytransrate_{self.runout}" / "assemblies.csv"
        strandeval_done = self.reports_dir / f"{self.runout}.strandeval.done"
        qualreport_done = self.reports_dir / f"qualreport.{self.runout}.done"
        cleanup_done = self.reports_dir / f"{self.runout}.cleanup.done"

        self.step("run_filtershort", short_fastas, assembly_fastas, self.run_filtershort)

        def pool_branch(cpu=None, mem=None):
            self.step("build_pool", [pool_fasta], short_fastas, self.build_pool)
            self.step("score_pool", [pool_csv], [pool_fasta, c1, c2],
                      partial(self.score_pool, cpu=cpu, mem=mem))

        def assembly_diamonds(cpu=None, mem=None):
            for a in self.diamond_priority:
                fasta, out = self.assembly_fasta(a), self.diamond_txt(a)
                self.step(
                    f"diamond_{a.diamond_label}", [out], [fasta],
                    partial(self.run_diamond_one, fasta, out, cpu=cpu),
                )

        # The pool is scored while any swissprot pass not already done runs
        # beside it -- under oyster.py only Trinity's, whose lane has just
        # finished; under chowder.py every assembly's. They are independent
        # until twotrack_select, which groups the pooled contigs by those
        # hits. The pool-side gate names everything score_pool reads,
        # corrected reads included, since it is checked before the inner
        # steps get a look.
        self.run_beside(
            ("pool_branch", [pool_csv], short_fastas + [c1, c2], pool_branch),
            ("assembly_diamonds", diamond_outs[1:], assembly_fastas, assembly_diamonds),
            *self.side_job_budget(max_cpu=self.cpu),
        )
        self.step(
            "twotrack_select", [good_list],
            [pool_fasta, pool_csv] + diamond_outs[1:], self.twotrack_select,
        )
        self.step("shuck", [shucked_fasta], [good_list, pool_fasta], self.shuck)
        self.after_pick(c1, c2, diamond_outs, diamond_shucked, shucked_fasta, uniq_outs,
                        list1, list2, list3, list5, list6, list7, newbies, working_shucked,
                        orp_intermediate, orp_diamond_txt, unique_orp_done, shucked_idx, quant_sf,
                        filter_done, low_txt, high_txt, orp_fasta, busco_done, pytransrate_csv,
                        strandeval_done, qualreport_done, cleanup_done, pipeline_start)

    def after_pick(self, c1, c2, diamond_outs, diamond_shucked, shucked_fasta, uniq_outs,
                   list1, list2, list3, list5, list6, list7, newbies, working_shucked,
                   orp_intermediate, orp_diamond_txt, unique_orp_done, shucked_idx, quant_sf,
                   filter_done, low_txt, high_txt, orp_fasta, busco_done, pytransrate_csv,
                   strandeval_done, qualreport_done, cleanup_done, pipeline_start):
        """Everything after good_list exists."""

        # Every assembly's diamond pass is done by now (twotrack_select needed
        # them). shucked's hits are looked up from them rather than searched
        # for again -- see hits_for.
        assembly_diamonds = diamond_outs[1:]
        self.step(
            "diamond_shucked", [diamond_shucked], [shucked_fasta] + assembly_diamonds,
            partial(self.hits_for, shucked_fasta, diamond_shucked),
        )
        self.step("diamond_uniq", uniq_outs, diamond_outs, self.diamond_uniq)
        self.step("make_list1", [list1], [diamond_shucked], self.make_list1)
        self.step("make_list2", [list2], [self.diamond_txt(a) for a in self.assemblies], self.make_list2)
        self.step("make_list3", [list3], [list1, list2], self.make_list3)
        self.step("make_list5", [list5], [list3] + [self.diamond_txt(a) for a in self.diamond_priority], self.make_list5)
        self.step("make_list6", [list6], [shucked_fasta], self.make_list6)
        self.step("make_list7", [list7], [list6, list5], self.make_list7)
        self.step("posthack", [newbies, working_shucked], [list7], self.posthack)
        self.step("cdhit", [orp_intermediate], [working_shucked], self.cdhit)

        # A lookup like diamond_shucked's, not a search.
        self.step("orp_diamond", [orp_diamond_txt], [orp_intermediate] + assembly_diamonds,
                  self.orp_diamond)
        self.step("orp_uniq", [unique_orp_done], [orp_diamond_txt], self.orp_uniq)
        salmon_stamp = self.stamp_tool_version("orp", "salmon", self.quants_dir / "salmon.version")
        self.step("salmon_index", [shucked_idx], [orp_intermediate, salmon_stamp], self.salmon_index)
        self.step("salmon", [quant_sf], [shucked_idx, c1, c2], self.salmon)
        self.step("filter", [filter_done], [orp_intermediate, quant_sf, orp_diamond_txt], self.filter_tpm)
        self.step(
            "secondfilter", [orp_fasta],
            [filter_done, low_txt, high_txt, orp_intermediate, quant_sf, orp_diamond_txt],
            self.secondfilter,
        )
        # BUSCO gets the whole of --cpu/--busco-threads to itself. strandeval
        # (a few minutes) runs beside the final pyTransRate scoring (much
        # longer) on a few threads of its own, rather than taking half its
        # cores for the whole of its run.
        self.step("busco", [busco_done], [orp_fasta], self.busco)
        self.adopt_pre_rename_reports()
        self.run_beside(
            ("pytransrate", [pytransrate_csv], [orp_fasta, c1, c2], self.pytransrate),
            ("strandeval", [strandeval_done], [orp_fasta, c1, c2], self.strandeval),
            *self.side_job_budget(),
        )
        # Everything the report quotes, so a re-scored pyTransRate or BUSCO
        # (say, against newer reads) rewrites the report rather than leaving
        # the old numbers in it.
        self.step(
            "reportgen", [qualreport_done],
            [unique_orp_done, orp_fasta, busco_done, pytransrate_csv, strandeval_done] + uniq_outs,
            self.reportgen,
        )
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
    p.add_argument("--trimmed-corrected-reads", dest="corrected_reads", action="store_true",
                   help="the reads have already been through trimmomatic and "
                        "rcorrector; skip both and assemble them as they are "
                        "(default: off)")
    p.add_argument("--spades1-kmer", type=parse_kmer_spec, default="auto",
                   help="rnaSPAdes k-mer(s) for the spadesauto assembly: 'auto' to let "
                        "rnaSPAdes pick its documented default pair from read length, or a "
                        "comma-separated list of odd sizes under 128 (default: auto)")
    p.add_argument("--spades2-kmer", type=parse_kmer_spec, default="60%,75%",
                   help="rnaSPAdes k-mer(s) for the spadeshigh assembly: percentages of "
                        "max read length ('60%%,75%%'), an explicit comma-separated list of "
                        "odd sizes under 128, or 'auto' (default: 60%%,75%%)")
    p.add_argument("--transabyss-kmer", type=int, default=32, help="Trans-ABySS k-mer (default: 32)")
    p.add_argument(
        "--max-parallel", type=int, default=2,
        help="2 or more (the default) runs a short independent job beside a long "
             "one on a few threads of its own: the remaining diamond passes beside "
             "score_pool, and strandeval beside pyTransRate. 1 runs them one after "
             "the other. The assemblers' two stage-pairings (see "
             "TRINITY_PHASE1_SHARE/TRINITY_PHASE2_SHARE) are unaffected (default: 2)",
    )
    p.add_argument(
        "--no-cleanup", "--keep-intermediates", dest="no_cleanup", action="store_true",
        help="keep every file a run produces, for debugging: skips the end-of-run "
             "cleanup (shuck/, quants/, diamond/, the working assemblies), "
             "the reclaim of the trimmed reads, Trinity's --full_cleanup, the "
             "removal of the rnaSPAdes and Trans-ABySS working directories and "
             "of strandeval's BAM and bwa index, and leaves the four assemblies "
             "and the corrected reads uncompressed. Re-running without it cleans "
             "up afterwards. --keep-intermediates is an older name for the same "
             "flag (default: off)",
    )
    p.add_argument(
        "--pytransrate-args", default="",
        help="extra arguments passed verbatim to both pyTransRate runs, as one "
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
