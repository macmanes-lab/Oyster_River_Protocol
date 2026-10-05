#!/usr/bin/env python3
"""Release check: oyster.py, chowder.py and strandeval.py on the sample data.

usage: python tests/release_check.py [--tier quick|preflight|full|all]
                                      [--workdir DIR] [--cpu N] [--mem GB]
                                      [--jobs N] [--only CASE ...] [--list]

The test of function to run before cutting a release or pushing a feature.
Every case runs the real entry points on sampledata/test.{1,2}.fq.gz (27k
pairs of 100 bp reads) and checks what they leave behind, not just the exit
status. Three tiers, each a superset of the cost of the one before:

  quick      no conda needed, seconds. Argument parsing and its rejections,
             the k-mer spec parser, chowder's label and order rules, contig
             renaming at ingest, and tests/test_twotrack_select.py.
  preflight  needs the conda envs, a few seconds each. The refusals that
             happen after the tool check and before any real work: missing
             reads, a k at or above read length, a missing BUSCO lineage,
             chowder given one, an empty, a missing or a duplicate-named
             assembly.
  full       needs the conda envs and the databases; minutes per run, most of
             an hour in all with --jobs 1. Complete runs through every
             option, then the resume, rerun and cleanup paths:

               oyster_default      defaults; a rerun is a no-op
               oyster_options      --strand RF --normalize-reads --tpm-filt 1
                                   explicit k-mers, --transabyss-kmer,
                                   --max-parallel 1, --busco-threads,
                                   --pytransrate-args
               oyster_resume       --keep-intermediates keeps everything; a
                                   step left half-done is re-run; rerunning
                                   without the flag cleans up
               oyster_corrected    --trimmed-corrected-reads --strand FR,
                                   multi-k / auto k-mers      (oyster_default)
               chowder_default     its four assemblies, gzipped; rerun no-op
                                                               (oyster_default)
               chowder_options     --labels --assembly-order given
                                   --corrected-reads --tpm-filt 1
                                   --max-parallel 1 --no-cleanup
                                                               (oyster_default)
               chowder_collide     two inputs with one filename, --seed
                                                               (oyster_default)
               strandeval_standalone  scripts/strandeval.py, with and without
                                   --no-cleanup                (oyster_default)

Each case works in <workdir>/<case>/ (cleared when the case starts) and logs
every command and its output to <workdir>/logs/<case>.log. A case is skipped
when its prerequisite fails; when the prerequisite is not selected, its
finished run from an earlier pass in the same workdir is used instead -- so
`--only chowder_default` after a full pass re-tests chowder without
reassembling. The summary goes to stdout
and <workdir>/summary.txt; the exit status is non-zero if anything failed.

Run it from a shell where `conda` is on PATH (the orp envs are reached with
`conda run`, as the pipeline does); the harness itself is stdlib Python and
runs the entry points with the same interpreter. On the cluster,
tests/release_check.sbatch wraps it.
"""

import argparse
import concurrent.futures
import contextlib
import gzip
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO))

OYSTER = REPO / "oyster.py"
CHOWDER = REPO / "chowder.py"
STRANDEVAL = REPO / "scripts" / "strandeval.py"
READ1 = REPO / "sampledata" / "test.1.fq.gz"
READ2 = REPO / "sampledata" / "test.2.fq.gz"
VERSION = (REPO / "version.txt").read_text().strip()

OYSTER_ASSEMBLIES = ("spadesauto.fasta", "spadeshigh.fasta", "transabyss.fasta",
                     "trinity.Trinity.fasta")


class CheckFailed(Exception):
    pass


class Skip(Exception):
    pass


def require(cond, message):
    if not cond:
        raise CheckFailed(message)


# -- the harness ---------------------------------------------------------------

class Context:
    def __init__(self, args):
        self.workdir = Path(args.workdir).resolve()
        self.cpu = args.cpu
        self.mem = args.mem
        self.logs = self.workdir / "logs"
        self.logs.mkdir(parents=True, exist_ok=True)
        self.conda = shutil.which("conda") is not None

    def case_dir(self, name, fresh=False):
        d = self.workdir / name
        if fresh and d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True, exist_ok=True)
        return d


class Case:
    def __init__(self, name, tier, func, needs=()):
        self.name, self.tier, self.func, self.needs = name, tier, func, tuple(needs)
        self.doc = (func.__doc__ or "").strip().splitlines()[0] if func.__doc__ else ""


CASES = []


def case(tier, needs=()):
    def register(func):
        CASES.append(Case(func.__name__.replace("case_", ""), tier, func, needs))
        return func
    return register


class Runner:
    """One case's commands, all logged to <workdir>/logs/<case>.log."""

    def __init__(self, ctx, name):
        self.ctx, self.name = ctx, name
        self.log = ctx.logs / f"{name}.log"
        self.log.write_text("")

    def __call__(self, cmd, expect_rc=0, cwd=None, timeout=None):
        """Run `cmd`; return its combined stdout+stderr. expect_rc=None takes
        any status, "nonzero" any failure."""
        cmd = [str(c) for c in cmd]
        start = time.time()
        with open(self.log, "a") as log:
            log.write(f"\n$ {' '.join(cmd)}\n")
            log.flush()
            proc = subprocess.run(cmd, cwd=str(cwd or self.ctx.workdir),
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  universal_newlines=True, timeout=timeout)
            log.write(proc.stdout)
            log.write(f"\n[exit {proc.returncode} after {int(time.time() - start)}s]\n")
        if expect_rc == "nonzero":
            ok = proc.returncode != 0
        else:
            ok = expect_rc is None or proc.returncode == expect_rc
        if not ok:
            tail = "\n".join(proc.stdout.rstrip().splitlines()[-15:])
            raise CheckFailed(f"exit {proc.returncode} (wanted {expect_rc}) from "
                              f"{Path(cmd[1]).name if len(cmd) > 1 else cmd[0]}; "
                              f"see {self.log}\n{tail}")
        return proc.stdout

    def oyster(self, d, *args, **kw):
        return self(self.pipeline_cmd(OYSTER, d, args), **kw)

    def chowder(self, d, *args, **kw):
        return self(self.pipeline_cmd(CHOWDER, d, args), **kw)

    def pipeline_cmd(self, script, d, args):
        base = [sys.executable, script, "--dir", d]
        if "--cpu" not in args:
            base += ["--cpu", self.ctx.cpu]
        if "--mem" not in args:
            base += ["--mem", self.ctx.mem]
        return base + list(args)


# -- output checks ---------------------------------------------------------------

def fasta_names(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def commands(log, tool):
    """The `+ ...` command lines oyster.py printed that run `tool`."""
    return [l for l in log.splitlines()
            if l.startswith("+ ") and re.search(rf"(^|[\s/]){re.escape(tool)}(\s|$)", l)]


def qualreport(d, runout):
    path = d / "reports" / f"qualreport.{runout}"
    require(path.is_file(), f"no quality report at {path}")
    text = path.read_text()
    fields = {}
    for line in text.splitlines():
        m = re.match(r"\*+\s+(.+?)\s*~*>\s*(.*)$", line)
        if m:
            fields[m.group(1).strip()] = m.group(2).strip()
    return text, fields


def check_finished(d, runout, kind, cleaned=True, corrected=False, labels=None):
    """What every finished run must have left behind.

    kind is "oyster" or "chowder"; labels, for chowder, the labels its
    contigs must be prefixed with.
    """
    asm = d / "assemblies"
    reports = d / "reports"
    orp = asm / f"{runout}.ORP.fasta"
    require(orp.is_file(), f"no final assembly at {orp}")
    names = fasta_names(orp)
    require(names, f"{orp} has no sequences")
    require(len(names) == len(set(names)), f"{orp} has duplicate contig names")
    if labels:
        stray = [n for n in names if not any(n.startswith(l + "_") for l in labels)]
        require(not stray, f"{len(stray)} contigs in {orp.name} carry no input label "
                           f"(e.g. {stray[0] if stray else ''})")

    text, fields = qualreport(d, runout)
    want = "ORP chowder (merge only" if kind == "chowder" else "using the ORP version"
    require(want in text.splitlines()[0], f"qualreport header lacks {want!r}")
    require(VERSION in text.splitlines()[0], f"qualreport header lacks version {VERSION}")
    busco = fields.get("BUSCO SCORE", "")
    require(re.match(r"C:\d", busco), f"no BUSCO score in the qualreport (got {busco!r})")
    try:
        score = float(fields.get("PYTRANSRATE SCORE", ""))
    except ValueError:
        raise CheckFailed(f"no pyTransRate score in the qualreport "
                          f"(got {fields.get('PYTRANSRATE SCORE')!r})")
    require(0 < score <= 1, f"pyTransRate score {score} out of (0, 1]")
    require(re.match(r"^\d+$", fields.get("UNIQUE GENES ORP", "")),
            f"no ORP unique-gene count (got {fields.get('UNIQUE GENES ORP')!r})")
    require(re.match(r"^[\d.]+%$", fields.get("READS MAPPED AS PROPER PAIRS", "")),
            "no proper-pair rate in the qualreport")
    require("STRAND EXAMINATION HISTOGRAM" in text, "no strand histogram in the qualreport")

    require((reports / f"run_{runout}.ORP").is_dir(), "BUSCO directory missing from reports/")
    require(not (reports / f"run_{runout}.ORP" / f"run_{runout}.ORP").exists(),
            "BUSCO output nested inside itself")
    require(list((reports / f"pytransrate_{runout}").rglob("assemblies.csv")),
            "pyTransRate assemblies.csv missing")
    timing = (reports / f"{runout}.timing.log").read_text()
    require(re.search(r"^TOTAL\s+\d\d:\d\d:\d\d", timing, re.M), "timing log has no TOTAL")
    running = list(reports.glob(f".{runout}.*.running"))
    require(not running, f"unfinished-step markers left behind: {[p.name for p in running]}")

    if kind == "oyster":
        inputs = [asm / f"{runout}.{a}" for a in OYSTER_ASSEMBLIES]
    else:
        inputs = sorted((asm / "ingested").glob(f"{runout}.*.fasta")) + \
            sorted((asm / "ingested").glob(f"{runout}.*.fasta.gz"))
        require(inputs, "no ingested assemblies under assemblies/ingested/")
        inputs = sorted({Path(str(p)[:-3]) if str(p).endswith(".gz") else p for p in inputs})
    cor = [d / "rcorr" / f"{runout}.TRIM_{i}P.cor.fq" for i in (1, 2)]

    if cleaned:
        require((reports / f"{runout}.cleanup.done").is_file(), "no cleanup.done")
        for gone in (d / "shuck", d / "quants", asm / "diamond", asm / "working",
                     asm / f"{runout}.ORP.intermediate.fasta", d / f"{runout}.sorted.bam"):
            require(not gone.exists(), f"cleanup left {gone.relative_to(d)} behind")
        for a in inputs:
            require(Path(f"{a}.gz").is_file(), f"cleanup did not keep {a.name}.gz")
            require(not a.exists(), f"cleanup kept {a.name} uncompressed as well")
        if not corrected:
            for c in cor:
                require(Path(f"{c}.gz").is_file(), f"corrected reads {c.name}.gz not kept")
        for t in ("TRIM_1P.fastq", "TRIM_2P.fastq"):
            require(not (d / "rcorr" / f"{runout}.{t}").exists(), f"trimmed {t} not reclaimed")
    else:
        require(not (reports / f"{runout}.cleanup.done").exists(),
                "cleanup.done written despite --no-cleanup")
        for kept in (d / "shuck" / runout, d / "quants", asm / "diamond", asm / "working",
                     asm / f"{runout}.ORP.intermediate.fasta", d / f"{runout}.sorted.bam"):
            require(kept.exists(), f"--no-cleanup lost {kept.relative_to(d)}")
        for a in inputs:
            require(a.is_file(), f"--no-cleanup did not keep {a.name} uncompressed")
    return names, fields


def log_has(log, text, why):
    require(text in log, f"{why}: {text!r} not in the log")


def log_lacks(log, text, why):
    require(text not in log, f"{why}: {text!r} in the log")


def finished_run(ctx, name, runout=None):
    """The directory of a prerequisite's finished, cleaned run, or Skip."""
    d = ctx.workdir / name
    if not (d / "reports" / f"{runout or name}.cleanup.done").is_file():
        raise Skip(f"needs a finished {name} run in {ctx.workdir}")
    return d


# -- quick tier: no conda --------------------------------------------------------

@case("quick")
def case_sample_data(ctx, run):
    """the sample read pair is present and paired"""
    counts = []
    for r in (READ1, READ2):
        require(r.is_file(), f"{r} missing")
        with gzip.open(r, "rt") as f:
            counts.append(sum(1 for _ in f) // 4)
    require(counts[0] == counts[1] and counts[0] > 0, f"read counts differ: {counts}")


@case("quick")
def case_cli_version_help(ctx, run):
    """--version and --help on all three entry points"""
    out = run([sys.executable, OYSTER, "--version"])
    require(VERSION in out and "chowder" not in out, f"oyster --version said {out!r}")
    out = run([sys.executable, CHOWDER, "--version"])
    require(VERSION in out and "(chowder)" in out, f"chowder --version said {out!r}")
    for script in (OYSTER, CHOWDER, STRANDEVAL):
        run([sys.executable, script, "--help"])


@case("quick")
def case_cli_rejections(ctx, run):
    """argparse refuses bad flags before any work"""
    r = ["--read1", READ1, "--read2", READ2]
    bad_oyster = [
        (["--read1", READ1], "--read2"),
        (r + ["--spades1-kmer", "32"], "odd"),
        (r + ["--spades1-kmer", "51,31"], "ascending"),
        (r + ["--spades2-kmer", "31,31"], "distinct"),
        (r + ["--spades2-kmer", "50%,61"], "mix"),
        (r + ["--spades2-kmer", "120%"], "below 100%"),
        (r + ["--spades2-kmer", "131"], "less than 128"),
        (r + ["--spades1-kmer", "big"], "expected"),
        (r + ["--strand", "XX"], "invalid choice"),
        (r + ["--cpu", "eight"], "invalid int"),
    ]
    for args, msg in bad_oyster:
        out = run([sys.executable, OYSTER] + args, expect_rc=2)
        require(msg in out, f"oyster {args[-2:]}: expected {msg!r} in {out[-300:]!r}")
    bad_chowder = [
        (r, "--assemblies"),
        (r + ["--assemblies", "a.fa", "b.fa", "--assembly-order", "random"], "invalid choice"),
        (r + ["--assemblies", "a.fa", "b.fa", "--strand", "RF"], "unrecognized"),
    ]
    for args, msg in bad_chowder:
        out = run([sys.executable, CHOWDER] + args, expect_rc=2)
        require(msg in out, f"chowder {args[-2:]}: expected {msg!r} in {out[-300:]!r}")
    # Caught in Chowder.__init__, before the tool check.
    out = run([sys.executable, CHOWDER] + r + ["--assemblies", "a.fa", "b.fa",
                                               "--labels", "one"], expect_rc="nonzero")
    require("--labels: got 1 for 2 assemblies" in out, f"--labels mismatch: {out[-300:]!r}")
    out = run([sys.executable, STRANDEVAL, "--assembly", "/nonexistent.fa",
               "--read1", READ1, "--read2", READ2], expect_rc="nonzero")
    require("no such file" in out, f"strandeval missing assembly: {out[-300:]!r}")


@case("quick")
def case_kmer_specs(ctx, run):
    """k-mer spec parsing and resolution against read length"""
    from oyster import parse_kmer_spec, resolve_kmers
    require(parse_kmer_spec("auto") is None, "auto")
    require(parse_kmer_spec(" AUTO ") is None, "AUTO")
    require(parse_kmer_spec("31,51") == [31, 51], "31,51")
    require(parse_kmer_spec("60%,75%") == [0.6, 0.75], "60%,75%")
    require(resolve_kmers([0.6, 0.75], 100) == [59, 75], "60/75% of 100")
    require(resolve_kmers([0.6, 0.75], 150) == [89, 113], "60/75% of 150")
    with contextlib.redirect_stdout(io.StringIO()):
        clamped = resolve_kmers([0.6, 0.75], 250)
    require(clamped == [127], "clamped to 127 and deduplicated")
    require(resolve_kmers([31, 51], 100) == [31, 51], "absolute passes through")
    require(resolve_kmers(None, 100) is None, "auto passes through")
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            resolve_kmers([0.1], 100)
        raise CheckFailed("10% of 100 (k=9) was accepted")
    except SystemExit:
        pass


@case("quick")
def case_chowder_labels_and_order(ctx, run):
    """chowder's labels, collisions, reserved names and seeded order"""
    from chowder import derive_labels, sanitise_label, shuffled_order
    require(sanitise_label("my-asm.Trinity.fasta.gz") == "my_asm_Trinity", "sanitise")
    require(sanitise_label("...fa") == "assembly", "empty label")
    labels = derive_labels([Path("a/trinity.fasta"), Path("b/trinity.fasta")])
    require(labels == ["trinity_1", "trinity_2"], f"collision: {labels}")
    labels = derive_labels([Path("b/trinity.fasta"), Path("a/trinity.fasta")])
    require(labels == ["trinity_2", "trinity_1"], f"collision numbered by path: {labels}")
    labels = derive_labels([Path("trinity.fasta"), Path("trinity.fa"), Path("trinity_1.fasta")])
    require(len(set(l.lower() for l in labels)) == 3 and "trinity_1" in labels,
            f"suffix must dodge taken labels: {labels}")
    labels = derive_labels([Path("orp.fasta"), Path("shucked.fa.gz")])
    require(labels == ["orp_input", "shucked_input"], f"reserved: {labels}")
    labels = derive_labels([Path("x.fa"), Path("y.fa")], ["My Asm", "other"])
    require(labels == ["My_Asm", "other"], f"explicit labels: {labels}")
    pairs = [(c, Path(c)) for c in "abcdef"]
    first = shuffled_order(pairs, 23894)
    require(first == shuffled_order(list(reversed(pairs)), 23894),
            "shuffled order depends on input order")
    require(first != shuffled_order(pairs, 7) or first != shuffled_order(pairs, 8),
            "--seed has no effect")


@case("quick")
def case_chowder_ingest_rename(ctx, run):
    """ingest prefixes contig names, keeps descriptions, reads .gz, refuses dupes"""
    from chowder import Chowder
    from oyster import Assembly
    a = Assembly("lab.fasta", "lab", "lab", "LAB")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = tmp / "in.fa.gz"
        with gzip.open(src, "wt") as f:
            f.write(">c1 len=4\nACGT\n> c2 desc here\nAC\nGT\n>\nTTTT\n")
        out = tmp / "out.fa"
        n = Chowder.ingest_one(None, a, src, out)
        require(n == 3, f"counted {n} contigs")
        got = out.read_text()
        want = ">lab_c1 len=4\nACGT\n>lab_c2 desc here\nAC\nGT\n>lab_contig3\nTTTT\n"
        require(got == want, f"renamed fasta:\n{got}")
        dup = tmp / "dup.fa"
        dup.write_text(">x\nA\n>x\nC\n")
        try:
            Chowder.ingest_one(None, a, dup, tmp / "o2.fa")
            raise CheckFailed("duplicate contig names were accepted")
        except SystemExit as e:
            require("more than once" in str(e), f"duplicate message: {e}")


@case("quick")
def case_version_compare(ctx, run):
    """version_below, as the pyTransRate preflight uses it"""
    from oyster import version_below
    require(version_below("2.2.1", "2.2.2"), "2.2.1 < 2.2.2")
    require(not version_below("2.2.2", "2.2.1"), "2.2.2 !< 2.2.1")
    require(not version_below("2.10.0", "2.2.1"), "2.10.0 !< 2.2.1")
    require(not version_below("2.2.2", "2.2.2"), "2.2.2 !< 2.2.2")


@case("quick")
def case_twotrack_select(ctx, run):
    """tests/test_twotrack_select.py (in the orp env when conda is here)"""
    script = HERE / "test_twotrack_select.py"
    if ctx.conda:
        out = run(["conda", "run", "--no-capture-output", "-n", "orp", "python", script])
    else:
        out = run([sys.executable, script])
    if "skip" in out.lower():
        raise Skip("its tools are not on PATH")


# -- preflight tier: conda, no real work -------------------------------------------

def tiny_fasta(path, names=("c1", "c2")):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{n}\n{'ACGT' * 60}\n" for n in names))
    return path


@case("preflight")
def case_preflight_oyster(ctx, run):
    """oyster.py refuses missing reads, oversized/undersized k, a missing lineage"""
    d = ctx.case_dir("preflight_oyster", fresh=True)
    r = ["--read1", READ1, "--read2", READ2]
    out = run.oyster(d, "--read1", d / "nope.fq.gz", "--read2", READ2, "--runout", "p1",
                     expect_rc="nonzero")
    log_has(out, "READ1 FILE DOES NOT EXIST", "missing read1")
    out = run.oyster(d, *r, "--runout", "p2", "--spades1-kmer", "101", expect_rc="nonzero")
    log_has(out, "NOT AT LEAST 101 BP LONG", "k above read length")
    out = run.oyster(d, *r, "--runout", "p3", "--spades2-kmer", "10%,20%", expect_rc="nonzero")
    log_has(out, "TOO SMALL TO ASSEMBLE WITH", "percentage k below 15")
    out = run.oyster(d, *r, "--runout", "p4", "--lineage", "nonexistent_odb99",
                     expect_rc="nonzero")
    log_has(out, "BUSCO lineage 'nonexistent_odb99'", "missing lineage")
    require(not list((d / "rcorr").glob("*.fq*")), "a refused run trimmed reads anyway")


@case("preflight")
def case_preflight_chowder(ctx, run):
    """chowder.py refuses one, empty, missing, duplicate-named or self-staged assemblies"""
    d = ctx.case_dir("preflight_chowder", fresh=True)
    r = ["--read1", READ1, "--read2", READ2]
    a, b = tiny_fasta(d / "in" / "a.fa"), tiny_fasta(d / "in" / "b.fa")
    out = run.chowder(d, *r, "--runout", "c1", "--assemblies", a, expect_rc="nonzero")
    log_has(out, "give it at least two", "single assembly")
    empty = d / "in" / "empty.fa"
    empty.write_text("")
    out = run.chowder(d, *r, "--runout", "c2", "--assemblies", a, empty, expect_rc="nonzero")
    log_has(out, "assembly is empty", "empty assembly")
    out = run.chowder(d, *r, "--runout", "c3", "--assemblies", a, d / "in" / "missing.fa",
                      expect_rc="nonzero")
    log_has(out, "assembly not found", "missing assembly")
    dup = tiny_fasta(d / "in" / "dup.fa", names=("x", "y", "x"))
    out = run.chowder(d, *r, "--runout", "c4", "--assemblies", a, dup, expect_rc="nonzero")
    log_has(out, "more than once", "duplicate contig names")
    staged = tiny_fasta(d / "assemblies" / "ingested" / "c5.a.fasta")
    out = run.chowder(d, *r, "--runout", "c5", "--assemblies", staged, b, expect_rc="nonzero")
    log_has(out, "where chowder writes", "input inside the staging directory")
    out = run.chowder(d, *r, "--runout", "c6", "--assemblies", a, b,
                      "--lineage", "nonexistent_odb99", expect_rc="nonzero")
    log_has(out, "BUSCO lineage", "missing lineage")
    out = run.chowder(d, "--read1", d / "nope.fq", "--read2", READ2, "--runout", "c7",
                      "--assemblies", a, b, expect_rc="nonzero")
    log_has(out, "READ1 FILE DOES NOT EXIST", "missing read1")


# -- full tier: complete runs ----------------------------------------------------

@case("full")
def case_oyster_default(ctx, run):
    """full oyster.py run, defaults; then a rerun is a no-op"""
    d = ctx.case_dir("oyster_default", fresh=True)
    out = run.oyster(d, "--read1", READ1, "--read2", READ2, "--runout", "oyster_default")
    check_finished(d, "oyster_default", "oyster")
    trinity = commands(out, "Trinity")
    require(len(trinity) == 2, f"expected Trinity's two phases, found {len(trinity)} runs")
    require(all("--no_normalize_reads" in c for c in trinity), "Trinity normalized by default")
    require("--full_cleanup" in trinity[-1], "Trinity phase 2 ran without --full_cleanup")
    log_lacks(out, "--SS_lib_type", "unstranded by default")
    log_has(out, "[kmer] spadesauto: k = auto", "default spades1 k")
    log_has(out, "[kmer] spadeshigh: k = 59,75", "default spades2 k (60%,75% of 100)")
    require(all("--kmer 32" in c for c in commands(out, "transabyss")), "default Trans-ABySS k")
    log_has(out, "with strandeval beside it", "--max-parallel 2 runs strandeval beside")
    log_has(out, "THERE IS NO LOW STUFF", "--tpm-filt 0 drops nothing")

    orp = d / "assemblies" / "oyster_default.ORP.fasta"
    mtime = orp.stat().st_mtime
    out = run.oyster(d, "--read1", READ1, "--read2", READ2, "--runout", "oyster_default")
    log_has(out, "already finished", "rerun of a finished run")
    require(orp.stat().st_mtime == mtime, "rerun rewrote the assembly")


@case("full")
def case_oyster_options(ctx, run):
    """full oyster.py run with every assembly/report option off its default"""
    d = ctx.case_dir("oyster_options", fresh=True)
    out = run.oyster(
        d, "--read1", READ1, "--read2", READ2, "--runout", "oyster_options",
        "--strand", "RF", "--normalize-reads", "--tpm-filt", "1",
        "--spades1-kmer", "31", "--spades2-kmer", "45,65", "--transabyss-kmer", "25",
        "--max-parallel", "1", "--busco-threads", "2",
        "--pytransrate-args", "--location-size 5",
    )
    check_finished(d, "oyster_options", "oyster")
    trinity = commands(out, "Trinity")
    require(trinity and all("--SS_lib_type RF" in c for c in trinity), "Trinity not run RF")
    require(not any("--no_normalize_reads" in c for c in trinity),
            "--normalize-reads did not reach Trinity")
    spades = commands(out, "rnaspades.py")
    require(len(spades) >= 2 and all("--ss-rf" in c for c in spades), "rnaSPAdes not run RF")
    require(any("-k 31 " in c for c in spades), "--spades1-kmer 31 not passed")
    require(any("-k 45,65 " in c for c in spades), "--spades2-kmer 45,65 not passed")
    tab = commands(out, "transabyss")
    require(tab and all("--SS" in c and "--kmer 25" in c for c in tab),
            "Trans-ABySS not run --SS --kmer 25")
    require(all("--cpu 2" in c for c in commands(out, "busco")), "--busco-threads 2 not passed")
    ptr = commands(out, "pytransrate")
    require(len(ptr) >= 2 and all("--location-size 5" in c for c in ptr),
            "--pytransrate-args not passed to both pyTransRate runs")
    log_lacks(out, "beside it", "--max-parallel 1 runs nothing side by side")
    log_has(out, "PART: TPM_FILT MAKE LOW AND HIGH", "the TPM filter ran")


@case("full")
def case_oyster_resume(ctx, run):
    """--keep-intermediates keeps all; a half-done step re-runs; cleanup follows"""
    d = ctx.case_dir("oyster_resume", fresh=True)
    args = ["--read1", READ1, "--read2", READ2, "--runout", "oyster_resume"]
    out = run.oyster(d, *args, "--keep-intermediates")
    check_finished(d, "oyster_resume", "oyster", cleaned=False)
    log_has(out, "--no-cleanup given", "cleanup skipped")
    require(not any("--full_cleanup" in c for c in commands(out, "Trinity")),
            "Trinity --full_cleanup used despite --keep-intermediates")
    require((d / "assemblies" / "oyster_resume.trinity").is_dir(), "Trinity work dir removed")
    require((d / "assemblies" / "oyster_resume.transabyss").is_dir(), "Trans-ABySS work dir removed")
    require(list((d / "assemblies").glob("oyster_resume.spades_k*")), "rnaSPAdes work dirs removed")

    # A salmon run that "died": its marker is still there, its output looks fine.
    marker = d / "reports" / ".oyster_resume.salmon.running"
    marker.touch()
    out = run.oyster(d, *args)
    log_has(out, "[salmon] an attempt started", "a left-behind marker forces a re-run")
    log_has(out, "=== salmon -- start", "salmon re-ran")
    log_has(out, "[run_rcorrector] up to date, skipping", "finished steps are skipped")
    log_lacks(out, "=== run_trimmomatic -- start", "trimming is not redone")
    log_lacks(out, "=== run_transabyss -- start", "assembling is not redone")
    check_finished(d, "oyster_resume", "oyster")


@case("full", needs=("oyster_default",))
def case_oyster_corrected(ctx, run):
    """--trimmed-corrected-reads (gzipped), --strand FR, multi-k and auto k"""
    src = finished_run(ctx, "oyster_default")
    c1, c2 = (src / "rcorr" / f"oyster_default.TRIM_{i}P.cor.fq.gz" for i in (1, 2))
    d = ctx.case_dir("oyster_corrected", fresh=True)
    out = run.oyster(d, "--read1", c1, "--read2", c2, "--runout", "oyster_corrected",
                     "--trimmed-corrected-reads", "--strand", "FR",
                     "--spades1-kmer", "35,55", "--spades2-kmer", "auto")
    check_finished(d, "oyster_corrected", "oyster", corrected=True)
    log_has(out, "trimmomatic and rcorrector skipped", "corrected reads used as given")
    log_lacks(out, "=== run_trimmomatic", "trimmomatic skipped")
    log_lacks(out, "=== run_rcorrector", "rcorrector skipped")
    require(all("--SS_lib_type FR" in c for c in commands(out, "Trinity")), "Trinity not run FR")
    spades = commands(out, "rnaspades.py")
    require(all("--ss-fr" in c for c in spades), "rnaSPAdes not run FR")
    require(any("-k 35,55 " in c for c in spades), "--spades1-kmer 35,55 not passed")
    require(sum(" -k " in c for c in spades) == 1, "--spades2-kmer auto still passed -k")
    require(c1.is_file() and c2.is_file(), "the user's corrected reads were removed")
    require(not (d / "rcorr" / "oyster_corrected.TRIM_1P.cor.fq").exists(),
            "the decompressed copy of the reads was not removed")


@case("full", needs=("oyster_default",))
def case_chowder_default(ctx, run):
    """chowder.py on oyster_default's four gzipped assemblies; rerun no-op"""
    src = finished_run(ctx, "oyster_default")
    inputs = [src / "assemblies" / f"oyster_default.{a}.gz" for a in OYSTER_ASSEMBLIES]
    d = ctx.case_dir("chowder_default", fresh=True)
    args = ["--read1", READ1, "--read2", READ2, "--runout", "chowder_default",
            "--assemblies", *inputs]
    out = run.chowder(d, *args)
    labels = ["oyster_default_spadesauto", "oyster_default_spadeshigh",
              "oyster_default_transabyss", "oyster_default_trinity_Trinity"]
    _, fields = check_finished(d, "chowder_default", "chowder", labels=labels)
    log_has(out, "This is chowder, NOT a full ORP run", "chowder banner")
    for tool in ("Trinity", "rnaspades.py", "transabyss"):
        require(not commands(out, tool), f"chowder ran {tool}")
    ingest = (d / "assemblies" / "chowder_default.ingest.done").read_text()
    require(ingest.startswith("# assembly order: shuffled, seed 23894"), "order not recorded")
    require(sorted(l.split("\t")[1] for l in ingest.splitlines()[1:]) == sorted(labels),
            f"ingest.done labels: {ingest}")
    for l in labels:
        require(f"UNIQUE GENES {l.upper()}" in fields, f"no qualreport line for {l}")

    # The same set listed in another order gives the same order.
    reorder = ctx.case_dir("chowder_default/reordered", fresh=True)
    out = run.chowder(reorder, "--read1", READ1, "--read2", READ2, "--runout", "reord",
                      "--assemblies", *reversed(inputs), "--lineage", "nonexistent_odb99",
                      expect_rc="nonzero")
    order = [l.split()[0] for l in out.split("Merging 4 assemblies", 1)[-1].splitlines()[2:6]] \
        if "Merging 4 assemblies" in out else []
    want = [l.split("\t")[1] for l in ingest.splitlines()[1:]]
    require(order == want, f"listing the inputs reversed changed the order: {order} vs {want}")

    orp = d / "assemblies" / "chowder_default.ORP.fasta"
    mtime = orp.stat().st_mtime
    out = run.chowder(d, *args)
    log_has(out, "already finished", "rerun of a finished merge")
    require(orp.stat().st_mtime == mtime, "rerun rewrote the assembly")


@case("full", needs=("oyster_default",))
def case_chowder_options(ctx, run):
    """chowder.py --labels, given order, corrected reads, TPM 1, no parallel, no cleanup"""
    src = finished_run(ctx, "oyster_default")
    a = src / "assemblies"
    inputs = [a / "oyster_default.trinity.Trinity.fasta.gz",
              a / "oyster_default.transabyss.fasta.gz",
              a / "oyster_default.spadeshigh.fasta.gz",
              a / "oyster_default.spadesauto.fasta.gz"]
    labels = ["tri", "tab", "sphi", "spauto"]
    c1, c2 = (src / "rcorr" / f"oyster_default.TRIM_{i}P.cor.fq.gz" for i in (1, 2))
    d = ctx.case_dir("chowder_options", fresh=True)
    out = run.chowder(d, "--read1", c1, "--read2", c2, "--runout", "chowder_options",
                      "--assemblies", *inputs, "--labels", *labels,
                      "--assembly-order", "given", "--corrected-reads", "--tpm-filt", "1",
                      "--max-parallel", "1", "--no-cleanup", "--busco-threads", "2")
    names, _ = check_finished(d, "chowder_options", "chowder", cleaned=False,
                              corrected=True, labels=labels)
    ingest = (d / "assemblies" / "chowder_options.ingest.done").read_text()
    require(ingest.startswith("# assembly order: as given"), "given order not recorded")
    order = [l.split("\t")[1] for l in ingest.splitlines()[1:]]
    require(order == labels, f"--assembly-order given not honoured: {order}")
    log_has(out, "trimmomatic and rcorrector skipped", "--corrected-reads")
    log_lacks(out, "beside it", "--max-parallel 1")
    require(all("--cpu 2" in c for c in commands(out, "busco")), "--busco-threads 2")

    # --tpm-filt 1: every contig kept is expressed, or kept for its swissprot hit.
    w = d / "assemblies" / "working"
    keep = set()
    for f in ("chowder_options.HIGHEXP.txt", "chowder_options.donotremove.list"):
        if (w / f).exists():
            keep |= {l.strip() for l in (w / f).read_text().splitlines() if l.strip()}
    low = (w / "chowder_options.LOWEXP.txt")
    require(low.exists(), "no LOWEXP list")
    stray = [n for n in names if n not in keep]
    require(not stray, f"{len(stray)} contigs below TPM 1 with no hit kept, e.g. {stray[:3]}")


@case("full", needs=("oyster_default",))
def case_chowder_collide(ctx, run):
    """two inputs with the same filename, --seed, the two-assembly minimum"""
    src = finished_run(ctx, "oyster_default")
    d = ctx.case_dir("chowder_collide", fresh=True)
    x, y = d / "inputs" / "x" / "asm.fasta.gz", d / "inputs" / "y" / "asm.fasta.gz"
    for dst, name in ((x, "trinity.Trinity.fasta.gz"), (y, "transabyss.fasta.gz")):
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(src / "assemblies" / f"oyster_default.{name}", dst)
    run.chowder(d, "--read1", READ1, "--read2", READ2, "--runout", "chowder_collide",
                "--assemblies", y, x, "--seed", "7")
    check_finished(d, "chowder_collide", "chowder", labels=["asm_1", "asm_2"])
    ingest = (d / "assemblies" / "chowder_collide.ingest.done").read_text()
    require(ingest.startswith("# assembly order: shuffled, seed 7"), "--seed not recorded")
    rows = {l.split("\t")[1]: l.split("\t")[2] for l in ingest.splitlines()[1:]}
    require(rows.get("asm_1") == str(x) and rows.get("asm_2") == str(y),
            f"colliding labels not numbered by path: {rows}")


@case("full", needs=("oyster_default",))
def case_strandeval_standalone(ctx, run):
    """scripts/strandeval.py on oyster_default's assembly, with and without --no-cleanup"""
    src = finished_run(ctx, "oyster_default")
    d = ctx.case_dir("strandeval_standalone", fresh=True)
    asm = src / "assemblies" / "oyster_default.ORP.fasta"
    base = [sys.executable, STRANDEVAL, "--assembly", asm, "--read1", READ1, "--read2", READ2,
            "--cpu", ctx.cpu, "--pairs", "5000", "--dir", d]
    out = run(base + ["--runout", "s1"])
    log_has(out, "STRAND EXAMINATION HISTOGRAM", "histogram printed")
    summary = d / "reports" / "s1.strandeval_summary.txt"
    require(summary.is_file() and "STRAND EXAMINATION" in summary.read_text(), "no summary")
    require("properly paired" in (d / "reports" / "s1.flagstat").read_text(), "no flagstat")
    require(not (d / "s1.sorted.bam").exists(), "BAM not cleaned up")
    run(base + ["--runout", "s2", "--no-cleanup"])
    require((d / "s2.sorted.bam").is_file() and (d / "s2.dat").is_file(),
            "--no-cleanup did not keep the BAM and table")
    out = run([sys.executable, STRANDEVAL, "--assembly", asm, "--read1", READ1, "--read2", READ2,
               "--pairs", "5000", "--cpu", ctx.cpu], cwd=d)
    require((d / "reports" / "oyster_default.strandeval_summary.txt").is_file(),
            "default --runout/--dir should be the assembly name and the cwd")


# -- driver ----------------------------------------------------------------------

TIERS = {"quick": ("quick",), "preflight": ("quick", "preflight"),
         "full": ("quick", "preflight", "full"), "all": ("quick", "preflight", "full")}


def run_case(ctx, c):
    start = time.time()
    try:
        if c.tier != "quick" and not ctx.conda:
            raise Skip("conda is not on PATH")
        c.func(ctx, Runner(ctx, c.name))
        status, detail = "PASS", ""
    except Skip as e:
        status, detail = "SKIP", str(e)
    except CheckFailed as e:
        status, detail = "FAIL", str(e)
    except Exception:
        status, detail = "FAIL", traceback.format_exc()
    return status, detail, time.time() - start


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__.split("\n\n", 1)[1])
    p.add_argument("--tier", choices=sorted(TIERS), default="all",
                   help="how far to go (default: all)")
    p.add_argument("--workdir", default="orp_release_check",
                   help="where the runs go (default: ./orp_release_check)")
    p.add_argument("--cpu", type=int, default=8, help="--cpu for each run (default: 8)")
    p.add_argument("--mem", type=int, default=32, help="--mem for each run, GB (default: 32)")
    p.add_argument("--jobs", type=int, default=1,
                   help="full-tier runs at once; budget --jobs x --cpu cores (default: 1)")
    p.add_argument("--only", nargs="+", metavar="CASE", help="run just these cases")
    p.add_argument("--list", action="store_true", help="list the cases and exit")
    args = p.parse_args()

    if args.list:
        for c in CASES:
            needs = f"  (needs {', '.join(c.needs)})" if c.needs else ""
            print(f"{c.tier:<10} {c.name:<24} {c.doc}{needs}")
        return
    names = {c.name for c in CASES}
    if args.only:
        unknown = set(args.only) - names
        if unknown:
            sys.exit(f"unknown case(s): {', '.join(sorted(unknown))}; see --list")
        selected = [c for c in CASES if c.name in args.only]
    else:
        selected = [c for c in CASES if c.tier in TIERS[args.tier]]

    ctx = Context(args)
    print(f"ORP {VERSION} release check -> {ctx.workdir}")
    print(f"  python {sys.version.split()[0]}, conda {'found' if ctx.conda else 'NOT found'}, "
          f"{len(selected)} case(s), --cpu {args.cpu} --mem {args.mem} --jobs {args.jobs}\n")

    results, lock = {}, threading.Lock()

    def report(c, status, detail, secs):
        with lock:
            results[c.name] = (status, detail, secs)
            print(f"[{status}] {c.name:<24} {int(secs // 60):3d}m{int(secs % 60):02d}s"
                  + (f"  {detail.splitlines()[0]}" if detail and status != "PASS" else ""),
                  flush=True)

    # Quick and preflight cases are cheap: in order, one at a time. Full runs
    # go through a small scheduler that honours `needs`.
    for c in [c for c in selected if c.tier != "full"]:
        report(c, *run_case(ctx, c))
    pending = [c for c in selected if c.tier == "full"]
    running = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        while pending or running:
            for c in list(pending):
                states = [results.get(n, (None,))[0] for n in c.needs]
                waiting = any(n in running or n in [q.name for q in pending] for n in c.needs)
                if waiting:
                    continue
                pending.remove(c)
                if any(s in ("FAIL", "SKIP") for s in states):
                    report(c, "SKIP", f"prerequisite {', '.join(c.needs)} did not pass", 0)
                    continue
                running[c.name] = pool.submit(run_case, ctx, c)
            if not running:
                continue
            done, _ = concurrent.futures.wait(list(running.values()),
                                              return_when=concurrent.futures.FIRST_COMPLETED)
            for name in [n for n, f in running.items() if f in done]:
                c = next(c for c in selected if c.name == name)
                report(c, *running.pop(name).result())

    counts = {s: sum(1 for r in results.values() if r[0] == s) for s in ("PASS", "FAIL", "SKIP")}
    lines = [f"ORP {VERSION} release check, {time.strftime('%Y-%m-%d %H:%M:%S')}",
             f"{counts['PASS']} passed, {counts['FAIL']} failed, {counts['SKIP']} skipped", ""]
    for c in selected:
        status, detail, secs = results[c.name]
        lines.append(f"[{status}] {c.name}  ({int(secs)}s)")
        if status != "PASS" and detail:
            lines += ["    " + l for l in detail.rstrip().splitlines()]
    (ctx.workdir / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines[1:2]) + f"\nsummary: {ctx.workdir / 'summary.txt'}"
          f"\nlogs:    {ctx.logs}/")
    sys.exit(1 if counts["FAIL"] else 0)


if __name__ == "__main__":
    main()
