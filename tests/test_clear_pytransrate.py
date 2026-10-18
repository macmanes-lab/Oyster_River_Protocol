#!/usr/bin/env python3
"""Test oyster.py's clear_pytransrate_outdir: whole and current are both required.

usage: python tests/test_clear_pytransrate.py

A pyTransRate -o with a finished snap index, a complete BAM, its salmon
quant and logs. Cleared for a retry against the same assembly, all of it is
kept. Cleared after the assembly was rewritten (a re-filter), the index, BAM
and quant go -- they were built from the old assembly -- and logs and the
index lock stay. Cleared after the reads changed under an unchanged assembly,
the BAM, its read count and quant go and the index stays. Exits non-zero on a
failure.
"""

import os
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from oyster import BGZF_EOF, Pipeline  # noqa: E402

failures = []


def check(cond, msg):
    if not cond:
        failures.append(msg)
        print(f"FAIL: {msg}")


def build(tmp):
    tmp = Path(tmp)
    for r in ("r1.fq", "r2.fq"):
        (tmp / r).write_text("@r\nACGT\n+\nIIII\n")
        os.utime(tmp / r, (time.time() - 7200,) * 2)
    assembly = tmp / "run.ORP.fasta"
    assembly.write_text(">a\nACGT\n>b\nACGT\n")
    old = time.time() - 3600
    os.utime(assembly, (old, old))
    out = tmp / "pytransrate_run"
    (out / "run.ORP").mkdir(parents=True)
    (out / "run.ORP" / "GenomeIndex").write_text("x")
    (out / "run.ORP.index.lock").write_text("1")
    (out / "reads.bam").write_bytes(b"\x1f\x8b" + b"\0" * 30 + BGZF_EOF)
    (out / "reads.bam.align.done").write_text("")
    (out / "r1-r2-read_count.txt").write_text("100\n")
    (out / "salmon").mkdir()
    (out / "salmon" / "quant.sf").write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\na\t4\t4\t1\t1\nb\t4\t4\t1\t1\n")
    (out / "logs").mkdir()
    (out / "logs" / "snap.log").write_text("log")
    (out / "assemblies.csv").write_text("old")
    # Dated explicitly, a minute after the assembly: filesystems with coarse
    # timestamps (GPFS) would otherwise give a just-rewritten assembly the
    # same mtime as files made a moment before it.
    for path in [out, *out.rglob("*")]:
        os.utime(path, (old + 60, old + 60))
    return assembly, out


with tempfile.TemporaryDirectory() as tmp:
    assembly, out = build(tmp)
    Pipeline.clear_pytransrate_outdir(out, assembly, (Path(tmp) / "r1.fq", Path(tmp) / "r2.fq"))
    for name in ("run.ORP", "reads.bam", "reads.bam.align.done", "salmon", "logs",
                 "run.ORP.index.lock", "r1-r2-read_count.txt"):
        check((out / name).exists(), f"retry against the same assembly keeps {name}")
    check(not (out / "assemblies.csv").exists(), "retry clears assemblies.csv")

with tempfile.TemporaryDirectory() as tmp:
    assembly, out = build(tmp)
    assembly.write_text(">a\nACGT\n>b\nACGT\n>c\nACGT\n")  # re-filtered: newer than everything
    Pipeline.clear_pytransrate_outdir(out, assembly)
    for name in ("run.ORP", "reads.bam", "reads.bam.align.done", "salmon", "r1-r2-read_count.txt"):
        check(not (out / name).exists(), f"a newer assembly clears {name}")
    for name in ("logs", "run.ORP.index.lock"):
        check((out / name).exists(), f"a newer assembly keeps {name}")

with tempfile.TemporaryDirectory() as tmp:
    assembly, out = build(tmp)
    # A resumed run given different reads: same assembly, same filenames.
    (Path(tmp) / "r1.fq").write_text("@r\nTTTT\n+\nIIII\n" * 2)
    Pipeline.clear_pytransrate_outdir(out, assembly, (Path(tmp) / "r1.fq", Path(tmp) / "r2.fq"))
    for name in ("reads.bam", "reads.bam.align.done", "salmon", "r1-r2-read_count.txt"):
        check(not (out / name).exists(), f"newer reads clear {name}")
    for name in ("run.ORP", "logs", "run.ORP.index.lock"):
        check((out / name).exists(), f"newer reads keep {name}")

if failures:
    print(f"{len(failures)} failure(s)")
    sys.exit(1)
print("all clear_pytransrate_outdir checks passed")
