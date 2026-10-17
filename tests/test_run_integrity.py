#!/usr/bin/env python3
"""Test that a run keeps to its own files and that a resume notices changed inputs.

usage: python tests/test_run_integrity.py

No conda needed. Each check builds a real Pipeline (or Chowder) in a
temporary --dir and calls the method under test directly:

  run_name         --runout refuses '.', so <runout>.* can't match another run
  cleanup          cleaning `sample` leaves `sample_rep2`'s files alone
  relative_reads   read paths typed relative to the launch directory still
                   name the same files from --dir, symlinks unresolved
  settings_stamp   dated to its outputs when new, untouched when unchanged,
                   newer than the outputs when a setting changes
  reads_stamp      a read file swapped for a different, older one is noticed
  corrected_relink --trimmed-corrected-reads pointed at another, older file
                   relinks, and what was built from the old link goes stale
  trinity_dir      Trinity's checkpoints are cleared when phase 1 is re-run
                   on other inputs, kept on a plain resume
  chowder_ids      label `a` + contig `b_c` and label `a_b` + contig `c`
                   are caught instead of both becoming a_b_c
  chowder_order    a different --seed re-ingests
  strand_exam      paths are quoted, bwa and samtools sort share --cpu, and
                   the wrap hook decides the env

Exits non-zero on a failure.
"""

import argparse
import contextlib
import io
import os
import shlex
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import oyster  # noqa: E402
import chowder  # noqa: E402
from oyster import Assembly, Pipeline, run_name, strand_exam  # noqa: E402

failures = []
FASTQ = "@r\nACGT\n+\nIIII\n"


def check(cond, msg):
    if not cond:
        failures.append(msg)
        print(f"FAIL: {msg}")


def age(path, seconds):
    """Set `path`'s mtime `seconds` into the past."""
    t = time.time() - seconds
    os.utime(path, (t, t))


def pipeline(work, *extra, runout="sample", reads=("r1.fq", "r2.fq")):
    sys.argv = ["oyster.py", "--read1", str(reads[0]), "--read2", str(reads[1]),
                "--runout", runout, "--dir", str(work), *extra]
    with contextlib.redirect_stdout(io.StringIO()):
        p = Pipeline(oyster.parse_args())
        p.setup()
    return p


def quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def check_run_name():
    for good in ("sample", "SRR1789336", "my_run-2"):
        check(run_name(good) == good, f"run_name refused {good!r}")
    for bad in ("sample.rep2", "a/b", "has space", ""):
        try:
            run_name(bad)
            check(False, f"run_name accepted {bad!r}")
        except argparse.ArgumentTypeError:
            pass


def check_cleanup(tmp):
    work = tmp / "cleanup"
    work.mkdir()
    for r in ("r1.fq", "r2.fq"):
        (work / r).write_text(FASTQ)
    p = pipeline(work, reads=(work / "r1.fq", work / "r2.fq"))
    other = pipeline(work, runout="sample_rep2", reads=(work / "r1.fq", work / "r2.fq"))
    mine, theirs = [], []
    for run, out in ((p, mine), (other, theirs)):
        for d in (run.assemblies_working, run.diamond_dir, run.quants_dir):
            f = d / f"{run.runout}.list1"
            f.write_text("x\n")
            out.append(f)
        f = run.assemblies_dir / f"{run.runout}.trinity.Trinity.fasta.gene_trans_map"
        f.write_text("x\n")
        out.append(f)
    quiet(p.cleanup)
    check(not any(f.exists() for f in mine), f"cleanup left its own files: {mine}")
    check(all(f.exists() for f in theirs),
          f"cleanup of 'sample' removed sample_rep2's: {[f for f in theirs if not f.exists()]}")


def check_relative_reads(tmp):
    launch = tmp / "launch"
    (launch / "reads").mkdir(parents=True)
    (launch / "reads" / "real1.fq").write_text(FASTQ)
    (launch / "reads" / "r1.fq").symlink_to("real1.fq")
    (launch / "reads" / "r2.fq").write_text(FASTQ)
    (launch / "lineage_odb12").mkdir()
    cwd = os.getcwd()
    os.chdir(launch)
    try:
        p = pipeline(tmp / "elsewhere", "--lineage", "lineage_odb12",
                     reads=("reads/r1.fq", "reads/r2.fq"))
    finally:
        os.chdir(cwd)
    check(p.read1 == launch / "reads" / "r1.fq", f"read1 is {p.read1}")
    check(p.read1.is_absolute() and p.read2.is_absolute(), "read paths not absolute")
    check(p.read1.name == "r1.fq", "read1's symlink was resolved")
    check(p.lineage == str(launch / "lineage_odb12"), f"lineage dir is {p.lineage}")
    p2 = pipeline(tmp / "elsewhere", "--lineage", "eukaryota_odb12.2")
    check(p2.lineage == "eukaryota_odb12.2", f"lineage name became {p2.lineage}")


def check_settings_stamp(tmp):
    work = tmp / "stamp"
    p = pipeline(work)
    out = work / "out.txt"
    out.write_text("x")
    age(out, 600)
    stamp = p.settings_stamp(work / "s.settings", "k=1\n", [out])
    check(not p.needs_run([out], [stamp]), "a new stamp made existing outputs stale")
    check(abs(stamp.stat().st_mtime - out.stat().st_mtime) < 1, "new stamp not dated to its anchor")
    before = stamp.stat().st_mtime
    p.settings_stamp(work / "s.settings", "k=1\n", [out])
    check(stamp.stat().st_mtime == before, "an unchanged stamp was rewritten")
    p.settings_stamp(work / "s.settings", "k=2\n", [out])
    check(p.needs_run([out], [stamp]), "a changed setting left the output up to date")
    fresh = p.settings_stamp(work / "t.settings", "k=1\n", [work / "missing"])
    check(fresh.read_text() == "k=1\n", "stamp without anchors")


def check_reads_stamp(tmp):
    work = tmp / "reads"
    work.mkdir()
    for r in ("r1.fq", "r2.fq", "older.fq"):
        (work / r).write_text(FASTQ if r != "older.fq" else FASTQ * 2)
    age(work / "older.fq", 7200)
    p = pipeline(work, reads=(work / "r1.fq", work / "r2.fq"))
    trimmed = work / "trimmed"
    trimmed.write_text("x")
    stamp = p.reads_stamp([trimmed])
    check(not p.needs_run([trimmed], [p.read1, p.read2, stamp]), "first stamp made outputs stale")
    p = pipeline(work, reads=(work / "older.fq", work / "r2.fq"))
    stamp = p.reads_stamp([trimmed])
    check(p.needs_run([trimmed], [p.read1, p.read2, stamp]),
          "read1 swapped for an older file went unnoticed")


def check_corrected_relink(tmp):
    work = tmp / "relink"
    work.mkdir()
    for r in ("a1.fq", "a2.fq", "b1.fq"):
        (work / r).write_text(FASTQ)
    age(work / "b1.fq", 7200)
    p = pipeline(work, "--trimmed-corrected-reads", reads=(work / "a1.fq", work / "a2.fq"))
    quiet(p.use_corrected_reads)
    check(p.cor1().resolve() == (work / "a1.fq").resolve(), "first link")
    built = work / "assembly.fasta"
    built.write_text(">x\nA\n")
    age(built, -5)  # built after the link was made
    p = pipeline(work, "--trimmed-corrected-reads", reads=(work / "b1.fq", work / "a2.fq"))
    quiet(p.use_corrected_reads)
    check(p.cor1().resolve() == (work / "b1.fq").resolve(),
          f"cor1 still points at {p.cor1().resolve()}")
    os.utime(p.cor1(), (time.time() + 10,) * 2, follow_symlinks=False)
    check(p.needs_run([built], [p.cor1(), p.cor2()]),
          "an assembly from the old reads stayed up to date after the relink")


def check_trinity_dir(tmp):
    work = tmp / "trinity"
    p = pipeline(work)
    for c in (p.cor1(), p.cor2()):
        c.write_text(FASTQ)
    out = p.trinity_out_dir()
    out.mkdir(parents=True)
    (out / "inchworm.fa.finished").write_text("")
    quiet(p.clear_stale_trinity_dir)
    check((out / "inchworm.fa.finished").exists(), "a directory without orp.inputs was cleared")
    quiet(p.clear_stale_trinity_dir)
    check((out / "inchworm.fa.finished").exists(), "a plain resume cleared Trinity's checkpoints")
    p = pipeline(work, "--strand", "RF")
    quiet(p.clear_stale_trinity_dir)
    check(not (out / "inchworm.fa.finished").exists(), "--strand changed but checkpoints were kept")
    check("strand=RF" in (out / "orp.inputs").read_text(), "orp.inputs not rewritten")


def check_chowder_ids(tmp):
    check(chowder.labels_can_collide(["a", "a_b"]), "a / a_b")
    check(not chowder.labels_can_collide(["trinity_1", "trinity_2", "spades"]), "unrelated labels")
    work = tmp / "ids"
    work.mkdir()
    (work / "one.fa").write_text(">b_c\nACGT\n")
    (work / "two.fa").write_text(">c\nACGT\n")
    taken = {}
    quiet(chowder.Chowder.ingest_one, None, Assembly("a.fasta", "a", "a", "A"),
          work / "one.fa", work / "o1.fa", taken)
    try:
        quiet(chowder.Chowder.ingest_one, None, Assembly("a_b.fasta", "a_b", "a_b", "A_B"),
              work / "two.fa", work / "o2.fa", taken)
        check(False, "a_b_c from two inputs was accepted")
    except SystemExit as e:
        check("both become 'a_b_c'" in str(e), f"collision message: {e}")


def make_chowder(work, seed):
    sys.argv = ["chowder.py", "--read1", str(work / "r1.fq"), "--read2", str(work / "r2.fq"),
                "--assemblies", *(str(work / f"{n}.fa") for n in "abcd"),
                "--seed", str(seed), "--runout", "ch", "--dir", str(work)]
    with contextlib.redirect_stdout(io.StringIO()):
        c = chowder.Chowder(chowder.parse_args())
        c.setup()
    return c


def check_chowder_order(tmp):
    work = tmp / "order"
    work.mkdir()
    for n in "abcd":
        (work / f"{n}.fa").write_text(f">{n}\nACGT\n")
    c = make_chowder(work, 1)
    outputs = c.assembly_fasta_paths() + [c.ingest_done]
    for o in outputs:
        o.parent.mkdir(parents=True, exist_ok=True)
        o.write_text("x")
    stamp = c.ingest_stamp()
    check(not c.needs_run(outputs, c.sources + [stamp]), "first ingest stamp made outputs stale")
    seeds = [s for s in range(2, 50)
             if [a.diamond_label for a in make_chowder(work, s).assemblies]
             != [a.diamond_label for a in c.assemblies]]
    c = make_chowder(work, seeds[0])
    stamp = c.ingest_stamp()
    check(c.needs_run(outputs, c.sources + [stamp]), "a different --seed kept the old ingest")


def check_strand_exam(tmp):
    work = tmp / "strand dir"
    work.mkdir()
    (work / "s.dat").write_text("header\n" + "t1 0 0 0 0.5\n")
    ran = []
    summary = strand_exam(
        lambda cmd, cwd=None: ran.append(cmd), work / "my asm.fa", work / "r 1.fq",
        work / "r2;x.fq", "s", work, 8, flagstat_out=work / "s.flagstat",
        wrap=lambda env, cmd: ["echo", env, *map(str, cmd)], perllib="/lib", strand_script="/x.pl")
    script = next(c[-1] for c in ran if c[:3] == ["bash", "-o", "pipefail"])
    check(shlex.quote(str(work / "r 1.fq")) in script, f"read path not quoted: {script}")
    check(shlex.quote(str(work / "r2;x.fq")) in script, f"metacharacter not quoted: {script}")
    check("samtools view" not in script, "samtools view still in the pipe")
    check("bwa mem -t 6 " in script, f"bwa threads: {script}")
    check("-@ 1 " in script and "-m 512M" in script, f"samtools sort threads/memory: {script}")
    check(ran[0][:2] == ["echo", "orp_trinity"], f"wrap not used for bwa index: {ran[0]}")
    check("STRAND EXAMINATION HISTOGRAM" in summary, "summary banner")


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp).resolve()
        for name, fn in (("run_name", check_run_name), ("cleanup", check_cleanup),
                         ("relative_reads", check_relative_reads),
                         ("settings_stamp", check_settings_stamp),
                         ("reads_stamp", check_reads_stamp),
                         ("corrected_relink", check_corrected_relink),
                         ("trinity_dir", check_trinity_dir), ("chowder_ids", check_chowder_ids),
                         ("chowder_order", check_chowder_order),
                         ("strand_exam", check_strand_exam)):
            before = len(failures)
            fn() if fn is check_run_name else fn(tmp)
            print(f"{'ok  ' if len(failures) == before else 'FAIL'}  {name}")
    if failures:
        print(f"\n{len(failures)} check(s) failed")
        return 1
    print("\nall run integrity checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
