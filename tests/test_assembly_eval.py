#!/usr/bin/env python3
"""Test that assembly_eval.py redoes a step whose inputs changed since it finished.

usage: python tests/test_assembly_eval.py

No conda needed: the three steps, preflight and the report are stubbed, so
what is checked is only which steps main() runs. A plain resume runs none; a
changed assembly redoes all three, changed reads pytransrate and strandeval,
a changed --lineage BUSCO alone. A sentinel from before they recorded
anything is trusted while newer than the inputs. Exits non-zero on a failure.
"""

import contextlib
import io
import os
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import assembly_eval  # noqa: E402

failures = []
STEPS = set(assembly_eval.STEPS)


def check(cond, msg):
    if not cond:
        failures.append(msg)
        print(f"FAIL: {msg}")


def age(path, seconds):
    t = time.time() - seconds
    os.utime(path, (t, t))


def run(work, *extra):
    """Run main() on the files in `work`; return the set of steps it ran."""
    sys.argv = ["assembly_eval.py", "--assembly", str(work / "a.fasta"),
                "--read1", str(work / "r1.fq"), "--read2", str(work / "r2.fq"),
                "--outdir", str(work / "eval"), "--no-conda", *extra]
    ev = assembly_eval.Eval(assembly_eval.parse_args())
    ran = set()
    for step in STEPS:
        setattr(ev, step, lambda cpu, mem, step=step: ran.add(step))
    ev.preflight = ev.reportgen = lambda: None
    with contextlib.redirect_stdout(io.StringIO()):
        ev.main()
    return ran


with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    for name in ("a.fasta", "r1.fq", "r2.fq"):
        (work / name).write_text("x\n")
        age(work / name, 7200)

    check(run(work) == STEPS, "a fresh run runs every step")
    check(run(work) == set(), "a plain resume reruns something")

    (work / "a.fasta").write_text("xy\n")
    age(work / "a.fasta", 3600)  # changed, but older than every sentinel
    check(run(work) == STEPS, "a changed assembly did not redo every step")

    (work / "r2.fq").write_text("xy\n")
    age(work / "r2.fq", 3600)
    check(run(work) == {"pytransrate", "strandeval"},
          "changed reads did not redo exactly pytransrate and strandeval")

    check(run(work, "--lineage", "metazoa_odb12") == {"busco"},
          "a changed --lineage did not redo exactly BUSCO")
    check(run(work, "--lineage", "metazoa_odb12", "--pytransrate-args", "--location-size 5")
          == {"pytransrate"}, "changed --pytransrate-args did not redo pytransrate")
    check(run(work, "--lineage", "metazoa_odb12", "--pytransrate-args", "--location-size 5",
              "--force") == STEPS, "--force did not redo every step")

    sentinels = work / "eval" / "reports"
    for step in STEPS:
        (sentinels / f"a.{step}.done").write_text("")  # from before sentinels recorded inputs
    check(run(work, "--lineage", "metazoa_odb12", "--pytransrate-args", "--location-size 5")
          == set(), "an empty sentinel newer than the inputs was not trusted")
    for step in STEPS:
        (sentinels / f"a.{step}.done").write_text("")
        age(sentinels / f"a.{step}.done", 1800)
    (work / "a.fasta").write_text("xyz\n")  # newer than the empty sentinels
    check(run(work, "--lineage", "metazoa_odb12", "--pytransrate-args", "--location-size 5")
          == STEPS, "an empty sentinel older than the assembly was trusted")

if failures:
    print(f"{len(failures)} failure(s)")
    sys.exit(1)
print("all assembly_eval resume checks passed")
