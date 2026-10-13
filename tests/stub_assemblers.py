#!/usr/bin/env python3
"""stub_assemblers.py - run_assemblers' lane logic with the assemblers stubbed.

No conda needed. Every assembler step is a sleep that writes its output,
except Trinity phase 2, which goes through the real run_trinity_phase2 --
session, grow, stop -- with a fake Trinity (this file, --fake-trinity) that
works through recursive_trinity.cmds the way ParaFly does: it skips the
lines already in recursive_trinity.cmds.completed and appends each one it
finishes, faster the more cores it is given. Like ParaFly's components, it
keeps a working dir (read_partitions/.../c<i>.trinity.reads.fa.out) while on
a command and removes it when done, and it notes any it finds left over --
a half-done component the real Trinity could not resume. It also starts a
grandchild (`sleep`), so a stop that missed the rest of the session would
show.

    python3 tests/stub_assemblers.py            # all cases
    python3 tests/stub_assemblers.py grow fail  # some

Times are in seconds, scaled down: phase 2 is COMMANDS * UNIT / cpu.
ParaFly's own resume behaviour is tested on the cluster instead
(experiments/transabyss_mpi/parafly_resume.sbatch).
"""

import contextlib
import io
import os
import re
import signal
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

COMMANDS = 200
UNIT = 0.2  # core-seconds per component: 2.0 s of phase 2 at 20 cpu


def fake_trinity(outdir, cpu, fasta, launches, fail_on_cpu):
    """Phase 2 as far as oyster.py can tell: ParaFly's files, then a fasta."""
    with open(launches, "a") as f:
        f.write(f"{cpu}\n")
    cmds = Path(outdir) / "recursive_trinity.cmds"
    done_f = cmds.with_name(cmds.name + ".completed")
    done = set(done_f.read_text().splitlines()) if done_f.exists() else set()
    todo = [c for c in cmds.read_text().splitlines() if c not in done]
    unit = float(os.environ.get("STUB_UNIT", UNIT))
    workdirs = [Path(re.search(r'--output "([^"]+)"', c).group(1)) for c in todo]
    stale = sum(d.is_dir() for d in workdirs)
    with open(launches, "a") as f:
        f.write(f"skipped {len(done)} ran {len(todo)} stale {stale}\n")
    # A grandchild, as ParaFly's components are: Trinity waits for them on
    # its way out, so only a stop can leave one behind -- and the stop must
    # reach it.
    grandchild = subprocess.Popen(["sleep", "1000"], stdin=subprocess.DEVNULL,
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        if cpu == fail_on_cpu:
            time.sleep(0.2)
            sys.exit(3)
        with open(done_f, "a") as out:
            for c, workdir in zip(todo, workdirs):
                workdir.mkdir(parents=True, exist_ok=True)
                time.sleep(unit / cpu)
                out.write(c + "\n")
                out.flush()
                workdir.rmdir()
        Path(fasta).write_text(">TRINITY_DN1_c0_g1_i1 len=4\nACGT\n")
    finally:
        grandchild.kill()
        grandchild.wait()


def make_pipeline(work, mpi, ta_secs, phase1_secs, launches, fail_on_cpu=0, cpu=40):
    import oyster
    sys.argv = ["oyster.py", "--read1", str(work / "r1.fq"), "--read2", str(work / "r2.fq"),
                "--runout", "stub", "--cpu", str(cpu), "--mem", "100", "--dir", str(work),
                "--normalize-reads"]
    (work / "r1.fq").write_text("@r\nACGT\n+\nIIII\n")
    (work / "r2.fq").write_text("@r\nACGT\n+\nIIII\n")
    p = oyster.Pipeline(oyster.parse_args())
    p.setup()
    p.timing_init()
    p.transabyss_use_mpi = mpi
    for c in (p.cor1(), p.cor2()):
        c.write_text("@r\nACGT\n+\nIIII\n")
    p.compress_async = lambda path: None

    def write_after(secs, path):
        time.sleep(secs)
        Path(path).write_text(">x\nACGT\n")

    def transabyss(cpu=None, mem=None):
        # ta_secs: seconds, or ("components", n) to finish once phase 2 has
        # finished n components -- the grow decision without timing races.
        p.budgets["transabyss"] = cpu
        if isinstance(ta_secs, tuple):
            done = p.trinity_out_dir() / "recursive_trinity.cmds.completed"
            deadline = time.time() + 60
            while time.time() < deadline:
                if done.exists() and oyster.count_lines(done) >= ta_secs[1]:
                    break
                time.sleep(0.02)
            write_after(0, p.assembly_fasta(oyster.TRANSABYSS))
        else:
            write_after(ta_secs, p.assembly_fasta(oyster.TRANSABYSS))

    def phase1(cpu=None, mem=None):
        p.budgets["phase1"] = cpu
        time.sleep(phase1_secs)
        out = p.trinity_out_dir()
        out.mkdir(parents=True, exist_ok=True)
        parts = out / "read_partitions" / "Fb_0"
        (out / "recursive_trinity.cmds").write_text("".join(
            f'Trinity --single "{parts}/CBin_{i // 100}/c{i}.trinity.reads.fa" '
            f'--output "{parts}/CBin_{i // 100}/c{i}.trinity.reads.fa.out" --CPU 1\n'
            for i in range(COMMANDS)))
        p.trinity_phase1_done().touch()

    def spades(outname):
        def run(cpu=None, mem=None):
            p.budgets["spades"] = cpu
            write_after(0.05, p.assembly_fasta(getattr(oyster, outname)))
        return run

    p.budgets = {}
    p.run_transabyss = transabyss
    p.run_trinity_phase1 = phase1
    p.run_spadesauto = spades("SPADES_AUTO")
    p.run_spadeshigh = spades("SPADES_HIGH")
    p.run_diamond_one = lambda query, out, cpu=None, mem=None: Path(out).write_text("")
    fasta = p.assemblies_dir / f"{p.runout}.trinity.Trinity.fasta"
    p._trinity_phase2_cmd = lambda cpu, mem: [
        sys.executable, __file__, "--fake-trinity", str(p.trinity_out_dir()), str(cpu),
        str(fasta), str(launches), str(fail_on_cpu)]
    p.sessions_seen = []
    start = p._start_trinity_phase2

    def recording_start(cpu, mem):
        proc = start(cpu, mem)
        p.sessions_seen.append(proc.pid)
        return proc
    p._start_trinity_phase2 = recording_start
    return p


def gone(pgid):
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return True
    return False


def run_case(name, mpi, ta_secs, phase1_secs, fail_on_cpu=0):
    work = Path(tempfile.mkdtemp(prefix=f"stub_{name}_"))
    launches = work / "launches.txt"
    p = make_pipeline(work, mpi, ta_secs, phase1_secs, launches, fail_on_cpu)
    error = None
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            p.run_assemblers()
    except Exception as e:  # noqa: BLE001 - the failure case expects one
        error = e
    p.log_text = buf.getvalue()
    print(p.log_text, end="")
    runs = launches.read_text().splitlines() if launches.exists() else []
    cpus = [int(r) for r in runs if r.isdigit()]
    skips = [r for r in runs if r.startswith("skipped")]
    timing = p.timing_log.read_text()
    marker = p.step_marker("run_trinity_phase2").exists()
    leftover = [pid for pid in p.sessions_seen if not gone(pid)]
    return p, cpus, skips, timing, marker, error, leftover


def check(label, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'} {label}{(': ' + detail) if detail else ''}")
    return ok


def case_early():
    """Trans-ABySS done before Stage A ends: phase 2 on all of --cpu, no restart."""
    p, cpus, skips, timing, marker, error, left = run_case("early", True, 0.1, 0.6)
    return all([
        check("no error", error is None, repr(error)),
        check("share 0.5 under MPI: 20 cpu", p.budgets.get("transabyss") == 20, str(p.budgets)),
        check("one launch, on 40", cpus == [40], str(cpus)),
        check("no grow entry", "run_trinity_phase2_grow" not in timing),
        check("no sessions left", not left, str(left)),
    ])


def case_grow():
    """Trans-ABySS done early in phase 2: one restart, on all cores, skipping what finished."""
    p, cpus, skips, timing, marker, error, left = run_case("grow", True, ("components", COMMANDS // 10), 0.2)
    ran_first = int(skips[1].split()[1]) if len(skips) > 1 else -1
    stale = int(skips[1].split()[5]) if len(skips) > 1 else -1
    return all([
        check("no error", error is None, repr(error)),
        check("launches 20 then 40", cpus == [20, 40], str(cpus)),
        check("relaunch skipped the finished ones", 0 < ran_first < COMMANDS, " | ".join(skips)),
        check("relaunch found no half-done component dirs", stale == 0, " | ".join(skips)),
        check("the stopped one's dir removal logged", "unfinished component dir(s)" in p.log_text),
        check("grow entry in timing log", "run_trinity_phase2_grow\t" in timing),
        check("grow logged", "run_trinity_phase2 grows to 40 cpu (Trans-ABySS done; " in p.log_text),
        check("phase 2 finished", "run_trinity_phase2\t" in timing and not marker),
        check("first session (and its grandchild) gone", not left, str(left)),
    ])


def case_late():
    """Trans-ABySS done near the end of phase 2: no restart."""
    p, cpus, skips, timing, marker, error, left = run_case("late", True, ("components", COMMANDS * 9 // 10), 0.2)
    return all([
        check("no error", error is None, repr(error)),
        check("one launch, on 20", cpus == [20], str(cpus)),
        check("no grow entry", "run_trinity_phase2_grow" not in timing),
        check("logged why not", "run_trinity_phase2 stays at 20 cpu" in p.log_text),
    ])


def case_fail():
    """The relaunched phase 2 fails: the step fails as before, marker left."""
    p, cpus, skips, timing, marker, error, left = run_case("fail", True, ("components", COMMANDS // 10), 0.2, fail_on_cpu=40)
    import subprocess as sp
    return all([
        check("CalledProcessError raised", isinstance(error, sp.CalledProcessError), repr(error)),
        check("launches 20 then 40", cpus == [20, 40], str(cpus)),
        check("run_trinity_phase2 marker left", marker),
        check("Trans-ABySS lane still finished", "diamond_transabyss\t" in timing),
        check("no sessions left", not left, str(left)),
    ])


def case_threaded():
    """Threaded: share 0.25; Trans-ABySS the long pole, so phase 2 never grows."""
    p, cpus, skips, timing, marker, error, left = run_case("threaded", False, 3.0, 0.2)
    return all([
        check("no error", error is None, repr(error)),
        check("share 0.25 threaded: 10 cpu", p.budgets.get("transabyss") == 10, str(p.budgets)),
        check("one launch, on 30", cpus == [30], str(cpus)),
        check("no grow entry", "run_trinity_phase2_grow" not in timing),
    ])


def case_sigterm():
    """SIGTERM to oyster.py alone while phase 2 runs: phase 2's session dies too."""
    work = Path(tempfile.mkdtemp(prefix="stub_sigterm_"))
    pidfile = work / "phase2.pid"
    with open(work / "child.log", "w") as log:
        child = subprocess.Popen([sys.executable, __file__, "--sigterm-child", str(work), str(pidfile)],
                                 stdout=log, stderr=subprocess.STDOUT)
    deadline = time.time() + 20
    while not pidfile.exists() and time.time() < deadline:
        time.sleep(0.1)
    if not pidfile.exists():
        child.kill()
        return check("phase 2 started", False)
    pgid = int(pidfile.read_text())
    time.sleep(0.5)
    child.send_signal(signal.SIGTERM)
    rc = child.wait(timeout=60)
    time.sleep(0.5)
    ok = gone(pgid)
    if not ok:
        os.killpg(pgid, signal.SIGKILL)
    return all([
        check("oyster.py died of SIGTERM", rc == -signal.SIGTERM, f"rc {rc}"),
        check("phase 2 session (and grandchild) gone", ok, f"pgid {pgid}"),
    ])


def sigterm_child(work, pidfile):
    import oyster
    oyster.forward_signals_to_sessions()
    p = make_pipeline(Path(work), True, 1000, 0.1, Path(work) / "launches.txt")
    start = p._start_trinity_phase2

    def write_pid(cpu, mem):
        proc = start(cpu, mem)
        Path(pidfile).write_text(str(proc.pid))
        return proc
    p._start_trinity_phase2 = write_pid
    os.environ["STUB_UNIT"] = "1.0"  # 10 s of phase 2: still running at the signal
    p.run_assemblers()


CASES = {"early": case_early, "grow": case_grow, "late": case_late, "fail": case_fail,
         "threaded": case_threaded, "sigterm": case_sigterm}


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--fake-trinity":
        outdir, cpu, fasta, launches, fail_on_cpu = sys.argv[2:7]
        fake_trinity(outdir, int(cpu), fasta, launches, int(fail_on_cpu))
        return
    if len(sys.argv) > 1 and sys.argv[1] == "--sigterm-child":
        sigterm_child(sys.argv[2], sys.argv[3])
        return
    names = sys.argv[1:] or list(CASES)
    unknown = [n for n in names if n not in CASES]
    if unknown:
        sys.exit(f"unknown case(s) {' '.join(unknown)}; cases: {' '.join(CASES)}\n{__doc__}")
    failed = []
    for name in names:
        print(f"{name}: {CASES[name].__doc__}")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ok = CASES[name]()
        out = buf.getvalue()
        checks = [l for l in out.splitlines() if l.startswith("  ok") or l.startswith("  FAIL")]
        print("\n".join(checks))
        if not ok:
            failed.append(name)
            print("  --- pipeline output ---\n" + "\n".join(
                "  | " + l for l in out.splitlines() if l not in checks))
    print(f"\n{len(names) - len(failed)}/{len(names)} cases passed" + (f"; failed: {' '.join(failed)}" if failed else ""))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
