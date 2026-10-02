#!/usr/bin/env python3
"""Arms side by side, per sample, for validate_arms.sbatch.

usage: validate_summary.py VALIDATE_DIR SAMPLES_FILE [--arms control candidate ...]
                           [--jobs JOBID ...]

One row per sample and arm (final contigs, BUSCO in genes, transrate). Then,
for every arm after the first (the control), the per-sample difference from
the control: median, range and total over the samples where both are
finished, plus how many samples moved each way. A sample counts as finished
for an arm once both BUSCO and transrate have reported.

With --jobs, adds each run's elapsed time and MaxRSS from sacct, matched to
its sample and arm by the first line of its log (validate_arms_<job>_<task>.log),
which names both.
"""
import argparse
import re
import statistics
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_arms import arm_row  # noqa: E402

COLS = [("final", "contigs"), ("BUSCO C", "C"), ("BUSCO S", "S"), ("BUSCO D", "D"),
        ("BUSCO F", "F"), ("BUSCO M", "M"), ("transrate score", "score"),
        ("transrate optimal_score", "optimal"), ("transrate p_good_mapping", "good"),
        ("unique genes", "uniq")]
N_BUSCO = 125


def value(row, key):
    v = row.get(key)
    if v is None:
        return None
    v = float(v)
    return round(v * N_BUSCO / 100) if key.startswith("BUSCO") else v


def sacct(jobids, log_dir):
    """{(run, arm): (elapsed, max GB)}, from sacct plus each task's log header."""
    out = {}
    for job in jobids:
        # stdout=PIPE rather than capture_output: Premise's system python3 is 3.6.
        res = subprocess.run(["sacct", "-j", job, "-n", "-P", "-o", "JobID,Elapsed,MaxRSS"],
                             stdout=subprocess.PIPE, universal_newlines=True).stdout
        per_task = {}
        for line in res.splitlines():
            jid, elapsed, rss = line.split("|")
            base = jid.split(".")[0]
            if "_" not in base:
                continue
            el, mx = per_task.get(base, ("", 0.0))
            if "." not in jid:
                el = elapsed
            if rss:
                mult = {"K": 1 / 2**20, "M": 1 / 2**10, "G": 1.0}.get(rss[-1], 1 / 2**30)
                mx = max(mx, float(rss.rstrip("KMG")) * mult)
            per_task[base] = (el, mx)
        for base, use in per_task.items():
            log = Path(log_dir) / f"validate_arms_{base}.log"
            if not log.is_file():
                continue
            m = re.match(r"(\S+) (\S+) on ", log.open().readline())
            if m:
                out[(m.group(1), m.group(2))] = use
    return out


def fmt(v):
    if v is None:
        return "-"
    if isinstance(v, float) and 0 < abs(v) < 1:
        return f"{v:.3f}"
    return f"{v:,.0f}"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("validate_dir")
    p.add_argument("samples_file")
    p.add_argument("--arms", nargs="+", default=["control", "candidate"])
    p.add_argument("--jobs", nargs="*", default=[])
    p.add_argument("--log-dir", default=".")
    args = p.parse_args()

    runs = [l.strip() for l in open(args.samples_file) if l.strip()]
    usage = sacct(args.jobs, args.log_dir) if args.jobs else {}
    w = max(len(a) for a in args.arms) + 2
    head = f"{'sample':12s}{'arm':{w}s}" + "".join(f"{c:>9s}" for _, c in COLS)
    if usage:
        head += f"{'elapsed':>10s}{'maxGB':>7s}"
    print(head)

    control, tests = args.arms[0], args.arms[1:]
    deltas = {t: {c: [] for _, c in COLS} for t in tests}
    for run in runs:
        rows = {}
        for arm in args.arms:
            r = arm_row(Path(args.validate_dir) / run / arm)
            # Finished means both reports are in; ORP.fasta alone appears
            # before BUSCO and transrate have run.
            if r and not ("BUSCO C" in r and "transrate score" in r):
                r = None
            rows[arm] = r
            line = f"{run:12s}{arm:{w}s}"
            line += "".join(f"{fmt(value(r, key)) if r else 'running':>9s}" for key, _ in COLS)
            if usage:
                el, mx = usage.get((run, arm), ("", 0))
                line += f"{el:>10s}{mx:7.0f}"
            print(line)
        for t in tests:
            if rows[control] and rows[t]:
                for key, c in COLS:
                    a, b = value(rows[control], key), value(rows[t], key)
                    if a is not None and b is not None:
                        deltas[t][c].append(b - a)

    for t in tests:
        d = deltas[t]
        n = len(d["C"])
        if not n:
            continue
        print(f"\n{t} - {control} over {n} samples (BUSCO in genes)")
        print(f"{'':22s}" + "".join(f"{c:>9s}" for _, c in COLS))
        for label, fn in (("median", statistics.median), ("min", min), ("max", max),
                          ("total", sum)):
            print(f"{label:22s}" + "".join(f"{fmt(fn(d[c])) if d[c] else '-':>9s}"
                                           for _, c in COLS))
        moved = lambda c, sign: sum((x * sign) > 0 for x in d[c])
        print(f"  duplicated: {moved('D', -1)} down, {moved('D', 1)} up; missing: "
              f"{moved('M', -1)} down, {moved('M', 1)} up; unique genes: "
              f"{moved('uniq', 1)} up, {moved('uniq', -1)} down; transrate score: "
              f"{moved('score', 1)} up, {moved('score', -1)} down")


if __name__ == "__main__":
    main()
