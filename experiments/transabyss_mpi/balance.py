#!/usr/bin/env python3
"""balance.py - the Trans-ABySS share at which a run's two assembler lanes
would have finished together, and the line fitted through those shares.

Usage:
    balance.py RUN [RUN ...]            one row per run, then the fit
    balance.py --points                 the fit over POINTS below only

RUN is TIMING_LOG[:SLURM_LOG], e.g.
    SRR1138704_dev24/reports/SRR1138704_dev24.timing.log:orp_ta_mpi_1321504.log

The timing log gives each lane's span: run_transabyss's start to the end of
diamond_transabyss, and to the end of run_trinity_phase2. The slurm log's
`=== Assemblers: transabyss (N cpu, mpi) || Trinity + SPAdes (M cpu)` and
`[assemblers] Trans-ABySS share ... R1 X GB` lines give the cores and the
corrected R1's size. Without a slurm log, pass --ta-cpu, --lane-cpu, --gb.

Balanced share = TA core-minutes / (TA + Trinity-lane core-minutes): the split
under which both lanes, each taken to scale linearly with its cores, finish
together. It is independent of node speed (both lanes slow alike), so runs on
slow nodes count. It assumes Trinity's lane is the whole lane for its span;
phase 2 taking all of --cpu (when Trans-ABySS finished before Stage A did)
breaks that, and the row says so. Same-dataset runs have spread +-0.02.

NOTES.md 2026-10-10 has the runs behind POINTS.
"""

import argparse
import re
import sys
from datetime import datetime, timedelta

#: (label, uncompressed corrected R1 GB, balanced share), --cpu 40,
#: --normalize-reads, MPI without OMPI_MCA_mpi_yield_when_idle (dev22, dev24).
#: dev23's runs had that setting, which slowed Trans-ABySS; they are left out.
POINTS = [
    ("SRR1138704 dev24", 3.64, 0.188),
    ("SRR866209 dev24", 7.35, 0.255),
    ("DRR031870 dev24", 17.47, 0.412),
    ("DRR031870 dev22", 17.47, 0.370),
]

PRETTY = re.compile(r"^(\S+)\s+(\d+):(\d\d):(\d\d)\s+\(started ([\d-]+ [\d:]+)\)")
ASSEMBLERS = re.compile(r"=== Assemblers: transabyss \((\d+) cpu[^)]*\) \|\| Trinity \+ SPAdes \((\d+) cpu\)")
SHARE = re.compile(r"Trans-ABySS share [\d.]+ of --cpu: .*R1 ([\d.]+) GB")
STAGE_B = re.compile(r"=== Stage B: run_trinity_phase2 \((\d+) cpu, ([^)]*)\)")


def read_timing(path):
    """{step: (start datetime, seconds)} from either form of the timing log."""
    steps = {}
    for line in open(path):
        line = line.rstrip("\n")
        if "\t" in line:
            parts = line.split("\t")
            if len(parts) >= 3 and parts[1].isdigit():
                steps[parts[0]] = (datetime.strptime(parts[2], "%Y-%m-%d %H:%M:%S"), int(parts[1]))
            continue
        m = PRETTY.match(line)
        if m:
            secs = int(m.group(2)) * 3600 + int(m.group(3)) * 60 + int(m.group(4))
            steps[m.group(1)] = (datetime.strptime(m.group(5), "%Y-%m-%d %H:%M:%S"), secs)
    return steps


def read_slurm(path):
    ta = lane = gb = None
    phase2 = None
    for line in open(path):
        m = ASSEMBLERS.search(line)
        if m:
            ta, lane = int(m.group(1)), int(m.group(2))
        m = SHARE.search(line)
        if m:
            gb = float(m.group(1))
        m = STAGE_B.search(line)
        if m:
            phase2 = (int(m.group(1)), m.group(2))
    return ta, lane, gb, phase2


def end(steps, name):
    start, secs = steps[name]
    return start + timedelta(seconds=secs)


def fmt(minutes):
    return f"{int(minutes // 60)}h{int(minutes % 60):02d}m"


def fit(points):
    xs = [p[1] for p in points]
    ys = [p[2] for p in points]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return my - b * mx, b


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*", help="TIMING_LOG[:SLURM_LOG]")
    ap.add_argument("--ta-cpu", type=int)
    ap.add_argument("--lane-cpu", type=int)
    ap.add_argument("--gb", type=float, help="uncompressed corrected R1, GB")
    ap.add_argument("--points", action="store_true", help="include POINTS in the fit (default when no RUN)")
    args = ap.parse_args()

    rows = []
    for spec in args.runs:
        timing, _, slurm = spec.partition(":")
        ta, lane, gb, phase2 = read_slurm(slurm) if slurm else (None, None, None, None)
        ta, lane, gb = ta or args.ta_cpu, lane or args.lane_cpu, gb or args.gb
        if not (ta and lane):
            sys.exit(f"{spec}: no core counts; give its slurm log or --ta-cpu/--lane-cpu")
        steps = read_timing(timing)
        missing = [s for s in ("run_transabyss", "diamond_transabyss", "run_trinity_phase2") if s not in steps]
        if missing:
            sys.exit(f"{timing}: no {', '.join(missing)} yet")
        t0 = steps["run_transabyss"][0]
        ta_min = (end(steps, "diamond_transabyss") - t0).total_seconds() / 60
        lane_min = (end(steps, "run_trinity_phase2") - t0).total_seconds() / 60
        share = ta * ta_min / (ta * ta_min + lane * lane_min)
        rc = steps.get("run_rcorrector", (None, 0))[1]
        note = ""
        if phase2 and phase2[0] != lane:
            note = f"  ! phase 2 ran on {phase2[0]} cpu ({phase2[1]}): share overstated"
        slower = "TA" if ta_min > lane_min else "Trinity"
        print(f"{timing}\n    cores TA/lane {ta}/{lane}  R1 {gb if gb else '?'} GB  rcorrector {rc // 60}m{rc % 60:02d}s\n"
              f"    TA lane {fmt(ta_min)}  Trinity lane {fmt(lane_min)}  "
              f"slower: {slower} by {fmt(abs(ta_min - lane_min))}\n"
              f"    balanced share {share:.3f}{note}")
        if gb:
            rows.append((timing, gb, share))

    points = (POINTS if (args.points or not args.runs) else []) + rows
    print(f"\nfit over {len(points)} point(s):")
    for label, gb, share in points:
        print(f"    {gb:6.2f} GB  {share:.3f}  {label}")
    line = fit(points)
    if line:
        a, b = line
        print(f"    share = {a:.3f} + {b:.4f} * GB   (oyster.py: TRANSABYSS_MPI_SHARE_BASE, _PER_GB)")
        for gb in (2, 5, 10, 15, 20, 25, 30):
            s = a + b * gb
            print(f"      {gb:3d} GB -> {s:.3f}, {min(0.45, max(0.2, s)):.3f} clamped, "
                  f"{round(40 * min(0.45, max(0.2, s)))} of 40 cpu")
    else:
        print("    need two or more read sizes to fit a line")


if __name__ == "__main__":
    main()
