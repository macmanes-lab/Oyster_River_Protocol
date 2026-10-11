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

Each row leads with the run's TOTAL wall time -- the number users see, and
the one every change here is for -- and what a balanced split would have
saved off it: both lanes finishing together at (TA core-min + lane core-min)
/ (TA + lane cores), against the later of the two as run.

Balanced share = TA core-minutes / (TA + Trinity-lane core-minutes): the split
under which both lanes, each taken to scale linearly with its cores, finish
together. It is independent of node speed (both lanes slow alike), so runs on
slow nodes count. It assumes Trinity's lane is the whole lane for its span;
phase 2 taking all of --cpu (when Trans-ABySS finished before Stage A did)
breaks that, and the row says so. Same-dataset runs have spread +-0.02.

NOTES.md 2026-10-10 has the runs behind POINTS.
"""

import argparse
import math
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
    ("DRR046632 dev24", 18.44, 0.399),
    ("DRR30368 dev24", 21.83, 0.447),
    ("DRR053698 dev24", 33.23, 0.511),
    ("SRR796576 dev24", 3.73, 0.168),
    ("SRR1138707 dev24", 6.24, 0.212),
    ("SRR1139191 dev24", 10.90, 0.217),
    ("SRR868662 dev24", 25.54, 0.323),
    ("SRR802670 dev24", 53.42, 0.497),
]
# With the SRR runs in, no shape fits (rms ~0.05): Trinity's work per GB of
# reads varies 5-fold between datasets while Trans-ABySS's is ~constant.
# NOTES.md 2026-10-11, "read size does not predict it".

#: Shapes tried for share against GB. Square root fit best over POINTS
#: (rms 0.014 vs 0.026 linear, 0.018 log): the share flattens with depth,
#: since normalization keeps more of a deep library's reads for Trinity too.
MODELS = (("linear", lambda gb: gb), ("sqrt", math.sqrt), ("log", math.log))
CLAMP = (0.2, 0.6)

TOTAL = re.compile(r"^TOTAL\s+(\d+):(\d\d):(\d\d)")
PRETTY = re.compile(r"^(\S+)\s+(\d+):(\d\d):(\d\d)\s+\(started ([\d-]+ [\d:]+)\)")
ASSEMBLERS = re.compile(r"=== Assemblers: transabyss \((\d+) cpu[^)]*\) \|\| Trinity \+ SPAdes \((\d+) cpu\)")
SHARE = re.compile(r"Trans-ABySS share [\d.]+ of --cpu: .*R1 ([\d.]+) GB")
STAGE_B = re.compile(r"=== Stage B: run_trinity_phase2 \((\d+) cpu, ([^)]*)\)")


def read_timing(path):
    """{step: (start datetime, seconds)} from either form of the timing log,
    with the run's wall-clock TOTAL in seconds under "TOTAL" once it has one."""
    steps = {}
    for line in open(path):
        line = line.rstrip("\n")
        m = TOTAL.match(line)
        if m:
            steps["TOTAL"] = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3))
            continue
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


def fit(xs, ys):
    """(a, b) of the least-squares line y = a + b*x, or None for one x."""
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
        balanced_min = (ta * ta_min + lane * lane_min) / (ta + lane)
        saving = max(ta_min, lane_min) - balanced_min
        if "TOTAL" in steps:
            total = steps["TOTAL"] / 60
            total_s = f"TOTAL {fmt(total)}, ~{fmt(total - saving)} if balanced (-{fmt(saving)})"
        else:
            total_s = f"TOTAL not yet (still running); balanced would save ~{fmt(saving)}"
        note = ""
        if phase2 and phase2[0] != lane:
            note = f"  ! phase 2 ran on {phase2[0]} cpu ({phase2[1]}): share overstated"
        slower = "TA" if ta_min > lane_min else "Trinity"
        print(f"{timing}\n    {total_s}\n    cores TA/lane {ta}/{lane}  R1 {gb if gb else '?'} GB  rcorrector {rc // 60}m{rc % 60:02d}s\n"
              f"    TA lane {fmt(ta_min)}  Trinity lane {fmt(lane_min)}  "
              f"slower: {slower} by {fmt(abs(ta_min - lane_min))}\n"
              f"    balanced share {share:.3f}{note}")
        if gb:
            rows.append((timing, gb, share))

    points = list(POINTS) if (args.points or not args.runs) else []
    for label, gb, share in rows:
        # A run already in POINTS (same R1 size and share) would count twice.
        known = next((k for k, g, sh in points if abs(g - gb) < 0.01 and abs(sh - share) < 0.0015), None)
        if known:
            print(f"    ({label} is POINTS' {known!r}: counted once)")
        else:
            points.append((label, gb, share))
    print(f"\nfit over {len(points)} point(s):")
    for label, gb, share in points:
        print(f"    {gb:6.2f} GB  {share:.3f}  {label}")
    if len({gb for _, gb, _ in points}) < 2:
        print("    need two or more read sizes to fit")
        return
    ys = [share for _, _, share in points]
    lo, hi = CLAMP
    for name, f in MODELS:
        a, b = fit([f(gb) for _, gb, _ in points], ys)
        resid = [y - (a + b * f(gb)) for (_, gb, _), y in zip(points, ys)]
        rms = math.sqrt(sum(e * e for e in resid) / len(resid))
        at = "  ".join(f"{gb}GB:{min(hi, max(lo, a + b * f(gb))):.2f}" for gb in (3, 7, 15, 20, 30, 50))
        print(f"    {name:6s} share = {a:.3f} + {b:.4f} * {name}(GB)  rms {rms:.3f}  "
              f"max {max(map(abs, resid)):.3f}\n           clamped {lo}-{hi}: {at}")


if __name__ == "__main__":
    main()
