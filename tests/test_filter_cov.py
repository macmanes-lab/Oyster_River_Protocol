#!/usr/bin/env python3
"""Test oyster.py's filter_tpm: --tpm-filt, --min-cov (a number or auto), and
what it reports; filter_settings_stamp; parse_min_cov; and the --min-cov
default.

usage: python tests/test_filter_cov.py

Runs filter_tpm on a hand-made quant.sf, with salmon's mean fragment length
set to 200, so coverage = NumReads x 200 / Length:

  contig  length  TPM   reads  cov    hit
  big     2000    10    100    10     no
  mid     1000     0.5    5     1     no
  short    200     8      2     2     no    high TPM from 2 reads
  weak    2000     0.2    4     0.4   yes
  zero    1000     0      0     0     no

Exits non-zero on a failure.
"""

import gzip
import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import argparse  # noqa: E402

from oyster import Pipeline, min_cov_auto, parse_min_cov  # noqa: E402

QUANT = [
    ("big", 2000, 10, 100),
    ("mid", 1000, 0.5, 5),
    ("short", 200, 8, 2),
    ("weak", 2000, 0.2, 4),
    ("zero", 1000, 0, 0),
]

failures = []


def check(cond, msg):
    if not cond:
        failures.append(msg)
        print(f"FAIL: {msg}")


def run(tmp, tpm_filt, min_cov, frag_len=200.0, pairs=None):
    """filter_tpm in a fresh run tree; returns (high, low, candidates)."""
    p = object.__new__(Pipeline)
    p.dir = Path(tmp)
    p.runout = "t"
    p.quants_dir = p.dir / "quants"
    p.assemblies_dir = p.dir / "assemblies"
    p.assemblies_working = p.assemblies_dir / "working"
    p.reports_dir = p.dir / "reports"
    p.tpm_filt, p.min_cov = tpm_filt, min_cov
    qd = p.quants_dir / "salmon_shucked_t"
    (qd / "aux_info").mkdir(parents=True, exist_ok=True)
    for d in (p.assemblies_working, p.reports_dir):
        d.mkdir(parents=True, exist_ok=True)
    with open(qd / "quant.sf", "w") as f:
        f.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        for name, length, tpm, reads in QUANT:
            f.write(f"{name}\t{length}\t{length - 150}\t{tpm}\t{reads}\n")
    meta = {"frag_length_mean": frag_len} if frag_len else {}
    if pairs is not None:
        meta["num_processed"] = pairs
    (qd / "aux_info" / "meta_info.json").write_text(json.dumps(meta))
    (p.assemblies_dir / "t.ORP.diamond.txt").write_text("weak\tsp|P1|X\t90\n")
    p.filter_tpm()
    read = lambda name: set((p.assemblies_working / f"t.{name}.txt").read_text().split())
    cand = {}
    for line in (p.reports_dir / "t.filter_candidates.tsv").read_text().splitlines():
        if line.startswith("#") or line.startswith("rule\t"):
            continue
        cols = line.split("\t")
        if cols[1] == "NA":
            cand[cols[0]] = None
            continue
        cand[cols[0]] = [int(c) for c in cols[1:5]] + [float(cols[5])]
    cand["#"] = {line.split("\t")[0][2:]: line.split("\t")[1:]
                 for line in (p.reports_dir / "t.filter_candidates.tsv").read_text().splitlines()
                 if line.startswith("# ")}
    with gzip.open(p.reports_dir / "t.quant.sf.gz", "rt") as f:
        check(len(f.read().splitlines()) == len(QUANT) + 1, "quant.sf.gz is not the whole quant.sf")
    return read("HIGHEXP"), read("LOWEXP"), cand


with tempfile.TemporaryDirectory() as tmp:
    high, low, _ = run(tmp, 0, 0)
    check(not low, f"defaults drop nothing, got LOWEXP {sorted(low)}")

with tempfile.TemporaryDirectory() as tmp:
    high, low, _ = run(tmp, 1, 0)
    check(high == {"big", "short"}, f"--tpm-filt 1 keeps big and short, got {sorted(high)}")

with tempfile.TemporaryDirectory() as tmp:
    high, low, cand = run(tmp, 0, 1)
    check(high == {"big", "mid", "short"}, f"--min-cov 1 keeps big, mid, short, got {sorted(high)}")
    # pass, nohit_pass, hit_rescued, final, frac_reads_in_final
    check(cand["applied:cov>=1"][:4] == [3, 3, 1, 4], f"applied row {cand['applied:cov>=1']}")
    check(cand["cov>=2,reads>=5"][:4] == [1, 1, 1, 2],
          f"cov>=2,reads>=5 keeps big plus the rescued hit, got {cand['cov>=2,reads>=5']}")
    check(abs(cand["tpm>=1"][4] - 106 / 111) < 1e-3, f"tpm>=1 read fraction {cand['tpm>=1'][4]}")

with tempfile.TemporaryDirectory() as tmp:
    high, low, _ = run(tmp, 1, 2)
    check(high == {"big", "short"}, f"both thresholds must pass, got {sorted(high)}")

with tempfile.TemporaryDirectory() as tmp:
    try:
        run(tmp, 0, 1, frag_len=None)
        check(False, "--min-cov without a fragment length should fail")
    except RuntimeError:
        pass
    high, low, cand = run(tmp, 0, 0, frag_len=None)
    check(not low, "no fragment length is fine when --min-cov is 0")
    check(cand["cov>=1"] is None and cand["tpm>=1"] is not None,
          "without a fragment length, coverage rows are NA and the rest are not")

# --min-cov auto: max(1, 2 x (pairs / 60M) ^ 0.85), from salmon's num_processed.
check(min_cov_auto(60_000_000) == 2, f"auto at 60M pairs is {min_cov_auto(60_000_000)}, not 2")
check(min_cov_auto(15_000_000) == 1, f"auto at 15M pairs is {min_cov_auto(15_000_000)}, not the floor of 1")
check(min_cov_auto(120_000_000) == 3.61, f"auto at 120M pairs is {min_cov_auto(120_000_000)}")
check(min_cov_auto(480_093_281) == 11.71, f"auto at tuco's 480M pairs is {min_cov_auto(480_093_281)}")

with tempfile.TemporaryDirectory() as tmp:
    high, low, cand = run(tmp, 0, "auto", pairs=60_000_000)
    check(high == {"big", "short"}, f"auto at 60M (cov 2) keeps big and short, got {sorted(high)}")
    check(cand["#"]["min_cov"] == ["2", "auto"], f"filter_candidates says min_cov {cand['#']['min_cov']}")
    check(cand["#"]["num_processed"] == ["60000000"], f"num_processed {cand['#']['num_processed']}")
    check("applied:cov>=2" in cand, f"no applied:cov>=2 row in {sorted(k for k in cand if k != '#')}")

with tempfile.TemporaryDirectory() as tmp:
    high, low, _ = run(tmp, 0, "auto", pairs=15_000_000)
    check(high == {"big", "mid", "short"}, f"auto at 15M (cov 1) keeps big, mid, short, got {sorted(high)}")

with tempfile.TemporaryDirectory() as tmp:
    high, low, _ = run(tmp, 0, "auto", pairs=480_000_000)
    check(not high, f"auto at 480M (cov 11.7) keeps no no-hit contig here, got {sorted(high)}")

with tempfile.TemporaryDirectory() as tmp:
    try:
        run(tmp, 0, "auto")
        check(False, "--min-cov auto without num_processed should fail")
    except RuntimeError:
        pass
    high, low, cand = run(tmp, 0, 2)
    check(cand["#"]["min_cov"] == ["2"], f"a numeric --min-cov is not marked auto: {cand['#']['min_cov']}")
    check(cand["#"]["num_processed"] == ["NA"], "a missing num_processed is not reported as NA")

check(parse_min_cov("auto") == "auto" and parse_min_cov("AUTO") == "auto", "parse_min_cov('auto')")
check(parse_min_cov("2.5") == 2.5 and parse_min_cov("0") == 0, "parse_min_cov of a number")
for bad in ("-1", "x", "nan", "inf"):
    try:
        parse_min_cov(bad)
        check(False, f"parse_min_cov accepted {bad!r}")
    except argparse.ArgumentTypeError:
        pass

# filter_settings_stamp: a changed threshold must make the filter step stale.
def stamper(tmp, tpm_filt, min_cov):
    p = object.__new__(Pipeline)
    p.runout = "t"
    p.assemblies_working = Path(tmp) / "assemblies" / "working"
    p.tpm_filt, p.min_cov = tpm_filt, min_cov
    return p


with tempfile.TemporaryDirectory() as tmp:
    done = Path(tmp) / "assemblies" / "t.filter.done"
    stamp = stamper(tmp, 0, 2).filter_settings_stamp(done)
    check(stamp.read_text() == "tpm_filt=0 min_cov=2\n", f"new stamp says {stamp.read_text()!r}")
    os.utime(stamp, (1000, 1000))
    stamper(tmp, 0, 2).filter_settings_stamp(done)
    check(stamp.stat().st_mtime == 1000, "unchanged settings rewrote the stamp")
    done.touch()
    os.utime(done, (2000, 2000))
    stamper(tmp, 0, 5).filter_settings_stamp(done)
    check(stamp.read_text() == "tpm_filt=0 min_cov=5\n", f"changed stamp says {stamp.read_text()!r}")
    check(stamp.stat().st_mtime > done.stat().st_mtime,
          "a changed --min-cov left the stamp older than filter.done, so no re-filter")
    # A dev30 run (min_cov=2) resumed under the auto default re-filters.
    os.utime(stamp, (1000, 1000))
    stamper(tmp, 0, "auto").filter_settings_stamp(done)
    check(stamp.read_text() == "tpm_filt=0 min_cov=auto\n", f"auto stamp says {stamp.read_text()!r}")
    check(stamp.stat().st_mtime > done.stat().st_mtime, "switching to auto did not re-filter")

with tempfile.TemporaryDirectory() as tmp:
    # A run filtered before the stamp existed: dated to its filter.done, so
    # resuming it doesn't redo everything after the filter.
    done = Path(tmp) / "assemblies" / "t.filter.done"
    done.parent.mkdir(parents=True)
    done.touch()
    os.utime(done, (3000, 3000))
    stamp = stamper(tmp, 0, 2).filter_settings_stamp(done)
    check(stamp.stat().st_mtime == 3000, "a first stamp on an old run is not dated to filter.done")

# The default is --min-cov auto on both entry points.
import chowder  # noqa: E402
import oyster  # noqa: E402

argv = sys.argv
try:
    sys.argv = ["oyster.py", "--read1", "a.fq", "--read2", "b.fq"]
    check(oyster.parse_args().min_cov == "auto", "oyster.py --min-cov default is not auto")
    sys.argv = ["oyster.py", "--read1", "a.fq", "--read2", "b.fq", "--min-cov", "5"]
    check(oyster.parse_args().min_cov == 5, "oyster.py --min-cov 5 is not 5")
    sys.argv = ["chowder.py", "--read1", "a.fq", "--read2", "b.fq", "--assemblies", "x.fa", "y.fa"]
    check(chowder.parse_args().min_cov == "auto", "chowder.py --min-cov default is not auto")
finally:
    sys.argv = argv

if failures:
    print(f"{len(failures)} failure(s)")
    sys.exit(1)
print("all filter_tpm checks passed")
