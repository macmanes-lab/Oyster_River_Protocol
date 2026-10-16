#!/usr/bin/env python3
"""Where no-hit contigs stop looking like noise, per sequencing depth.

usage: conda run -n orp python depth_curve.py ROOT NAME [NAME ...] [--tsv PREFIX]

For ROOT/<NAME>_cov (subsample.sbatch's tuco_d<NNN>_cov, tuco_cov, or a
sweep.sbatch dataset), bins the intermediate's contigs without a swissprot
hit by coverage (salmon NumReads x mean fragment length / length, as
--min-cov computes it) and gives each bin's coding signal: the share with an
ORF over 149 codons, against the same share for their bases shuffled
(collect.py's ORF columns). The excess of one over the other is the signal;
noise has none. (collect.py's other column, open over >= 80% of the
length, is no use here: short contigs are open end to end by chance, and on
the sweep and tuco quants its excess is highest in the lowest bins.)

Why not just match removed sets across depths: every expressed contig's
coverage grows in proportion to depth, so the threshold that removes the
same contigs at every depth always scales linearly -- that is a TPM cutoff
again. What should set the threshold is the coverage where the signal
appears, and how that point moves with depth is the question:

  floor      the lowest excess of any bin at coverage >= 1: the noise.
             On the sweep and tuco quants it sits at coverage 2-3.
  crossover  the lower edge of the first bin at or above the floor's whose
             excess, and that of the next two bins, is at least twice the
             floor: where the signal clears the noise for good. (Searching
             from coverage 1 instead found 1 on SRR868662, whose floor is
             at 5 and whose bins below it are already over twice that.)
  Bins under MIN_BIN contigs are skipped. Bins under coverage 1 are left
  out of both: their excess is erratic (a few thousand contigs, some of
  them copies of real genes whose reads salmon gave to the other copy).
  top5_excess, the excess of the 5% of no-hit contigs with the most
  coverage, is printed for scale.

With three or more datasets with a known depth, fits
crossover = a x (pairs / 100M) ^ alpha. alpha near 0 says a fixed --min-cov
is right at any depth; near 1, a depth-scaled one (TPM-like); in between, a
rule like a x (pairs / 100M) ^ alpha.

Depth is salmon's num_processed (read pairs given to it): ROOT/<NAME>_cov/
depth.txt, else quants/*/aux_info/meta_info.json, else the
'# num_processed' line of the filter_candidates file.

Prints the per-dataset summary; --tsv PREFIX also writes PREFIX.bins.tsv
and PREFIX.summary.tsv.
"""

import argparse
import gzip
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import collect  # noqa: E402

EDGES = [0, 0.25, 0.5, 1, 1.5, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100, math.inf]
REMOVED_AT = [1, 2, 5, 10, 20]
TOP_SHARE = 0.05
MIN_BIN = 200


def depth(ds, cand):
    f = ds / "depth.txt"
    if f.exists():
        return int(f.read_text().split()[0])
    for meta in ds.glob("quants/*/aux_info/meta_info.json"):
        with open(meta) as fh:
            return int(json.load(fh)["num_processed"])
    for line in open(cand):
        if line.startswith("# num_processed"):
            return int(line.split("\t")[1])
    return None


def coverage(ds, orfs):
    cand = next(ds.glob("reports/*.filter_candidates.tsv"))
    fl = float(next(l.split("\t")[1] for l in open(cand) if l.startswith("# frag_length_mean")))
    cov = {}
    with gzip.open(next(ds.glob("reports/*.quant.sf.gz")), "rt") as f:
        next(f)
        for line in f:
            n, length, _, _, reads = line.rstrip("\n").split("\t")
            if n in orfs:
                cov[n] = float(reads) * fl / float(length)
    return cov, depth(ds, cand)


def signal(names, orfs):
    """(n, long-ORF share, shuffled long-ORF share, open-80% share) of a set of no-hit contigs."""
    if not names:
        return 0, None, None, None
    is_long = sum(1 for n in names if orfs[n][1] > collect.LONG_ORF_CODONS)
    is_long_shuf = sum(1 for n in names if orfs[n][2] > collect.LONG_ORF_CODONS)
    is_open = sum(1 for n in names if orfs[n][1] * 3 >= collect.ORF_SPAN * orfs[n][0])
    k = len(names)
    return k, is_long / k, is_long_shuf / k, is_open / k


def fmt(x, nd=3):
    return "" if x is None else f"{x:.{nd}f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root", type=Path)
    ap.add_argument("names", nargs="+")
    ap.add_argument("--tsv", metavar="PREFIX")
    args = ap.parse_args()

    bin_rows, summary = [], []
    for name in args.names:
        ds = args.root / f"{name}_cov"
        orfs = collect.nohit_orfs(ds)
        if orfs is None:
            sys.exit(f"{ds}: no nohit_orfs.tsv and no intermediate/diamond to build it")
        cov, pairs = coverage(ds, orfs)

        ranked = sorted(cov, key=cov.get, reverse=True)
        _, o, s, _ = signal(ranked[:max(1, int(len(ranked) * TOP_SHARE))], orfs)
        ref = o - s

        excess = []
        for lo, hi in zip(EDGES, EDGES[1:]):
            k, o, s, opn = signal([n for n, c in cov.items() if lo <= c < hi], orfs)
            ex = None if k < MIN_BIN else o - s
            excess.append((lo, ex))
            bin_rows.append([name, pairs, lo, hi, k, fmt(o), fmt(s), fmt(None if o is None else o - s),
                             fmt(opn)])

        upper = [(lo, ex) for lo, ex in excess if lo >= 1]
        floor = min((ex for _, ex in upper if ex is not None), default=None)
        crossover = None
        if floor is not None:
            start = next(i for i, (_, ex) in enumerate(upper) if ex == floor)
            for i, (lo, _) in enumerate(upper[start:], start):
                run = [ex for _, ex in upper[i:i + 3]]
                if len(run) == 3 and all(ex is not None and ex >= 2 * floor for ex in run):
                    crossover = lo
                    break

        removed = [sum(1 for c in cov.values() if c < t) / len(cov) for t in REMOVED_AT]
        summary.append([name, pairs, len(cov), fmt(ref), fmt(floor), crossover]
                       + [fmt(r) for r in removed])

    head = (["dataset", "pairs", "nohit", "top5_excess", "floor", "crossover_cov"]
            + [f"removed_cov{t}" for t in REMOVED_AT])
    print("\t".join(head))
    for row in summary:
        print("\t".join("" if v is None else str(v) for v in row))

    pts = [(math.log(p / 1e8), math.log(c)) for _, p, _, _, _, c, *_ in summary if p and c]
    if len(pts) >= 3:
        mx = sum(x for x, _ in pts) / len(pts)
        my = sum(y for _, y in pts) / len(pts)
        alpha = (sum((x - mx) * (y - my) for x, y in pts)
                 / sum((x - mx) ** 2 for x, _ in pts))
        a = math.exp(my - alpha * mx)
        print(f"\nfit over {len(pts)} datasets: crossover = {a:.2f} x (pairs / 100M) ^ {alpha:.2f}")
    else:
        print(f"\nno fit: {len(pts)} dataset(s) with both a depth and a crossover (need 3)")

    if args.tsv:
        with open(f"{args.tsv}.bins.tsv", "w") as f:
            f.write("dataset\tpairs\tcov_lo\tcov_hi\tn\torf150\torf150_shuf\texcess\topen80\n")
            for row in bin_rows:
                f.write("\t".join(str(v) for v in row) + "\n")
        with open(f"{args.tsv}.summary.tsv", "w") as f:
            f.write("\t".join(head) + "\n")
            for row in summary:
                f.write("\t".join("" if v is None else str(v) for v in row) + "\n")


if __name__ == "__main__":
    main()
