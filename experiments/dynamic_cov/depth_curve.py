#!/usr/bin/env python3
"""Where no-hit contigs stop looking like noise, per sequencing depth.

usage: conda run -n orp python depth_curve.py ROOT NAME [NAME ...] [--tsv PREFIX]

For ROOT/<NAME>_cov (subsample.sbatch's tuco_d<NNN>_cov, tuco_cov, or a
sweep.sbatch dataset), bins the intermediate's contigs without a swissprot
hit by coverage (salmon NumReads x mean fragment length / length, as
--min-cov computes it) and gives each bin's coding signal: the share open
over >= 80% of their length, against the same share for their bases
shuffled (collect.py's ORF columns). The excess of one over the other is
the signal; noise has none.

Why not just match removed sets across depths: every expressed contig's
coverage grows in proportion to depth, so the threshold that removes the
same contigs at every depth always scales linearly -- that is a TPM cutoff
again. What should set the threshold is the coverage where the signal
appears, and how that point moves with depth is the question:

  crossover  the lower edge of the first bin, from the bottom, whose excess
             reaches half the reference excess and stays there in the next
             bin. The reference is the excess of the 5% of no-hit contigs
             with the most coverage, which is at every depth the same
             mostly-real set.

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
    """(n, open share, shuffled open share, long-ORF share) of a set of no-hit contigs."""
    if not names:
        return 0, None, None, None
    is_open = sum(1 for n in names if orfs[n][1] * 3 >= collect.ORF_SPAN * orfs[n][0])
    is_open_shuf = sum(1 for n in names if orfs[n][2] * 3 >= collect.ORF_SPAN * orfs[n][0])
    is_long = sum(1 for n in names if orfs[n][1] > collect.LONG_ORF_CODONS)
    k = len(names)
    return k, is_open / k, is_open_shuf / k, is_long / k


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
            k, o, s, lng = signal([n for n, c in cov.items() if lo <= c < hi], orfs)
            ex = None if o is None else o - s
            excess.append((lo, ex))
            bin_rows.append([name, pairs, lo, hi, k, fmt(o), fmt(s), fmt(ex), fmt(lng)])

        crossover = None
        for (lo, ex), (_, nxt) in zip(excess, excess[1:]):
            if ex is not None and nxt is not None and ex >= ref / 2 and nxt >= ref / 2:
                crossover = lo
                break

        removed = [sum(1 for c in cov.values() if c < t) / len(cov) for t in REMOVED_AT]
        summary.append([name, pairs, len(cov), fmt(ref), crossover] + [fmt(r) for r in removed])

    head = (["dataset", "pairs", "nohit", "top5_excess", "crossover_cov"]
            + [f"removed_cov{t}" for t in REMOVED_AT])
    print("\t".join(head))
    for row in summary:
        print("\t".join("" if v is None else str(v) for v in row))

    pts = [(math.log(p / 1e8), math.log(c)) for _, p, _, _, c, *_ in summary if p and c]
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
            f.write("dataset\tpairs\tcov_lo\tcov_hi\tn\topen80\topen80_shuf\texcess\torf150\n")
            for row in bin_rows:
                f.write("\t".join(str(v) for v in row) + "\n")
        with open(f"{args.tsv}.summary.tsv", "w") as f:
            f.write("\t".join(head) + "\n")
            for row in summary:
                f.write("\t".join("" if v is None else str(v) for v in row) + "\n")


if __name__ == "__main__":
    main()
