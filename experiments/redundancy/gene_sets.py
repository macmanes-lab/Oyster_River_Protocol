#!/usr/bin/env python3
"""Swissprot genes lost and gained by each arm, against the control.

usage: gene_sets.py VALIDATE_DIR SAMPLES_FILE --arms control ARM [ARM ...]
                    [--show N]

qualreport's UNIQUE GENES ORP is a count: distinct swissprot genes hit by
diamond blastx of ORP.intermediate.fasta, a gene being the swissprot name
with the species dropped (ACTB_HUMAN and ACTB_MOUSE are both ACTB) --
extract_gene_ids in oyster.py. A count can rise while genes are lost, so
this compares the sets: for each sample and arm, how many of the control's
genes are missing (lost) and how many are new (gained).

Not every "lost" name is a lost gene. Two arms can keep near-identical
contigs whose best swissprot hit falls on different members of one family
(4CLL5 in one, 4CLL4 in the other), and a weak hit can drop out of diamond's
--top 0.1 window without anything biological changing. --min-bits restricts
both sides to strong hits: "lost" to genes the control hits at best bitscore
>= N and the arm misses entirely, "gained" to the reverse. Strong-hit losses
are the ones to worry about.

Reads assemblies/<run>.ORP.diamond.txt, or the same member of run.tar for
arms that validate_arms.sbatch ran in scratch and archived.
"""
import argparse
import io
import tarfile
from pathlib import Path


def gene_ids(lines):
    """{gene: best bitscore}; the gene key is extract_gene_ids's."""
    ids = {}
    for line in lines:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 2:
            continue
        parts = cols[1].split("|")
        if len(parts) >= 3:
            g = parts[2].split("_")[0]
            bits = float(cols[11]) if len(cols) > 11 else 0.0
            if bits > ids.get(g, -1.0):
                ids[g] = bits
    return ids


def arm_genes(arm_dir):
    arm_dir = Path(arm_dir)
    hits = list(arm_dir.glob("assemblies/*.ORP.diamond.txt"))
    if hits:
        with open(hits[0]) as f:
            return gene_ids(f)
    tar = arm_dir / "run.tar"
    if tar.is_file():
        with tarfile.open(tar) as t:
            for m in t:
                if m.name.endswith(".ORP.diamond.txt") and "/assemblies/" in "/" + m.name.lstrip("./"):
                    with t.extractfile(m) as f:
                        return gene_ids(io.TextIOWrapper(f))
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("validate_dir")
    p.add_argument("samples_file")
    p.add_argument("--arms", nargs="+", required=True)
    p.add_argument("--show", type=int, default=0, help="list up to N lost genes per sample")
    p.add_argument("--min-bits", type=float, default=0,
                   help="count a lost gene only if the control's best hit to it scored "
                        "at least this bitscore")
    args = p.parse_args()

    control, tests = args.arms[0], args.arms[1:]
    runs = [l.strip() for l in open(args.samples_file) if l.strip()]
    totals = {t: [0, 0, 0] for t in tests}
    print(f"{'sample':12s}{'control':>9s}" +
          "".join(f"{t + ' lost':>18s}{'gained':>8s}{'net':>7s}" for t in tests))
    for run in runs:
        base = arm_genes(Path(args.validate_dir) / run / control)
        if base is None:
            continue
        line = f"{run:12s}{len(base):9d}"
        lost_lists = {}
        for t in tests:
            g = arm_genes(Path(args.validate_dir) / run / t)
            if g is None:
                line += f"{'-':>18s}{'-':>8s}{'-':>7s}"
                continue
            lost = {x for x in base.keys() - g.keys() if base[x] >= args.min_bits}
            gained = {x for x in g.keys() - base.keys() if g[x] >= args.min_bits}
            lost_lists[t] = sorted(lost)
            totals[t][0] += len(lost)
            totals[t][1] += len(gained)
            totals[t][2] += 1
            line += f"{len(lost):18d}{len(gained):8d}{len(gained) - len(lost):+7d}"
        print(line)
        if args.show:
            for t, lost in lost_lists.items():
                if lost:
                    print(f"    {t} lost: {', '.join(lost[:args.show])}"
                          + (" ..." if len(lost) > args.show else ""))
    if args.min_bits:
        print(f"\n(lost = genes the {control} hits at bitscore >= {args.min_bits:g} and the arm "
              f"does not hit at all; gained = the reverse, also at >= {args.min_bits:g})")
    for t in tests:
        lost, gained, n = totals[t]
        if n:
            print(f"\n{t}: over {n} samples, {lost} genes lost and {gained} gained "
                  f"against {control} (net {gained - lost:+d})")


if __name__ == "__main__":
    main()
