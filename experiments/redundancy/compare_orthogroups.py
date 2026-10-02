#!/usr/bin/env python3
"""Compare two OrthoFinder runs over the same contigs, e.g. diamond vs blastn.

usage: compare_orthogroups.py --groups NAME=Orthogroups.txt [--groups ...]
                              [--pairs SAMPLE.orp.pairs.tsv]

For each run: how many orthogroups, how many are singletons, and how big the
rest are. makeorthout keeps at most one contig per orthogroup, so the group
count is roughly the size of orthomerged.fasta before the rescue and
cd-hit-est steps.

With --pairs (dup_busco_pairs.py output for the same sample's ORP.fasta),
it also reports, for each pair class, the share of duplicated-BUSCO contig
pairs that the run puts in one orthogroup. A pair in one orthogroup is a
pair makeorthout would have cut to one contig. Contigs missing from every
run (under 201 bp, so rescued by posthack rather than searched) are counted
separately.

With --busco (BUSCO's full_table.tsv over the pooled search input, from
busco_pool.sbatch), it reports how well each run's orthogroups line up with
genes, in both directions:

  split genes     BUSCO genes whose contigs fall in more than one
                  orthogroup. makeorthout keeps one contig per orthogroup,
                  so each extra orthogroup is a potential extra copy.
  merged groups   orthogroups holding contigs of two or more BUSCO genes.
                  makeorthout keeps one contig, so all but one of those
                  genes are lost there.
  genes at risk   genes whose contigs sit only in merged groups, so
                  makeorthout could drop the gene altogether (the posthack
                  rescue may bring some back).

Lower inflation should reduce split genes and raise merged groups. The
useful setting is the lowest one that leaves genes at risk near zero.
"""
import argparse
import collections
import csv
import re
import statistics


def read_groups(path):
    group_of, sizes = {}, []
    with open(path) as f:
        for line in f:
            toks = line.split()
            if not toks:
                continue
            sizes.append(len(toks) - 1)
            for member in toks[1:]:
                group_of[member] = toks[0]
    return group_of, sizes


def read_busco(path):
    """{busco_id: {contig, ...}} over Complete, Duplicated and Fragmented rows."""
    genes = collections.defaultdict(set)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3 and cols[1] in ("Complete", "Duplicated", "Fragmented"):
                genes[cols[0]].add(re.sub(r":\d+-\d+$", "", cols[2]))
    return genes


def busco_report(runs, path):
    genes = read_busco(path)
    print(f"\nBUSCO genes against orthogroups ({len(genes)} genes found in the pool, "
          f"{sum(map(len, genes.values()))} contigs)")
    print(f"{'':28s}" + "".join(f"{n:>14s}" for n in runs))
    stats = {}
    for name, (group_of, _) in runs.items():
        gene_groups = {g: {group_of[c] for c in cs if c in group_of}
                       for g, cs in genes.items()}
        genes_in = collections.defaultdict(set)
        for g, gs in gene_groups.items():
            for grp in gs:
                genes_in[grp].add(g)
        merged = {grp for grp, gs in genes_in.items() if len(gs) > 1}
        stats[name] = [
            sum(len(gs) > 1 for gs in gene_groups.values()),
            sum(max(len(gs) - 1, 0) for gs in gene_groups.values()),
            len(merged),
            sum(bool(gs) and gs <= merged for gs in gene_groups.values()),
        ]
    labels = ["split genes", "  extra orthogroups", "merged groups", "genes at risk"]
    for i, label in enumerate(labels):
        print(f"{label:28s}" + "".join(f"{stats[n][i]:>14}" for n in runs))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--groups", action="append", required=True, metavar="NAME=PATH")
    p.add_argument("--pairs")
    p.add_argument("--busco", help="full_table.tsv from BUSCO over the pooled contigs")
    args = p.parse_args()

    runs = {}
    for spec in args.groups:
        name, path = spec.split("=", 1)
        runs[name] = read_groups(path)

    print(f"{'':28s}" + "".join(f"{n:>14s}" for n in runs))
    rows = [
        ("contigs clustered", lambda g, s: sum(s)),
        ("orthogroups", lambda g, s: len(s)),
        ("  singletons", lambda g, s: sum(x == 1 for x in s)),
        ("  size 2-4", lambda g, s: sum(2 <= x <= 4 for x in s)),
        ("  size 5+", lambda g, s: sum(x >= 5 for x in s)),
        ("mean group size", lambda g, s: round(statistics.mean(s), 2)),
    ]
    for label, fn in rows:
        print(f"{label:28s}" + "".join(f"{fn(*runs[n]):>14}" for n in runs))

    if args.busco:
        busco_report(runs, args.busco)
    if not args.pairs:
        return
    pairs = list(csv.DictReader(open(args.pairs), delimiter="\t"))
    print(f"\nduplicated-BUSCO pairs in one orthogroup ({args.pairs})")
    print(f"{'class':28s}{'pairs':>7s}{'unsearched':>11s}" + "".join(f"{n:>14s}" for n in runs))
    by_cls = collections.defaultdict(list)
    for r in pairs:
        by_cls[r["cls"]].append(r)
        by_cls["all"].append(r)
    for cls in ("redundant", "ends", "isoform", "variant", "divergent", "all"):
        rs = by_cls.get(cls, [])
        if not rs:
            continue
        searched = [r for r in rs if all(r[k] in g for g, _ in runs.values()
                                         for k in ("query", "subject"))]
        cells = []
        for g, _ in runs.values():
            together = sum(g[r["query"]] == g[r["subject"]] for r in searched)
            cells.append(f"{together / len(searched):.2f}" if searched else "-")
        print(f"{cls:28s}{len(rs):7d}{len(rs) - len(searched):11d}"
              + "".join(f"{c:>14s}" for c in cells))


if __name__ == "__main__":
    main()
