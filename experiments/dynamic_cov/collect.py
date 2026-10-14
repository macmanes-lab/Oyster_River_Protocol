#!/usr/bin/env python3
"""One table of the filter sweep: every dataset x level that sweep.sbatch finished.

usage: collect.py ROOT [--tsv OUT]

ROOT holds one directory per dataset (sweep.sbatch's $DIR), each with
levels/<level>/. Prints, per dataset and level: contigs and bases in the
final assembly, the share of contigs with no swissprot hit kept by the
filter, BUSCO, unique genes, pyTransRate's score, optimal score and its
mapping metrics, and proper-pair mapping.
"""

import argparse
import csv
import re
from pathlib import Path

LEVEL_ORDER = ["none", "tpm1", "cov0.5", "cov1", "cov2", "cov5", "cov10"]

QUAL = {
    "busco": r"BUSCO SCORE ~+>\s+(\S+)",
    "score": r"PYTRANSRATE SCORE ~+>\s+(\S+)",
    "optimal": r"PYTRANSRATE OPTIMAL SCORE ~+>\s+(\S+)",
    "unique_genes": r"UNIQUE GENES ORP ~+>\s+(\S+)",
    "proper_pairs": r"READS MAPPED AS PROPER PAIRS ~+>\s+(\S+)",
}
PYTR = ["n_seqs", "n_bases", "n50", "p_fragments_mapped", "p_good_mapping",
        "p_contigs_uncovered", "p_contigs_lowcovered", "p_contigs_segmented"]
COLUMNS = (["dataset", "level", "nohit_kept", "C", "S", "D", "F", "M", "unique_genes",
            "score", "optimal"] + PYTR + ["proper_pairs"])


def level_key(name):
    return (LEVEL_ORDER.index(name) if name in LEVEL_ORDER else len(LEVEL_ORDER), name)


def read_level(d):
    row = {}
    qual = next(d.glob("qualreport.*"), None)
    if qual:
        text = qual.read_text()
        for key, pat in QUAL.items():
            m = re.search(pat, text)
            row[key] = m.group(1) if m else ""
        m = re.match(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", row.pop("busco"))
        if m:
            row.update(zip("CSDFM", m.groups()))
    pytr = d / "pytransrate.assemblies.csv"
    if pytr.exists():
        with open(pytr) as f:
            rec = next(csv.DictReader(f))
        row.update({k: rec.get(k, "") for k in PYTR})
    cand = next(d.glob("*.filter_candidates.tsv"), None)
    if cand:
        text = cand.read_text()
        total = re.search(r"^# contigs\t(\d+)", text, re.M)
        hits = re.search(r"^# with_swissprot_hit\t(\d+)", text, re.M)
        applied = re.search(r"^applied:\S+\t\d+\t(\d+)\t", text, re.M)
        if total and hits and applied:
            nohit = int(total.group(1)) - int(hits.group(1))
            row["nohit_kept"] = f"{applied.group(1)}/{nohit}"
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root", type=Path)
    ap.add_argument("--tsv", type=Path)
    args = ap.parse_args()
    rows = []
    for ds in sorted(p for p in args.root.iterdir() if (p / "levels").is_dir()):
        for lv in sorted((p for p in (ds / "levels").iterdir() if p.is_dir()), key=lambda p: level_key(p.name)):
            rows.append({"dataset": ds.name, "level": lv.name, **read_level(lv)})
    lines = ["\t".join(COLUMNS)] + ["\t".join(str(r.get(c, "")) for c in COLUMNS) for r in rows]
    if args.tsv:
        args.tsv.write_text("\n".join(lines) + "\n")
    widths = [max(len(l.split("\t")[i]) for l in lines) for i in range(len(COLUMNS))]
    for l in lines:
        print("  ".join(c.ljust(w) for c, w in zip(l.split("\t"), widths)))


if __name__ == "__main__":
    main()
