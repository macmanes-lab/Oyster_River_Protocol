#!/usr/bin/env python3
"""One csv row per component assembly scored by components_score.sbatch.

For each line of components/manifest.tsv (SRR, assembler, assembly, reads,
lineage), reads what assembly_eval.py left in components/<SRR>/<assembler>/:

    reports/run_<run>.ORP/**/short_summary*.txt    BUSCO
    reports/pytransrate_<run>/**/assemblies.csv    pyTransRate (every column)
    reports/<run>.flagstat                         proper pairs (strandeval's mapping)

where <run> is <SRR>_<assembler>. The per-assembler unique-gene count is
not recomputed: ORP's own run already counted it (qualreport UNIQUE GENES
<ASSEMBLER>), and collect_lineage410.py carries it as unique_genes_<assembler>,
so the packaging step joins it from lineage410_metrics.csv.

A component still running gets a row with status "not finished" and blank
metrics, so the table always has one row per manifest line (200 for the
first 50 runs, 236 with the nine orp_ additions).

    ./collect_components.py -o components_metrics.csv     # run from lineage_410/
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

BUSCO_RE = re.compile(
    r"C:(?P<complete>[\d.]+)%\[S:(?P<single>[\d.]+)%,"
    r"D:(?P<duplicated>[\d.]+)%\],F:(?P<fragmented>[\d.]+)%,"
    r"M:(?P<missing>[\d.]+)%,n:(?P<n>\d+)"
)
BUSCO_FIELDS = ["complete", "single", "duplicated", "fragmented", "missing", "n"]


def read_busco(d):
    line = None
    if d.is_dir():
        for p in sorted(d.rglob("short_summary*.txt")):
            for text in p.read_text(errors="replace").splitlines():
                if BUSCO_RE.search(text):
                    line = text
    return BUSCO_RE.search(line).groupdict() if line else {}


def read_pytransrate(d):
    hits = sorted(d.rglob("assemblies.csv")) if d.is_dir() else []
    if not hits:
        return {}
    rows = list(csv.DictReader(hits[0].open()))
    row = dict(rows[0]) if rows else {}
    row.pop("assembly", None)
    return row


def read_proper_pairs(p):
    if not p.is_file():
        return ""
    for line in p.read_text(errors="replace").splitlines():
        if "properly paired" in line:
            m = re.search(r"\(([\d.]+)%", line)
            if m:
                return m.group(1)
    return ""


def main():
    runs = Path(os.environ.get("RUNS", "/mnt/home/macmaneslab/macmanes/compare/lineage_410"))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=runs / "components" / "manifest.tsv")
    ap.add_argument("--components", type=Path, default=runs / "components")
    ap.add_argument("-o", "--out", type=Path, default=Path("components_metrics.csv"))
    args = ap.parse_args()

    records, tr_cols = [], []
    for line in args.manifest.read_text().splitlines():
        if not line.strip():
            continue
        srr, asm, fasta, r1, r2, lineage = line.split("\t")
        run = "{}_{}".format(srr, asm)
        rep = args.components / srr / asm / "reports"
        rec = {"srr": srr, "assembler": asm, "lineage": lineage}
        busco = read_busco(rep / "run_{}.ORP".format(run))
        tr = read_pytransrate(rep / "pytransrate_{}".format(run))
        rec.update({"busco_" + k: v for k, v in busco.items()})
        rec.update({"tr_" + k: v for k, v in tr.items()})
        rec["proper_pairs_pct"] = read_proper_pairs(rep / "{}.flagstat".format(run))
        for c in tr:
            if c not in tr_cols:
                tr_cols.append(c)
        done = (rep / "qualreport.{}".format(run)).is_file()
        rec["status"] = "complete" if done and busco and tr else "not finished"
        records.append(rec)

    header = (["srr", "assembler", "lineage", "status"]
              + ["busco_" + f for f in BUSCO_FIELDS] + ["proper_pairs_pct"]
              + ["tr_" + c for c in tr_cols])
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)
    done = sum(r["status"] == "complete" for r in records)
    print("wrote {} rows to {} ({} complete)".format(len(records), args.out, done), file=sys.stderr)


if __name__ == "__main__":
    main()
