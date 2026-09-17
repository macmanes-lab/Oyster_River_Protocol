#!/usr/bin/env python3
"""Collect pytransrate and BUSCO metrics from the 108 ORP run directories.

Each run made by orp_array.sbatch lives in <runs>/<SRR>/ with the shape oyster.py
gives it, so everything here is read from where the pipeline actually writes it:

    reports/transrate_<SRR>/**/assemblies.csv   pytransrate, final ORP assembly
    reports/run_<SRR>.ORP/**/short*.txt         BUSCO
    assemblies/working/<SRR>.unique.ORP.txt     unique SwissProt genes
    assemblies/<SRR>.flagstat                   reads mapped as proper pairs

Every pytransrate column is carried through under its own header name rather
than by position, so this does not have to be touched when that csv changes.
Runs that are missing a file get a blank cell and a note in `status`, so an
incomplete sample shows up as a gap in the table instead of dropping out of it.

    ./collect_metrics.py --runs /path/to/orp_runs -o metrics.csv
"""

import argparse
import csv
import re
import sys
from pathlib import Path

BUSCO_RE = re.compile(
    r"C:(?P<busco_complete>[\d.]+)%\[S:(?P<busco_single>[\d.]+)%,"
    r"D:(?P<busco_duplicated>[\d.]+)%\],F:(?P<busco_fragmented>[\d.]+)%,"
    r"M:(?P<busco_missing>[\d.]+)%,n:(?P<busco_n>\d+)"
)
BUSCO_FIELDS = [
    "busco_complete", "busco_single", "busco_duplicated",
    "busco_fragmented", "busco_missing", "busco_n",
]


def first(paths):
    """The first match in sorted order, or None. rglob's order is arbitrary."""
    return next(iter(sorted(paths)), None)


def read_busco(reports, srr):
    """Parse the C:/S:/D:/F:/M:/n: line out of BUSCO's short summary.

    Matches reportgen's own reading of it (oyster.py qualreport): last such line
    in the file wins, because the summary repeats it inside the plain-text block.
    """
    busco_dir = reports / f"run_{srr}.ORP"
    if not busco_dir.is_dir():
        return {}, "no busco dir"
    line = None
    for p in sorted(busco_dir.rglob("short*.txt")):
        for text in p.read_text(errors="replace").splitlines():
            if BUSCO_RE.search(text):
                line = text
    if line is None:
        return {}, "no busco score line"
    return BUSCO_RE.search(line).groupdict(), ""


def read_transrate(reports, srr):
    csv_path = first((reports / f"transrate_{srr}").rglob("assemblies.csv")) \
        if (reports / f"transrate_{srr}").is_dir() else None
    if csv_path is None:
        return {}, "no transrate csv"
    rows = list(csv.DictReader(csv_path.open()))
    if not rows:
        return {}, "empty transrate csv"
    row = dict(rows[0])
    # The full path of the assembly that was scored: same information as the SRR
    # column and long enough to make the table unreadable.
    row.pop("assembly", None)
    return row, ""


def read_unique_genes(rundir, srr):
    p = rundir / "assemblies" / "working" / f"{srr}.unique.ORP.txt"
    return p.read_text().strip() if p.is_file() else ""


def read_proper_pairs(rundir, srr):
    p = rundir / "assemblies" / f"{srr}.flagstat"
    if not p.is_file():
        return ""
    for line in p.read_text(errors="replace").splitlines():
        if "properly paired" in line:
            fields = line.split()
            if len(fields) > 5:
                return fields[5].lstrip("(")
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True, type=Path,
                    help="directory holding one subdirectory per run")
    ap.add_argument("--manifest", type=Path, default=None,
                    help="manifest.tsv, to carry the tsa code across and to "
                         "report samples that produced no run directory at all")
    ap.add_argument("-o", "--out", type=Path, default=Path("metrics.csv"))
    args = ap.parse_args()

    tsa_of, order = {}, []
    if args.manifest and args.manifest.is_file():
        # Manifest columns: tsa_code, TSA assembly, R1, R2, run name. Column 5
        # is what the job passed to --runout and so what the run directory is
        # called; taking it verbatim rather than re-deriving it from the R1
        # filename is what makes the two agree when an accession was shared by
        # two samples and had to be disambiguated.
        for line in args.manifest.read_text().splitlines():
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 5 or not f[4].strip():
                print(f"manifest line has no run name in column 5: {line[:60]}...",
                      file=sys.stderr)
                continue
            run = f[4].strip()
            order.append(run)
            tsa_of[run] = f[0].strip()
    else:
        order = sorted(d.name for d in args.runs.iterdir()
                       if d.is_dir() and d.name != "logs")

    records, transrate_cols = [], []
    for srr in order:
        rundir = args.runs / srr
        rec = {"run": srr, "tsa": tsa_of.get(srr, "")}
        if not rundir.is_dir():
            rec["status"] = "no run directory"
            records.append(rec)
            continue

        reports = rundir / "reports"
        busco, busco_note = read_busco(reports, srr)
        tr, tr_note = read_transrate(reports, srr)
        rec.update(busco)
        rec.update(tr)
        rec["unique_genes_ORP"] = read_unique_genes(rundir, srr)
        rec["proper_pairs"] = read_proper_pairs(rundir, srr)

        for col in tr:
            if col not in transrate_cols:
                transrate_cols.append(col)

        notes = [n for n in (busco_note, tr_note) if n]
        if (reports / f"qualreport.{srr}.done").is_file() and not notes:
            rec["status"] = "complete"
        else:
            rec["status"] = "; ".join(notes) or "incomplete"
        records.append(rec)

    header = (["run", "tsa", "status"] + BUSCO_FIELDS
              + ["unique_genes_ORP", "proper_pairs"] + transrate_cols)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        for rec in records:
            w.writerow(rec)

    done = sum(1 for r in records if r.get("status") == "complete")
    print(f"{done}/{len(records)} complete -> {args.out}", file=sys.stderr)
    for r in records:
        if r.get("status") != "complete":
            print(f"  {r['run']}: {r.get('status')}", file=sys.stderr)


if __name__ == "__main__":
    main()
