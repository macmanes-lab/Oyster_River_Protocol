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

    ./collect_metrics.py -o metrics.csv        # paths default off $COMPARE
"""

import argparse
import csv
import os
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


QUALREPORT_RE = {
    "unique_genes_ORP": re.compile(r"UNIQUE GENES ORP\s*~*>\s*(\S+)"),
    "proper_pairs": re.compile(r"READS MAPPED AS PROPER PAIRS\s*~*>\s*(\S+)"),
}


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


def read_from_qualreport(reports, srr):
    """Unique genes and proper pairs as reportgen recorded them.

    Both are read above from the files the pipeline computes them in --
    assemblies/working/<run>.unique.ORP.txt and <run>.flagstat -- and ORP's
    end-of-run cleanup deletes both (oyster.py's cleanup list). So for exactly
    the runs that finished, the direct reads return nothing, and the only
    surviving copy is the text report reportgen wrote before the cleanup ran.
    """
    p = reports / f"qualreport.{srr}"
    if not p.is_file():
        return {}
    text = p.read_text(errors="replace")
    out = {}
    for field, pattern in QUALREPORT_RE.items():
        m = pattern.search(text)
        if m:
            out[field] = m.group(1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    compare = Path(os.environ.get("COMPARE", "/mnt/home/macmaneslab/macmanes/compare"))
    ap.add_argument("--runs", type=Path, default=compare / "orp_runs",
                    help="directory holding one subdirectory per run "
                         "(default: $COMPARE/orp_runs)")
    ap.add_argument("--manifest", type=Path, default=compare / "manifest.tsv",
                    help="manifest.tsv, for the tsa code and accession of each "
                         "run and to report samples that produced no run "
                         "directory at all (default: $COMPARE/manifest.tsv)")
    ap.add_argument("-o", "--out", type=Path, default=Path("metrics.csv"))
    ap.add_argument("--include-partial", action="store_true",
                    help="also write a row for every run that has not produced "
                         "both metric sources yet, with its cells left empty "
                         "(default: list those on stderr and leave them out)")
    args = ap.parse_args()

    meta, order = {}, []
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
            tsa = f[0].strip()
            # The accession comes from the R1 filename rather than from the run
            # name, which carries a _tsa_XXXX suffix when two samples shared an
            # accession -- the SRR itself is what that suffix was added to.
            m = re.search(r"[SED]RR\d{4,}", Path(f[2].strip()).name)
            order.append(run)
            meta[run] = {
                "tsa": tsa,
                # tsa_GADU -> GADU, orp_ABCD -> ABCD. The prefix says which
                # collection a sample came from and is kept in `tsa`; `code` is
                # the bare identifier, for joining against tables that carry it
                # without one.
                "code": re.sub(r"^(tsa|orp)_", "", tsa),
                "srr": m.group(0) if m else "",
            }
    else:
        order = sorted(d.name for d in args.runs.iterdir()
                       if d.is_dir() and d.name != "logs")
        for run in order:
            m = re.search(r"[SED]RR\d{4,}", run)
            meta[run] = {"tsa": "", "code": "", "srr": m.group(0) if m else ""}

    records, skipped, transrate_cols = [], [], []
    for srr in order:
        rundir = args.runs / srr
        info = meta.get(srr, {})
        rec = {"run": srr, "srr": info.get("srr", ""),
               "code": info.get("code", ""), "tsa": info.get("tsa", "")}
        if not rundir.is_dir():
            if args.include_partial:
                rec["status"] = "no run directory"
                records.append(rec)
            else:
                skipped.append((srr, "not started"))
            continue

        reports = rundir / "reports"
        busco, busco_note = read_busco(reports, srr)
        tr, tr_note = read_transrate(reports, srr)
        rec.update(busco)
        rec.update(tr)
        rec["unique_genes_ORP"] = read_unique_genes(rundir, srr)
        rec["proper_pairs"] = read_proper_pairs(rundir, srr)
        # Only where cleanup has taken the sources away; a live run's own files
        # stay authoritative.
        for field, value in read_from_qualreport(reports, srr).items():
            if not rec.get(field):
                rec[field] = value

        for col in tr:
            if col not in transrate_cols:
                transrate_cols.append(col)

        # A row is worth a line in the csv when both metric sources are there.
        # Anything less is a run still working its way up to them, and printing
        # it as a line of empty commas puts a placeholder in the table that has
        # to be filtered out of every later analysis.
        notes = [n for n in (busco_note, tr_note) if n]
        if notes and not args.include_partial:
            skipped.append((srr, "; ".join(notes)))
            continue
        # Past that, status distinguishes a run that reached the end from one
        # that has the metrics but has not finished: the last steps after
        # pytransrate and BUSCO are what write qualreport.<run>.done.
        if notes:
            rec["status"] = "; ".join(notes)
        elif (reports / f"qualreport.{srr}.done").is_file():
            rec["status"] = "complete"
        else:
            rec["status"] = "still running"
        records.append(rec)

    header = (["run", "srr", "code", "tsa", "status"] + BUSCO_FIELDS
              + ["unique_genes_ORP", "proper_pairs"] + transrate_cols)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        for rec in records:
            w.writerow(rec)

    done = sum(1 for r in records if r.get("status") == "complete")
    running = sum(1 for r in records if r.get("status") == "still running")
    parts = [f"{done} finished"]
    if running:
        parts.append(f"{running} with metrics but still running")
    # Only reachable with --include-partial; without it these were skipped.
    empty = len(records) - done - running
    if empty:
        parts.append(f"{empty} with no metrics")
    print(f"wrote {len(records)} rows to {args.out} ({', '.join(parts)})",
          file=sys.stderr)
    if skipped:
        print(f"left out {len(skipped)} of {len(order)} runs, no metrics yet:",
              file=sys.stderr)
        for run, why in skipped:
            print(f"  {run}: {why}", file=sys.stderr)


if __name__ == "__main__":
    main()
