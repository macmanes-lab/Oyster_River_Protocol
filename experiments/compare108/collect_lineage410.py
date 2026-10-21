#!/usr/bin/env python3
"""One csv row per lineage sample: the ORP run against the original assembly.

Run after orp_lineage410.sbatch finishes. Reads the 7-column manifest that
build_lineage410_manifest.sh wrote (pair, original assembly, R1, R2, SRR,
group, lineage) and, for each sample, everything ORP and the original-assembly
scoring left on disk:

    <runs>/<SRR>/reports/qualreport.<runout>        unique genes (ORP and each
                                                     assembler), proper pairs
    <runs>/<SRR>/reports/run_<runout>.ORP/          BUSCO, ORP assembly
    <runs>/<SRR>/reports/pytransrate_<runout>/      pyTransRate, ORP assembly
    <runs>/<SRR>/reports/<runout>.timing.log        command line, step times, TOTAL
    <runs>/<SRR>/reports/<runout>.filter_candidates.tsv
                                                     min-cov applied, contigs before
                                                     the filter, salmon read counts
    <runs>/logs/<SRR>.log                           read pairs in/after trimming,
                                                     read length, spades k, slurm job
    $COMPARE/busco_original/<pair>/short_summary*   BUSCO, original assembly
    $COMPARE/pytransrate/results/<pair>/assemblies.csv
                                                     pyTransRate, original assembly
    $COMPARE/orig_unique_genes/<pair>.unique.txt    unique genes, original assembly
                                                     (orig_unique_genes.sbatch)
    $COMPARE/pairs/<pair>/pair.json                 organism, original assembler
    sacct                                           elapsed, CPU time, peak memory

ORP columns are prefixed orp_, original-assembly columns orig_, and d_ columns
are ORP minus original. A sample missing a source gets blank cells and a note
in `status` rather than dropping out of the table.

    ./collect_lineage410.py -o lineage410_metrics.csv
    # the dev30 runs, for a before/after table:
    ./collect_lineage410.py --runs '$COMPARE/orp_{group}' --runout-tag v4_1_0-dev30 \\
        -o lineage_dev30_metrics.csv
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

BUSCO_RE = re.compile(
    r"C:(?P<complete>[\d.]+)%\[S:(?P<single>[\d.]+)%,"
    r"D:(?P<duplicated>[\d.]+)%\],F:(?P<fragmented>[\d.]+)%,"
    r"M:(?P<missing>[\d.]+)%,n:(?P<n>\d+)"
)
BUSCO_FIELDS = ["complete", "single", "duplicated", "fragmented", "missing", "n"]
ASSEMBLERS = ["TRINITY", "SPADESAUTO", "SPADESHIGH", "TRANSABYSS"]

# Steps worth a column of their own: the four assemblers, the merge, and the
# two scorers. Every step stays in the timing log; these are the ones a
# runtime figure is built from. Trinity's two phases are summed.
TIMED_STEPS = {
    "trinity": ("run_trinity_phase1", "run_trinity_phase2"),
    "spadesauto": ("run_spadesauto",),
    "spadeshigh": ("run_spadeshigh",),
    "transabyss": ("run_transabyss",),
    "rcorrector": ("run_rcorrector",),
    "merge": ("pool_branch", "twotrack_select", "cdhit"),
    "busco": ("busco",),
    "pytransrate": ("pytransrate",),
}


def hms_to_hours(text):
    h, m, s = (int(x) for x in text.split(":"))
    return round((h * 3600 + m * 60 + s) / 3600, 3)


def slurm_to_hours(text):
    """sacct's [D-]HH:MM:SS or MM:SS.mmm, as hours."""
    if not text:
        return ""
    days = 0
    if "-" in text:
        d, text = text.split("-", 1)
        days = int(d)
    parts = [float(x) for x in text.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    h, m, s = parts
    return round(days * 24 + h + m / 60 + s / 3600, 3)


def rss_to_gb(text):
    if not text:
        return ""
    unit = text[-1].upper()
    scale = {"K": 1 / 1024**2, "M": 1 / 1024, "G": 1, "T": 1024}.get(unit)
    if scale is None:
        return round(float(text) / 1024**3, 2)
    return round(float(text[:-1]) * scale, 2)


def read_busco(directory):
    """C/S/D/F/M/n from BUSCO's short summary; the last score line wins, as in reportgen."""
    if not directory.is_dir():
        return {}
    line = None
    for p in sorted(directory.rglob("short_summary*.txt")):
        for text in p.read_text(errors="replace").splitlines():
            if BUSCO_RE.search(text):
                line = text
    return BUSCO_RE.search(line).groupdict() if line else {}


def read_pytransrate(csv_path):
    if not csv_path or not csv_path.is_file():
        return {}
    rows = list(csv.DictReader(csv_path.open()))
    if not rows:
        return {}
    row = dict(rows[0])
    row.pop("assembly", None)
    return row


def read_qualreport(path):
    out = {}
    if not path.is_file():
        return out
    text = path.read_text(errors="replace")
    for name in ["ORP"] + ASSEMBLERS:
        m = re.search(rf"UNIQUE GENES {name}\s*~*>\s*(\S+)", text)
        if m:
            out[f"unique_genes_{name.lower()}"] = m.group(1)
    m = re.search(r"READS MAPPED AS PROPER PAIRS\s*~*>\s*(\S+)", text)
    if m:
        out["proper_pairs_pct"] = m.group(1).rstrip("%")
    m = re.search(r"version (\S+) \*", text)
    if m:
        out["orp_version"] = m.group(1)
    m = re.search(r"mean\s*:\s*(\S+)", text)
    if m:
        out["strand_mean"] = m.group(1)
    return out


def read_timing(path):
    out = {}
    if not path.is_file():
        return out
    steps = {}
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("Command:"):
            out["orp_command"] = line[len("Command:"):].strip()
            continue
        m = re.match(r"(\S+)\s+(\d+:\d\d:\d\d)\b", line)
        if m:
            steps[m.group(1)] = hms_to_hours(m.group(2))
    if "TOTAL" in steps:
        out["orp_total_hours"] = steps["TOTAL"]
    for col, names in TIMED_STEPS.items():
        if any(n in steps for n in names):
            out[f"hours_{col}"] = round(sum(steps.get(n, 0) for n in names), 3)
    return out


def read_filter_candidates(path):
    """min-cov applied and the salmon numbers it came from."""
    out = {}
    if not path.is_file():
        return out
    for line in path.read_text(errors="replace").splitlines():
        f = line.split("\t")
        if line.startswith("# contigs"):
            out["contigs_before_filter"] = f[1]
        elif line.startswith("# reads_assigned"):
            out["salmon_reads_assigned"] = f[1]
        elif line.startswith("# frag_length_mean"):
            out["frag_length_mean"] = round(float(f[1]), 1)
        elif line.startswith("# min_cov"):
            out["min_cov_applied"] = f[1]
            out["min_cov_mode"] = f[2] if len(f) > 2 else "fixed"
        elif line.startswith("applied:"):
            # rule, pass, nohit_pass, hit_rescued, final, frac_reads_in_final
            out["contigs_after_filter"] = f[4]
            out["frac_reads_kept_by_filter"] = f[5]
            # Files from before the "# min_cov" header carry the threshold
            # only in the rule name, applied:cov>=2.
            m = re.search(r"cov>=([\d.]+)", f[0])
            if m and "min_cov_applied" not in out:
                out["min_cov_applied"] = m.group(1)
                out["min_cov_mode"] = "fixed"
    return out


def read_task_log(path):
    out = {}
    if not path.exists():
        return out
    text = path.read_text(errors="replace")
    m = re.search(r"Input Read Pairs: (\d+) Both Surviving: (\d+)", text)
    if m:
        out["read_pairs_in"] = m.group(1)
        out["read_pairs_trimmed"] = m.group(2)
        out["pct_pairs_surviving_trim"] = round(100 * int(m.group(2)) / int(m.group(1)), 2)
    m = re.search(r"\[kmer\] spadeshigh: k = (\S+)\s+\(reads: mean (\d+) bp, max (\d+) bp\)", text)
    if m:
        out["spadeshigh_k"] = m.group(1)
        out["read_len_mean"] = m.group(2)
        out["read_len_max"] = m.group(3)
    m = re.search(r"\[provenance\] slurm job (\S+)", text)
    if m:
        out["slurm_job"] = m.group(1)
    m = re.search(r"\[provenance\] host (\S+?),", text)
    if m:
        out["node"] = m.group(1)
    return out


def read_sacct(job):
    """Elapsed, CPU time and peak memory of one array task, from its batch step."""
    if not job:
        return {}
    try:
        res = subprocess.run(
            ["sacct", "-j", job, "-n", "-P", "--units=G",
             "--format=JobID,Elapsed,TotalCPU,MaxRSS,State"],
            # not capture_output/text: Premise's anaconda python is 3.6
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return {}
    out = {}
    for line in res.stdout.splitlines():
        jobid, elapsed, cpu, rss, state = (line.split("|") + [""] * 5)[:5]
        if jobid == job:
            out["slurm_state"] = state.split()[0] if state else ""
            out["slurm_elapsed_hours"] = slurm_to_hours(elapsed)
            out["slurm_cpu_hours"] = slurm_to_hours(cpu)
        if rss:
            gb = rss_to_gb(rss)
            if gb != "" and gb > (out.get("slurm_maxrss_gb") or 0):
                out["slurm_maxrss_gb"] = gb
    return out


def number(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    compare = Path(os.environ.get("COMPARE", "/mnt/home/macmaneslab/macmanes/compare"))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=compare / "lineage_410" / "manifest.tsv")
    ap.add_argument("--runs", default=str(compare / "lineage_410"),
                    help="run tree; may contain {group} and $COMPARE "
                         "(default: $COMPARE/lineage_410)")
    ap.add_argument("--runout-tag", default=None,
                    help="the part of the runout after <SRR>_, e.g. v4_1_0 or "
                         "v4_1_0-dev30 (default: the only qualreport in the run, "
                         "or the newest if there are several)")
    ap.add_argument("--busco-original", type=Path, default=compare / "busco_original")
    ap.add_argument("--pytransrate-original", type=Path,
                    default=compare / "pytransrate" / "results")
    ap.add_argument("--orig-unique-genes", type=Path, default=compare / "orig_unique_genes")
    ap.add_argument("--pairs", type=Path, default=compare / "pairs")
    ap.add_argument("--no-sacct", action="store_true", help="skip the sacct lookups")
    ap.add_argument("-o", "--out", type=Path, default=Path("lineage410_metrics.csv"))
    args = ap.parse_args()

    runs_template = args.runs.replace("$COMPARE", str(compare))
    records, tr_cols, problems = [], [], []

    for line in args.manifest.read_text().splitlines():
        if not line.strip():
            continue
        pair, orig_fasta, r1, r2, srr, group, lineage = line.split("\t")
        runs = Path(runs_template.format(group=group))
        rundir = runs / srr
        reports = rundir / "reports"
        rec = {"pair_id": pair, "srr": srr, "group": group, "lineage": lineage}

        pj = args.pairs / pair / "pair.json"
        if pj.is_file():
            meta = json.loads(pj.read_text())
            rec["organism"] = meta.get("organism", "")
            rec["orig_assembler"] = meta.get("assembler", "") or (
                "ORP 2018 orthomerged" if pair.startswith("orp_") else "")
            rec["orig_n_contigs_published"] = meta.get("n_contigs_published", "")
        rec["read_gb"] = round(sum(Path(p).stat().st_size for p in (r1, r2)
                                   if Path(p).is_file()) / 1e9, 2)

        # Which runout: the tag if given, otherwise whatever qualreport is there.
        if args.runout_tag:
            runout = f"{srr}_{args.runout_tag}"
        else:
            found = sorted((p for p in reports.glob(f"qualreport.{srr}_*")
                            if not p.name.endswith(".done")),
                           key=lambda p: p.stat().st_mtime) if reports.is_dir() else []
            runout = found[-1].name[len("qualreport."):] if found else f"{srr}_?"
        rec["runout"] = runout

        notes = []
        rec.update(read_task_log(runs / "logs" / f"{srr}.log"))
        rec.update(read_qualreport(reports / f"qualreport.{runout}"))
        rec.update(read_timing(reports / f"{runout}.timing.log"))
        rec.update(read_filter_candidates(reports / f"{runout}.filter_candidates.tsv"))
        mode = reports / f"{runout}.transabyss.mode"
        if mode.is_file():
            rec["transabyss_mode"] = mode.read_text().strip()
        if not args.no_sacct:
            rec.update(read_sacct(rec.get("slurm_job")))

        busco = read_busco(reports / f"run_{runout}.ORP")
        if not busco:
            notes.append("no ORP busco")
        rec.update({f"orp_busco_{k}": v for k, v in busco.items()})
        obusco = read_busco(args.busco_original / pair)
        if not obusco:
            notes.append("no original busco")
        rec.update({f"orig_busco_{k}": v for k, v in obusco.items()})

        tr = read_pytransrate(next(iter(sorted(
            (reports / f"pytransrate_{runout}").rglob("assemblies.csv"))), None)
            if (reports / f"pytransrate_{runout}").is_dir() else None)
        if not tr:
            notes.append("no ORP pytransrate")
        otr = read_pytransrate(args.pytransrate_original / pair / "assemblies.csv")
        if not otr:
            notes.append("no original pytransrate")
        for col in list(tr) + list(otr):
            if col not in tr_cols:
                tr_cols.append(col)
        rec.update({f"orp_tr_{k}": v for k, v in tr.items()})
        rec.update({f"orig_tr_{k}": v for k, v in otr.items()})

        # The headline differences, ORP minus original.
        for col in ("complete", "single", "duplicated", "fragmented", "missing"):
            a, b = number(busco.get(col)), number(obusco.get(col))
            if a is not None and b is not None:
                rec[f"d_busco_{col}"] = round(a - b, 2)
        for col in ("score", "optimal_score", "p_good_mapping", "n_seqs", "n50"):
            a, b = number(tr.get(col)), number(otr.get(col))
            if a is not None and b is not None:
                rec[f"d_tr_{col}"] = round(a - b, 5)
        ou = args.orig_unique_genes / f"{pair}.unique.txt"
        rec["orig_unique_genes"] = ou.read_text().strip() if ou.is_file() else ""
        if not rec["orig_unique_genes"]:
            notes.append("no original unique genes")
        a, b = number(rec.get("unique_genes_orp")), number(rec["orig_unique_genes"])
        if a is not None and b is not None:
            rec["d_unique_genes"] = int(a - b)
        singles = [number(rec.get(f"unique_genes_{a.lower()}")) for a in ASSEMBLERS]
        singles = [s for s in singles if s is not None]
        orp_genes = number(rec.get("unique_genes_orp"))
        if singles and orp_genes is not None:
            rec["unique_genes_best_single"] = int(max(singles))
            rec["unique_genes_gain_over_best_single"] = int(orp_genes - max(singles))

        if not (reports / f"qualreport.{runout}.done").is_file():
            notes.insert(0, "not finished" if rundir.is_dir() else "no run directory")
        rec["status"] = "; ".join(notes) or "complete"
        if notes:
            problems.append((pair, srr, rec["status"]))
        records.append(rec)

    lead = ["pair_id", "srr", "organism", "group", "lineage", "status", "orp_version",
            "runout", "orig_assembler", "orig_n_contigs_published",
            "read_gb", "read_pairs_in", "read_pairs_trimmed", "pct_pairs_surviving_trim",
            "read_len_mean", "read_len_max", "spadeshigh_k",
            "min_cov_applied", "min_cov_mode", "salmon_reads_assigned", "frag_length_mean",
            "contigs_before_filter", "contigs_after_filter", "frac_reads_kept_by_filter",
            "transabyss_mode"]
    busco_cols = ([f"orp_busco_{c}" for c in BUSCO_FIELDS]
                  + [f"orig_busco_{c}" for c in BUSCO_FIELDS]
                  + [f"d_busco_{c}" for c in BUSCO_FIELDS[:5]])
    gene_cols = (["unique_genes_orp", "orig_unique_genes", "d_unique_genes"]
                 + [f"unique_genes_{a.lower()}" for a in ASSEMBLERS]
                 + ["unique_genes_best_single", "unique_genes_gain_over_best_single",
                    "proper_pairs_pct", "strand_mean"])
    time_cols = (["orp_total_hours"] + [f"hours_{c}" for c in TIMED_STEPS]
                 + ["slurm_elapsed_hours", "slurm_cpu_hours", "slurm_maxrss_gb",
                    "slurm_state", "slurm_job", "node"])
    tr_out = ([f"orp_tr_{c}" for c in tr_cols] + [f"orig_tr_{c}" for c in tr_cols]
              + [f"d_tr_{c}" for c in ("score", "optimal_score", "p_good_mapping",
                                        "n_seqs", "n50")])
    header = lead + busco_cols + gene_cols + time_cols + tr_out + ["orp_command"]

    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)

    done = sum(r["status"] == "complete" for r in records)
    print(f"wrote {len(records)} rows, {len(header)} columns to {args.out} "
          f"({done} complete)", file=sys.stderr)
    for pair, srr, why in problems:
        print(f"  {pair} ({srr}): {why}", file=sys.stderr)


if __name__ == "__main__":
    main()
