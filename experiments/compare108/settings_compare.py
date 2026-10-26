#!/usr/bin/env python3
"""Section 4.4 table: dev30 settings vs the 4.0.0 release runs, for all 59 libraries.

dev30 settings = Trinity read normalization on and a fixed --min-cov 2; the
release runs have normalization off and --min-cov auto. The first 50 libraries
have a dev30-build run (compare/orp_<group>/); the nine orp_ additions have a
run of the release with the dev30 settings (compare/dev30_settings/,
orp_dev30set.sbatch). Each row says which.

Inputs are collect_lineage410.py tables:
  --release   lineage_410/lineage410_metrics.csv            (59 rows)
  --dev30     dev30_metrics.csv     (--runs '$COMPARE/orp_{group}' --runout-tag v4_1_0-dev30)
  --dev30set  dev30set_metrics.csv  (--runs '$COMPARE/dev30_settings')

Writes settings_dev30_vs_400.csv (the Settings_dev30_vs_400 sheet of
Supplementary Data S1) and prints the numbers Section 4.4 quotes.
"""

import argparse
import csv
import statistics as st
import sys

METRICS = [  # label, column, scale
    ("BUSCO complete (%)", "orp_busco_complete", 1),
    ("BUSCO duplicated (%)", "orp_busco_duplicated", 1),
    ("Unique Swiss-Prot genes", "unique_genes_orp", 1),
    ("pyTransRate score", "orp_tr_score", 1),
    ("Good mappings (%)", "orp_tr_p_good_mapping", 100),
    ("Contigs", "orp_tr_n_seqs", 1),
    ("Total wall time (h)", "orp_total_hours", 1),
    ("Trinity time (h)", "hours_trinity", 1),
    ("CPU-hours", "slurm_cpu_hours", 1),
]


def load(path):
    with open(path) as fh:
        return {r["pair_id"]: r for r in csv.DictReader(fh)}


def finished(r):
    """True when the ORP run finished; notes about the original (published)
    assembly, such as 'no original pytransrate', do not matter here."""
    notes = (r or {}).get("status", "")
    return notes == "complete" or (bool(notes) and all(
        n.strip().startswith("no original") for n in notes.split(";")))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", required=True)
    ap.add_argument("--dev30", required=True)
    ap.add_argument("--dev30set", required=True)
    ap.add_argument("-o", "--out", default="settings_dev30_vs_400.csv")
    args = ap.parse_args()
    rel, d30, d30s = load(args.release), load(args.dev30), load(args.dev30set)

    rows, missing = [], []
    for pair, r in rel.items():
        if finished(d30.get(pair)):
            dev, src = d30[pair], "dev30 build"
        elif finished(d30s.get(pair)):
            dev, src = d30s[pair], "4.0.0 with dev30 settings"
        else:
            missing.append(pair)
            continue
        row = {"Pair": pair, "Organism": r.get("organism", ""), "Group": r.get("group", ""),
               "Read pairs": r.get("read_pairs_in", ""), "--min-cov applied (4.0.0)": r.get("min_cov_applied", ""),
               "dev30-settings run": src}
        for label, col, scale in METRICS:
            a, b = num(dev.get(col)), num(r.get(col))
            a = a * scale if a is not None else None
            b = b * scale if b is not None else None
            row["%s: dev30" % label] = a if a is not None else ""
            row["%s: 4.0.0" % label] = b if b is not None else ""
            row["%s: 4.0.0 − dev30" % label] = round(b - a, 5) if a is not None and b is not None else ""
        rows.append(row)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote %d rows to %s" % (len(rows), args.out), file=sys.stderr)
    if missing:
        print("no finished dev30-settings run for: " + ", ".join(missing), file=sys.stderr)

    def d(label):
        return [x["%s: 4.0.0 − dev30" % label] for x in rows if x["%s: 4.0.0 − dev30" % label] != ""]

    def col(label, side):
        return [x["%s: %s" % (label, side)] for x in rows if x["%s: %s" % (label, side)] != ""]

    n = len(rows)
    c, g, s = d("BUSCO complete (%)"), d("Unique Swiss-Prot genes"), d("pyTransRate score")
    print("\nSection 4.4 numbers (n = %d)" % n)
    print("  BUSCO complete: max |diff| %.1f, median %.1f" % (max(abs(x) for x in c), st.median(c)))
    print("  unique genes: %+d to %+d" % (min(g), max(g)))
    print("  pyTransRate score: %+.3f to %+.3f, median %+.3f, lower in %d of %d"
          % (min(s), max(s), st.median(s), sum(1 for x in s if x < 0), len(s)))
    lo = [x for x in rows if num(x["--min-cov applied (4.0.0)"]) is not None and num(x["--min-cov applied (4.0.0)"]) < 2]
    hi = [x for x in rows if num(x["--min-cov applied (4.0.0)"]) is not None and num(x["--min-cov applied (4.0.0)"]) >= 2]
    for name, grp in (("threshold < 2", lo), ("threshold >= 2", hi)):
        pct = [100.0 * x["Contigs: 4.0.0 − dev30"] / x["Contigs: dev30"] for x in grp if x["Contigs: 4.0.0 − dev30"] != ""]
        bc = [abs(x["BUSCO complete (%): 4.0.0 − dev30"]) for x in grp if x["BUSCO complete (%): 4.0.0 − dev30"] != ""]
        if pct:
            print("  %s: %d libraries, contigs %+.0f%% to %+.0f%% (median %+.1f%%), max |BUSCO diff| %.1f, thresholds %s-%s"
                  % (name, len(grp), min(pct), max(pct), st.median(pct), max(bc),
                     min(x["--min-cov applied (4.0.0)"] for x in grp), max(x["--min-cov applied (4.0.0)"] for x in grp)))
    for label in ("Trinity time (h)", "Total wall time (h)", "CPU-hours"):
        a, b = col(label, "dev30"), col(label, "4.0.0")
        if a and b:
            dd = d(label)
            print("  %s: total %.0f -> %.0f h (%+.0f%%), longer in %d of %d"
                  % (label, sum(a), sum(b), 100.0 * (sum(b) - sum(a)) / sum(a), sum(1 for x in dd if x > 0), len(dd)))
    print("  rows by dev30-settings run: " + ", ".join(
        "%s %d" % (k, sum(1 for x in rows if x["dev30-settings run"] == k))
        for k in ("dev30 build", "4.0.0 with dev30 settings")))


if __name__ == "__main__":
    main()
