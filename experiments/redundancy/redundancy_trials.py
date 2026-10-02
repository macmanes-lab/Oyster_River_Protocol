#!/usr/bin/env python3
"""One table of every trial: BUSCO, unique genes, transrate, good mappings.

usage: redundancy_trials.py REDUNDANCY_TESTS_DIR [--out redundancy_trials.csv]

Walks the run directories this experiment produced and writes one CSV row per
sample and trial, for runs that finished (BUSCO and transrate both reported):

  validate_arms/<sample>/<trial>            the 10-sample grid
  chowder_arms/SRR1139197/<arm>             SRR1139197's first three arms
  repick_arms/SRR1139197/<arm>.<rule>       SRR1139197's pick-rule re-picks

Columns: sample, trial, search, inflation, pick, contigs, BUSCO complete /
single / duplicated / fragmented / missing (genes of 125), unique genes
(qualreport's UNIQUE GENES ORP), transrate score, optimal score, good
mappings (pytransrate's p_good_mapping, the share of read pairs mapping
properly to the assembly) and proper pairs (qualreport's READS MAPPED AS
PROPER PAIRS, from strandeval's bwa alignment of a 400k-read subsample). Then prints, per trial, the median
of each metric and the per-sample median change against the control.
"""
import argparse
import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_arms import arm_row  # noqa: E402

# trial -> (search, inflation, pick)
TRIALS = {
    "control": ("diamond", "12", "score"),
    "blastn_I12": ("blastn", "12", "score"),
    "candidate": ("blastn", "3", "score_orf"),
    "control_protein": ("diamond", "12", "protein"),
    "blastn_I12_protein": ("blastn", "12", "protein"),
    "candidate_protein": ("blastn", "3", "protein"),
    "twotrack": ("swissprot gene + cd-hit", "-", "two-track"),
}
# SRR1139197-only arms: directory -> (trial name, search, inflation, pick)
EXTRA = {
    "chowder_arms/SRR1139197/blastn_I2": ("blastn_I2", "blastn", "2", "score"),
    "chowder_arms/SRR1139197/blastn_I3": ("blastn_I3", "blastn", "3", "score"),
}
for arm, (search, infl) in {"diamond_I12": ("diamond", "12"), "blastn_I2": ("blastn", "2"),
                            "blastn_I3": ("blastn", "3")}.items():
    for rule in ("score_len", "score_orf", "near_best"):
        if arm == "blastn_I3" and rule == "score_orf":
            continue  # identical to candidate for SRR1139197
        EXTRA[f"repick_arms/SRR1139197/{arm}.{rule}"] = (f"{arm}_{rule}", search, infl, rule)

COLS = ["sample", "trial", "search", "inflation", "pick", "contigs", "busco_complete",
        "busco_single", "busco_duplicated", "busco_fragmented", "busco_missing",
        "unique_genes", "transrate_score", "transrate_optimal", "good_mappings", "proper_pairs"]
METRICS = COLS[5:]


def row(d, sample, trial, search, infl, pick):
    r = arm_row(Path(d))
    if not r or "BUSCO C" not in r or "transrate score" not in r:
        return None
    g = lambda k: round(float(r[k]) * 125 / 100)
    return dict(sample=sample, trial=trial, search=search, inflation=infl, pick=pick,
                contigs=r["final"], busco_complete=g("BUSCO C"), busco_single=g("BUSCO S"),
                busco_duplicated=g("BUSCO D"), busco_fragmented=g("BUSCO F"),
                busco_missing=g("BUSCO M"), unique_genes=r.get("unique genes"),
                transrate_score=round(float(r["transrate score"]), 4),
                transrate_optimal=round(float(r["transrate optimal_score"]), 4),
                good_mappings=round(float(r["transrate p_good_mapping"]), 4),
                proper_pairs=round(r["proper pairs"], 4) if "proper pairs" in r else None)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("root")
    p.add_argument("--out", default="redundancy_trials.csv")
    args = p.parse_args()
    root = Path(args.root)

    rows = []
    for sdir in sorted((root / "validate_arms").iterdir()):
        for trial, (search, infl, pick) in TRIALS.items():
            if (sdir / trial).exists():
                x = row(sdir / trial, sdir.name, trial, search, infl, pick)
                if x:
                    rows.append(x)
    for d, (trial, search, infl, pick) in EXTRA.items():
        if (root / d).exists():
            x = row(root / d, "SRR1139197", trial, search, infl, pick)
            if x:
                rows.append(x)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.out}\n")

    ctrl = {r["sample"]: r for r in rows if r["trial"] == "control"}
    short = {"contigs": "contigs", "busco_complete": "C", "busco_single": "S",
             "busco_duplicated": "D", "busco_fragmented": "F", "busco_missing": "M",
             "unique_genes": "uniq", "transrate_score": "score",
             "transrate_optimal": "optimal", "good_mappings": "good", "proper_pairs": "proper"}
    print(f"{'trial (median over samples)':30s}{'n':>3s}" + "".join(f"{short[m]:>9s}" for m in METRICS))
    for trial in list(TRIALS) + [v[0] for v in EXTRA.values()]:
        sub = [r for r in rows if r["trial"] == trial]
        if not sub:
            continue
        med = lambda m: statistics.median(r[m] for r in sub)
        print(f"{trial:30s}{len(sub):3d}" + "".join(
            f"{med(m):9.3f}" if isinstance(sub[0][m], float) else f"{med(m):9,.0f}" for m in METRICS))
    print(f"\n{'change vs control (median)':30s}{'n':>3s}" + "".join(f"{short[m]:>9s}" for m in METRICS))
    for trial in TRIALS:
        if trial == "control":
            continue
        sub = [r for r in rows if r["trial"] == trial and r["sample"] in ctrl]
        if not sub:
            continue
        d = lambda m: statistics.median(r[m] - ctrl[r["sample"]][m] for r in sub)
        print(f"{trial:30s}{len(sub):3d}" + "".join(
            f"{d(m):+9.3f}" if isinstance(sub[0][m], float) else f"{d(m):+9,.0f}" for m in METRICS))


if __name__ == "__main__":
    main()
