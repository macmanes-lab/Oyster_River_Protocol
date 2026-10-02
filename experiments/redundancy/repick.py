#!/usr/bin/env python3
"""Re-pick one contig per orthogroup under an alternative rule.

usage: repick.py --contigs-csv merged/contigs.csv --orthogroups Orthogroups.txt
                 --rule RULE --out good.list [--pick-script scripts/pick_best_contigs.py]
       repick.py --contigs-csv ... --orthogroups ... --compare [--watch CONTIG ...]

makeorthout keeps the member with the highest pytransrate contig score.
That score does not depend on length, so in an orthogroup that correctly
holds one gene from several assemblers, a short fragment with clean read
support can beat the full-length contig. The chowder_arms test lost 5 BUSCOs
that way under blastn. The rules:

  score        today's rule: highest score, which must be above 0. Kept as
               a sanity check -- its output must match good.<run>.list
               byte for byte.
  score_len    highest score x length.
  score_orf    highest score x ORF length (pytransrate's orf_length).
  near_best    the longest member whose score is at least --near (0.8)
               of the group's best.

Every rule keeps the score > 0 floor, breaks ties by first member seen, and
writes groups in the same lexicographic order as pick_best_contigs.py
(group order reaches cd-hit-est's tie-breaks).

--compare prints, for every rule, how many picks differ from `score`, the
total and median length picked, and whether each --watch contig is picked.
"""
import argparse
import csv
import importlib.util
import statistics
import sys
from pathlib import Path

RULES = ("score", "score_len", "score_orf", "near_best")


def load_picker(path):
    spec = importlib.util.spec_from_file_location("pick_best_contigs", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_metrics(path):
    """{contig: (score, length, orf_length)}, the row with the highest score
    when a contig appears twice -- the same row pick_best_contigs.py uses."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                rec = (float(row["score"]), int(row["length"]), int(float(row["orf_length"])))
            except (KeyError, ValueError):
                continue
            name = row["contig_name"]
            if name not in out or rec[0] > out[name][0]:
                out[name] = rec
    return out


def pick(members, metrics, rule, near):
    scored = [(m, metrics[m]) for m in members if m and m in metrics and metrics[m][0] > 0]
    if not scored:
        return None
    if rule == "score":
        key = lambda mr: mr[1][0]
    elif rule == "score_len":
        key = lambda mr: mr[1][0] * mr[1][1]
    elif rule == "score_orf":
        key = lambda mr: mr[1][0] * mr[1][2]
    elif rule == "near_best":
        best = max(r[0] for _, r in scored)
        scored = [mr for mr in scored if mr[1][0] >= near * best]
        key = lambda mr: mr[1][1]
    else:
        raise ValueError(rule)
    want, top = None, None
    for m, r in scored:  # strict > keeps the first of equals
        k = key((m, r))
        if top is None or k > top:
            want, top = m, k
    return want


def ordered_groups(picker, orthogroups):
    groups = list(picker.read_orthogroups(orthogroups))
    groups.sort(key=lambda item: f"{item[0]}.groups")
    return groups


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--contigs-csv", required=True)
    p.add_argument("--orthogroups", required=True)
    p.add_argument("--rule", choices=RULES)
    p.add_argument("--near", type=float, default=0.8)
    p.add_argument("--out")
    p.add_argument("--compare", action="store_true")
    p.add_argument("--watch", nargs="*", default=[])
    p.add_argument("--pick-script",
                   default=str(Path(__file__).resolve().parents[2] / "scripts" / "pick_best_contigs.py"))
    args = p.parse_args()

    picker = load_picker(args.pick_script)
    metrics = load_metrics(args.contigs_csv)
    groups = ordered_groups(picker, args.orthogroups)

    if args.compare:
        base = [pick(m, metrics, "score", args.near) for _, m in groups]
        print(f"{'rule':12s}{'picks':>9s}{'changed':>9s}{'total bp':>14s}{'median bp':>11s}"
              + "".join(f"  {w[:28]:>28s}" for w in args.watch))
        for rule in RULES:
            got = [pick(m, metrics, rule, args.near) for _, m in groups]
            chosen = [g for g in got if g]
            lens = [metrics[g][1] for g in chosen]
            changed = sum(a != b for a, b in zip(got, base))
            picked = set(chosen)
            print(f"{rule:12s}{len(chosen):9d}{changed:9d}{sum(lens):14,d}{statistics.median(lens):11.0f}"
                  + "".join(f"  {('picked' if w in picked else '-'):>28s}" for w in args.watch))
        return

    if not args.rule or not args.out:
        sys.exit("--rule and --out are required unless --compare")
    with open(args.out, "w") as out:
        for _, members in groups:
            want = pick(members, metrics, args.rule, args.near)
            if want is not None:
                out.write(want + "\n")


if __name__ == "__main__":
    main()
