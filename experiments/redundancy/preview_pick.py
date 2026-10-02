#!/usr/bin/env python3
"""Offline preview of pick rules on finished arms: which lost BUSCOs come back?

usage: preview_pick.py VALIDATE_DIR SAMPLES_FILE --pick-script pick_best_contigs.py
                       [--arms candidate blastn_I12] [--rules score score_orf protein ...]

For each sample and arm, re-picks the arm's saved Orthogroups.txt with every
rule (scripts/pick_best_contigs.py's own functions, so this is the picker
itself) and reports:

  changed    picks that differ from today's rule (score)
  bp         total length picked
  recovered  of the BUSCOs this arm lost against the control (present there,
             missing here), how many have the control's BUSCO-carrying contig
             picked under the rule

Recovered is a preview, not a result: a picked contig still has to survive
the rescue step, cd-hit-est and the TPM filter, and BUSCO has to call it.
The per-assembly diamond files are read from the control arm, since every
arm of a sample ingests the same assemblies under the same names.
"""
import argparse
import collections
import glob
import importlib.util
import re
from pathlib import Path


def load_picker(path):
    spec = importlib.util.spec_from_file_location("pick_best_contigs", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def busco(arm_dir):
    path = glob.glob(f"{arm_dir}/reports/run_*.ORP/run_*/full_table.tsv")[0]
    status, contigs = {}, collections.defaultdict(set)
    for line in open(path):
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        status[cols[0]] = cols[1]
        if len(cols) > 2 and cols[2]:
            contigs[cols[0]].add(re.sub(r":\d+-\d+$", "", cols[2]))
    return status, contigs


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("validate_dir")
    p.add_argument("samples_file")
    p.add_argument("--pick-script", required=True)
    p.add_argument("--arms", nargs="+", default=["candidate", "blastn_I12"])
    p.add_argument("--rules", nargs="+",
                   default=["score", "score_orf", "protein", "protein_len"])
    args = p.parse_args()
    pk = load_picker(args.pick_script)

    runs = [l.strip() for l in open(args.samples_file) if l.strip()]
    tot = {(a, r): [0, 0] for a in args.arms for r in args.rules}
    print(f"{'sample':12s}{'arm':12s}{'lost':>5s}" +
          "".join(f"{r:>22s}" for r in args.rules))
    print(f"{'':29s}" + "".join(f"{'changed  Mbp  rec':>22s}" for _ in args.rules))
    for run in runs:
        ctrl = Path(args.validate_dir) / run / "control"
        c_status, c_contigs = busco(ctrl)
        diamond = sorted(p for p in glob.glob(f"{ctrl}/assemblies/diamond/*.diamond.txt")
                         if "orthomerged" not in p)
        prot = {"protein": pk.load_protein(diamond, 11), "protein_len": pk.load_protein(diamond, 3)}
        for arm in args.arms:
            d = Path(args.validate_dir) / run / arm
            groups_txt = glob.glob(f"{d}/orthofuse/*/search/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt")
            csvs = glob.glob(f"{d}/orthofuse/*/merged/contigs.csv")
            if not groups_txt or not csvs:
                continue
            status, _ = busco(d)
            lost = [b for b, s in c_status.items() if s != "Missing" and status.get(b) == "Missing"]
            scores = pk.load_scores(csvs[0])
            metrics = pk.load_metrics(csvs[0])
            groups = list(pk.read_orthogroups(groups_txt[0]))
            base = [pk.best_in_group(m, scores) for _, m in groups]
            line = f"{run:12s}{arm:12s}{len(lost):5d}"
            for rule in args.rules:
                if rule == "score":
                    picks = base
                elif rule in pk.PROTEIN_RULES:
                    picks = [pk.best_by_protein(m, metrics, prot[rule]) for _, m in groups]
                else:
                    picks = [pk.best_by_rule(m, metrics, rule) for _, m in groups]
                chosen = set(x for x in picks if x)
                changed = sum(a != b for a, b in zip(picks, base))
                bp = sum(metrics[x][1] for x in chosen) / 1e6
                rec = sum(bool(c_contigs[b] & chosen) for b in lost)
                tot[(arm, rule)][0] += rec
                tot[(arm, rule)][1] += len(lost)
                line += f"{changed:>10d}{bp:6.1f}{rec:4d}/{len(lost):<2d}"
            print(line)
    print()
    for arm in args.arms:
        print(f"{arm}: lost BUSCOs whose control contig is picked -- " + ", ".join(
            f"{r} {tot[(arm, r)][0]}/{tot[(arm, r)][1]}" for r in args.rules))


if __name__ == "__main__":
    main()
