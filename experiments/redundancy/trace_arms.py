#!/usr/bin/env python3
"""Trace BUSCO changes between a control arm and test arms back through the merge.

usage: trace_arms.py ARMS_DIR RUN CONTROL ARM [ARM ...]

Needs --no-cleanup chowder.py runs (chowder_arms.sbatch). For each test arm:

  lost      BUSCOs complete or fragmented in the control, missing here. For
            each, the control's contig for it: which orthogroup it fell in
            under this arm, how big that group is, and which contig
            makeorthout picked from it instead.
  gained    the reverse.
  dups      for each duplicated BUSCO left in this arm, where its copies came
            from: makeorthout picks from separate orthogroups, or the diamond
            rescue (list7), which adds contigs back after the pick.
"""
import collections
import glob
import re
import sys


def full_table(arm_dir, run):
    path = glob.glob(f"{arm_dir}/reports/run_{run}.ORP/run_*/full_table.tsv")[0]
    status, contigs = {}, collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        status[cols[0]] = cols[1]
        if len(cols) > 2 and cols[2]:
            contigs[cols[0]].append(re.sub(r":\d+-\d+$", "", cols[2]))
    return status, contigs


def orthogroups(arm_dir, run):
    path = glob.glob(f"{arm_dir}/orthofuse/{run}/search/OrthoFinder/Results_*/"
                     "Orthogroups/Orthogroups.txt")[0]
    group_of, members = {}, {}
    for line in open(path):
        toks = line.split()
        if toks:
            members[toks[0]] = toks[1:]
            for m in toks[1:]:
                group_of[m] = toks[0]
    return group_of, members


def id_set(path):
    return {line.strip() for line in open(path) if line.strip()}


def label(contig):
    return contig.split("_", 1)[0]


def main():
    arms_dir, run, control, *arms = sys.argv[1:]
    c_status, c_contigs = full_table(f"{arms_dir}/{control}", f"{run}_{control}")
    for arm in arms:
        d, r = f"{arms_dir}/{arm}", f"{run}_{arm}"
        status, contigs = full_table(d, r)
        group_of, members = orthogroups(d, r)
        picked = id_set(f"{d}/orthofuse/{r}/good.{r}.list")
        rescued = id_set(f"{d}/assemblies/diamond/{r}.list7")
        lost = [b for b, s in c_status.items() if s != "Missing" and status.get(b) == "Missing"]
        gained = [b for b, s in c_status.items() if s == "Missing" and status.get(b) != "Missing"]
        print(f"=== {arm} vs {control}: {len(lost)} BUSCOs lost, {len(gained)} gained")
        for b in lost:
            for c in c_contigs[b][:2]:
                grp = group_of.get(c)
                m = members.get(grp, [])
                pick = next((x for x in m if x in picked), None)
                mix = ", ".join(f"{k}:{v}" for k, v in collections.Counter(map(label, m)).most_common())
                print(f"  lost {b} ({c_status[b]} in control) via {c}")
                print(f"      orthogroup {grp}, {len(m)} contigs [{mix}]; picked {pick}")
        for b in gained:
            print(f"  gained {b} ({status[b]} here) via {contigs[b][0] if contigs[b] else '?'}")
        why = collections.Counter()
        for b, s in status.items():
            if s != "Duplicated":
                continue
            src = ["rescued" if c in rescued else "picked" if c in picked else "other"
                   for c in contigs[b]]
            groups = {group_of.get(c) for c, t in zip(contigs[b], src) if t == "picked"}
            if "rescued" in src:
                why["has a rescued copy"] += 1
            elif all(t == "picked" for t in src) and len(groups) == len(src):
                why["picks from separate orthogroups"] += 1
            else:
                why["other: " + ",".join(src)] += 1
        print(f"  duplicated BUSCOs here: {sum(why.values())} -> {dict(why)}")


if __name__ == "__main__":
    main()
