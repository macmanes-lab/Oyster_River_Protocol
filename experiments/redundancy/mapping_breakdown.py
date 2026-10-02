#!/usr/bin/env python3
"""Mapped vs good-mapping rates per arm, from pytransrate's assemblies.csv.

usage: mapping_breakdown.py VALIDATE_DIR SAMPLES_FILE ARM [ARM ...]

pytransrate aligns every read pair to the final assembly. A fragment is
"mapped" if it aligns at all, and a "good" mapping only if both mates align,
in the right orientation, to the same contig, at a plausible distance -- so
a drop in good mappings can be reads that no longer align, or pairs that
still align but are split or broken (bad_mappings). potential_bridges counts
pairs whose mates land on two different contigs.
"""
import csv
import glob
import sys


def main():
    vdir, samples_file, *arms = sys.argv[1:]
    samples = [l.strip() for l in open(samples_file) if l.strip()]
    print(f"{'sample':11s}{'arm':11s}{'mapped':>8s}{'good':>8s}{'bad':>8s}{'bridges':>10s}")
    for s in samples:
        for arm in arms:
            f = glob.glob(f"{vdir}/{s}/{arm}/reports/transrate_*/assemblies.csv")
            if not f:
                continue
            r = next(csv.DictReader(open(f[0])))
            frag = float(r["fragments"])
            print(f"{s:11s}{arm:11s}{float(r['p_fragments_mapped']):8.3f}"
                  f"{float(r['p_good_mapping']):8.3f}{float(r['bad_mappings']) / frag:8.3f}"
                  f"{int(float(r['potential_bridges'])):10,d}")


if __name__ == "__main__":
    main()
