#!/usr/bin/env python3

#usage: python pick_best_contigs.py contigs.csv orthofuse_dir good.list.txt
#
#contigs.csv is transrate's per-contig metrics file (comma-delimited, contig
#ID in column 1, score in column 9). orthofuse_dir contains one *.groups
#file per orthogroup (one member contig ID per line, produced by groups.done
#upstream). For each orthogroup, picks the highest-scoring member contig
#(score must be > 0, ties keep the first one seen - same rule the old
#`awk -v max=0 '{if($9>max){want=$1;max=$9}}'` used) and writes that contig
#ID to good.list.txt, one per line.

import csv
import glob
import os
import sys


def load_scores(contigs_csv_path):
    if not os.path.isfile(contigs_csv_path):
        sys.exit(f"pick_best_contigs.py: contigs.csv not found at '{contigs_csv_path}'")
    scores = {}
    with open(contigs_csv_path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)  # header
        for row in reader:
            if len(row) < 9:
                continue
            try:
                score = float(row[8])
            except ValueError:
                continue
            contig_id = row[0]
            if contig_id not in scores or score > scores[contig_id]:
                scores[contig_id] = score
    return scores


def best_in_group(groups_path, scores):
    max_score = 0.0
    want = None
    with open(groups_path) as handle:
        for line in handle:
            contig_id = line.strip()
            if not contig_id:
                continue
            score = scores.get(contig_id)
            if score is not None and score > max_score:
                max_score = score
                want = contig_id
    return want


def main():
    contigs_csv_path, orthofuse_dir, out_path = sys.argv[1:4]
    scores = load_scores(contigs_csv_path)

    with open(out_path, "w") as out:
        for groups_path in sorted(glob.glob(os.path.join(orthofuse_dir, "*groups"))):
            want = best_in_group(groups_path, scores)
            if want is not None:
                out.write(want + "\n")


if __name__ == "__main__":
    main()
