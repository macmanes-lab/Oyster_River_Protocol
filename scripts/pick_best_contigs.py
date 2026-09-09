#!/usr/bin/env python3

#usage: python pick_best_contigs.py contigs.csv Orthogroups.txt good.list.txt
#
#contigs.csv is transrate's per-contig metrics file (comma-delimited, contig
#ID in column 1, score in column 9). Orthogroups.txt is OrthoFinder's output,
#one orthogroup per line: a label token followed by its member contig IDs.
#For each orthogroup, picks the highest-scoring member contig (score must be
#> 0, ties keep the first one seen - same rule the old
#`awk -v max=0 '{if($9>max){want=$1;max=$9}}'` used) and writes that contig
#ID to good.list.txt, one per line.
#
#This used to read one `<i>.groups` file per orthogroup, written by a
#`makegroups` step in oyster.py that split Orthogroups.txt line by line and
#deleted again once this script had run - of order 1e5 small files written,
#globbed back in and unlinked, all to move data between two Python processes.
#Reading Orthogroups.txt directly is the same work without the round-trip.
#
#Group ORDER is preserved exactly, and deliberately: the old glob-and-sort
#ordered groups by *filename* ("1.groups", "10.groups", "100.groups",
#"2.groups", ...), which is lexicographic, not numeric. That order carries
#into good.list.txt, from there into contig order in orthomerged.fasta, and
#so into cd-hit-est, where input order breaks length ties and can change
#which representative survives into the final assembly. Emitting groups in
#Orthogroups.txt's own line order instead would be an assembly-changing
#change dressed up as a cleanup.

import csv
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


def read_orthogroups(orthogroups_path):
    """Yield (index, members) per line of Orthogroups.txt, index 1-based.

    The leading token on each line is OrthoFinder's group label and is
    dropped. Every line is numbered, blank ones included, so an index here
    is the same index the `<i>.groups` file for that line used to carry.
    """
    if not os.path.isfile(orthogroups_path):
        sys.exit(f"pick_best_contigs.py: Orthogroups.txt not found at '{orthogroups_path}'")
    with open(orthogroups_path) as handle:
        for index, line in enumerate(handle, start=1):
            yield index, line.split()[1:]


def best_in_group(members, scores):
    max_score = 0.0
    want = None
    for contig_id in members:
        if not contig_id:
            continue
        score = scores.get(contig_id)
        if score is not None and score > max_score:
            max_score = score
            want = contig_id
    return want


def main():
    contigs_csv_path, orthogroups_path, out_path = sys.argv[1:4]
    scores = load_scores(contigs_csv_path)

    groups = list(read_orthogroups(orthogroups_path))
    # Lexicographic on the filename the old implementation would have
    # globbed - see the note at the top of this file. Sorting the indices as
    # bare strings gives the same order (the '.' of the suffix sorts below
    # every digit), but keeping the suffix here makes that not something a
    # reader has to work out.
    groups.sort(key=lambda item: f"{item[0]}.groups")

    with open(out_path, "w") as out:
        for _index, members in groups:
            want = best_in_group(members, scores)
            if want is not None:
                out.write(want + "\n")


if __name__ == "__main__":
    main()
