#!/usr/bin/env python3

#usage: python pick_best_contigs.py contigs.csv Orthogroups.txt good.list.txt
#           [--rule RULE] [--diamond blastx.txt ...]
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
#--rule picks by something other than the score alone (experimental; the
#default, `score`, is the rule above and its output is unchanged). The
#contig score does not depend on length, so a short fragment with clean read
#support can beat the full-length contig in an orthogroup that holds one gene
#from several assemblers (experiments/redundancy):
#  score_len   highest score x length
#  score_orf   highest score x ORF length (contigs.csv column 4)
#  near_best   the longest member scoring at least 0.8 of the group's best
#  protein     the member with the strongest swissprot hit (best bitscore in
#              the --diamond files: ORP's per-assembly diamond blastx output);
#              ties, and groups where no member has a hit, fall back to score
#  protein_len the same, by the longest aligned protein stretch (diamond's
#              length column) rather than bitscore
#All keep the score > 0 floor and the first-seen tie-break.
#
#The protein rules exist because every BUSCO the length rules still lost was
#a group that held the gene's full-length contig and kept a member without
#the gene's protein match.
#
#Group ORDER is preserved exactly, and deliberately: the old glob-and-sort
#ordered groups by *filename* ("1.groups", "10.groups", "100.groups",
#"2.groups", ...), which is lexicographic, not numeric. That order carries
#into good.list.txt, from there into contig order in shucked.fasta, and
#so into cd-hit-est, where input order breaks length ties and can change
#which representative survives into the final assembly. Emitting groups in
#Orthogroups.txt's own line order instead would be an assembly-changing
#change dressed up as a cleanup.

import csv
import os
import sys

RULES = ("score", "score_len", "score_orf", "near_best", "protein", "protein_len")
PROTEIN_RULES = ("protein", "protein_len")
NEAR_BEST = 0.8


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


def load_metrics(contigs_csv_path):
    """{contig: (score, length, orf_length)} from the highest-scoring row,
    the same row load_scores keeps."""
    if not os.path.isfile(contigs_csv_path):
        sys.exit(f"pick_best_contigs.py: contigs.csv not found at '{contigs_csv_path}'")
    metrics = {}
    with open(contigs_csv_path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)  # header
        for row in reader:
            if len(row) < 9:
                continue
            try:
                rec = (float(row[8]), int(row[1]), int(float(row[3])))
            except ValueError:
                continue
            if row[0] not in metrics or rec[0] > metrics[row[0]][0]:
                metrics[row[0]] = rec
    return metrics


def load_protein(diamond_paths, column):
    """{contig: best value of `column` over its hits}. Column 11 is bitscore
    and 3 is alignment length, in diamond's default outfmt 6."""
    best = {}
    for path in diamond_paths:
        if not os.path.isfile(path):
            sys.exit(f"pick_best_contigs.py: diamond output not found at '{path}'")
        with open(path) as handle:
            for line in handle:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue
                try:
                    v = float(cols[column])
                except ValueError:
                    continue
                if v > best.get(cols[0], 0.0):
                    best[cols[0]] = v
    return best


def best_by_protein(members, metrics, protein):
    """Strongest protein evidence among members with score > 0; ties, and
    groups with no hit at all, go to the highest transrate score."""
    scored = [(m, metrics[m]) for m in members if m and m in metrics and metrics[m][0] > 0]
    if not scored:
        return None
    want, top = None, None
    for contig_id, rec in scored:
        key = (protein.get(contig_id, 0.0), rec[0])
        if top is None or key > top:
            want, top = contig_id, key
    return want


def best_by_rule(members, metrics, rule):
    scored = [(m, metrics[m]) for m in members if m and m in metrics and metrics[m][0] > 0]
    if not scored:
        return None
    if rule == "score_len":
        key = lambda rec: rec[0] * rec[1]
    elif rule == "score_orf":
        key = lambda rec: rec[0] * rec[2]
    else:  # near_best
        floor = NEAR_BEST * max(rec[0] for _, rec in scored)
        scored = [(m, rec) for m, rec in scored if rec[0] >= floor]
        key = lambda rec: rec[1]
    want, top = None, None
    for contig_id, rec in scored:
        if top is None or key(rec) > top:
            want, top = contig_id, key(rec)
    return want


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
    args = sys.argv[1:]
    rule = "score"
    diamond_paths = []
    if "--diamond" in args:
        i = args.index("--diamond")
        j = i + 1
        while j < len(args) and not args[j].startswith("--"):
            j += 1
        diamond_paths = args[i + 1:j]
        del args[i:j]
    if "--rule" in args:
        i = args.index("--rule")
        rule = args[i + 1]
        del args[i:i + 2]
        if rule not in RULES:
            sys.exit(f"pick_best_contigs.py: --rule must be one of {', '.join(RULES)}")
    if rule in PROTEIN_RULES and not diamond_paths:
        sys.exit(f"pick_best_contigs.py: --rule {rule} needs --diamond files")
    contigs_csv_path, orthogroups_path, out_path = args[:3]
    if rule == "score":
        scores = load_scores(contigs_csv_path)
        choose = lambda members: best_in_group(members, scores)
    elif rule in PROTEIN_RULES:
        metrics = load_metrics(contigs_csv_path)
        protein = load_protein(diamond_paths, 11 if rule == "protein" else 3)
        choose = lambda members: best_by_protein(members, metrics, protein)
    else:
        metrics = load_metrics(contigs_csv_path)
        choose = lambda members: best_by_rule(members, metrics, rule)

    groups = list(read_orthogroups(orthogroups_path))
    # Lexicographic on the filename the old implementation would have
    # globbed - see the note at the top of this file. Sorting the indices as
    # bare strings gives the same order (the '.' of the suffix sorts below
    # every digit), but keeping the suffix here makes that not something a
    # reader has to work out.
    groups.sort(key=lambda item: f"{item[0]}.groups")

    with open(out_path, "w") as out:
        for _index, members in groups:
            want = choose(members)
            if want is not None:
                out.write(want + "\n")


if __name__ == "__main__":
    main()
