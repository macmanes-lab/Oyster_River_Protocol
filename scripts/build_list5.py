#!/usr/bin/env python3

#usage: python build_list5.py list3.txt list5.txt transabyss.diamond.txt spades75.diamond.txt spades55.diamond.txt trinity.diamond.txt
#
#for each gene ID in list3, finds the first diamond hit (searching the
#diamond.txt files in the order given) whose subject gene ID matches, and
#writes the sorted, deduplicated set of contig IDs for those hits to list5.

import sys


def gene_id(diamond_line):
    subject = diamond_line.split("\t", 2)[1]
    accession_field = subject.split("|")
    gene_species = accession_field[2] if len(accession_field) > 2 else subject
    return gene_species.split("_")[0]


def main():
    list3_path, list5_path = sys.argv[1], sys.argv[2]
    diamond_paths = sys.argv[3:]

    with open(list3_path) as handle:
        wanted = {line.strip() for line in handle if line.strip()}

    contig_for_gene = {}
    for diamond_path in diamond_paths:
        with open(diamond_path) as handle:
            for line in handle:
                if "\t" not in line:
                    continue
                gid = gene_id(line)
                if gid in wanted and gid not in contig_for_gene:
                    contig_for_gene[gid] = line.split("\t", 1)[0]

    with open(list5_path, "w") as out:
        for contig in sorted(set(contig_for_gene.values())):
            out.write(contig + "\n")


if __name__ == "__main__":
    main()
