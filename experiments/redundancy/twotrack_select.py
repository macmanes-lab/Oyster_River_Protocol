#!/usr/bin/env python3
"""Two-track redundancy removal: group by swissprot gene, not by orthogroup.

usage: twotrack_select.py --merged merged.fasta --contigs-csv contigs.csv
                          --diamond a.diamond.txt [b ...] --sprot uniprot_sprot.fasta
                          --out good.list [--cdhit cd-hit-est] [--threads N]
                          [--nohit-id 0.95] [--nohit-cov 0.90]

Replaces OrthoFinder plus the per-orthogroup pick. Writes the list of contigs
from the pooled assembly (merged.fasta) to carry forward, in the format
makeorthout's good.<run>.list uses, so everything downstream -- the diamond
rescue, cd-hit-est, salmon and the TPM filter -- runs unchanged.

Track 1, contigs with a swissprot hit. Each contig belongs to the gene of its
best hit (highest bitscore), the gene being the swissprot name with the
species dropped -- extract_gene_ids in oyster.py, so exactly what qualreport's
UNIQUE GENES counts. One contig is kept per gene: the one whose best hit
covers the largest fraction of the swissprot protein, ties to the higher
transrate contig score, then the higher TPM (pytransrate's salmon estimate on
the pool), then the name. Nothing is excluded, so every gene the pool hits
keeps a representative.

Track 2, contigs with no hit. Nucleotide redundancy is removed with
cd-hit-est, both strands, local identity: a contig is dropped when it is at
least --nohit-id identical to a longer one over at least --nohit-cov of its
own length. The cluster representatives are kept here; ORP's TPM filter
downstream (--tpm-filt 1) then drops those below 1 TPM, measured on the
reduced assembly rather than split across redundant copies, while keeping
every contig with a swissprot hit whatever its expression.

Contigs that are in merged.fasta but absent from contigs.csv (pytransrate
did not score them) still count; they rank last within a gene.
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict


def fasta_lengths(path, key=lambda header: header.split()[0]):
    lens, name, n = {}, None, 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name, n = key(line[1:].rstrip("\n")), 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


def best_hits(paths):
    """{contig: (bitscore, subject id, sstart, send)} for each contig's best hit."""
    best = {}
    for p in paths:
        with open(p) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 12:
                    continue
                bits = float(c[11])
                if c[0] not in best or bits > best[c[0]][0]:
                    best[c[0]] = (bits, c[1], int(c[8]), int(c[9]))
    return best


def metrics(path):
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                out[row["contig_name"]] = (float(row["score"]), float(row["tpm"]))
            except (KeyError, ValueError):
                continue
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--merged", required=True)
    p.add_argument("--contigs-csv", required=True)
    p.add_argument("--diamond", nargs="+", required=True)
    p.add_argument("--sprot", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cdhit", default="cd-hit-est")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--nohit-id", type=float, default=0.95)
    p.add_argument("--nohit-cov", type=float, default=0.90)
    args = p.parse_args()

    pool = fasta_lengths(args.merged)
    slen = fasta_lengths(args.sprot)  # keyed by the full sp|ACC|NAME token
    hits = best_hits(args.diamond)
    m = metrics(args.contigs_csv)

    # Track 1: one contig per swissprot gene.
    genes = defaultdict(list)
    for contig, (bits, sid, ss, se) in hits.items():
        if contig not in pool:
            continue  # under 201 bp, so not in the pool; the rescue step covers those
        parts = sid.split("|")
        if len(parts) < 3:
            continue
        gene = parts[2].split("_")[0]
        cov = (abs(se - ss) + 1) / slen[sid] if slen.get(sid) else 0.0
        score, tpm = m.get(contig, (-1.0, -1.0))
        genes[gene].append(((cov, score, tpm), contig))
    keep1 = set()
    for gene, members in genes.items():
        # Highest coverage, then score, then TPM; the name breaks exact ties
        # deterministically (smallest name wins).
        members.sort(key=lambda x: (-x[0][0], -x[0][1], -x[0][2], x[1]))
        keep1.add(members[0][1])
    hit_contigs = {c for g in genes.values() for _, c in g}

    # Track 2: contigs with no hit, nucleotide-deduplicated on both strands.
    nohit = [c for c in pool if c not in hit_contigs]
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(args.out))) as tmp:
        fa, out = os.path.join(tmp, "nohit.fa"), os.path.join(tmp, "nohit.cdhit.fa")
        want, keep = set(nohit), False
        with open(args.merged) as f, open(fa, "w") as o:
            for line in f:
                if line.startswith(">"):
                    keep = line[1:].split()[0] in want
                if keep:
                    o.write(line)
        word = 10 if args.nohit_id >= 0.95 else 8
        subprocess.run([args.cdhit, "-i", fa, "-o", out, "-c", str(args.nohit_id),
                        "-G", "0", "-aS", str(args.nohit_cov), "-r", "1", "-n", str(word),
                        "-M", "0", "-T", str(args.threads), "-d", "0"],
                       check=True, stdout=subprocess.DEVNULL)
        keep2 = set(fasta_lengths(out))

    with open(args.out, "w") as o:
        for c in pool:  # merged.fasta order, so the list reads like the pool
            if c in keep1 or c in keep2:
                o.write(c + "\n")
    print(f"pool {len(pool):,} contigs: {len(hit_contigs):,} with a swissprot hit in "
          f"{len(genes):,} genes -> {len(keep1):,} kept (track 1); "
          f"{len(nohit):,} without -> {len(keep2):,} after cd-hit-est "
          f"-c {args.nohit_id} -aS {args.nohit_cov} (track 2); "
          f"{len(keep1) + len(keep2):,} carried forward")


if __name__ == "__main__":
    main()
