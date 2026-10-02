#!/usr/bin/env python3
"""Does ORP's orthogroup step group a transcript with its own reverse complement?

usage: orthofinder_strand_test.py --fasta asm.fasta[.gz] --workdir DIR
                                  [--n 200] [--inflation 12] [--seed 1]

ORP clusters contigs from its four assemblers with `orthofinder -d`, and
this OrthoFinder's -d does not change the search: it still runs diamond
*blastp*, reading each base as an amino acid. A protein search has one
strand. Unstranded assemblers emit a transcript in whichever orientation
they happen to build it, so if the same transcript comes out forward from
one assembler and reverse-complemented from another, the two copies share
no hit, land in different orthogroups, each wins its own group, and both
reach the final assembly -- where BUSCO counts them as a duplicate.

This test builds four fake "assemblies" from N real contigs, with the
answer known for every pair:

  fwd       the contigs as given
  fwd_mut   the same, with 1% random substitutions   (same strand)
  rc        the reverse complement of fwd            (opposite strand)
  rc_mut    the reverse complement of fwd_mut        (opposite strand)

runs OrthoFinder with ORP's own flags on them, and reports, for each kind
of partner, the fraction of contigs whose partner shares its orthogroup.
All four copies of a contig are one transcript, so a correct clustering
scores 1.0 on every row. A same-strand row near 1.0 with opposite-strand
rows near 0 is the strand blindness described above.

Run it in an environment with orthofinder on PATH (orp_orthofinder).
"""
import argparse
import gzip
import os
import random
import shutil
import subprocess
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


STOPS = {"TAA", "TAG", "TGA"}


def longest_orf(seq):
    """Longest stop-free stretch, in codons, over the three forward frames."""
    best = 0
    for frame in range(3):
        run = 0
        for i in range(frame, len(seq) - 2, 3):
            if seq[i:i + 3] in STOPS:
                best, run = max(best, run), 0
            else:
                run += 1
        best = max(best, run)
    return best


def orient(seq):
    """seq, or its reverse complement if that carries the longer ORF.

    Symmetric by construction -- orient(s) and orient(revcomp(s)) come out
    identical unless the two strands tie -- so it gives every contig a
    canonical orientation whether or not it codes for anything.
    """
    rc = seq.translate(COMP)[::-1]
    return rc if longest_orf(rc) > longest_orf(seq) else seq


def mutate(seq, rate, rng):
    out = list(seq)
    for i, b in enumerate(out):
        if b in "ACGT" and rng.random() < rate:
            out[i] = rng.choice([x for x in "ACGT" if x != b])
    return "".join(out)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--fasta", required=True, help="source of real contigs, e.g. a Trinity assembly")
    p.add_argument("--workdir", required=True)
    p.add_argument("--n", type=int, default=200)
    p.add_argument("--min-len", type=int, default=800)
    p.add_argument("--inflation", default="12", help="MCL -I; ORP uses 12")
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--search", default=None,
                   help="OrthoFinder -S program (default: OrthoFinder's own, diamond "
                        "blastp; ORP passes a diamond_orp_<threads> copy of it)")
    p.add_argument("--orient", action="store_true",
                   help="orient every contig by its longest ORF before the search, "
                        "the candidate fix")
    args = p.parse_args()

    rng = random.Random(args.seed)
    # One contig per Trinity gene, so no two sources are isoforms of each
    # other and every "should group" pair below is exactly one transcript.
    seen_genes, pool = set(), []
    for name, seq in read_fasta(args.fasta):
        gene = name.rsplit("_i", 1)[0]
        if len(seq) >= args.min_len and gene not in seen_genes and "N" not in seq.upper():
            seen_genes.add(gene)
            pool.append((name, seq.upper()))
    if len(pool) < args.n:
        sys.exit(f"only {len(pool)} usable contigs in {args.fasta}")
    picked = rng.sample(pool, args.n)

    fa_dir = os.path.join(args.workdir, "fasta")
    shutil.rmtree(args.workdir, ignore_errors=True)
    os.makedirs(fa_dir)
    seqs = {}
    for i, (_, s) in enumerate(picked):
        seqs[("fwd", i)] = s
        seqs[("fwd_mut", i)] = mutate(s, 0.01, rng)
        seqs[("rc", i)] = s.translate(COMP)[::-1]
        seqs[("rc_mut", i)] = seqs[("fwd_mut", i)].translate(COMP)[::-1]
    for kind in ("fwd", "fwd_mut", "rc", "rc_mut"):
        with open(os.path.join(fa_dir, f"{kind}.fa"), "w") as f:
            for i in range(args.n):
                seq = orient(seqs[(kind, i)]) if args.orient else seqs[(kind, i)]
                f.write(f">{kind}_{i}\n{seq}\n")

    subprocess.run(["orthofinder", "-d", "-I", str(args.inflation), "-f", fa_dir,
                    "-og", "-t", str(args.threads), "-a", "1",
                    *(["-S", args.search] if args.search else [])],
                   check=True, stdout=subprocess.DEVNULL)
    groups_txt = None
    for root, _, files in os.walk(fa_dir):
        if "Orthogroups.txt" in files:
            groups_txt = os.path.join(root, "Orthogroups.txt")
    if groups_txt is None:
        sys.exit("OrthoFinder wrote no Orthogroups.txt")

    group_of = {}
    n_groups = 0
    with open(groups_txt) as f:
        for line in f:
            toks = line.split()
            if not toks:
                continue
            n_groups += 1
            for member in toks[1:]:
                group_of[member] = toks[0]

    def together(a, b):
        hits = sum(group_of.get(f"{a}_{i}") is not None and
                   group_of.get(f"{a}_{i}") == group_of.get(f"{b}_{i}")
                   for i in range(args.n))
        return hits / args.n

    all4 = sum(len({group_of.get(f"{k}_{i}") for k in ("fwd", "fwd_mut", "rc", "rc_mut")}) == 1
               for i in range(args.n)) / args.n
    print(f"OrthoFinder -d -I {args.inflation}{' -S ' + args.search if args.search else ''}{' (oriented by longest ORF)' if args.orient else ''}, {args.n} transcripts x 4 copies, "
          f"{n_groups} orthogroups (ideal: {args.n})")
    print(f"  fwd   ~ fwd_mut  same strand, 1% subs    {together('fwd', 'fwd_mut'):.2f}")
    print(f"  rc    ~ rc_mut   same strand, 1% subs    {together('rc', 'rc_mut'):.2f}")
    print(f"  fwd   ~ rc       opposite strand, exact  {together('fwd', 'rc'):.2f}")
    print(f"  fwd   ~ rc_mut   opposite strand, 1% subs {together('fwd', 'rc_mut'):.2f}")
    print(f"  all four copies in one orthogroup        {all4:.2f}")


if __name__ == "__main__":
    main()
