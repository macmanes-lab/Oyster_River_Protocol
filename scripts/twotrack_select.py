#!/usr/bin/env python3

#usage: python twotrack_select.py --pool pool.fasta --contigs-csv contigs.csv
#           --diamond a.diamond.txt [b ...] --out good.list
#           [--sprot uniprot_sprot.fasta] [--table selection.tsv] [--threads N]
#
#The two-track merge, which replaces OrthoFinder plus the per-orthogroup pick
#(ORP 4.1.0). It writes the list of pooled contigs (pool.fasta) to carry
#forward, in the format of makeorthout's good.<run>.list, so everything
#downstream -- the diamond rescue, cd-hit-est, salmon and the TPM filter --
#runs unchanged.
#
#Why: OrthoFinder's orthogroups do not line up with genes. Some split one gene
#(each part keeps a contig: the duplicated BUSCOs), others mix several genes
#(keeping one member drops the rest: the lost genes). Grouping by swissprot
#gene makes the groups the very thing BUSCO and UNIQUE GENES count. Tested on
#10 samples against ORP 4.0 (experiments/redundancy/README.md).
#
#Track 1, contigs with a swissprot hit (ORP's per-assembly diamond blastx). A
#contig belongs to the gene of its best hit (highest bitscore), the gene being
#the swissprot name without its species -- extract_gene_ids in oyster.py, so
#exactly what qualreport's UNIQUE GENES counts.
#  Step 1, the representative: among the gene's contigs whose best hit covers
#    at least --near (0.9) of the best protein coverage any of them reaches,
#    the longest; ties to the higher transrate contig score, then TPM, then
#    name. Protein coverage alone favoured coding-only contigs and lost the
#    UTRs and ends where many reads land; length among the near-best restores
#    them without adding a contig.
#  Step 2, distinct copies: a dropped contig of a kept gene comes back when
#    (a) the representative covers less than --distinct-cov (0.5) of it at
#    >= 95% identity -- another isoform, a paralog under the same swissprot
#    name, a different region -- and (b) it is expressed at >= --rescue-tpm
#    (1) TPM. Such copies are first deduplicated among themselves
#    (cd-hit-est), then added in order of TPM up to --max-per-gene (3)
#    contigs per gene. This puts back sequence the representative lacks, at
#    some cost in duplication.
#Every gene in the pool keeps at least its representative.
#
#Track 2, contigs with no hit: cd-hit-est on both strands with local identity
#(-c 0.95 -G 0 -aS 0.9 -r 1), so a contig goes when it is >= 95% identical to
#a longer one over >= 90% of its own length. ORP's TPM filter downstream
#(--tpm-filt) then drops those below threshold, measured on the reduced
#assembly, while keeping every contig with a swissprot hit.
#
#Protein coverage needs each swissprot protein's length, read from --sprot
#(the uniprot_sprot.fasta the diamond database was built from). Without it,
#contigs are ranked by aligned length instead.
#
#TPM is pytransrate's salmon estimate on the pool (contigs.csv), split across
#redundant copies, so it understates the expression of any one of them.
#
#--table writes one row per pooled contig: contig, track, gene, role
#(representative / rescued / dropped / nohit_kept / nohit_dropped), protein
#coverage, length, transrate score, TPM.

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict


def fasta_lengths(path):
    lens, name, n = {}, None, 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


def write_subset(src, ids, dst):
    keep = False
    with open(src) as inf, open(dst, "w") as out:
        for line in inf:
            if line.startswith(">"):
                keep = line[1:].split()[0] in ids
            if keep:
                out.write(line)


def best_hits(paths):
    """{contig: (bitscore, subject, sstart, send, aligned length)}"""
    best = {}
    for path in paths:
        if not os.path.isfile(path):
            sys.exit(f"twotrack_select.py: diamond output not found at '{path}'")
        with open(path) as handle:
            for line in handle:
                c = line.rstrip("\n").split("\t")
                if len(c) < 12:
                    continue
                bits = float(c[11])
                if c[0] not in best or bits > best[c[0]][0]:
                    best[c[0]] = (bits, c[1], int(c[8]), int(c[9]), int(c[3]))
    return best


def contig_metrics(path):
    out = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                out[row["contig_name"]] = (float(row["score"]), float(row["tpm"]))
            except (KeyError, ValueError):
                continue
    return out


def cdhit(cdhit_bin, src, dst, ident, cov, threads, mem_mb=0):
    word = 10 if ident >= 0.95 else 8
    subprocess.run([cdhit_bin, "-i", src, "-o", dst, "-c", str(ident), "-G", "0",
                    "-aS", str(cov), "-r", "1", "-n", str(word), "-M", str(mem_mb),
                    "-T", str(threads), "-d", "0"],
                   check=True, stdout=subprocess.DEVNULL)
    return set(fasta_lengths(dst))


def coverage_by_rep(blastn_bin, makeblastdb_bin, query_fa, rep_fa, rep_of, tmp, threads):
    """{query contig: fraction of it covered by its own gene's representative},
    megablast at >= 95% identity, union of HSPs on that one subject."""
    db = os.path.join(tmp, "reps")
    subprocess.run([makeblastdb_bin, "-in", rep_fa, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL)
    out = subprocess.run(
        [blastn_bin, "-task", "megablast", "-query", query_fa, "-db", db,
         "-perc_identity", "95", "-evalue", "1e-10", "-max_target_seqs", "20",
         "-num_threads", str(threads), "-outfmt", "6 qseqid sseqid qstart qend"],
        check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout
    spans = defaultdict(list)
    for line in out.splitlines():
        q, s, qs, qe = line.split("\t")
        if rep_of.get(q) == s:
            a, b = sorted((int(qs), int(qe)))
            spans[q].append((a, b))
    covered = {}
    for q, iv in spans.items():
        iv.sort()
        total, cur_a, cur_b = 0, None, None
        for a, b in iv:
            if cur_b is None or a > cur_b + 1:
                if cur_b is not None:
                    total += cur_b - cur_a + 1
                cur_a, cur_b = a, b
            else:
                cur_b = max(cur_b, b)
        total += cur_b - cur_a + 1
        covered[q] = total
    return covered


def main():
    p = argparse.ArgumentParser(description="ORP's two-track merge selection")
    p.add_argument("--pool", required=True)
    p.add_argument("--contigs-csv", required=True)
    p.add_argument("--diamond", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--sprot", default=None)
    p.add_argument("--table", default=None)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--mem-mb", type=int, default=0,
                   help="cd-hit-est -M, in MB; 0 is cd-hit's 'no limit' (default: 0)")
    p.add_argument("--near", type=float, default=0.9)
    p.add_argument("--distinct-cov", type=float, default=0.5)
    p.add_argument("--rescue-tpm", type=float, default=1.0)
    p.add_argument("--max-per-gene", type=int, default=3)
    p.add_argument("--nohit-id", type=float, default=0.95)
    p.add_argument("--nohit-cov", type=float, default=0.90)
    p.add_argument("--cdhit", default="cd-hit-est")
    p.add_argument("--blastn", default="blastn")
    p.add_argument("--makeblastdb", default="makeblastdb")
    args = p.parse_args()

    pool = fasta_lengths(args.pool)
    hits = best_hits(args.diamond)
    metrics = contig_metrics(args.contigs_csv)
    slen = fasta_lengths(args.sprot) if args.sprot and os.path.isfile(args.sprot) else {}
    if not slen:
        print("twotrack_select.py: no swissprot fasta; ranking by aligned length, not coverage")

    # Track 1, step 1: one representative per gene.
    genes = defaultdict(list)
    info = {}
    for contig, (bits, sid, ss, se, alen) in hits.items():
        if contig not in pool:
            continue  # under the pool's length floor; the diamond rescue covers those
        parts = sid.split("|")
        if len(parts) < 3:
            continue
        gene = parts[2].split("_")[0]
        cov = (abs(se - ss) + 1) / slen[sid] if slen.get(sid) else float(alen)
        score, tpm = metrics.get(contig, (-1.0, -1.0))
        info[contig] = (gene, cov, pool[contig], score, tpm)
        genes[gene].append(contig)
    rep_of_gene = {}
    for gene, members in genes.items():
        best_cov = max(info[c][1] for c in members)
        near = [c for c in members if info[c][1] >= args.near * best_cov]
        near.sort(key=lambda c: (-info[c][2], -info[c][3], -info[c][4], c))
        rep_of_gene[gene] = near[0]
    reps = set(rep_of_gene.values())

    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(args.out))) as tmp:
        # Track 1, step 2: distinct, expressed copies of kept genes.
        others = [c for c in info if c not in reps and info[c][4] >= args.rescue_tpm]
        rep_of = {c: rep_of_gene[info[c][0]] for c in others}
        rescued = set()
        if others:
            q_fa, r_fa = os.path.join(tmp, "others.fa"), os.path.join(tmp, "reps.fa")
            write_subset(args.pool, set(others), q_fa)
            write_subset(args.pool, reps, r_fa)
            covered = coverage_by_rep(args.blastn, args.makeblastdb, q_fa, r_fa, rep_of,
                                      tmp, args.threads)
            distinct = {c for c in others
                        if covered.get(c, 0) / pool[c] < args.distinct_cov}
            if distinct:
                d_fa, d_out = os.path.join(tmp, "distinct.fa"), os.path.join(tmp, "distinct.cd.fa")
                write_subset(args.pool, distinct, d_fa)
                distinct = cdhit(args.cdhit, d_fa, d_out, 0.95, 0.9, args.threads, args.mem_mb)
            by_gene = defaultdict(list)
            for c in distinct:
                by_gene[info[c][0]].append(c)
            for gene, cs in by_gene.items():
                cs.sort(key=lambda c: (-info[c][4], -info[c][2], c))
                rescued.update(cs[:max(0, args.max_per_gene - 1)])

        # Track 2: contigs with no hit, nucleotide-deduplicated on both strands.
        nohit = [c for c in pool if c not in info]
        n_fa, n_out = os.path.join(tmp, "nohit.fa"), os.path.join(tmp, "nohit.cd.fa")
        write_subset(args.pool, set(nohit), n_fa)
        keep2 = set()
        if nohit:
            keep2 = cdhit(args.cdhit, n_fa, n_out, args.nohit_id, args.nohit_cov,
                          args.threads, args.mem_mb)

    keep = reps | rescued | keep2
    with open(args.out, "w") as out:
        for c in pool:  # pool.fasta order
            if c in keep:
                out.write(c + "\n")
    if args.table:
        with open(args.table, "w") as out:
            out.write("contig\ttrack\tgene\trole\tprotein_cov\tlength\tscore\ttpm\n")
            for c in pool:
                if c in info:
                    g, cov, ln, sc, tpm = info[c]
                    role = "representative" if c in reps else "rescued" if c in rescued else "dropped"
                    out.write(f"{c}\thit\t{g}\t{role}\t{cov:.3f}\t{ln}\t{sc:.4f}\t{tpm:.3f}\n")
                else:
                    m = metrics.get(c, (-1.0, -1.0))
                    role = "nohit_kept" if c in keep2 else "nohit_dropped"
                    out.write(f"{c}\tnohit\t-\t{role}\t-\t{pool[c]}\t{m[0]:.4f}\t{m[1]:.3f}\n")
    print(f"two-track: {len(pool):,} pooled contigs; {len(info):,} with a swissprot hit in "
          f"{len(genes):,} genes -> {len(reps):,} representatives + {len(rescued):,} distinct "
          f"copies rescued; {len(nohit):,} without a hit -> {len(keep2):,} after cd-hit-est; "
          f"{len(keep):,} carried forward")


if __name__ == "__main__":
    main()
