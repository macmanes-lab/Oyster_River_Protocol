#!/usr/bin/env python3
"""Where did the reads go? Control contigs missing from a trial, weighted by reads.

usage: reads_lost.py CONTROL_DIR TRIAL_DIR OUT_PREFIX [--blastn DIR] [--threads N]

Two-track (and any trial) maps fewer reads properly than the control. This
attributes the difference to the control contigs that the trial no longer
carries, weighted by how many reads each carried in the control (NumReads
from the control's salmon run on ORP.intermediate.fasta, restricted to the
contigs in its final ORP.fasta).

Each control contig is searched against the trial's final ORP.fasta with
megablast (>= 95% identity) and called:
  represented  its best trial hit covers >= 90% of it
  partial      50-90% covered (e.g. the trial kept a shorter version)
  absent       < 50% covered
Contigs not fully represented are then classed by the control's swissprot
hits (ORP.diamond.txt) and the trial's gene set:
  copy of kept gene   it has a hit and the trial keeps that gene (another
                      isoform, a UTR-extended version, a duplicate)
  gene lost           it has a hit and the trial has no contig for that gene
  no hit              no swissprot hit
Writes OUT_PREFIX.contigs.tsv (one row per control contig) and prints the
read share of each class.
"""
import argparse
import glob
import os
import subprocess
import tempfile
from collections import defaultdict


def gene_of_hits(path):
    """{contig: gene of its best hit}, gene as in oyster.py's extract_gene_ids."""
    best = {}
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            parts = c[1].split("|")
            if len(parts) < 3:
                continue
            b = float(c[11])
            if c[0] not in best or b > best[c[0]][0]:
                best[c[0]] = (b, parts[2].split("_")[0])
    return {k: v[1] for k, v in best.items()}


def fasta_lengths(path):
    lens, name = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lens[name] = 0
            else:
                lens[name] += len(line.strip())
    return lens


def one(pattern):
    hits = glob.glob(pattern)
    if not hits:
        raise SystemExit(f"nothing matches {pattern}")
    return hits[0]


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("control_dir")
    p.add_argument("trial_dir")
    p.add_argument("out_prefix")
    p.add_argument("--blastn", default="/mnt/home/macmaneslab/macmanes/orp_envs/orp/bin")
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args()
    c, t = args.control_dir, args.trial_dir

    c_fa = one(f"{c}/assemblies/*.ORP.fasta")
    t_fa = one(f"{t}/assemblies/*.ORP.fasta")
    c_len = fasta_lengths(c_fa)
    reads = {}
    with open(one(f"{c}/quants/salmon_orthomerged_*/quant.sf")) as f:
        next(f)
        for line in f:
            x = line.rstrip("\n").split("\t")
            reads[x[0]] = float(x[4])
    c_gene = gene_of_hits(one(f"{c}/assemblies/*.ORP.diamond.txt"))
    t_gene = set(gene_of_hits(one(f"{t}/assemblies/*.ORP.diamond.txt")).values())

    # Coverage of each control contig by its best trial contig (megablast).
    cov = defaultdict(float)
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(args.out_prefix))) as tmp:
        db = os.path.join(tmp, "trial")
        subprocess.run([f"{args.blastn}/makeblastdb", "-in", t_fa, "-dbtype", "nucl", "-out", db],
                       check=True, stdout=subprocess.DEVNULL)
        res = subprocess.run(
            [f"{args.blastn}/blastn", "-task", "megablast", "-query", c_fa, "-db", db,
             "-perc_identity", "95", "-evalue", "1e-20", "-max_target_seqs", "5",
             "-num_threads", str(args.threads), "-outfmt", "6 qseqid sseqid qstart qend"],
            check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout  # Python 3.6 on Premise
    spans = defaultdict(lambda: defaultdict(list))
    for line in res.splitlines():
        q, s, qs, qe = line.split("\t")
        a, b = sorted((int(qs), int(qe)))
        spans[q][s].append((a, b))
    for q, subj in spans.items():
        best = 0
        for s, iv in subj.items():
            iv.sort()
            covered, cur_a, cur_b = 0, None, None
            for a, b in iv:
                if cur_b is None or a > cur_b + 1:
                    if cur_b is not None:
                        covered += cur_b - cur_a + 1
                    cur_a, cur_b = a, b
                else:
                    cur_b = max(cur_b, b)
            covered += cur_b - cur_a + 1
            best = max(best, covered)
        cov[q] = best / c_len[q]

    total = sum(reads.get(x, 0.0) for x in c_len)
    share = defaultdict(float)
    count = defaultdict(int)
    with open(args.out_prefix + ".contigs.tsv", "w") as out:
        out.write("contig\tlength\treads\tcoverage\tstatus\tclass\tgene\n")
        for x, ln in c_len.items():
            cv = cov.get(x, 0.0)
            status = "represented" if cv >= 0.9 else "partial" if cv >= 0.5 else "absent"
            g = c_gene.get(x)
            cls = ("-" if status == "represented" else
                   "no hit" if g is None else
                   "copy of kept gene" if g in t_gene else "gene lost")
            key = status if status == "represented" else f"{status}: {cls}"
            share[key] += reads.get(x, 0.0)
            count[key] += 1
            out.write(f"{x}\t{ln}\t{reads.get(x, 0.0):.1f}\t{cv:.3f}\t{status}\t{cls}\t{g or '-'}\n")
    name = os.path.basename(args.out_prefix)
    print(f"{name}: {len(c_len):,} control contigs, {total:,.0f} reads assigned")
    for key in sorted(share, key=lambda k: -share[k]):
        print(f"  {key:32s} {count[key]:8,d} contigs  {share[key] / total:7.2%} of reads")


if __name__ == "__main__":
    main()
