#!/usr/bin/env python3
"""Classify the contigs behind each duplicated BUSCO: isoform, or redundant copy?

usage: dup_busco_pairs.py --full-table full_table.tsv --fasta asm.fasta[.gz]
                          --label NAME --out-prefix PREFIX

A BUSCO is "Duplicated" whenever two or more contigs each carry a complete
copy. That is what an assembly that resolves isoforms should produce, and
also what an assembly that keeps the same transcript two or three times
produces, so the duplicated percentage alone cannot tell good from bad.
This script looks at the contigs themselves. For every duplicated BUSCO it
aligns each pair of its contigs with blastn (shorter as query, longer as
subject, both strands) and puts the pair in one class:

  isoform    the shared sequence is near-identical, and one contig carries
             an internal block the other lacks (an indel of --min-indel bp
             or more between two colinear HSPs, or inside one). That is the
             signature of an alternative exon, or a retained intron.
  redundant  near-identical over nearly all of the shorter contig and no
             such block: the same transcript, kept twice. Nothing is gained
             by keeping both.
  ends       near-identical where they overlap, no internal block, but the
             overlap covers too little of the shorter contig to call it
             redundant: staggered fragments, or different first/last exons.
  variant    they align over most of the shorter contig, but at identity
             below the near-identical cutoff: alleles, recent paralogs, or
             a copy carrying assembly errors.
  divergent  no usable alignment at the DNA level: BUSCO's protein search
             matched two contigs that blastn cannot align. Paralogs, or
             hits to different parts of the protein.

Pairs are also labelled with the assembler each contig came from (by name),
whether two Trinity or rnaSPAdes contigs share a gene ID (the assembler's
own claim that they are isoforms), and the strand they align on.

"minus" strand matters for ORP specifically: the orthogroup step runs
OrthoFinder, which searches DNA with diamond *blastp* -- one strand only --
so two copies of a transcript assembled in opposite orientations can never
be grouped, and both go on to the final assembly. A redundant pair on the
minus strand is a pair the orthogroup step could not have caught.

Writes PREFIX.pairs.tsv (one row per contig pair) and PREFIX.summary.tsv
(one row for this assembly). Needs blastn on PATH (the orp env has it).
"""
import argparse
import collections
import csv
import gzip
import itertools
import os
import re
import subprocess
import sys
import tempfile

TRINITY_GENE = [
    re.compile(r"(TRINITY_DN\d+_c\d+_g\d+)_i\d+"),
    re.compile(r"\b(c\d+_g\d+)_i\d+"),
    re.compile(r"\b(comp\d+_c\d+)_seq\d+"),
]
SPADES_GENE = re.compile(r"^NODE_\d+_length_\d+_cov_[\d.]+_(g\d+)_i\d+")


# chowder.py prefixes every contig with its input's label; --strip-label
# removes it before the name is read for assembler and gene ID.
STRIP_LABEL = None


def unlabel(name):
    return STRIP_LABEL.sub("", name, count=1) if STRIP_LABEL else name


def source(name, header):
    name, header = unlabel(name), unlabel(header)
    if name.startswith("TRINITY_") or any(p.search(header) for p in TRINITY_GENE[1:]):
        return "trinity"
    if name.startswith("NODE_"):
        return "spades"
    if re.match(r"^[RS]\d+$", name):
        return "transabyss"
    return "other"


def gene_id(name, header):
    # Gene IDs restart in every assembly -- rnaSPAdes' g7 at k=55 is not its
    # g7 at k=75 -- so keep the label in the ID when there is one. Without
    # labels (oyster.py's own names) the two spades runs cannot be told
    # apart, and same_gene can over-call between them.
    m = STRIP_LABEL.match(name) if STRIP_LABEL else None
    gid = _gene_id(unlabel(name), unlabel(header))
    return f"{m.group(0)}{gid}" if m and gid else gid


def _gene_id(name, header):
    m = SPADES_GENE.match(name)
    if m:
        return "spades:" + m.group(1)
    for pat in TRINITY_GENE:
        m = pat.search(header)
        if m:
            return "trinity:" + m.group(1)
    return None


def read_full_table(path):
    """{busco_id: [contig, ...]} for the Duplicated rows."""
    dups = collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3 or cols[1] != "Duplicated":
                continue
            # BUSCO 6 transcriptome mode appends the hit's span, name:start-end.
            contig = re.sub(r":\d+-\d+$", "", cols[2])
            if contig not in dups[cols[0]]:
                dups[cols[0]].append(contig)
    return dups


def read_counts(path):
    counts = collections.Counter()
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 2:
                counts[cols[0]] = cols[1]
    return collections.Counter(counts.values())


def read_fasta(path, wanted):
    """{name: (header, seq)} for the contigs in wanted."""
    opener = gzip.open if path.endswith(".gz") else open
    out = {}
    name = header = None
    seq = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name in wanted:
                    out[name] = (header, "".join(seq))
                header = line[1:].strip()
                name = header.split()[0]
                seq = []
            elif name in wanted:
                seq.append(line.strip())
    if name in wanted:
        out[name] = (header, "".join(seq))
    return out


def read_contig_metrics(path):
    """{contig: {"tpm": float, "score": float}} from pytransrate's contigs.csv."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                out[row["contig_name"]] = dict(tpm=float(row["tpm"]),
                                               score=float(row["score"]))
            except (KeyError, ValueError):
                continue
    return out


def blast_pairs(seqs, pairs):
    """Run blastn over each pair; {(q, s): [hsp, ...]} with q the shorter."""
    hits = collections.defaultdict(list)
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "dups.fa")
        with open(fa, "w") as f:
            for name, (_, seq) in seqs.items():
                f.write(f">{name}\n{seq}\n")
        # All-vs-all over only the duplicated contigs -- a few hundred
        # sequences -- is cheaper than one blastn per pair.
        fmt = ("6 qseqid sseqid qstart qend sstart send nident length gaps "
               "bitscore qseq sseq")
        res = subprocess.run(
            ["blastn", "-task", "blastn", "-query", fa, "-subject", fa,
             "-evalue", "1e-10", "-dust", "no", "-outfmt", fmt,
             "-max_hsps", "50"],
            check=True, capture_output=True, text=True)
    wanted = set(pairs)
    for line in res.stdout.splitlines():
        q, s, qs, qe, ss, se, nident, length, gaps, bits, qaln, saln = \
            line.split("\t")
        if (q, s) not in wanted:
            continue
        # blastn's pident divides by the gapped length, so one clean 60 bp
        # exon indel drags an otherwise identical pair below any identity
        # cutoff. Keep identity and indels apart: identity over the
        # ungapped columns, and the longest single gap run separately.
        gap_runs = [len(r) for r in re.findall(r"-+", qaln + " " + saln)]
        hits[(q, s)].append(dict(qs=int(qs), qe=int(qe), ss=int(ss), se=int(se),
                                 nident=int(nident),
                                 ungapped=int(length) - int(gaps),
                                 length=int(length),
                                 max_gap=max(gap_runs, default=0),
                                 bits=float(bits)))
    return hits


def chain(hsps):
    """Best-scoring strand's HSPs, colinear in both sequences, in q order."""
    if not hsps:
        return None, []
    strand = "plus" if max(hsps, key=lambda h: h["bits"])["se"] > \
        max(hsps, key=lambda h: h["bits"])["ss"] else "minus"
    same = [h for h in hsps if (h["se"] > h["ss"]) == (strand == "plus")]
    # Greedy by score: keep an HSP only if it sits colinearly beside every
    # HSP already kept, which drops repeats and tandem hits.
    kept = []
    for h in sorted(same, key=lambda h: -h["bits"]):
        s_lo, s_hi = sorted((h["ss"], h["se"]))
        ok = True
        for k in kept:
            k_lo, k_hi = sorted((k["ss"], k["se"]))
            q_before = h["qe"] <= k["qs"] + 10
            q_after = h["qs"] >= k["qe"] - 10
            if strand == "plus":
                s_before, s_after = s_hi <= k_lo + 10, s_lo >= k_hi - 10
            else:
                s_before, s_after = s_lo >= k_hi - 10, s_hi <= k_lo + 10
            if not ((q_before and s_before) or (q_after and s_after)):
                ok = False
                break
        if ok:
            kept.append(h)
    return strand, sorted(kept, key=lambda h: h["qs"])


def classify(qlen, hsps, args):
    strand, kept = chain([h for h in hsps if h["length"] >= args.min_hsp])
    if not kept:
        return dict(strand="-", qcov=0.0, pid=0.0, max_indel=0, n_hsp=0,
                    cls="divergent")
    covered = set()
    for h in kept:
        covered.update(range(h["qs"], h["qe"] + 1))
    qcov = len(covered) / qlen
    pid = 100 * sum(h["nident"] for h in kept) / sum(h["ungapped"] for h in kept)
    # An internal block present in one contig and not the other shows up
    # as a gap between neighbouring HSPs that is longer on one side, or as
    # one long gap run inside a single HSP.
    max_indel = max(h["max_gap"] for h in kept)
    for a, b in zip(kept, kept[1:]):
        dq = b["qs"] - a["qe"]
        ds = (min(b["ss"], b["se"]) - max(a["ss"], a["se"]) if strand == "plus"
              else min(a["ss"], a["se"]) - max(b["ss"], b["se"]))
        max_indel = max(max_indel, abs(dq - ds))
    if pid >= args.near_identical and max_indel >= args.min_indel:
        cls = "isoform"
    elif pid >= args.near_identical and qcov >= args.redundant_cov:
        cls = "redundant"
    elif pid >= args.near_identical:
        cls = "ends"
    elif qcov >= 0.5:
        cls = "variant"
    else:
        cls = "divergent"
    return dict(strand=strand, qcov=round(qcov, 3), pid=round(pid, 2),
                max_indel=max_indel, n_hsp=len(kept), cls=cls)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--full-table", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--near-identical", type=float, default=98.0,
                   help="percent identity at or above which a pair counts as the "
                        "same sequence (default 98, cd-hit-est's -c in ORP)")
    p.add_argument("--redundant-cov", type=float, default=0.9,
                   help="fraction of the shorter contig a near-identical, "
                        "indel-free alignment must cover to call it redundant")
    p.add_argument("--min-indel", type=int, default=20,
                   help="internal indel, in bp, that marks an isoform")
    p.add_argument("--min-hsp", type=int, default=50,
                   help="ignore HSPs shorter than this")
    p.add_argument("--strip-label", nargs="+", metavar="LABEL",
                   help="chowder.py input labels to strip from the front of contig "
                        "names (e.g. spades55 spades75 transabyss trinity)")
    p.add_argument("--contigs-csv",
                   help="pytransrate's per-contig contigs.csv for this assembly; "
                        "adds each contig's TPM and score, and the minor "
                        "contig's share of the pair's TPM")
    p.add_argument("--minor-share", type=float, default=0.05,
                   help="an isoform pair whose minor contig gets less than this "
                        "share of the pair's TPM is counted as weakly supported")
    args = p.parse_args()
    if args.strip_label:
        global STRIP_LABEL
        STRIP_LABEL = re.compile(
            "^(?:" + "|".join(re.escape(l) for l in args.strip_label) + ")_")

    dups = read_full_table(args.full_table)
    status = read_counts(args.full_table)
    wanted = {c for contigs in dups.values() for c in contigs}
    seqs = read_fasta(args.fasta, wanted)
    missing = wanted - seqs.keys()
    if missing:
        sys.exit(f"{len(missing)} duplicated-BUSCO contigs not in {args.fasta}, "
                 f"e.g. {sorted(missing)[:3]}")

    pairs = {}
    for busco, contigs in dups.items():
        for a, b in itertools.combinations(contigs, 2):
            q, s = (a, b) if len(seqs[a][1]) <= len(seqs[b][1]) else (b, a)
            pairs[(q, s)] = busco
    hits = blast_pairs(seqs, pairs) if pairs else {}

    metrics = read_contig_metrics(args.contigs_csv) if args.contigs_csv else {}

    rows = []
    for (q, s), busco in pairs.items():
        qh, qseq = seqs[q]
        sh, sseq = seqs[s]
        c = classify(len(qseq), hits.get((q, s), []), args)
        gq, gs = gene_id(q, qh), gene_id(s, sh)
        rows.append(dict(
            label=args.label, busco=busco, query=q, subject=s,
            qlen=len(qseq), slen=len(sseq),
            src_q=source(q, qh), src_s=source(s, sh),
            same_gene="yes" if gq and gq == gs else "no", **c))
        mq, ms = metrics.get(q), metrics.get(s)
        if mq and ms:
            total = mq["tpm"] + ms["tpm"]
            rows[-1].update(
                tpm_q=mq["tpm"], tpm_s=ms["tpm"],
                score_q=mq["score"], score_s=ms["score"],
                minor_share=round(min(mq["tpm"], ms["tpm"]) / total, 4) if total else 0.0)
        else:
            rows[-1].update(tpm_q="NA", tpm_s="NA", score_q="NA", score_s="NA",
                            minor_share="NA")

    cols = ["label", "busco", "query", "subject", "qlen", "slen", "src_q",
            "src_s", "same_gene", "strand", "qcov", "pid", "max_indel",
            "n_hsp", "cls", "tpm_q", "tpm_s", "score_q", "score_s", "minor_share"]
    with open(args.out_prefix + ".pairs.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    # Per BUSCO: how many copies are surplus? Count the contigs that are
    # redundant with (or an end-trimmed version of) another kept contig --
    # a greedy lower bound on what a better dedup could drop with nothing
    # lost, since those pairs carry no sequence the other lacks.
    classes = collections.Counter(r["cls"] for r in rows)
    by_busco = collections.defaultdict(list)
    for r in rows:
        by_busco[r["busco"]].append(r)
    only_redundant = 0
    surplus = 0
    for busco, contigs in dups.items():
        rs = by_busco[busco]
        drop = set()
        for r in sorted(rs, key=lambda r: r["qlen"]):
            if r["cls"] in ("redundant", "ends") and r["subject"] not in drop:
                drop.add(r["query"])
        surplus += len(drop)
        if len(contigs) - len(drop) <= 1:
            only_redundant += 1
    redundant = [r for r in rows if r["cls"] == "redundant"]
    isoform = [r for r in rows if r["cls"] == "isoform"]
    summary = collections.OrderedDict(
        label=args.label,
        complete_single=status.get("Complete", 0),
        duplicated=len(dups),
        fragmented=status.get("Fragmented", 0),
        missing=status.get("Missing", 0),
        dup_contigs=len(wanted),
        pairs=len(rows),
        **{f"pairs_{k}": classes.get(k, 0)
           for k in ("isoform", "redundant", "ends", "variant", "divergent")},
        pairs_same_assembler=sum(r["src_q"] == r["src_s"] for r in rows),
        pairs_same_gene_id=sum(r["same_gene"] == "yes" for r in rows),
        isoform_exon_sized=sum(r["max_indel"] >= 50 for r in isoform),
        isoform_minor_weak=(sum(r["minor_share"] < args.minor_share for r in isoform)
                            if metrics else "NA"),
        redundant_minus_strand=sum(r["strand"] == "minus" for r in redundant),
        redundant_cross_assembler=sum(r["src_q"] != r["src_s"] for r in redundant),
        surplus_contigs=surplus,
        dups_explained_by_redundancy=only_redundant,
    )
    with open(args.out_prefix + ".summary.tsv", "w") as f:
        f.write("\t".join(summary) + "\n")
        f.write("\t".join(str(v) for v in summary.values()) + "\n")
    print("\t".join(f"{k}={v}" for k, v in summary.items()))


if __name__ == "__main__":
    main()
