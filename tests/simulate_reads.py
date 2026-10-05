#!/usr/bin/env python3
"""Simulate a small stranded RNA-seq library whose answer is known.

usage: python tests/simulate_reads.py --proteins BUSCO_ancestral --out DIR [--seed N]

The release check's sampledata pair is real reads, so a run on it can only
be judged on whether it finished and wrote sane files. This library is built
so a run can be judged on what it got right:

  busco_high   16 BUSCO genes (back-translated from the lineage's ancestral
               proteins), 80x       -> recovered; BUSCO finds them
  busco_low     4 BUSCO genes, 5x   -> below --tpm-filt, but kept for their
                                       swissprot hit
  nohit_high    5 random transcripts, 80x  -> kept by the no-hit track
  nohit_low     5 random transcripts, 5x   -> removed by --tpm-filt

Reads are 2x100 bp phred+33, dUTP-stranded (RF: read 1 is antisense), with
0.2% substitutions and a few Ns (both at Q2), quality falling off toward 3', and ~5% of fragments shorter than a read so their
reads run into the TruSeq adapter -- which must not survive into the
assembly. Everything is drawn from one seeded RNG, so a given --seed and
protein file give byte-identical output.

Writes DIR/sim_1.fq.gz, DIR/sim_2.fq.gz, DIR/transcripts.fa and
DIR/truth.tsv (id, class, length, pairs, expected TPM), plus
DIR/tpm_filt.txt: a --tpm-filt halfway (geometrically) between the low and
high classes' expected TPM.

The BUSCO proteins come from the lineage dataset already installed for the
pipeline (busco_dbs/<lineage>/ancestral, or refseq_db.faa.gz), so nothing
needs downloading and the set is fixed for a given lineage version.
"""

import argparse
import gzip
import math
import random
import re
from pathlib import Path

READ_LEN = 100
FRAG_MEAN, FRAG_SD = 220, 40
SHORT_FRAG_SHARE = 0.05
SUB_RATE = 0.002
N_READ_SHARE = 0.02
HIGH_COV, LOW_COV = 80, 5
CLASSES = (("busco_high", 16, HIGH_COV), ("busco_low", 4, LOW_COV),
           ("nohit_high", 5, HIGH_COV), ("nohit_low", 5, LOW_COV))
# Read-through adapters, as in barcodes/barcodes.fa (PE2_rc / PE1_rc).
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"

CODONS = {
    "A": ["GCT", "GCC", "GCA", "GCG"], "R": ["CGT", "CGC", "AGA", "AGG"],
    "N": ["AAT", "AAC"], "D": ["GAT", "GAC"], "C": ["TGT", "TGC"],
    "Q": ["CAA", "CAG"], "E": ["GAA", "GAG"], "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"], "I": ["ATT", "ATC", "ATA"], "L": ["CTT", "CTC", "CTG", "TTG"],
    "K": ["AAA", "AAG"], "M": ["ATG"], "F": ["TTT", "TTC"], "P": ["CCT", "CCC", "CCA", "CCG"],
    "S": ["TCT", "TCC", "AGC", "AGT"], "T": ["ACT", "ACC", "ACA", "ACG"], "W": ["TGG"],
    "Y": ["TAT", "TAC"], "V": ["GTT", "GTC", "GTA", "GTG"],
}
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs, name = {}, None
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def find_busco_proteins(busco_root, lineage):
    """The lineage's ancestral protein file under busco_root, or None."""
    base = re.sub(r"\.\d+$", "", Path(lineage).name)
    for depth in ("*", "*/*", "*/*/*"):
        for d in sorted(Path(busco_root).glob(depth)):
            if d.is_dir() and d.name.startswith(base):
                for name in ("ancestral", "ancestral_variants", "refseq_db.faa.gz"):
                    if (d / name).is_file():
                        return d / name
    return None


def random_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def back_translate(rng, protein):
    aas = [a for a in protein.upper() if a != "*"]
    body = "".join(rng.choice(CODONS.get(a) or CODONS[rng.choice("ACDEFGHIKLMNPQRSTVWY")])
                   for a in aas[1:])
    return "ATG" + body + rng.choice(["TAA", "TAG", "TGA"])


LONG_ORF = re.compile(r"ATG(?:(?!TAA|TAG|TGA)...){100}")


def no_long_orf(rng, n):
    """Random sequence, re-drawn until neither strand carries a 100-codon
    ORF, so a 'no hit' transcript can't hit swissprot by chance."""
    while True:
        s = random_seq(rng, n)
        if not (LONG_ORF.search(s) or LONG_ORF.search(revcomp(s))):
            return s


def base_quality(i):
    """Phred+33, falling off toward the 3' end the way Illumina's does. The
    tail has to reach below ';' (Q26): a file of nothing but 'I' is valid
    phred+33 and phred+64 alike, and trimmomatic refuses to guess."""
    return "I" if i < 70 else "F" if i < 85 else ":" if i < 95 else "5"


def mutate(rng, read):
    """Substitutions and Ns, each given quality '#' (Q2): (seq, qual)."""
    bases = list(read)
    qual = [base_quality(i) for i in range(len(bases))]
    for i, b in enumerate(bases):
        if rng.random() < SUB_RATE:
            bases[i], qual[i] = rng.choice([c for c in "ACGT" if c != b]), "#"
    if rng.random() < N_READ_SHARE:
        i = rng.randrange(len(bases))
        bases[i], qual[i] = "N", "#"
    return "".join(bases), "".join(qual)


def to_read(rng, frag, adapter):
    """The first READ_LEN bases off one end, running into adapter (then A)
    when the fragment is shorter than a read."""
    seq = (frag + adapter + "A" * READ_LEN)[:READ_LEN]
    return mutate(rng, seq)  # (seq, qual)


def simulate(proteins_path, out, seed=11):
    rng = random.Random(seed)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)

    proteins = read_fasta(proteins_path)
    # One protein per BUSCO id (refseq_db carries several per id), 250-700 aa.
    by_id = {}
    for name in sorted(proteins):
        busco_id = name.split("_")[0]
        if busco_id not in by_id and 250 <= len(proteins[name].rstrip("*")) <= 700:
            by_id[busco_id] = proteins[name]
    ids = sorted(by_id)
    rng.shuffle(ids)
    need = sum(n for c, n, _ in CLASSES if c.startswith("busco"))
    if len(ids) < need:
        raise SystemExit(f"{proteins_path}: only {len(ids)} usable proteins, need {need}")

    transcripts = []  # (id, class, seq, coverage)
    for cls, n, cov in CLASSES:
        for i in range(n):
            if cls.startswith("busco"):
                bid = ids.pop()
                seq = random_seq(rng, 120) + back_translate(rng, by_id[bid]) + random_seq(rng, 180)
                tid = f"{cls}_{bid}"
            else:
                seq = no_long_orf(rng, rng.randint(1200, 2000))
                tid = f"{cls}_{i + 1}"
            transcripts.append((tid, cls, seq, cov))

    rows, reads = [], []
    for tid, cls, seq, cov in transcripts:
        pairs = max(1, int(cov * len(seq) / (2 * READ_LEN)))
        for _ in range(pairs):
            if rng.random() < SHORT_FRAG_SHARE:
                flen = rng.randint(50, READ_LEN - 5)
            else:
                flen = int(min(max(rng.gauss(FRAG_MEAN, FRAG_SD), READ_LEN + 10), len(seq)))
            start = rng.randint(0, len(seq) - flen)
            frag = seq[start:start + flen]
            # dUTP / RF: read 1 is the antisense strand, read 2 the sense.
            reads.append((to_read(rng, revcomp(frag), ADAPTER_R1),
                          to_read(rng, frag, ADAPTER_R2)))
        rows.append((tid, cls, len(seq), pairs))
    # Library order, not transcript order.
    rng.shuffle(reads)

    gz = lambda p: gzip.GzipFile(filename="", mode="wb", fileobj=open(p, "wb"), mtime=0)
    with gz(out / "sim_1.fq.gz") as f1, gz(out / "sim_2.fq.gz") as f2:
        for n, ((s1, q1), (s2, q2)) in enumerate(reads, start=1):
            f1.write(f"@sim{n}/1\n{s1}\n+\n{q1}\n".encode())
            f2.write(f"@sim{n}/2\n{s2}\n+\n{q2}\n".encode())
    n_pair = len(reads)

    # Expected TPM: pairs per effective base, normalised to a million.
    rates = {tid: pairs / max(1, length - FRAG_MEAN + 1) for tid, _, length, pairs in rows}
    total = sum(rates.values())
    tpm = {tid: r / total * 1e6 for tid, r in rates.items()}
    low = [tpm[t] for t, c, _, _ in rows if c.endswith("_low")]
    high = [tpm[t] for t, c, _, _ in rows if c.endswith("_high")]
    tpm_filt = round(math.sqrt(max(low) * min(high)))

    with open(out / "transcripts.fa", "w") as f:
        for tid, _, seq, _ in transcripts:
            f.write(f">{tid}\n{seq}\n")
    with open(out / "truth.tsv", "w") as f:
        f.write("id\tclass\tlength\tpairs\texpected_tpm\n")
        for tid, cls, length, pairs in rows:
            f.write(f"{tid}\t{cls}\t{length}\t{pairs}\t{tpm[tid]:.1f}\n")
    (out / "tpm_filt.txt").write_text(f"{tpm_filt}\n")
    return n_pair, tpm_filt


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--proteins", required=True, help="BUSCO ancestral protein fasta")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--seed", type=int, default=11)
    args = p.parse_args()
    pairs, tpm_filt = simulate(args.proteins, args.out, args.seed)
    print(f"{pairs} read pairs in {args.out}; --tpm-filt {tpm_filt}")


if __name__ == "__main__":
    main()
