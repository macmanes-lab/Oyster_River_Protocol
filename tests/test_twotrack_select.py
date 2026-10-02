#!/usr/bin/env python3
"""Test scripts/twotrack_select.py on a small synthetic pool with a known answer.

usage: python tests/test_twotrack_select.py

Needs blastn, makeblastdb and cd-hit-est on PATH (all three are in the orp
conda env): `conda run -n orp python tests/test_twotrack_select.py`. Skips,
exit 0, when they are missing. Exits non-zero on a failure.

The pool, with the rule each contig tests:

  gene A (protein 300 aa)
    A1  1500 bp, hit covers 97% of the protein   step 1: longest near-best -> representative
    A2  1000 bp, part of A1, 97% coverage        near-best but shorter -> not the representative;
                                                 and contained in A1 -> not rescued
    A3   900 bp, unrelated sequence, 50%, 5 TPM  step 2: distinct + expressed -> rescued
    A4   800 bp, part of A1, 10 TPM              step 2: contained in A1 -> dropped
    A5   700 bp, unrelated, 0.2 TPM              step 2: below 1 TPM -> dropped
    A6   600 bp, unrelated, 3 TPM                step 2: distinct + expressed, but the cap of
    A7   600 bp, unrelated, 4 TPM                3 per gene keeps A1, A3 (5 TPM), A7 (4 TPM)
  gene B (protein 200 aa)
    B1   600 bp, 90% coverage                    step 1: the only near-best (>= 0.9 x 0.9)
    B2  1200 bp, 60% coverage, 0.5 TPM           longer, but not near-best -> not the
                                                 representative; below 1 TPM -> not rescued
  no swissprot hit
    N1  1000 bp                                  track 2: kept
    N2   950 bp, part of N1                      track 2: >= 95% id over >= 90% -> collapsed
    N3   500 bp, unrelated                       track 2: kept

Expected kept: A1, A3, A7, B1, N1, N3.
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "twotrack_select.py")
EXPECTED = {"A1", "A3", "A7", "B1", "N1", "N3"}


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, s in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")


def build(tmp):
    rng = random.Random(7)
    a1 = rand_seq(rng, 1500)
    n1 = rand_seq(rng, 1000)
    pool = {
        "A1": a1, "A2": a1[200:1200], "A3": rand_seq(rng, 900), "A4": a1[600:1400],
        "A5": rand_seq(rng, 700), "A6": rand_seq(rng, 600), "A7": rand_seq(rng, 600),
        "B1": rand_seq(rng, 600), "B2": rand_seq(rng, 1200),
        "N1": n1, "N2": n1[25:975], "N3": rand_seq(rng, 500),
    }
    write_fasta(os.path.join(tmp, "merged.fasta"), pool)
    write_fasta(os.path.join(tmp, "sprot.fasta"),
                {"sp|P00001|GENEA_HUMAN": "M" * 300, "sp|P00002|GENEB_HUMAN": "M" * 200})
    # contig -> (subject, sstart, send, bitscore), sstart/send in protein coordinates
    hits = {
        "A1": ("sp|P00001|GENEA_HUMAN", 1, 291, 500), "A2": ("sp|P00001|GENEA_HUMAN", 5, 295, 480),
        "A3": ("sp|P00001|GENEA_HUMAN", 1, 150, 200), "A4": ("sp|P00001|GENEA_HUMAN", 1, 200, 300),
        "A5": ("sp|P00001|GENEA_HUMAN", 1, 100, 150), "A6": ("sp|P00001|GENEA_HUMAN", 1, 100, 150),
        "A7": ("sp|P00001|GENEA_HUMAN", 1, 100, 150),
        "B1": ("sp|P00002|GENEB_HUMAN", 1, 180, 300), "B2": ("sp|P00002|GENEB_HUMAN", 1, 120, 200),
    }
    with open(os.path.join(tmp, "a.diamond.txt"), "w") as f:
        for c, (s, ss, se, b) in hits.items():
            f.write("\t".join(map(str, [c, s, 90.0, se - ss + 1, 0, 0, 1, 3 * (se - ss + 1),
                                        ss, se, 1e-50, b])) + "\n")
    tpm = {"A1": 20, "A2": 10, "A3": 5, "A4": 10, "A5": 0.2, "A6": 3, "A7": 4,
           "B1": 8, "B2": 0.5, "N1": 6, "N2": 2, "N3": 3}
    with open(os.path.join(tmp, "contigs.csv"), "w") as f:
        f.write("contig_name,score,tpm\n")
        for c in pool:
            f.write(f"{c},0.5,{tpm[c]}\n")


def run(tmp, with_sprot=True):
    out = os.path.join(tmp, "good.list")
    cmd = [sys.executable, SCRIPT, "--merged", os.path.join(tmp, "merged.fasta"),
           "--contigs-csv", os.path.join(tmp, "contigs.csv"),
           "--diamond", os.path.join(tmp, "a.diamond.txt"),
           "--table", os.path.join(tmp, "table.tsv"), "--out", out]
    if with_sprot:
        cmd += ["--sprot", os.path.join(tmp, "sprot.fasta")]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE)
    return {l.strip() for l in open(out) if l.strip()}


def main():
    missing = [t for t in ("blastn", "makeblastdb", "cd-hit-est") if not shutil.which(t)]
    if missing:
        print(f"SKIP: {', '.join(missing)} not on PATH (try: conda run -n orp python {sys.argv[0]})")
        return 0
    failures = 0
    with tempfile.TemporaryDirectory() as tmp:
        build(tmp)
        got = run(tmp)
        if got == EXPECTED:
            print(f"PASS  kept {sorted(got)}")
        else:
            failures += 1
            print(f"FAIL  kept {sorted(got)}\n      expected {sorted(EXPECTED)}\n"
                  f"      extra {sorted(got - EXPECTED)}, missing {sorted(EXPECTED - got)}")
        # Without the swissprot fasta the selector ranks by aligned length; it
        # must still keep one representative per gene and run to the end.
        got = run(tmp, with_sprot=False)
        reps = {l.split("\t")[0] for l in open(os.path.join(tmp, "table.tsv"))
                if "\trepresentative\t" in l}
        if len(reps) == 2 and reps <= got:
            print(f"PASS  no-sprot fallback: representatives {sorted(reps)}")
        else:
            failures += 1
            print(f"FAIL  no-sprot fallback: representatives {sorted(reps)}, kept {sorted(got)}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
