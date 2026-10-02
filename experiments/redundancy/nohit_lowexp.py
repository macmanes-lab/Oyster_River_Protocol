#!/usr/bin/env python3
"""How many contigs have no swissprot hit and low expression, at each stage?

usage: nohit_lowexp.py VALIDATE_DIR SAMPLES_FILE [--arm control]

For each sample's arm (a --no-cleanup chowder run), splits contigs by whether
they have a swissprot hit and by expression:

  picks   the orthogroup picks (good.<run>.list). Hit = a hit in the
          per-assembly diamond blastx output. TPM = pytransrate's salmon
          estimate on the whole pool, split across redundant copies, so it
          understates expression here.
  final   ORP.fasta, after cd-hit-est and ORP's TPM filter. Hit = a hit in
          ORP.diamond.txt (blastx of ORP.intermediate.fasta). TPM = ORP's own
          salmon run on ORP.intermediate.fasta, the number the TPM filter
          used.
"""
import argparse
import csv
import glob
from pathlib import Path


def ids_with_hits(paths):
    out = set()
    for p in paths:
        with open(p) as f:
            for line in f:
                out.add(line.split("\t", 1)[0])
    return out


def fasta_ids(path):
    with open(path) as f:
        return [l[1:].split()[0] for l in f if l.startswith(">")]


def bins(tpms):
    edges = [(0, 0.1), (0.1, 1), (1, 5), (5, float("inf"))]
    return [sum(1 for t in tpms if lo <= t < hi) for lo, hi in edges]


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("validate_dir")
    p.add_argument("samples_file")
    p.add_argument("--arm", default="control")
    args = p.parse_args()

    print(f"{'sample':11s}{'stage':7s}{'contigs':>9s}{'no hit':>9s}"
          f"{'nohit TPM<0.1':>15s}{'0.1-1':>8s}{'1-5':>8s}{'>=5':>8s}{'hit TPM<1':>11s}")
    for run in [l.strip() for l in open(args.samples_file) if l.strip()]:
        d = Path(args.validate_dir) / run / args.arm
        good = glob.glob(f"{d}/orthofuse/*/good.*.list")
        if not good:
            continue
        r = Path(good[0]).name[len("good."):-len(".list")]
        # Picks, against the pool's diamond hits and pytransrate TPM.
        hits = ids_with_hits(p for p in glob.glob(f"{d}/assemblies/diamond/{r}.*.diamond.txt")
                             if "orthomerged" not in p)
        pool_tpm = {}
        with open(f"{d}/orthofuse/{r}/merged/contigs.csv", newline="") as f:
            for row in csv.DictReader(f):
                pool_tpm[row["contig_name"]] = float(row["tpm"])
        picks = [l.strip() for l in open(good[0]) if l.strip()]
        nohit = [c for c in picks if c not in hits]
        b = bins([pool_tpm.get(c, 0.0) for c in nohit])
        hl = sum(1 for c in picks if c in hits and pool_tpm.get(c, 0.0) < 1)
        print(f"{run:11s}{'picks':7s}{len(picks):9,d}{len(nohit):9,d}{b[0]:15,d}{b[1]:8,d}{b[2]:8,d}{b[3]:8,d}{hl:11,d}")
        # Final assembly, against ORP's own diamond and salmon.
        orp = glob.glob(f"{d}/assemblies/{r}.ORP.fasta")
        quant = glob.glob(f"{d}/quants/salmon_orthomerged_{r}/quant.sf")
        dtxt = glob.glob(f"{d}/assemblies/{r}.ORP.diamond.txt")
        if not (orp and quant and dtxt):
            continue
        fhits = ids_with_hits(dtxt)
        tpm = {}
        with open(quant[0]) as f:
            next(f)
            for line in f:
                c = line.split("\t")
                tpm[c[0]] = float(c[3])
        final = fasta_ids(orp[0])
        fnohit = [c for c in final if c not in fhits]
        b = bins([tpm.get(c, 0.0) for c in fnohit])
        hl = sum(1 for c in final if c in fhits and tpm.get(c, 0.0) < 1)
        print(f"{'':11s}{'final':7s}{len(final):9,d}{len(fnohit):9,d}{b[0]:15,d}{b[1]:8,d}{b[2]:8,d}{b[3]:8,d}{hl:11,d}")


if __name__ == "__main__":
    main()
