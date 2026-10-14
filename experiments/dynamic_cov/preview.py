#!/usr/bin/env python3
"""What each filter would remove, before running BUSCO or pyTransRate on it.

usage: conda run -n orp python preview.py ROOT NAME [NAME ...]

For ROOT/<NAME>_cov (a sweep.sbatch run directory, or anything laid out like
one): of the intermediate's contigs without a swissprot hit, how many each
rule removes, how many of those have an ORF over 149 codons, and the
ORF-span shares collect.py reports, observed and shuffled. Needs only
assemblies/*.ORP.intermediate.fasta, assemblies/*.ORP.diamond.txt,
reports/*.quant.sf.gz and the '# frag_length_mean' line of
reports/*.filter_candidates.tsv.
"""

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import collect  # noqa: E402

RULES = [("tpm0.5", "tpm", 0.5), ("tpm1", "tpm", 1), ("cov0.5", "cov", 0.5), ("cov1", "cov", 1),
         ("cov2", "cov", 2), ("cov5", "cov", 5)]


def main():
    root, names = Path(sys.argv[1]), sys.argv[2:]
    print("dataset\trule\tnohit\tnohit_removed\tremoved_orf150_n\tremoved_orf80\tremoved_orf80_shuf\t"
          "kept_orf80\tkept_orf80_shuf")
    for name in names:
        ds = root / f"{name}_cov"
        orfs = collect.nohit_orfs(ds)
        cand = next(ds.glob("reports/*.filter_candidates.tsv"))
        fl = float(next(l.split("\t")[1] for l in open(cand) if l.startswith("# frag_length_mean")))
        q = {}
        with gzip.open(next(ds.glob("reports/*.quant.sf.gz")), "rt") as f:
            next(f)
            for line in f:
                n, length, _, tpm, reads = line.rstrip("\n").split("\t")
                if n in orfs:
                    q[n] = (float(tpm), float(reads) * fl / float(length))
        for label, kind, x in RULES:
            kept = {n for n, (tpm, cov) in q.items() if (tpm if kind == "tpm" else cov) >= x}
            c = collect.orf_columns(orfs, kept)
            long_removed = sum(1 for n in orfs if n not in kept and orfs[n][1] > collect.LONG_ORF_CODONS)
            print(f"{name}\t{label}\t{len(orfs)}\t{c['nohit_removed']}\t{long_removed}\t"
                  f"{c['nohit_removed_orf80']}\t{c['nohit_removed_orf80_shuf']}\t"
                  f"{c['nohit_kept_orf80']}\t{c['nohit_kept_orf80_shuf']}")


if __name__ == "__main__":
    main()
