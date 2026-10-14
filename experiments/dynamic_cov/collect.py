#!/usr/bin/env python3
"""One table of the filter sweep: every dataset x level that sweep.sbatch finished.

usage: conda run -n orp python collect.py ROOT [--tsv OUT]

ROOT holds one directory per dataset (sweep.sbatch's $DIR), each with
levels/<level>/. Prints, per dataset and level: contigs and bases in the
final assembly, the share of contigs with no swissprot hit kept by the
filter, BUSCO, unique genes, pyTransRate's score, optimal score and its
mapping metrics, and proper-pair mapping.

ORF columns. BUSCO and unique genes only see contigs with a swissprot hit,
which no filter removes, so they cannot tell the filters apart. These look
at the contigs without a hit instead:

  n_with_orf, mean_orf_percent  pyTransRate's, over the whole final
      assembly. n_with_orf counts ORFs over 149 codons, which no contig
      under ~450 bp can have, so it mostly follows length.
  nohit_kept_orf80, nohit_removed_orf80  of the no-hit contigs the filter
      kept / removed, the share whose longest ORF spans >= 80% of the
      contig. A fragment of coding sequence is open nearly end to end at any
      length; noise of a few hundred bases usually hits a stop.
  nohit_removed_orf150  of those removed, the share with an ORF over 149
      codons (pyTransRate's n_with_orf rule): candidate real genes that are
      just not in swissprot.

ORF lengths are pyTransRate's own (pytransrate.sequence.longest_orf, so run
this in the orp env), computed once per dataset over the intermediate's
no-hit contigs and cached in <dataset>/nohit_orfs.tsv.
"""

import argparse
import csv
import re
from pathlib import Path

LEVEL_ORDER = ["none", "tpm1", "cov0.5", "cov1", "cov2", "cov5", "cov10"]

QUAL = {
    "busco": r"BUSCO SCORE ~+>\s+(\S+)",
    "score": r"PYTRANSRATE SCORE ~+>\s+(\S+)",
    "optimal": r"PYTRANSRATE OPTIMAL SCORE ~+>\s+(\S+)",
    "unique_genes": r"UNIQUE GENES ORP ~+>\s+(\S+)",
    "proper_pairs": r"READS MAPPED AS PROPER PAIRS ~+>\s+(\S+)",
}
PYTR = ["n_seqs", "n_bases", "n50", "p_fragments_mapped", "p_good_mapping",
        "p_contigs_uncovered", "p_contigs_lowcovered", "p_contigs_segmented"]
ORF_COLUMNS = ["nohit_kept_orf80", "nohit_removed", "nohit_removed_orf80", "nohit_removed_orf150"]
COLUMNS = (["dataset", "level", "nohit_kept", "C", "S", "D", "F", "M", "unique_genes",
            "score", "optimal"] + PYTR + ["n_with_orf", "mean_orf_percent", "proper_pairs"]
           + ORF_COLUMNS)

#: The share of a contig its longest ORF must span to count as open.
ORF_SPAN = 0.8
#: pyTransRate's n_with_orf rule: an ORF over this many codons.
LONG_ORF_CODONS = 149


def level_key(name):
    return (LEVEL_ORDER.index(name) if name in LEVEL_ORDER else len(LEVEL_ORDER), name)


def read_level(d):
    row = {}
    qual = next(d.glob("qualreport.*"), None)
    if qual:
        text = qual.read_text()
        for key, pat in QUAL.items():
            m = re.search(pat, text)
            row[key] = m.group(1) if m else ""
        m = re.match(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", row.pop("busco"))
        if m:
            row.update(zip("CSDFM", m.groups()))
    pytr = d / "pytransrate.assemblies.csv"
    if pytr.exists():
        with open(pytr) as f:
            rec = next(csv.DictReader(f))
        row.update({k: rec.get(k, "") for k in PYTR + ["n_with_orf", "mean_orf_percent"]})
    cand = next(d.glob("*.filter_candidates.tsv"), None)
    if cand:
        text = cand.read_text()
        total = re.search(r"^# contigs\t(\d+)", text, re.M)
        hits = re.search(r"^# with_swissprot_hit\t(\d+)", text, re.M)
        applied = re.search(r"^applied:\S+\t\d+\t(\d+)\t", text, re.M)
        if total and hits and applied:
            nohit = int(total.group(1)) - int(hits.group(1))
            row["nohit_kept"] = f"{applied.group(1)}/{nohit}"
    return row


def fasta_records(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def fasta_ids(path):
    with open(path) as f:
        return {line[1:].split()[0] for line in f if line.startswith(">")}


def nohit_orfs(ds):
    """{contig: (length, longest ORF in codons)} for the intermediate's no-hit contigs."""
    cache = ds / "nohit_orfs.tsv"
    if cache.exists():
        with open(cache) as f:
            return {n: (int(l), int(o)) for n, l, o in (line.split() for line in f)}
    inter = next((ds / "assemblies").glob("*.ORP.intermediate.fasta"), None)
    diamond = next((ds / "assemblies").glob("*.ORP.diamond.txt"), None)
    if inter is None or diamond is None:
        return None
    from pytransrate.sequence import longest_orf
    with open(diamond) as f:
        hits = {line.split("\t", 1)[0] for line in f if line.strip()}
    orfs = {n: (len(s), longest_orf(s.upper())) for n, s in fasta_records(inter) if n not in hits}
    part = cache.with_suffix(".part")
    with open(part, "w") as f:
        for n, (length, orf) in orfs.items():
            f.write(f"{n}\t{length}\t{orf}\n")
    part.rename(cache)
    return orfs


def orf_columns(orfs, kept_ids):
    def share(names, test):
        return f"{sum(1 for n in names if test(*orfs[n])) / len(names):.3f}" if names else ""
    kept = [n for n in orfs if n in kept_ids]
    removed = [n for n in orfs if n not in kept_ids]
    is_open = lambda length, orf: orf * 3 >= ORF_SPAN * length
    is_long = lambda length, orf: orf > LONG_ORF_CODONS
    return {
        "nohit_kept_orf80": share(kept, is_open),
        "nohit_removed": len(removed),
        "nohit_removed_orf80": share(removed, is_open),
        "nohit_removed_orf150": share(removed, is_long),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root", type=Path)
    ap.add_argument("--tsv", type=Path)
    args = ap.parse_args()
    rows = []
    for ds in sorted(p for p in args.root.iterdir() if (p / "levels").is_dir()):
        levels = sorted((p for p in (ds / "levels").iterdir() if p.is_dir()), key=lambda p: level_key(p.name))
        orfs = nohit_orfs(ds) if levels else None
        for lv in levels:
            row = {"dataset": ds.name, "level": lv.name, **read_level(lv)}
            fasta = next(lv.glob("*.ORP.fasta"), None)
            if orfs is not None and fasta is not None:
                row.update(orf_columns(orfs, fasta_ids(fasta)))
            rows.append(row)
    lines = ["\t".join(COLUMNS)] + ["\t".join(str(r.get(c, "")) for c in COLUMNS) for r in rows]
    if args.tsv:
        args.tsv.write_text("\n".join(lines) + "\n")
    widths = [max(len(l.split("\t")[i]) for l in lines) for i in range(len(COLUMNS))]
    for l in lines:
        print("  ".join(c.ljust(w) for c, w in zip(l.split("\t"), widths)))


if __name__ == "__main__":
    main()
