#!/usr/bin/env python3
"""Side-by-side numbers for the chowder_arms.sbatch runs of one sample.

usage: compare_arms.py ARMS_DIR [ARM ...]

ARMS_DIR holds one run directory per arm (diamond_I12, blastn_I2, ...),
each a --no-cleanup chowder.py run with --runout <RUN>_<arm>. For each arm
this prints the contig count at every stage of the merge, so a change in the
final assembly can be traced to the stage that made it:

  pooled         contigs over 200 bp from all four assemblies (merged.fasta)
  orthogroups    groups in Orthogroups.txt
  picked         makeorthout's picks (good.<run>.list): one per group, score > 0
  rescued        contigs the diamond rescue added back (list7)
  after cd-hit   ORP.intermediate.fasta
  final          ORP.fasta, after the TPM filter

then BUSCO and transrate on the final assembly, and qualreport's UNIQUE
GENES ORP (swissprot genes hit, counted on ORP.intermediate.fasta).
"""
import csv
import json
import re
import sys
from pathlib import Path


def count_fasta(path):
    if not path.is_file():
        return None
    with open(path, "rb") as f:
        return sum(1 for line in f if line.startswith(b">"))


def count_lines(path):
    if not path or not path.is_file():
        return None
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def newest(paths):
    paths = [p for p in paths if p.is_file()]
    return max(paths, key=lambda p: p.stat().st_mtime) if paths else None


def busco(run_dir, run):
    summary = newest(run_dir.glob(f"reports/run_{run}.ORP/short_summary*.json"))
    if not summary:
        return {}
    data = json.load(open(summary))
    res = data.get("results", data)
    out = {}
    for key, label in (("Complete percentage", "C"), ("Single copy percentage", "S"),
                       ("Multi copy percentage", "D"), ("Fragmented percentage", "F"),
                       ("Missing percentage", "M")):
        if key in res:
            out[label] = float(res[key])
    if not out:  # older layout: parse the one-line string
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%",
                      res.get("one_line_summary", ""))
        if m:
            out = dict(zip("CSDFM", map(float, m.groups())))
    return out


def transrate(run_dir, run):
    path = run_dir / "reports" / f"transrate_{run}" / "assemblies.csv"
    if not path.is_file():
        return {}
    row = next(csv.DictReader(open(path)), {})
    keep = ("score", "optimal_score", "n_seqs", "n_bases", "mean_len", "n50",
            "p_good_mapping", "p_contigs_lowcovered", "mean_orf_percent")
    return {k: row[k] for k in keep if k in row}


def unique_genes(run_dir, run):
    """From qualreport: UNIQUE GENES ORP (distinct swissprot genes hit by
    diamond blastx of ORP.intermediate.fasta, after cd-hit and before the TPM
    filter) and READS MAPPED AS PROPER PAIRS (strandeval's bwa alignment of a
    400k-read subsample)."""
    path = run_dir / "reports" / f"qualreport.{run}"
    if not path.is_file():
        return {}
    text = path.read_text(errors="replace")
    out = {}
    m = re.search(r"UNIQUE GENES ORP\s*~*>\s*(\d+)", text)
    if m:
        out["unique genes"] = int(m.group(1))
    # strandeval's bwa alignment of a 400k-read subsample, samtools flagstat
    m = re.search(r"READS MAPPED AS PROPER PAIRS\s*~*>\s*([\d.]+)%", text)
    if m:
        out["proper pairs"] = float(m.group(1)) / 100
    return out


def arm_row(run_dir):
    run = next((p.name[:-len(".ORP.fasta")]
                for p in (run_dir / "assemblies").glob("*.ORP.fasta")), None)
    if run is None:
        return None
    of = run_dir / "orthofuse" / run
    groups = newest(of.glob("search/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt"))
    row = {
        "pooled": count_fasta(of / "merged.fasta"),
        "orthogroups": count_lines(groups),
        "picked": count_lines(of / f"good.{run}.list"),
        "rescued": count_lines(run_dir / "assemblies" / "diamond" / f"{run}.list7"),
        "after cd-hit": count_fasta(run_dir / "assemblies" / f"{run}.ORP.intermediate.fasta"),
        "final": count_fasta(run_dir / "assemblies" / f"{run}.ORP.fasta"),
    }
    row.update({f"BUSCO {k}": v for k, v in busco(run_dir, run).items()})
    row.update({f"transrate {k}": v for k, v in transrate(run_dir, run).items()})
    row.update(unique_genes(run_dir, run))
    return row


def main():
    arms_dir = Path(sys.argv[1])
    arms = sys.argv[2:] or sorted(p.name for p in arms_dir.iterdir() if p.is_dir())
    rows = {a: arm_row(arms_dir / a) for a in arms}
    rows = {a: r for a, r in rows.items() if r}
    if not rows:
        sys.exit(f"no finished arms under {arms_dir}")
    keys = list(dict.fromkeys(k for r in rows.values() for k in r))
    print(f"{'':26s}" + "".join(f"{a:>14s}" for a in rows))
    for k in keys:
        cells = []
        for r in rows.values():
            v = r.get(k)
            if isinstance(v, str):
                try:
                    v = float(v)
                except ValueError:
                    pass
            cells.append("-" if v is None else
                         f"{v:,.0f}" if isinstance(v, (int, float)) and abs(v) >= 100 else
                         f"{v:.3f}" if isinstance(v, float) and abs(v) < 1 else f"{v}")
        print(f"{k:26s}" + "".join(f"{c:>14s}" for c in cells))


if __name__ == "__main__":
    main()
