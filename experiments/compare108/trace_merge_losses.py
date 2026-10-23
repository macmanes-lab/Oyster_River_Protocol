#!/usr/bin/env python3
"""Trace what the ORP merge loses from its four component assemblies.

For each lineage_410 run, a BUSCO that is complete in at least one component
assembly (Trinity, rnaSPAdes auto-k, rnaSPAdes high-k, Trans-ABySS) but only
fragmented or missing in the merged ORP assembly is a merge loss; so is a
Swiss-Prot gene hit in a component but not in ORP. The merge is meant to be
lossless for gene content, so every loss is listed with the evidence for why
it happened:

  component contig(s)   the contigs BUSCO called complete in the component
  best_hit_gene         that contig's best DIAMOND hit (or "nohit"): a contig
                        with no hit goes down the no-hit path (CD-HIT-EST and
                        the --min-cov coverage filter), one with a hit is
                        grouped under that gene (representative + copy cap)
  orp_contigs_for_gene  how many ORP contigs have that gene as best hit
                        (3 = the per-gene cap was reached)
  in_orp                best match of the component contig in the ORP assembly
                        (minimap2 or blastn): identity and query coverage
  merge_files           ORP intermediate files that mention the contig id

and a cause label (nohit_path, copy_selection, gene_dropped, seq_present).
Lost genes get the same treatment plus whether the gene was ever any
component contig's BEST hit or only a secondary hit inside --top 0.1
(secondary_only: the merge groups by best hit, so such genes are unprotected).

DIAMOND is rerun on each of the five assemblies exactly as ORP runs it
(blastx -e 1e-8 --top 0.1, ORP's swissprot) and cached in the trace directory
(--use-run-hits takes the run's own tables instead, when it can find them). The
unique-gene counts recomputed from those hits are checked against the counts
in the run's qualreport, so the trace is known to describe the same numbers
the paper reports.

  # one run (array task):  trace_merge_losses.py --task N
  # all runs, after the array: trace_merge_losses.py --collate

Inputs (defaults under $COMPARE/lineage_410):
  manifest.tsv                                    pair ... SRR group lineage
  <SRR>/reports/qualreport.<run>[.done]           run name, unique genes
  <SRR>/reports/run_<run>.ORP/**/full_table.tsv   BUSCO, ORP
  <SRR>/assemblies/<run>.ORP.fasta[.gz]           ORP assembly
  <SRR>/assemblies/<run>.<asm>.fasta.gz           components (trinity: .trinity.Trinity.fasta.gz)
  components/<SRR>/<asm>/reports/run_<SRR>_<asm>.ORP/**/full_table.tsv
Outputs in merge_trace/<SRR>/: lost_buscos.tsv, lost_genes.tsv, summary.json
and, with --collate, merge_trace/merge_losses_summary.csv and
merge_losses_detail.tsv.
"""

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ASM = ["trinity", "spadesauto", "spadeshigh", "transabyss"]
COMPLETE = {"Complete", "Duplicated"}
TEXT_SUFFIXES = {".tsv", ".txt", ".clstr", ".log", ".csv", ".m8", ".out", ".list", ".ids"}
MAX_GREP_BYTES = 3 * 1024 ** 3


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def opener(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def run_name(rundir, srr):
    done = sorted((rundir / "reports").glob("qualreport.%s_v*.done" % srr))
    if len(done) != 1:
        return None
    return done[0].name[len("qualreport."):-len(".done")]


def qualreport_genes(path):
    out = {}
    if not path.is_file():
        return out
    text = path.read_text(errors="replace")
    for name in ["ORP"] + [a.upper() for a in ASM]:
        m = re.search(r"UNIQUE GENES %s\s*~*>\s*(\S+)" % name, text)
        if m:
            out[name.lower()] = m.group(1)
    return out


def busco_table(directory):
    """busco id -> (status, set of sequence ids) from BUSCO's full_table.tsv."""
    hits = sorted(directory.rglob("full_table.tsv")) if directory.is_dir() else []
    if not hits:
        return None, None
    table = {}
    for line in hits[0].open(errors="replace"):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        bid, status = f[0], f[1]
        seq = f[2] if len(f) > 2 else ""
        old_status, seqs = table.get(bid, (status, set()))
        if seq:
            seqs.add(seq)
        table[bid] = (old_status, seqs)
    return table, hits[0]


def fasta_ids(path):
    ids = set()
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def fasta_extract(path, wanted, out):
    """Write the records in `wanted` to `out`; return how many were written."""
    n, keep = 0, False
    with opener(path) as fh, open(str(out), "w") as oh:
        for line in fh:
            if line.startswith(">"):
                keep = line[1:].split()[0] in wanted
                n += keep
            if keep:
                oh.write(line)
    return n


def resolve(seqid, names):
    """BUSCO may decorate the contig id (coordinates, |-suffixes); map it back."""
    if seqid in names:
        return seqid
    for cut in (":", "|", " "):
        base = seqid.split(cut)[0]
        if base in names:
            return base
    base = re.sub(r"[:_]\d+-\d+$", "", seqid)
    return base if base in names else seqid


def gene_of(sseqid):
    # oyster.py's rule: third |-field, part before the "_" (sp|P12345|ACTB_HUMAN -> ACTB)
    f = sseqid.split("|")
    name = f[2] if len(f) >= 3 else sseqid
    return name.rsplit("_", 1)[0]


def read_hits(path):
    """contig -> (best gene, best bitscore, set of all hit genes) from a blast6 table."""
    best, genes = {}, defaultdict(set)
    with opener(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, g, bits = f[0], gene_of(f[1]), float(f[11])
            genes[q].add(g)
            if q not in best or bits > best[q][1]:
                best[q] = (g, bits)
    return {q: (best[q][0], best[q][1], genes[q]) for q in genes}


def looks_like_blast6(path):
    try:
        with opener(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                return len(f) >= 12 and "|" in f[1]
    except (OSError, UnicodeDecodeError, EOFError):
        return False
    return False


def find_hit_tables(rundir, run):
    """Existing DIAMOND tables under assemblies/, keyed by assembly label."""
    found = defaultdict(list)
    adir = rundir / "assemblies"
    if not adir.is_dir():
        return {}
    for p in adir.rglob("*"):
        n = p.name.lower()
        if not p.is_file() or n.endswith(".dmnd") or "diamond" not in n and not n.endswith((".m8", ".blastx")):
            continue
        if not looks_like_blast6(p):
            continue
        for label in ["orp"] + ASM:
            tag = ".orp." if label == "orp" else label
            if tag in n:
                found[label].append(p)
    # one table per label; several means the naming is not what we expect
    return {k: sorted(v, key=lambda x: x.stat().st_size)[-1] for k, v in found.items() if len(v) >= 1}


def diamond(fasta, db, out, cpu):
    """DIAMOND as ORP runs it, cached: out + '.done' marks a finished search (an
    empty table is a valid result). DIAMOND's own log is kept beside the table,
    and a run that exits 0 without writing output is reported with that log."""
    done = Path(str(out) + ".done")
    if out.is_file() and done.is_file():
        return out
    part = Path(str(out) + ".part.tsv")
    dlog = Path(str(out) + ".log")
    tmp = out.parent / ("tmp_" + out.stem)
    tmp.mkdir(exist_ok=True)
    cmd = ["diamond", "blastx", "-p", str(cpu), "-e", "1e-8", "--top", "0.1",
           "-q", str(fasta), "-d", str(db), "-o", str(part), "--tmpdir", str(tmp)]
    log("    diamond " + fasta.name)
    with open(str(dlog), "w") as lh:
        res = subprocess.run(cmd, stdout=lh, stderr=subprocess.STDOUT)
    tail = "".join(dlog.read_text(errors="replace").splitlines(True)[-15:])
    if res.returncode != 0:
        raise RuntimeError("diamond exited %d on %s:\n%s" % (res.returncode, fasta.name, tail))
    if not part.is_file():
        if re.search(r"Reported 0 pairwise alignments|0 queries aligned", tail):
            part.write_text("")
        else:
            raise RuntimeError("diamond exited 0 but wrote no output for %s:\n%s" % (fasta.name, tail))
    part.rename(out)
    done.write_text("")
    shutil.rmtree(str(tmp), ignore_errors=True)
    return out


def best_match_in_orp(queries, orp_fa, work, cpu):
    """query id -> (orp contig, identity %, query coverage %), via minimap2 or blastn."""
    out = {}
    if not queries.is_file() or queries.stat().st_size == 0:
        return out, "none"
    if shutil.which("minimap2"):
        paf = work / "lost_vs_orp.paf"
        with open(str(paf), "w") as fh:
            subprocess.run(["minimap2", "-c", "-x", "asm20", "-t", str(cpu), str(orp_fa), str(queries)],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        for line in paf.open():
            f = line.split("\t")
            q, qlen, qs, qe, t, match, alen = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[9]), int(f[10])
            cov, ident = 100.0 * (qe - qs) / qlen, 100.0 * match / max(alen, 1)
            if q not in out or cov * ident > out[q][1] * out[q][2]:
                out[q] = (t, round(ident, 2), round(cov, 2))
        return out, "minimap2"
    if shutil.which("blastn") and shutil.which("makeblastdb"):
        db = work / "orp_db"
        plain = work / "orp.fasta"
        if not Path(str(db) + ".nsq").exists() and not Path(str(db) + ".00.nsq").exists():
            with opener(orp_fa) as fh, open(str(plain), "w") as oh:
                shutil.copyfileobj(fh, oh)
            subprocess.run(["makeblastdb", "-in", str(plain), "-dbtype", "nucl", "-out", str(db)],
                           stdout=subprocess.DEVNULL, check=True)
        res = subprocess.run(["blastn", "-query", str(queries), "-db", str(db), "-num_threads", str(cpu),
                              "-max_target_seqs", "5", "-evalue", "1e-20",
                              "-outfmt", "6 qseqid sseqid pident length qlen"],
                             stdout=subprocess.PIPE, universal_newlines=True, check=True)
        for line in res.stdout.splitlines():
            q, t, pid, length, qlen = line.split("\t")
            cov = min(100.0, 100.0 * int(length) / int(qlen))
            if q not in out or cov * float(pid) > out[q][1] * out[q][2]:
                out[q] = (t, float(pid), round(cov, 2))
        if plain.exists():
            plain.unlink()
        return out, "blastn"
    return out, "none"


def merge_mentions(rundir, ids, work):
    """id -> list of ORP intermediate files (relative) that mention it."""
    mentions = defaultdict(list)
    if not ids:
        return mentions
    idfile = work / "lost_contig_ids.txt"
    idfile.write_text("\n".join(sorted(ids)) + "\n")
    for sub in ("assemblies", "reports"):
        root = rundir / sub
        if not root.is_dir():
            continue
        for p in root.rglob("*"):
            if not p.is_file() or p.suffix.lower() not in TEXT_SUFFIXES or p.stat().st_size > MAX_GREP_BYTES:
                continue
            res = subprocess.run(["grep", "-F", "-o", "-f", str(idfile), str(p)],
                                 stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
            for i in set(res.stdout.split()):
                if i in ids:
                    mentions[i].append(str(p.relative_to(rundir)))
    return mentions


def trace_one(line, args):
    pair, _orig, _r1, _r2, srr, group, lineage = line.split("\t")
    rundir = args.runs / srr
    run = run_name(rundir, srr)
    work = args.out / srr
    work.mkdir(parents=True, exist_ok=True)
    summary = {"pair_id": pair, "srr": srr, "group": group, "lineage": lineage, "run": run, "status": "ok"}
    if run is None:
        summary["status"] = "no finished run"
        (work / "summary.json").write_text(json.dumps(summary, indent=1))
        return summary
    log("=== %s (%s), run %s" % (pair, srr, run))

    orp_fa = next(iter(sorted((rundir / "assemblies").glob(run + ".ORP.fasta*"))), None)
    comp_fa = {a: rundir / "assemblies" / ("%s.trinity.Trinity.fasta.gz" % run if a == "trinity"
                                           else "%s.%s.fasta.gz" % (run, a)) for a in ASM}
    orp_busco, orp_tab = busco_table(rundir / "reports" / ("run_%s.ORP" % run))
    comp_busco = {}
    for a in ASM:
        t, _ = busco_table(args.runs / "components" / srr / a / "reports" / ("run_%s_%s.ORP" % (srr, a)))
        comp_busco[a] = t
    missing = [k for k, v in [("ORP assembly", orp_fa), ("ORP BUSCO table", orp_busco)] if v is None]
    missing += ["%s BUSCO table" % a for a in ASM if comp_busco[a] is None]
    missing += ["%s assembly" % a for a in ASM if not comp_fa[a].is_file()]
    if missing:
        summary["status"] = "missing: " + ", ".join(missing)
        log("    " + summary["status"])
        (work / "summary.json").write_text(json.dumps(summary, indent=1))
        return summary

    # DIAMOND hits: the run's own tables if they can be found, else rerun as ORP does.
    tables = find_hit_tables(rundir, run)
    hits, hit_source = {}, {}
    for label, fa in [("orp", orp_fa)] + [(a, comp_fa[a]) for a in ASM]:
        if label in tables and args.use_run_hits:
            path, hit_source[label] = tables[label], "run"
        else:
            path, hit_source[label] = diamond(fa, args.db, work / ("%s.diamond.txt" % label), args.cpu), "rerun"
        hits[label] = read_hits(path)
    genes = {k: set(g for v in h.values() for g in v[2]) for k, h in hits.items()}
    reported = qualreport_genes(rundir / "reports" / ("qualreport." + run))
    summary["hit_source"] = hit_source
    summary["unique_genes_recomputed"] = {k: len(v) for k, v in genes.items()}
    summary["unique_genes_qualreport"] = reported
    summary["gene_counts_match_qualreport"] = all(
        str(len(genes[k])) == reported.get(k, "") for k in genes if k in reported)

    orp_best_count = defaultdict(int)
    for q, (g, _b, _all) in hits["orp"].items():
        orp_best_count[g] += 1

    # --- BUSCO losses
    lost = []
    for bid, (ost, oseqs) in orp_busco.items():
        if ost in COMPLETE:
            continue
        donors = [(a, comp_busco[a][bid]) for a in ASM if comp_busco[a].get(bid, ("Missing",))[0] in COMPLETE]
        if donors:
            lost.append((bid, ost, oseqs, donors))
    gained = sum(1 for bid, (ost, _s) in orp_busco.items()
                 if ost in COMPLETE and not any(comp_busco[a].get(bid, ("Missing",))[0] in COMPLETE for a in ASM))
    union_complete = sum(1 for bid in orp_busco
                         if any(comp_busco[a].get(bid, ("Missing",))[0] in COMPLETE for a in ASM))

    # component contigs to look up, by assembler
    names = {a: None for a in ASM}
    wanted = defaultdict(set)
    for _bid, _ost, _oseqs, donors in lost:
        for a, (_st, seqs) in donors:
            if names[a] is None:
                names[a] = fasta_ids(comp_fa[a])
            for s in seqs:
                wanted[a].add(resolve(s, names[a]))

    # --- gene losses
    lost_genes = sorted(set().union(*[genes[a] for a in ASM]) - genes["orp"])
    gene_rows = []
    for g in lost_genes:
        where, best_in, contigs = [], [], []
        for a in ASM:
            cs = [q for q, (_bg, _b, allg) in hits[a].items() if g in allg]
            if not cs:
                continue
            where.append(a)
            as_best = [q for q in cs if hits[a][q][0] == g]
            if as_best:
                best_in.append(a)
            pick = sorted(as_best or cs, key=lambda q: -hits[a][q][1])[0]
            contigs.append("%s:%s" % (a, pick))
            if names[a] is None:
                names[a] = fasta_ids(comp_fa[a])
            wanted[a].add(pick)
        gene_rows.append({"gene": g, "in_components": ",".join(where),
                          "best_hit_in": ",".join(best_in) or "-",
                          "cause": "best_hit_dropped" if best_in else "secondary_only",
                          "example_contigs": ";".join(contigs)})

    # sequences of every component contig involved, matched against ORP
    qfa = work / "lost_contigs.fasta"
    with open(str(qfa), "w") as oh:
        for a in ASM:
            if wanted[a]:
                tmp = work / ("%s.tmp.fa" % a)
                fasta_extract(comp_fa[a], wanted[a], tmp)
                for line in tmp.open():
                    oh.write(">" + a + "::" + line[1:] if line.startswith(">") else line)
                tmp.unlink()
    match, matcher = best_match_in_orp(qfa, orp_fa, work, args.cpu)
    mentions = merge_mentions(rundir, set(i for a in ASM for i in wanted[a]), work)

    def in_orp(a, c):
        m = match.get("%s::%s" % (a, c))
        return ("%s %.1f%%id %.0f%%cov" % m) if m else "no match"

    def seq_present(a, c):
        m = match.get("%s::%s" % (a, c))
        return bool(m) and m[1] >= 95 and m[2] >= 90

    rows = []
    for bid, ost, oseqs, donors in lost:
        for a, (st, seqs) in donors:
            for s in sorted(seqs):
                c = resolve(s, names[a])
                bh = hits[a].get(c)
                gene = bh[0] if bh else "nohit"
                n_orp = orp_best_count.get(gene, 0) if bh else ""
                if seq_present(a, c):
                    cause = "seq_present"
                elif not bh:
                    cause = "nohit_path"
                elif n_orp == 0:
                    cause = "gene_dropped"
                else:
                    cause = "copy_selection"
                rows.append({"busco_id": bid, "orp_status": ost, "orp_seqs": ";".join(sorted(oseqs)) or "-",
                             "assembler": a, "comp_status": st, "comp_contig": c,
                             "best_hit_gene": gene, "orp_contigs_for_gene": n_orp,
                             "in_orp": in_orp(a, c), "cause": cause,
                             "merge_files": ";".join(sorted(set(mentions.get(c, [])))) or "-"})
    for r in gene_rows:
        ex = r["example_contigs"].split(";")[0]
        a, c = ex.split(":", 1)
        r["in_orp"] = in_orp(a, c)
        r["merge_files"] = ";".join(sorted(set(mentions.get(c, [])))) or "-"

    for name, data in (("lost_buscos.tsv", rows), ("lost_genes.tsv", gene_rows)):
        with open(str(work / name), "w", newline="") as fh:
            if data:
                w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t")
                w.writeheader()
                w.writerows(data)

    by_busco = {}
    for r in rows:  # one cause per lost BUSCO: the most informative across donors
        rank = ["gene_dropped", "nohit_path", "copy_selection", "seq_present"]
        old = by_busco.get(r["busco_id"])
        if old is None or rank.index(r["cause"]) < rank.index(old):
            by_busco[r["busco_id"]] = r["cause"]
    summary.update({
        "busco_n": len(orp_busco),
        "orp_complete": sum(1 for v in orp_busco.values() if v[0] in COMPLETE),
        "union_components_complete": union_complete,
        "lost_buscos": len(lost),
        "lost_buscos_fragmented_in_orp": sum(1 for x in lost if x[1] == "Fragmented"),
        "lost_buscos_missing_in_orp": sum(1 for x in lost if x[1] == "Missing"),
        "gained_buscos": gained,
        "lost_busco_causes": {k: sum(1 for v in by_busco.values() if v == k)
                              for k in ("gene_dropped", "nohit_path", "copy_selection", "seq_present")},
        "lost_genes": len(gene_rows),
        "lost_genes_secondary_only": sum(1 for r in gene_rows if r["cause"] == "secondary_only"),
        "lost_genes_best_hit_dropped": sum(1 for r in gene_rows if r["cause"] == "best_hit_dropped"),
        "matcher": matcher,
        "orp_busco_table": str(orp_tab),
    })
    (work / "summary.json").write_text(json.dumps(summary, indent=1))
    log("    lost %d BUSCOs (%s), %d genes (%d secondary-only); gene counts match qualreport: %s"
        % (summary["lost_buscos"], summary["lost_busco_causes"], summary["lost_genes"],
           summary["lost_genes_secondary_only"], summary["gene_counts_match_qualreport"]))
    return summary


def collate(args):
    rows, detail = [], []
    for line in args.manifest.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        pair, srr = f[0], f[4]
        sj = args.out / srr / "summary.json"
        if not sj.is_file():
            rows.append({"pair_id": pair, "srr": srr, "status": "not traced"})
            continue
        s = json.loads(sj.read_text())
        flat = {k: v for k, v in s.items() if not isinstance(v, dict)}
        for k in ("lost_busco_causes", "unique_genes_recomputed", "unique_genes_qualreport", "hit_source"):
            for kk, vv in (s.get(k) or {}).items():
                flat["%s_%s" % (k, kk)] = vv
        rows.append(flat)
        lb = args.out / srr / "lost_buscos.tsv"
        if lb.is_file() and lb.stat().st_size:
            for r in csv.DictReader(lb.open(), delimiter="\t"):
                r = dict(r)
                r.update({"pair_id": pair, "srr": srr, "kind": "busco"})
                detail.append(r)
        lg = args.out / srr / "lost_genes.tsv"
        if lg.is_file() and lg.stat().st_size:
            for r in csv.DictReader(lg.open(), delimiter="\t"):
                r = dict(r)
                r.update({"pair_id": pair, "srr": srr, "kind": "gene"})
                detail.append(r)
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(str(args.out / "merge_losses_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    dcols = []
    for r in detail:
        dcols += [k for k in r if k not in dcols]
    with open(str(args.out / "merge_losses_detail.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=dcols or ["pair_id"], delimiter="\t")
        w.writeheader()
        w.writerows(detail)
    traced = [r for r in rows if r.get("status") == "ok"]
    log("collated %d runs (%d traced): %d with a lost BUSCO, %d with a lost gene; %d lost BUSCOs, %d lost genes"
        % (len(rows), len(traced), sum(1 for r in traced if r.get("lost_buscos")),
           sum(1 for r in traced if r.get("lost_genes")),
           sum(r.get("lost_buscos", 0) for r in traced), sum(r.get("lost_genes", 0) for r in traced)))


def main():
    compare = Path(os.environ.get("COMPARE", "/mnt/home/macmaneslab/macmanes/compare"))
    orpdir = Path(os.environ.get("ORPDIR", str(Path.home() / "Oyster_River_Protocol")))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=Path, default=compare / "lineage_410")
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--db", type=Path, default=orpdir / "software" / "diamond" / "swissprot")
    ap.add_argument("--cpu", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "16")))
    ap.add_argument("--task", type=int, help="manifest line to trace (1-based)")
    ap.add_argument("--use-run-hits", action="store_true",
                    help="use DIAMOND tables found under the run's assemblies/ instead of rerunning "
                         "DIAMOND (faster, but only valid if their query ids are the component contig ids)")
    ap.add_argument("--collate", action="store_true")
    args = ap.parse_args()
    args.manifest = args.manifest or args.runs / "manifest.tsv"
    args.out = args.out or args.runs / "merge_trace"
    args.out.mkdir(parents=True, exist_ok=True)
    if args.collate:
        collate(args)
        return
    lines = [l for l in args.manifest.read_text().splitlines() if l.strip()]
    todo = [lines[args.task - 1]] if args.task else lines
    for line in todo:
        trace_one(line, args)


if __name__ == "__main__":
    main()
