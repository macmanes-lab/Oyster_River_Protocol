"""Why genes found among a component's hits are missing from ORP's gene count.

Reads trace_merge_losses.py's collated output (merge_trace/merge_losses_detail.tsv)
and the DIAMOND tables it cached per run. For each lost gene, the first example
component contig is classed by its blastn match in the ORP assembly:
  masked      inside an ORP contig (>=95% id, >=90% of its length) whose best
              bitscore is more than 1/0.9 times the gene's: --top 0.1 hides it
  evalue      inside an ORP contig with no hit at all
  near_tie    inside an ORP contig, within 10%, but not reported
  seq_absent  no ORP match at that threshold (upper bound on sequence lost)
Writes merge_trace/gene_masking_by_run.csv (S1 sheet Gene_masking). Python 3.6.

  python3 explain_gene_masking.py      # on Premise, any directory
"""
import csv, re, sys, os, collections
T = "/mnt/home/macmaneslab/macmanes/compare/lineage_410/merge_trace"
def gene(s):
    p = s.split("|"); return p[2].split("_")[0] if len(p) > 2 else None
def load(path):
    d = collections.defaultdict(dict)  # contig -> {gene: best bitscore}
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if len(c) < 12: continue
        g = gene(c[1]); b = float(c[11])
        if g and b > d[c[0]].get(g, 0): d[c[0]][g] = b
    return d
rows = collections.defaultdict(list)
with open(T + "/merge_losses_detail.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["kind"] == "gene": rows[r["srr"]].append(r)
tot = collections.Counter(); ex = collections.defaultdict(list); per = collections.defaultdict(collections.Counter); pair = {}
for srr, rs in sorted(rows.items()):
    tabs = {}
    def tab(l):
        if l not in tabs: tabs[l] = load(os.path.join(T, srr, l + ".diamond.txt"))
        return tabs[l]
    for r in rs:
        m = re.match(r"(\S+) ([\d.]+)%id (\d+)%cov", r["in_orp"])
        if not m or float(m.group(2)) < 95 or int(m.group(3)) < 90:
            tot["sequence not in ORP"] += 1; per[srr]["seq_absent"] += 1; pair[srr] = r["pair_id"]; continue
        oc = m.group(1)
        a, cc = r["example_contigs"].split(";")[0].split(":", 1)
        cg = tab(a).get(cc, {}).get(r["gene"])
        oh = tab("orp").get(oc, {})
        if not oh: k = "ORP contig has no hit at all"
        elif cg is None: k = "example contig lacks gene (?)"
        elif cg < 0.9 * max(oh.values()): k = "masked: gene bitscore < 90% of ORP contig's best (--top 0.1)"
        else: k = "other"
        tot[k] += 1; pair[srr] = r["pair_id"]
        per[srr][{"ORP contig has no hit at all": "evalue", "other": "near_tie"}.get(k, "masked" if k.startswith("masked") else "other")] += 1
        if len(ex[k]) < 3: ex[k].append((srr, r["gene"], cc, oc, cg, max(oh.values()) if oh else None, max(oh, key=oh.get) if oh else None))
for k, n in tot.most_common(): print(n, k)
for k, v in ex.items():
    print("--", k)
    for e in v: print("  ", e)

with open(T + "/gene_masking_by_run.csv", "w") as out:
    out.write("pair_id,srr,lost_gene_labels,masked_top01,evalue_longer_query,near_tie_unexplained,sequence_absent\n")
    for srr in sorted(per):
        c = per[srr]
        out.write("%s,%s,%d,%d,%d,%d,%d\n" % (pair[srr], srr, sum(c.values()), c["masked"], c["evalue"], c["near_tie"], c["seq_absent"]))
