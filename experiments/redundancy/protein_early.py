#!/usr/bin/env python3
"""Protein pick rule vs each source arm and vs the control, on named clean runs.

usage: protein_early.py VALIDATE_DIR SAMPLE:ARM[,ARM...] [...]

For each sample and source arm (control, candidate, blastn_I12) whose
<arm>_protein re-pick is listed as clean, prints both, then totals: protein
against its own source arm (same orthogroups, old pick rule) and against the
control, in BUSCO genes, transrate score, unique genes and strong-hit
(bitscore >= 200) swissprot genes lost/gained against the control.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_arms import arm_row  # noqa: E402
from gene_sets import arm_genes  # noqa: E402

LABEL = {"control": "diamond I12", "candidate": "blastn I3", "blastn_I12": "blastn I12"}
OLD = {"control": "score", "candidate": "score_orf", "blastn_I12": "score"}


def get(d):
    r = arm_row(d)
    if not r or "BUSCO C" not in r or "transrate score" not in r:
        return None
    b = lambda k: round(float(r[k]) * 1.25)
    return dict(C=b("BUSCO C"), D=b("BUSCO D"), M=b("BUSCO M"),
                score=float(r["transrate score"]), uniq=r.get("unique genes", 0),
                contigs=r["final"])


def main():
    vdir = Path(sys.argv[1])
    clean = {}
    for spec in sys.argv[2:]:
        s, arms = spec.split(":")
        clean[s] = arms.split(",")
    print(f"{'sample':11s}{'orthogroups':13s}{'rule':11s}{'C':>5s}{'D':>5s}{'M':>5s}"
          f"{'score':>8s}{'uniq':>8s}{'contigs':>9s}")
    tot = {}
    for s, arms in clean.items():
        ctrl = get(vdir / s / "control")
        gc = arm_genes(vdir / s / "control")
        for a in arms:
            base, prot = get(vdir / s / a), get(vdir / s / f"{a}_protein")
            if not (ctrl and base and prot):
                print(f"{s:11s}{LABEL[a]:13s}incomplete, skipped")
                continue
            for rule, r in ((OLD[a], base), ("protein", prot)):
                print(f"{s:11s}{LABEL[a]:13s}{rule:11s}{r['C']:5d}{r['D']:5d}{r['M']:5d}"
                      f"{r['score']:8.3f}{r['uniq']:8d}{r['contigs']:9,d}")
            gp = arm_genes(vdir / s / f"{a}_protein")
            t = tot.setdefault(a, {k: 0 for k in
                                   ("n", "dC", "dD", "dM", "dS", "dU", "cC", "cD", "cM", "cS",
                                    "cU", "sl", "sg")})
            t["n"] += 1
            for k in "CDM":
                t["d" + k] += prot[k] - base[k]
                t["c" + k] += prot[k] - ctrl[k]
            t["dS"] += prot["score"] - base["score"]
            t["cS"] += prot["score"] - ctrl["score"]
            t["dU"] += prot["uniq"] - base["uniq"]
            t["cU"] += prot["uniq"] - ctrl["uniq"]
            t["sl"] += sum(1 for g in gc if g not in gp and gc[g] >= 200)
            t["sg"] += sum(1 for g in gp if g not in gc and gp[g] >= 200)
    print()
    for a, t in tot.items():
        n = t["n"]
        print(f"{LABEL[a]} orthogroups, {n} samples")
        print(f"  protein vs {OLD[a]:9s}: complete {t['dC']:+d}, duplicated {t['dD']:+d}, "
              f"missing {t['dM']:+d}, transrate {t['dS'] / n:+.3f} mean, unique genes {t['dU']:+d}")
        print(f"  protein vs control  : complete {t['cC']:+d}, duplicated {t['cD']:+d}, "
              f"missing {t['cM']:+d}, transrate {t['cS'] / n:+.3f} mean, unique genes {t['cU']:+d}, "
              f"strong genes lost {t['sl']} / gained {t['sg']}")


if __name__ == "__main__":
    main()
