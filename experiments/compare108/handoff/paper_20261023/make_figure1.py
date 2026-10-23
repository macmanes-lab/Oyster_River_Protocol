"""Supplementary File S2: draw Figure 1 from sheet Datasets of Supplementary Data S1.

Usage: python make_figure1.py ORP_4.0.0_lineage_results.xlsx Figure1
Writes Figure1.pdf and Figure1.png. Requires pandas, openpyxl and matplotlib.
"""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

xlsx, stem = sys.argv[1], sys.argv[2]
d = pd.read_excel(xlsx, sheet_name="Datasets")
V = "ORP 4.0.0"
panels = [("BUSCO complete (%)", "A", (0, 100)), ("BUSCO duplicated (%)", "B", (0, 100)),
          ("pyTransRate score", "C", (0, 0.7)), ("Good mappings (%)", "D", (0, 100))]
styles = {"TSA": dict(marker="o", color="#2b6cb0", label="TSA deposited (n = 44)"),
          "ORP 2018 (orthomerged)": dict(marker="s", color="#dd6b20", label="ORP 2018 (n = 15)")}

plt.rcParams.update({"font.family": "Arial", "font.size": 8})
fig, axes = plt.subplots(2, 2, figsize=(6.6, 6.4))
for ax, (metric, letter, lim) in zip(axes.flat, panels):
    for src, st in styles.items():
        s = d[d["Previous assembly"] == src]
        ax.scatter(s[f"{metric}: previous"], s[f"{metric}: {V}"], s=22, alpha=0.85, edgecolor="white",
                   linewidth=0.4, marker=st["marker"], color=st["color"], label=st["label"], zorder=3)
    ax.plot(lim, lim, ls="--", lw=0.8, color="0.5", zorder=1)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel(f"{metric}, previous assembly"); ax.set_ylabel(f"{metric}, {V}")
    ax.text(-0.2, 1.04, letter, transform=ax.transAxes, fontsize=11, fontweight="bold")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
axes[0, 0].legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(f"{stem}.pdf")
fig.savefig(f"{stem}.png", dpi=300)
