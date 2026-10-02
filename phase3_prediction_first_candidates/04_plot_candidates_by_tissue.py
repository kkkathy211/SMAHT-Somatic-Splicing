# -*- coding: utf-8 -*-
# ============================================================================
# Phase 3 - Prediction-first candidate selection
#
# Candidate counts per tissue. Hatched bars = no LeafCutter perind, i.e. not testable.
#
# Original location in the analysis project:
#     0720/code/20_candidates_by_tissue_barplot.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
20_candidates_by_tissue_barplot.py
===================================
Bar plot of VEP-splice + SpliceAI≥0.2 candidate counts per tissue.

Design:
  - X: tissues sorted by count desc
  - Y: candidate count
  - Bar color: official SMaHT tissue color
  - Bar hatch: solid = LeafCutter available, hatched = not
  - Above bar: total count + IMPACT breakdown (H/M/L)
  - Below x-axis: ✓ / ✗ for LeafCutter perind availability
"""

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ── Official SMaHT tissue color palette ──────────────────────────────────────
TISSUE = {
    "3A":   {"smaht":"BLOO","full":"Whole blood",          "hex":"#FF382E"},
    "3B":   {"smaht":"BUCC","full":"Buccal swab",          "hex":"#75001F"},
    "3C":   {"smaht":"ESOP","full":"Esophagus",            "hex":"#FFB161"},
    "3E":   {"smaht":"COAS","full":"Colon-Ascending",      "hex":"#C0B5E7"},
    "3G":   {"smaht":"CODS","full":"Colon-Descending",     "hex":"#9D4BBF"},
    "3I":   {"smaht":"LIVR","full":"Liver",                "hex":"#FF4AA9"},
    "3K":   {"smaht":"ADGL","full":"Adrenal gland, L",     "hex":"#EA7A00"},
    "3M":   {"smaht":"ADGR","full":"Adrenal gland, R",     "hex":"#EA7A00"},
    "3O":   {"smaht":"AORT","full":"Aorta",                "hex":"#F99B9B"},
    "3Q":   {"smaht":"LUNG","full":"Lung",                 "hex":"#1DDD30"},
    "3S":   {"smaht":"HART","full":"Heart",                "hex":"#DDD824"},
    "3U":   {"smaht":"TESL","full":"Testis, L",            "hex":"#AAF8B4"},
    "3W":   {"smaht":"TESR","full":"Testis, R",            "hex":"#AAF8B4"},
    "3Y":   {"smaht":"OVAL","full":"Ovary, L",             "hex":"#FFAFD4"},
    "3AA":  {"smaht":"OVAR","full":"Ovary, R",             "hex":"#FFAFD4"},
    "3AC":  {"smaht":"FBRO","full":"Fibroblast",           "hex":"#FFE98D"},
    "3AD":  {"smaht":"SKSE","full":"Skin-Calf",            "hex":"#0EC1B8"},
    "3AF":  {"smaht":"SKNE","full":"Skin-Abdomen",         "hex":"#91E5DB"},
    "3AH":  {"smaht":"MUSC","full":"Skeletal muscle",      "hex":"#AA9A41"},
    "3AK":  {"smaht":"BRFL","full":"Brain-Frontal lobe",   "hex":"#BBE7FF"},
    "3AL":  {"smaht":"BRTL","full":"Brain-Temporal lobe",  "hex":"#76AEFF"},
    "3AM":  {"smaht":"BRCE","full":"Brain-Cerebellum",     "hex":"#1655C4"},
    "3AN":  {"smaht":"BRHL","full":"Brain-Hippocampus, L", "hex":"#002A66"},
    "3AO":  {"smaht":"BRHR","full":"Brain-Hippocampus, R", "hex":"#002A66"},
}
# Tissues with LeafCutter perind files (Laura's re-clustered set)
LEAFCUTTER_TISSUES = {"3AD","3AF","3AH","3AK","3AL","3AM","3C","3G","3Q","3S","3U"}

INPUT_CSV = "/Users/spectremac/Desktop/Differential_splicing/0720/output/vep_splice_candidates_v47_enriched.csv"
OUT_PNG   = "/Users/spectremac/Desktop/Differential_splicing/0720/output/candidates_by_tissue_barplot_SA0p5.png"
OUT_PDF   = OUT_PNG.replace(".png", ".pdf")

SPLICEAI_MIN = 0.5   # tightened from 0.2

plt.rcParams.update({
    "font.family":  "DejaVu Sans",
    "pdf.fonttype": 42,
})

# ── Load & count per tissue × IMPACT ─────────────────────────────────────────
df = pd.read_csv(INPUT_CSV)
df = df[df["SpliceAI"] >= SPLICEAI_MIN].copy()
counts = df.groupby(["TissueID", "Impact"]).size().unstack(fill_value=0)
# ensure H/M/L columns exist
for c in ["HIGH","MODERATE","LOW","MODIFIER"]:
    if c not in counts.columns: counts[c] = 0
counts["total"] = counts[["HIGH","MODERATE","LOW","MODIFIER"]].sum(axis=1)
counts = counts.sort_values("total", ascending=False)

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(max(9, 0.55*len(counts) + 2), 6.5))

x = list(range(len(counts)))
tids = counts.index.tolist()

# Base bar (total count)
for i, tid in enumerate(tids):
    info = TISSUE.get(tid, {"hex":"#888888","full":tid})
    has_lc = tid in LEAFCUTTER_TISSUES
    total = int(counts.loc[tid, "total"])

    # Bar: solid if LeafCutter available, hatched if not
    ax.bar(i, total,
           color=info["hex"],
           edgecolor="#1F2937",
           linewidth=0.8,
           hatch="" if has_lc else "///",
           alpha=1.0 if has_lc else 0.55)

    # Count above bar (no H/M/L breakdown)
    ax.text(i, total + 0.8, str(total),
            ha="center", va="bottom", fontsize=11, fontweight="bold",
            color="#1F2937")

# Y-axis
ymax_val = counts["total"].max()
ax.set_ylim(-ymax_val*0.15, ymax_val*1.25)
ax.set_ylabel(f"# candidates (VEP splice + SpliceAI ≥ {SPLICEAI_MIN})", fontsize=11)

# X-axis labels: tissue full name (rotated)
ax.set_xticks(x)
labels = [f"{TISSUE.get(t,{'full':t})['full']}\n({t})" for t in tids]
ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=9.5)

# LeafCutter row below x-axis
for i, tid in enumerate(tids):
    has_lc = tid in LEAFCUTTER_TISSUES
    ax.text(i, -ymax_val*0.11,
            "✓" if has_lc else "✗",
            ha="center", va="center", fontsize=14, fontweight="bold",
            color="#16A34A" if has_lc else "#DC2626")

# Row header for LeafCutter annotation
ax.text(-0.9, -ymax_val*0.11,
        "LeafCutter\navailable:",
        ha="right", va="center", fontsize=9, color="#374151", fontweight="bold")

# Legend for solid vs hatched
legend_elements = [
    Patch(facecolor="#888888", edgecolor="#1F2937", label="Has LeafCutter perind (testable)"),
    Patch(facecolor="#888888", edgecolor="#1F2937", hatch="///", alpha=0.55,
          label="No LeafCutter perind (not testable yet)"),
]
ax.legend(handles=legend_elements, loc="upper right", fontsize=9,
          frameon=True, edgecolor="lightgray", facecolor="white")

# Cosmetic
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_title(
    f"VEP-splice + SpliceAI ≥ {SPLICEAI_MIN} candidates by tissue  "
    f"(total: {counts['total'].sum()})",
    fontsize=13, fontweight="bold", loc="left", color="#111827", pad=12)

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=180, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_PDF,               bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved -> {OUT_PNG}")
print(f"Saved -> {OUT_PDF}")

# Print a text summary
print(f"\n=== Summary ===")
print(f"Total candidates: {counts['total'].sum()}")
print(f"Tissues represented: {len(counts)}")
print(f"With LeafCutter data: {counts[counts.index.isin(LEAFCUTTER_TISSUES)]['total'].sum()}")
print(f"WITHOUT LeafCutter data: {counts[~counts.index.isin(LEAFCUTTER_TISSUES)]['total'].sum()}")
print()
print(counts[["HIGH","MODERATE","LOW","MODIFIER","total"]].to_string())
