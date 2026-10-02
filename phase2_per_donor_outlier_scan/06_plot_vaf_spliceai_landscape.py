# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# SpliceAI x VAF landscape with the detectability zones shaded.
#
# Original location in the analysis project:
#     0619/code/04_vaf_spliceai_distribution.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
04_vaf_spliceai_distribution.py
================================

Research question
-----------------
Where are the V2 somatic variants in the (SpliceAI × VAF) plane, and which
zones are useful for the per-donor outlier validation work?

Output
------
- output/figures/04_vaf_spliceai_landscape.png
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from _utils import load_v2, TISSUE_PREFIX

OUT = ROOT / "output" / "figures" / "04_vaf_spliceai_landscape.png"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

# ---- Load everything (full V2, all SpliceAI tiers) ----
v2 = load_v2(min_spliceai=0.0)
print(f"Total V2 variants: {len(v2)}")

# Mark testable (donor has perind in matched tissue)
v2["tissue_full"] = v2["TissueID"].map(TISSUE_PREFIX)
v2["testable"]    = v2["tissue_full"].notna()

# Heroes from step 02 (annotated by name on plot)
heroes = pd.read_csv(ROOT / "output" / "tables" / "02_outlier_scan_all.csv")
heroes["abs_z"] = heroes["z"].abs()
TOP = heroes.query("abs_z >= 2 and shift_over_vaf < 1 and SpliceAI >= 0.2") \
              .sort_values("abs_z", ascending=False)

# ---- Figure: 2-panel ----
fig, axes = plt.subplots(1, 2, figsize=(14, 6.5),
                          gridspec_kw={"width_ratios": [1.4, 1]})

# === Panel 1: full landscape with zones ===
ax = axes[0]
# all variants (gray faint)
ax.scatter(v2["score"], v2["SR_VAF"], s=4, color="#C0C7D0", alpha=0.3,
            label=f"all V2 variants (n={len(v2)})")
# testable variants (darker)
test = v2[v2["testable"]]
ax.scatter(test["score"], test["SR_VAF"], s=8, color="#6B7889", alpha=0.55,
            label=f"testable (carrier has RNA, n={len(test)})")

# shade the zones we use
zones = [
    # (x0, x1, y0, y1, color, label)
    (0.20, 0.50, 0.05, 0.30, "#FBBF24", "LOOSE × mid-VAF"),
    (0.50, 0.80, 0.05, 0.30, "#F97316", "STRICT × mid-VAF"),
    (0.80, 1.01, 0.05, 0.30, "#DC2626", "high-conf × mid-VAF"),
    (0.20, 0.50, 0.10, 0.30, None,       None),
    (0.50, 0.80, 0.10, 0.30, None,       None),
    (0.80, 1.01, 0.10, 0.30, None,       None),
]
for x0, x1, y0, y1, color, label in zones:
    if color is None: continue
    ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, facecolor=color,
                            alpha=0.18, edgecolor=color, lw=1.2,
                            label=label))

# Physical-limit guide line: VAF = expected max |Δψ| if SpliceAI is perfect
ax.axhline(0.05, color="black", ls=":", lw=0.8, alpha=0.5)
ax.text(0.99, 0.052, "VAF=0.05 (signal below noise below this line)",
         ha="right", va="bottom", fontsize=9, color="#374151", style="italic")
ax.axhline(0.10, color="black", ls="--", lw=0.8, alpha=0.5)
ax.text(0.99, 0.102, "VAF=0.10 (clean case-study territory above this line)",
         ha="right", va="bottom", fontsize=9, color="#374151", style="italic")

# Annotate the heroes
for _, h in TOP.iterrows():
    ax.scatter(h["SpliceAI"], h["SR_VAF"], s=80,
                facecolor="#10B981", edgecolor="black", lw=1.0, zorder=10)
    ax.annotate(f"{h['gene']}\n({h['donor']})",
                 (h["SpliceAI"], h["SR_VAF"]),
                 xytext=(6, 6), textcoords="offset points",
                 fontsize=8.5, color="#065F46", fontweight="bold")

ax.set_xlim(0, 1.02); ax.set_ylim(0, 0.32)
ax.set_xlabel("SpliceAI max score", fontsize=11)
ax.set_ylabel("SR_VAF (somatic VAF, mean across cores)", fontsize=11)
ax.set_title("Where do V2 somatic variants sit?\n"
              "Shaded boxes = the zones we use; green dots = hero candidates",
              fontsize=12, fontweight="bold", loc="left")
ax.legend(loc="upper right", fontsize=8.5, frameon=False)
ax.grid(alpha=0.2)

# === Panel 2: counts table as text ===
ax2 = axes[1]
ax2.axis("off")
ax2.set_title("How many variants in each zone?",
                fontsize=12, fontweight="bold", loc="left")

# Build the count table
def count(d, score_range, vaf_range):
    s_lo, s_hi = score_range; v_lo, v_hi = vaf_range
    return ((d["score"]>=s_lo) & (d["score"]<s_hi) &
             (d["SR_VAF"]>=v_lo) & (d["SR_VAF"]<v_hi)).sum()

rows = []
for s_label, s_lo, s_hi in [("LOOSE  0.2–0.5",  0.2, 0.5),
                              ("STRICT 0.5–0.8",  0.5, 0.8),
                              ("high   ≥0.8",    0.8, 1.01)]:
    r_all  = count(v2,   (s_lo, s_hi), (0, 0.05))
    r_lowt = count(test, (s_lo, s_hi), (0, 0.05))
    r_mid  = count(v2,   (s_lo, s_hi), (0.05, 0.10))
    r_midt = count(test, (s_lo, s_hi), (0.05, 0.10))
    r_hi   = count(v2,   (s_lo, s_hi), (0.10, 1.01))
    r_hit  = count(test, (s_lo, s_hi), (0.10, 1.01))
    rows.append([s_label,
                 f"{r_lowt}/{r_all}",
                 f"{r_midt}/{r_mid}",
                 f"{r_hit}/{r_hi}"])

# Pretty text table
y0 = 0.95; dy = 0.10
ax2.text(0.0, y0, f"{'':<18}{'VAF<0.05':>12}{'0.05–0.10':>14}{'VAF≥0.10':>12}",
          family="monospace", fontsize=11, fontweight="bold")
for i, r in enumerate(rows):
    y = y0 - (i+1)*dy
    ax2.text(0.0, y, f"{r[0]:<18}{r[1]:>12}{r[2]:>14}{r[3]:>12}",
              family="monospace", fontsize=11)
ax2.text(0.0, y0 - 4.5*dy,
          "(format: testable / all)\n"
          "testable = donor has matched-tissue perind RNA",
          fontsize=9, color="#6B7280", style="italic")

# Recommendation box
ax2.text(0.0, 0.30,
          "Where to spend effort:",
          fontsize=11, fontweight="bold", color="#1F2937")
recs = [
    "•  Case study (lightning-talk material):",
    "    SpliceAI ≥ 0.5  AND  VAF ≥ 0.10",
    "    (very few — 4 total; only PLBD2 made the hero list)",
    "",
    "•  Lenient case study (Yilin's relaxation):",
    "    SpliceAI ≥ 0.2  AND  VAF ≥ 0.05",
    "    (22 candidates → 8 heroes after | z | + dist filter)",
    "",
    "•  Aggregate trend (population-level claim):",
    "    SpliceAI ≥ 0.2 (all 717 LOOSE+) regardless of VAF",
    "    (test: does |z| trend with SpliceAI × VAF?)",
]
for i, t in enumerate(recs):
    ax2.text(0.0, 0.24 - i*0.025, t, fontsize=9.5, color="#1F2937")

ax2.set_xlim(0, 1); ax2.set_ylim(0, 1)

fig.suptitle("V2 SpliceAI × VAF landscape — what we can actually test",
              fontsize=13, fontweight="bold", y=1.00)
fig.tight_layout()
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print(f"Saved: {OUT}")
