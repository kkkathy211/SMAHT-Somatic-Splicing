# -*- coding: utf-8 -*-
# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# Per-donor psi strip plot. Works at low VAF where pileups do not. Matches junctions on chr:start:end only.
#
# Original location in the analysis project:
#     0630/code/07_psi_stripplot.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
07_psi_stripplot.py
===================
Show carrier-as-outlier visually: for each case study's outlier junction,
plot every donor's psi value in the cohort as a strip/dot plot.

This is the visualization that ACTUALLY shows what LeafCutter's z-score
captures: carrier sits far from cohort's tight distribution.

Compared to sashimi/pileup which fails to show low-VAF signal, this
strip plot dramatically shows carrier vs cohort separation.
"""

import gzip
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT     = Path("/Users/spectremac/Desktop/Differential_splicing")
PERIND   = ROOT / "leafcutter_results" / "3AD_Skin_Calf" / "3AD_Skin_Calf_perind.counts.gz"
SCAN_CSV = ROOT / "0619" / "output" / "tables" / "02_outlier_scan_all.csv"
OUT_DIR  = ROOT / "0630" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Cases to plot (matching outlier scan)
CASES = [
    {"gene":"MCM7",     "carrier":"SMHT005", "chr":"chr7",  "cluster":"clu_19626_-",
     "jct_start":100100093, "jct_end":100100213,
     "spliceai":0.94, "vaf":0.021, "abs_z":2.78, "delta_psi":0.213},
    {"gene":"RSRC2",    "carrier":"SMHT005", "chr":"chr12", "cluster":"clu_4996_-",
     "jct_start":122519051, "jct_end":122520504,
     "spliceai":1.00, "vaf":0.014, "abs_z":0.77, "delta_psi":0.019},
    {"gene":"PLBD2",    "carrier":"SMHT016", "chr":"chr12", "cluster":"clu_5571_+",
     "jct_start":113384946, "jct_end":113385212,
     "spliceai":1.00, "vaf":0.058, "abs_z":3.55, "delta_psi":0.038},
    {"gene":"ISCU",     "carrier":"SMHT039", "chr":"chr12", "cluster":"clu_5515_+",
     "jct_start":108567757, "jct_end":108568831,
     "spliceai":0.99, "vaf":0.019, "abs_z":4.31, "delta_psi":0.011},
    {"gene":"C11orf54", "carrier":"SMHT039", "chr":"chr11", "cluster":"clu_4223_+",
     "jct_start":93754037, "jct_end":93755210,
     "spliceai":0.38, "vaf":0.056, "abs_z":3.84, "delta_psi":0.044},
]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size":   11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})


# ── Load perind file ──────────────────────────────────────────────────────
with gzip.open(str(PERIND), "rt") as f:
    header = f.readline().strip().split()
    # First column is junction id; remaining are sample paths ending in .regtools_junc.txt.gz
    sample_cols = header[1:]
    # extract donor from sample name
    donors = [s.split("-")[0] for s in sample_cols]
    print("Donors in perind:", donors)

    # index rows by junction id
    lines = f.readlines()

print(f"Total junctions in perind: {len(lines)}")

# Parse "K/N" -> psi
def psi_from_kn(s):
    try:
        k, n = s.split("/")
        k, n = float(k), float(n)
        return k / n if n > 0 else np.nan
    except:
        return np.nan


# ── Find junction row for each case + compute psi per donor ─────────────
def find_junction_psi(case):
    """Return dict {donor: psi} for the outlier junction of this case.
    Match by (chr, start, end), ignoring cluster id (leafcutter renumbers
    clusters between runs)."""
    prefix = f"{case['chr']}:{case['jct_start']}:{case['jct_end']}:"
    for line in lines:
        parts = line.strip().split()
        if parts[0].startswith(prefix):
            print(f"  Matched junction: {parts[0]}")
            kn = parts[1:]
            psi = [psi_from_kn(v) for v in kn]
            return dict(zip(donors, psi))
    print(f"  ! Junction {prefix}* NOT FOUND in perind")
    return None


# ── Plot: one panel per case, all in one figure ───────────────────────────
n_cases = len(CASES)
fig, axes = plt.subplots(1, n_cases, figsize=(2.5 * n_cases + 1, 5.5), sharey=False)
if n_cases == 1:
    axes = [axes]

for ax, case in zip(axes, CASES):
    print(f"\n=== {case['gene']} (carrier {case['carrier']}) ===")
    psi_map = find_junction_psi(case)
    if psi_map is None:
        ax.text(0.5, 0.5, "Junction not found",
                ha="center", va="center", transform=ax.transAxes)
        continue

    carrier = case["carrier"]
    # Get non-carrier psi values (excluding NaN)
    nc_psi = np.array([v for d, v in psi_map.items()
                       if d != carrier and not np.isnan(v)])
    c_psi = psi_map.get(carrier, np.nan)

    print(f"  Non-carrier psi (n={len(nc_psi)}): mean={nc_psi.mean():.3f}, std={nc_psi.std():.3f}")
    print(f"  Carrier psi:              {c_psi:.3f}")

    # Jitter for strip plot
    jitter = np.random.RandomState(42).uniform(-0.12, 0.12, size=len(nc_psi))

    # Cohort dots (non-carriers)
    ax.scatter(jitter, nc_psi,
               s=80, c="#4A90D9", edgecolor="black", lw=0.7, alpha=0.7,
               label=f"Cohort (n={len(nc_psi)})", zorder=3)

    # Carrier as red star, larger
    ax.scatter([0], [c_psi],
               s=250, marker="*", c="#DC2626", edgecolor="black", lw=1.0,
               label=f"Carrier {carrier}", zorder=5)

    # Reference line at cohort median
    median = np.median(nc_psi)
    ax.axhline(median, ls="--", color="gray", lw=0.8, alpha=0.6,
               zorder=1)
    ax.text(0.5, median, "cohort median", fontsize=8, color="gray",
            ha="right", va="bottom")

    # Formatting
    ax.set_xlim(-0.5, 0.5)
    ax.set_xticks([])
    ax.set_ylabel("ψ (junction usage ratio)", fontsize=10)

    # Title with all stats
    title = (f"{case['gene']}\n"
             f"SpliceAI={case['spliceai']:.2f}, VAF={case['vaf']:.3f}\n"
             f"|z|={case['abs_z']:.2f}, |Δψ|={case['delta_psi']:.3f}")
    ax.set_title(title, fontsize=10, fontweight="bold", pad=10)

    # Set y range to nicely include both cohort and carrier
    ymin = min(nc_psi.min(), c_psi) - 0.05
    ymax = max(nc_psi.max(), c_psi) + 0.05
    ax.set_ylim(ymin, ymax)

    # Legend only on first panel
    if ax == axes[0]:
        ax.legend(fontsize=9, frameon=True, loc="best")

fig.suptitle("ψ distribution: carrier vs cohort at each case's outlier junction",
             fontsize=13, fontweight="bold", x=0.02, ha="left", y=1.02)
plt.tight_layout()
out = OUT_DIR / "psi_stripplot_all_cases.png"
fig.savefig(str(out), dpi=170, bbox_inches="tight")
fig.savefig(str(OUT_DIR / "psi_stripplot_all_cases.pdf"), bbox_inches="tight")
plt.close(fig)
print(f"\nSaved: {out}")
