# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# Two-panel case-study figure: junction arcs + per-donor psi scatter.
#
# Original location in the analysis project:
#     0619/code/03_case_study_figures.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
03_case_study_figures.py
========================

Research question
-----------------
For the small set of "clean" SpliceAI validation hits identified by step 02
(high SpliceAI + |z|≥2 + shift_over_vaf<1, i.e. observed ψ shift within VAF
physical limit), produce one publication-quality figure per case showing the
carrier's intron usage vs the cohort.

Design (per Yilin's "previous fig too busy" feedback)
-----------------------------------------------------
2 panels per case (simpler than the old 3-panel arc+z+scatter):

  TOP   — arc plot: each junction in the cluster as a half-circle arc
          - carrier shown in red, arc thickness ∝ carrier's ψ
          - cohort median shown in gray dashed, thickness ∝ median ψ
          - variant position marked with vertical black line + label
          - SpliceAI + VAF annotated as title

  BOTTOM — per-donor scatter on the TOP-|z| junction:
          - one dot per donor, sorted by ψ
          - carrier highlighted (red, larger), all others gray
          - cohort median as horizontal dashed line
          - y-axis = "intron usage ψ"

Output
------
- output/figures/case_studies/<gene>_<donor>_<tissue>.png  (one per hero)
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from _utils import (TISSUES, TISSUE_PREFIX, load_perind, load_v2,
                     per_donor_psi, find_clusters_overlapping,
                     junctions_of_cluster, robust_zscore)

SCAN_CSV = ROOT / "output" / "tables" / "02_outlier_scan_all.csv"
OUT_DIR  = ROOT / "output" / "figures" / "case_studies"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# heroes filter: |z| ≥ 2, shift_over_vaf < 1, SpliceAI ≥ 0.2
HERO_FILTER = "abs_z >= 2 and shift_over_vaf < 1 and SpliceAI >= 0.20"

COL_CARRIER = "#D1495B"
COL_OTHER   = "#9CA3AF"
COL_MEDIAN  = "#374151"
COL_VARIANT = "#000000"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})


def plot_one_case(row, perind, psi_donor, donors):
    """Plot a single case-study figure. Returns the figure object."""
    gene, donor, tissue   = row["gene"], row["donor"], row["tissue"]
    var_chr, var_pos      = row["chr"], int(row["pos"])
    sai, vaf              = row["SpliceAI"], row["SR_VAF"]
    cluster_id            = row["cluster"]
    top_jct_start         = int(row["jct_start"])
    top_jct_end           = int(row["jct_end"])

    # All junctions in the cluster
    jcts = junctions_of_cluster(perind, cluster_id)
    j_idx = [j["idx"] for j in jcts]
    psi_clu = psi_donor[j_idx, :]                       # n_jct × n_donors
    d_i     = donors.index(donor)
    carrier = psi_clu[:, d_i]
    others_mask = np.array([i != d_i for i in range(len(donors))])
    others_med  = np.nanmedian(psi_clu[:, others_mask], axis=1)

    # Find top-|z| junction (the one driving the signal)
    top_k = None
    for k, j in enumerate(jcts):
        if j["start"] == top_jct_start and j["end"] == top_jct_end:
            top_k = k; break

    # === figure ===
    fig = plt.figure(figsize=(11, 6.5))
    gs  = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.4], hspace=0.4)

    # ===== TOP: arc plot =====
    ax1 = fig.add_subplot(gs[0])
    clu_start = min(j["start"] for j in jcts)
    clu_end   = max(j["end"]   for j in jcts)
    # extend x-limits to include the variant position (if outside cluster span)
    x_lo = min(clu_start, var_pos)
    x_hi = max(clu_end,   var_pos)
    pad  = max(500, (x_hi - x_lo) * 0.05)
    ax1.set_xlim(x_lo - pad, x_hi + pad)
    # flag the variant-cluster distance if variant is outside the cluster span
    var_dist = 0
    if   var_pos < clu_start: var_dist = clu_start - var_pos
    elif var_pos > clu_end:   var_dist = var_pos - clu_end

    # variant marker
    ax1.axvline(var_pos, color=COL_VARIANT, ls="--", lw=1.4, alpha=0.85, zorder=5)
    ax1.text(var_pos, 1.18, "⚡ variant", ha="center", va="bottom",
              fontsize=10, color=COL_VARIANT, fontweight="bold")

    # arcs — height fixed (so arcs at different widths are comparable),
    # but linewidth ∝ psi (visually salient)
    for k, j in enumerate(jcts):
        x_mid = (j["start"] + j["end"]) / 2
        w     = j["end"] - j["start"]
        # cohort median (gray, dashed underlay)
        lw_o = max(0.5, 8 * others_med[k])
        arc_o = Arc((x_mid, 0), w, 1.0, angle=0, theta1=0, theta2=180,
                     color=COL_OTHER, lw=lw_o, ls="--", alpha=0.55, zorder=2)
        ax1.add_patch(arc_o)
        # carrier (red, solid overlay)
        lw_c = max(0.8, 10 * carrier[k])
        arc_c = Arc((x_mid, 0), w, 1.0, angle=0, theta1=0, theta2=180,
                     color=COL_CARRIER, lw=lw_c, alpha=0.9, zorder=3)
        ax1.add_patch(arc_c)
        # annotation on the top-|z| junction
        if k == top_k:
            ax1.text(x_mid, 0.55,
                      f"ψ carrier={carrier[k]:.2f}\nψ others={others_med[k]:.2f}",
                      ha="center", va="center", fontsize=9, color=COL_CARRIER,
                      fontweight="bold",
                      bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                  edgecolor=COL_CARRIER, lw=0.8))

    ax1.set_ylim(-0.05, 1.4)
    ax1.set_yticks([])
    ax1.set_xlabel(f"{var_chr} genomic position", fontsize=10)
    ax1.spines[["top", "right", "left"]].set_visible(False)
    title = (f"{gene}  ·  {tissue}  ·  carrier {donor}    "
              f"SpliceAI = {sai:.2f}  |  VAF = {vaf:.3f}")
    if var_dist > 0:
        title += f"\n⚠  variant is {var_dist:,} bp outside the cluster span — check this is not coincidental"
    ax1.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=8)

    # legend
    from matplotlib.lines import Line2D
    leg = [
        Line2D([0],[0], color=COL_CARRIER, lw=3, label=f"carrier {donor}"),
        Line2D([0],[0], color=COL_OTHER,   lw=2, ls="--",
                label=f"cohort median (n={int(row['n_donors_cohort'])})"),
    ]
    ax1.legend(handles=leg, loc="upper right", fontsize=9, frameon=False)

    # ===== BOTTOM: per-donor scatter on top-|z| junction =====
    ax2 = fig.add_subplot(gs[1])
    vals = psi_clu[top_k, :]                       # one ψ per donor
    order = np.argsort(np.where(np.isnan(vals), -np.inf, vals))
    xs = np.arange(len(donors))
    for x_i, o in enumerate(order):
        is_carrier = (donors[o] == donor)
        color = COL_CARRIER if is_carrier else COL_OTHER
        size  = 130 if is_carrier else 70
        zo    = 5 if is_carrier else 2
        ax2.scatter(x_i, vals[o], s=size, color=color, edgecolor="black",
                     lw=0.5, zorder=zo)
        if is_carrier:
            ax2.annotate(f"{donors[o]}\nψ={vals[o]:.3f}",
                          (x_i, vals[o]),
                          xytext=(8, 12), textcoords="offset points",
                          fontsize=10, color=COL_CARRIER, fontweight="bold")
    med = np.nanmedian(vals[others_mask])
    ax2.axhline(med, color=COL_MEDIAN, ls="--", lw=1, alpha=0.7)
    ax2.text(len(donors)-0.5, med, f"  cohort median = {med:.3f}",
              ha="left", va="center", fontsize=9, color=COL_MEDIAN)

    ax2.set_xticks(xs)
    ax2.set_xticklabels([donors[o] for o in order], rotation=45,
                         ha="right", fontsize=9)
    ax2.set_ylabel("intron usage ψ\n(reads on junction / cluster reads)",
                    fontsize=10)
    ax2.set_ylim(-0.05, max(1.05, np.nanmax(vals)*1.1))
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.set_title(
        f"Per-donor ψ on top-|z| junction  {top_jct_start:,}–{top_jct_end:,}  "
        f"|  |z| = {abs(row['z']):.2f}  |  observed Δψ = {row['psi_shift']:+.3f}  "
        f"(VAF physical max: ±{vaf:.3f})",
        loc="left", fontsize=10)

    fig.tight_layout()
    return fig


def main():
    df = pd.read_csv(SCAN_CSV)
    df["abs_z"] = df["z"].abs()
    heroes = df.query(HERO_FILTER).sort_values("abs_z", ascending=False)
    print(f"  Heroes (|z|≥2, shift_over_vaf<1, SpliceAI≥0.20): {len(heroes)}")
    print(heroes[["gene","donor","tissue","SpliceAI","SR_VAF",
                    "z","shift_over_vaf"]].to_string(index=False))

    # cache perind loads (one per tissue we need)
    needed_tissues = heroes["tissue"].unique()
    print(f"\n  Loading perind for {len(needed_tissues)} tissues...")
    cache = {}
    for t in needed_tissues:
        p = load_perind(t)
        psi_d, donors_d = per_donor_psi(p)
        cache[t] = (p, psi_d, donors_d)

    print(f"\n  Generating figures to {OUT_DIR}/")
    for _, row in heroes.iterrows():
        p, psi_d, donors_d = cache[row["tissue"]]
        try:
            fig = plot_one_case(row, p, psi_d, donors_d)
            out = OUT_DIR / f"{row['gene']}_{row['donor']}_{row['tissue']}.png"
            fig.savefig(out, dpi=140, bbox_inches="tight")
            plt.close(fig)
            print(f"    ✓ {out.name}")
        except Exception as e:
            print(f"    ✗ {row['gene']}/{row['donor']}/{row['tissue']}: {e}")


if __name__ == "__main__":
    main()
