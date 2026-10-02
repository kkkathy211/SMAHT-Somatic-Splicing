# -*- coding: utf-8 -*-
# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# Cummings et al. 2017 Fig.2 layout: stacked carrier/control coverage + arcs + splice-site motif.
#
# Original location in the analysis project:
#     0625/13_pileup_cummings_fig2.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
13_pileup_cummings_fig2.py
==========================
Three figures (PLBD2, ISCU, C11orf54) in the EXACT layout of
Cummings et al. 2017 (Science Trans Med, PMC5548421) Figure 2.

Layout per gene:
  1. Top panel    : Control (mean of non-carriers) — BLUE filled coverage,
                    junction arcs with READ-COUNT labels
  2. Bottom panel : Carrier — RED filled coverage,
                    junction arcs with READ-COUNT labels
                    Same Y-axis as Control panel
  3. Gene track   : dark exon boxes + arrowed intron lines
                    Gene name to the right
  4. Sequence     : splice-site motifs (WT vs Variant)
                    variant nucleotide in RED

Uses `samtools depth` (via subprocess) for coverage so introns drop to 0,
and pysam for junction extraction.

Carriers:
  PLBD2    -> SMHT016
  ISCU     -> SMHT039
  C11orf54 -> SMHT039

Run on Broad:
    python3 13_pileup_cummings_fig2.py
"""

import subprocess
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
import pysam

BAM_DIR = Path("/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin")
OUT_DIR = Path("/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/pileup_from_Laura")

LIBSIZES = {
    "SMHT001": 137_854_051, "SMHT005":  50_585_123, "SMHT006": 104_900_557,
    "SMHT007": 175_064_688, "SMHT009": 116_122_104, "SMHT016": 183_010_215,
    "SMHT022": 103_429_558, "SMHT027": 136_061_708, "SMHT029": 119_051_166,
    "SMHT039": 191_231_191,
}

# Cummings-style colors
COL_CTRL = "#1F4E8C"      # control blue
COL_PAT  = "#C8102E"      # patient red
COL_EXON = "#3F3F3F"      # gene track dark gray
COL_RED  = "#C8102E"      # variant nucleotide

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})

# Canonical transcript exons (one per gene; gene track shows just the main one)
EXONS = {
    "PLBD2":    [(113358587,113358890),(113369116,113369209),(113372649,113372807),
                 (113374474,113374574),(113374793,113375007),(113380745,113380842),
                 (113384105,113384265),(113384851,113384946),(113385212,113385283),
                 (113386937,113387089),(113387744,113387906),(113388459,113391629)],
    "ISCU":     [(108562596,108562736),(108564279,108564392),(108565321,108565431),
                 (108567190,108567268),(108568831,108569368)],
    "C11orf54": [(93741672,93741728),(93747297,93747448),(93750346,93750444),
                 (93753682,93753755),(93753936,93754037),(93755210,93755386),
                 (93757316,93757465),(93759742,93759858),(93761515,93764749)],
}

# Per-case info: variant + zoom + sequence motifs
CASES = [
    {
        "gene":      "PLBD2",
        "carrier":   "SMHT016",
        "chrom":     "chr12",
        "strand":    "+",
        "var_pos":   113_385_204,
        "zoom_lo":   113_384_000,
        "zoom_hi":   113_387_300,
        "spliceai":  1.00,
        "vaf":       0.058,
        "abs_z":     3.55,
        "var_note":  "splice acceptor -8",
        "seq_wt":    [("polypyrimidine -8 ...", "AG"), ("(exon 9 start)", "")],
        "seq_var":   [("polypyrimidine -8 ...", "AG"), ("(variant at -8 in PPT)", "")],
        # which side of variant to put a dashed box on (highlights affected exon)
        "highlight": (113385212, 113385283),  # exon 9
    },
    {
        "gene":      "ISCU",
        "carrier":   "SMHT039",
        "chrom":     "chr12",
        "strand":    "+",
        "var_pos":   108_565_432,
        "zoom_lo":   108_563_900,
        "zoom_hi":   108_569_400,
        "spliceai":  0.99,
        "vaf":       0.019,
        "abs_z":     4.31,
        "var_note":  "splice donor +1 (G>A)",
        "seq_wt":    [("exon 3 end", "GT"), ("intron 3", "")],
        "seq_var":   [("exon 3 end", "AT"), ("variant +1: G->A", "")],
        "highlight": (108565321, 108565431),  # exon 3
    },
    {
        "gene":      "C11orf54",
        "carrier":   "SMHT039",
        "chrom":     "chr11",
        "strand":    "+",
        "var_pos":   93_757_419,
        "zoom_lo":   93_754_500,
        "zoom_hi":   93_760_300,
        "spliceai":  0.38,
        "vaf":       0.056,
        "abs_z":     3.84,
        "var_note":  "exonic (within exon 7)",
        "seq_wt":    [("exon 7", "intact"), ("", "")],
        "seq_var":   [("exon 7", "variant inside"), ("", "")],
        "highlight": (93757316, 93757465),  # exon 7
    },
]


# ── Pick one BAM per donor ───────────────────────────────────────────────────
bams_all = sorted(f for f in BAM_DIR.glob("*.bam") if not str(f).endswith(".bai"))
donor_bam = {}
for bam in bams_all:
    donor = bam.name.split("-")[0]
    if donor not in donor_bam or "A101" in bam.name:
        donor_bam[donor] = bam
donors_all = sorted(donor_bam.keys())
print("Donors: " + str(donors_all))


# ── Coverage via samtools depth (introns -> 0) ───────────────────────────────
def samtools_depth(bam_path, chrom, start, end):
    n = end - start
    cov = np.zeros(n, dtype=np.float64)
    region = chrom + ":" + str(start + 1) + "-" + str(end)
    cmd = ["samtools", "depth", "-a", "-r", region, str(bam_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        pos   = int(parts[1])
        depth = int(parts[2])
        idx   = pos - 1 - start
        if 0 <= idx < n:
            cov[idx] = depth
    return cov


# ── Junctions via pysam CIGAR-N (filtered) ───────────────────────────────────
def junctions(bam_path, chrom, start, end):
    jcts = defaultdict(int)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if (read.is_unmapped or read.is_secondary or
                read.is_supplementary or read.is_duplicate or read.is_qcfail):
                continue
            if read.mapping_quality < 1:
                continue
            rpos = read.reference_start
            for op, length in (read.cigartuples or []):
                if op == 0:
                    rpos += length
                elif op == 3:
                    if start <= rpos and (rpos + length) <= end:
                        jcts[(rpos, rpos + length)] += 1
                    rpos += length
                elif op in (2, 7, 8):
                    rpos += length
    return jcts


# ── Drawing helpers ──────────────────────────────────────────────────────────
def draw_arc_with_label(ax, x0, x1, value, color, apex_y, lw=1.6):
    xs = np.linspace(x0, x1, 100)
    mid = (x0 + x1) / 2
    halfw = (x1 - x0) / 2
    ys = apex_y * (1 - ((xs - mid) / halfw) ** 2)
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="round", zorder=5)
    ax.text(mid, apex_y * 1.02, "{:d}".format(int(round(value))),
            ha="center", va="bottom", fontsize=10,
            color="black", zorder=6)


def select_top(jdict, lo, hi, max_n=4):
    items = [(k, v) for k, v in jdict.items()
             if v > 0 and (k[1] - k[0]) > 30 and k[0] >= lo and k[1] <= hi]
    items.sort(key=lambda kv: -kv[1])
    return items[:max_n]


def draw_gene_track(ax, exon_list, lo, hi, strand, gene_name, highlight=None):
    visible = [(s, e) for (s, e) in exon_list if e >= lo and s <= hi]
    line_y = 0.5
    box_h  = 0.55
    # intron line + arrowheads
    line_lo = lo + 30
    line_hi = hi - 30
    ax.plot([line_lo, line_hi], [line_y, line_y],
            color=COL_EXON, lw=0.9, zorder=2)
    n_arrows = 16
    for k in range(1, n_arrows):
        x = line_lo + (line_hi - line_lo) * k / n_arrows
        dx = (line_hi - line_lo) * 0.012 * (1 if strand == "+" else -1)
        ax.annotate("", xy=(x + dx, line_y), xytext=(x, line_y),
                    arrowprops=dict(arrowstyle="->", color=COL_EXON,
                                    lw=0.7, mutation_scale=8), zorder=2)
    # exon boxes
    for s, e in visible:
        s_c = max(s, lo); e_c = min(e, hi)
        ax.add_patch(Rectangle((s_c, line_y - box_h / 2),
                               e_c - s_c, box_h,
                               facecolor=COL_EXON, edgecolor="black",
                               lw=0.8, zorder=3))
    # dashed highlight box around affected exon
    if highlight is not None:
        hs, he = highlight
        if he >= lo and hs <= hi:
            pad = (hi - lo) * 0.005
            ax.add_patch(Rectangle((hs - pad, line_y - box_h / 2 - 0.10),
                                   (he - hs) + 2 * pad, box_h + 0.20,
                                   facecolor="none", edgecolor="black",
                                   lw=1.0, linestyle="--", zorder=4))
    # gene name on the right
    ax.text(hi - 60, line_y + 0.05, gene_name,
            ha="right", va="center", fontsize=11,
            fontweight="bold", fontstyle="italic", color="black")
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.2, 1.3)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(False)


def draw_seq_panel(ax, var_pos, hs, he, seq_wt, seq_var):
    """
    Two-line text panel:  WT  :   ...left...  | ...right...
                          Var :   ...left...  | ...right...  (variant in red)
    Positions are aligned to the splice site (left edge = acceptor, right edge = donor)
    """
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(False)

    # row labels
    ax.text(0.05, 0.72, "WT",       ha="left", va="center",
            fontsize=10, color="#6B7280")
    ax.text(0.05, 0.28, "Variant",  ha="left", va="center",
            fontsize=10, color=COL_RED, fontweight="bold")

    # WT row content
    ax.text(0.30, 0.72, seq_wt[0][0], ha="right", va="center",
            fontsize=10, family="monospace", color="black")
    ax.text(0.32, 0.72, seq_wt[0][1], ha="left", va="center",
            fontsize=11, family="monospace", color="black", fontweight="bold")
    ax.text(0.55, 0.72, seq_wt[1][0], ha="left", va="center",
            fontsize=10, family="monospace", color="black")

    # Variant row content (red highlight on changed motif if present)
    ax.text(0.30, 0.28, seq_var[0][0], ha="right", va="center",
            fontsize=10, family="monospace", color="black")
    ax.text(0.32, 0.28, seq_var[0][1], ha="left", va="center",
            fontsize=11, family="monospace", color=COL_RED, fontweight="bold")
    ax.text(0.55, 0.28, seq_var[1][0], ha="left", va="center",
            fontsize=10, family="monospace", color=COL_RED)


# ── Main loop ────────────────────────────────────────────────────────────────
for case in CASES:
    gene    = case["gene"]
    carrier = case["carrier"]
    chrom   = case["chrom"]
    lo      = case["zoom_lo"]
    hi      = case["zoom_hi"]
    var_pos = case["var_pos"]

    print("\n=== " + gene + " (carrier=" + carrier + ") ===")

    # Scale everything to a common library size (100M reads)
    # so all donors are directly comparable, and arc numbers are
    # in the 100-500 range like Cummings (instead of 0-7 with RPM).
    SCALE = 1e8   # "per 100 million reads"

    cov_rpm  = {}
    jcts_rpm = {}
    for donor in donors_all:
        print("  " + donor + " ...", end=" ")
        cov  = samtools_depth(donor_bam[donor], chrom, lo, hi)
        jc   = junctions     (donor_bam[donor], chrom, lo, hi)
        lib  = LIBSIZES.get(donor, 1e8)
        cov_rpm[donor]  = cov / lib * SCALE
        jcts_rpm[donor] = {k: v / lib * SCALE for k, v in jc.items()}
        print("max=" + str(round(float(cov_rpm[donor].max()), 1)))

    noncarriers = [d for d in donors_all if d != carrier]
    ctrl_cov = np.stack([cov_rpm[d] for d in noncarriers]).mean(axis=0)
    carr_cov = cov_rpm[carrier]

    all_keys = set()
    for d in donors_all:
        all_keys.update(jcts_rpm[d].keys())
    ctrl_jcts = {k: float(np.mean([jcts_rpm[d].get(k, 0.0) for d in noncarriers]))
                 for k in all_keys}
    carr_jcts = {k: float(jcts_rpm[carrier].get(k, 0.0)) for k in all_keys}

    # ── Figure: 4 stacked panels (ctrl, carrier, gene track, seq) ──
    fig = plt.figure(figsize=(12, 7.8))
    gs  = fig.add_gridspec(4, 1,
                           height_ratios=[3.2, 3.2, 0.7, 0.8],
                           hspace=0.10)
    ax_ctrl = fig.add_subplot(gs[0, 0])
    ax_pat  = fig.add_subplot(gs[1, 0], sharex=ax_ctrl)
    ax_gene = fig.add_subplot(gs[2, 0], sharex=ax_ctrl)
    ax_seq  = fig.add_subplot(gs[3, 0])

    x = np.arange(lo, hi)
    ymax = max(float(ctrl_cov.max()), float(carr_cov.max())) * 1.55

    # Control panel
    ax_ctrl.fill_between(x, ctrl_cov, color=COL_CTRL, lw=0)
    ax_ctrl.set_xlim(lo, hi); ax_ctrl.set_ylim(0, ymax)
    ax_ctrl.set_ylabel("Coverage", fontsize=11)
    ax_ctrl.text(0.985, 0.92, "Control", transform=ax_ctrl.transAxes,
                 ha="right", va="top", fontsize=12, fontweight="bold",
                 color="black")
    ax_ctrl.tick_params(axis="x", labelbottom=False)
    for (j0, j1), v in select_top(ctrl_jcts, lo, hi):
        max_v = max(v for _, v in select_top(ctrl_jcts, lo, hi)) or 1
        apex = ymax * (0.35 + 0.55 * (v / max_v))
        draw_arc_with_label(ax_ctrl, j0, j1, v, COL_CTRL, apex)

    # Patient/Carrier panel
    ax_pat.fill_between(x, carr_cov, color=COL_PAT, lw=0)
    ax_pat.set_xlim(lo, hi); ax_pat.set_ylim(0, ymax)
    ax_pat.set_ylabel("Coverage", fontsize=11)
    ax_pat.text(0.985, 0.92,
                "Carrier  " + carrier + "  (" + gene + ")",
                transform=ax_pat.transAxes,
                ha="right", va="top", fontsize=12, fontweight="bold",
                color="black")
    ax_pat.xaxis.set_major_formatter(
        FuncFormatter(lambda v, p: "{:,.0f}".format(v)))
    ax_pat.tick_params(axis="x", labelbottom=True)
    for (j0, j1), v in select_top(carr_jcts, lo, hi):
        max_v = max(v for _, v in select_top(carr_jcts, lo, hi)) or 1
        apex = ymax * (0.35 + 0.55 * (v / max_v))
        draw_arc_with_label(ax_pat, j0, j1, v, COL_PAT, apex)

    # Gene track
    draw_gene_track(ax_gene, EXONS[gene], lo, hi,
                    case["strand"], gene, highlight=case["highlight"])
    ax_gene.axvline(var_pos, color=COL_PAT, lw=1.2, alpha=0.9, zorder=5)

    # Sequence panel
    hs, he = case["highlight"]
    draw_seq_panel(ax_seq, var_pos, hs, he,
                   case["seq_wt"], case["seq_var"])

    fig.suptitle(
        gene + "  |  " + chrom + ":" + "{:,}".format(var_pos) +
        "  (" + case["var_note"] + ")  |  " +
        "SpliceAI=" + "{:.2f}".format(case["spliceai"]) +
        "  |  VAF=" + "{:.3f}".format(case["vaf"]) +
        "  |  |z|=" + "{:.2f}".format(case["abs_z"]),
        fontsize=11, fontweight="bold", x=0.05, y=0.985,
        ha="left", color="#1F2937")

    out_png = OUT_DIR / (gene + "_cummings_fig2.png")
    out_pdf = OUT_DIR / (gene + "_cummings_fig2.pdf")
    fig.savefig(str(out_png), dpi=170, bbox_inches="tight")
    fig.savefig(str(out_pdf),               bbox_inches="tight")
    plt.close(fig)
    print("Saved -> " + str(out_png))
    print("Saved -> " + str(out_pdf))

print("\nAll done.")
