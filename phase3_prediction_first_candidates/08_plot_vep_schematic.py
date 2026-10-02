# -*- coding: utf-8 -*-
# ============================================================================
# Phase 3 - Prediction-first candidate selection
#
# Publication-style gene schematic showing VEP consequence and IMPACT per case.
#
# Original location in the analysis project:
#     0716/code/13_vep_schematic_3cases.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
13_vep_schematic_3cases.py  (v2 — publication style, Cummings et al. 2017)
==========================================================================
Gene schematic + VEP annotation for 3 hero case studies.

Design principles adopted from scientific journal figures (Cummings 2017,
Nature/Science convention):
  - Muted, mostly-grayscale palette; single accent color for the variant
  - Sans-serif but plain (no bolded decorative badges)
  - Gene name in italic (biology convention)
  - Compact intron/exon glyphs, thin lines
  - Panel letters (a, b, c) in top-left
  - Structured text annotation, no interpretive banners

Output: single PNG/PDF for insertion into a slide.
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec

OUT_DIR = Path("/Users/spectremac/Desktop/Differential_splicing/0716/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family":  "DejaVu Sans",
    "font.size":     9,
    "pdf.fonttype": 42,
    "axes.edgecolor":  "#1F2937",
    "axes.linewidth":  0.6,
})

# ── Palette — muted, publication style ───────────────────────────────────────
CLR_EXON     = "#1F2937"      # dark navy — exons
CLR_INTRON   = "#6B7280"      # medium gray — intron lines
CLR_TEXT     = "#111827"      # near-black — primary text
CLR_MUTED    = "#6B7280"      # secondary
CLR_LABEL    = "#374151"      # slightly lighter than primary for coord line
IMPACT_ACCENT = {              # single-hue accent used only for variant marker
    "HIGH":     "#B91C1C",    # muted red
    "MODERATE": "#B45309",    # muted amber
    "LOW":      "#A16207",    # dark ochre
    "MODIFIER": "#6B7280",    # gray
}

# ── Cases ────────────────────────────────────────────────────────────────────
CASES = [
    {
        "letter":      "a",
        "gene":        "ISCU",
        "chrom":       "chr12",
        "strand":      "+",
        "exons":       [(108562596, 108562736), (108564279, 108564392),
                        (108565321, 108565431), (108567190, 108567268),
                        (108568831, 108569368)],
        "var_pos":     108565432,
        "var_context": "intron 3/4 (canonical AG+1)",
        "spliceai":    0.99,
        "consequence": "splice_donor_variant",
        "impact":      "HIGH",
    },
    {
        "letter":      "b",
        "gene":        "PLBD2",
        "chrom":       "chr12",
        "strand":      "+",
        "exons":       [(113358587, 113358890), (113369116, 113369209),
                        (113372649, 113372807), (113374474, 113374574),
                        (113374793, 113375007), (113380745, 113380842),
                        (113384105, 113384265), (113384851, 113384946),
                        (113385212, 113385283), (113386937, 113387089),
                        (113387744, 113387906), (113388459, 113391629)],
        "var_pos":     113385204,
        "var_context": "intron 8/11 (polypyrimidine tract)",
        "spliceai":    1.00,
        "consequence": "splice_polypyrimidine_tract_variant",
        "impact":      "LOW",
    },
    {
        "letter":      "c",
        "gene":        "C11orf54",
        "chrom":       "chr11",
        "strand":      "+",
        "exons":       [(93741672, 93741728), (93747297, 93747448),
                        (93750346, 93750444), (93753682, 93753755),
                        (93753936, 93754037), (93755210, 93755386),
                        (93757316, 93757465), (93759742, 93759858),
                        (93761515, 93764749)],
        "var_pos":     93757419,
        "var_context": "exon 7/9 (p.G204E)",
        "spliceai":    0.38,
        "consequence": "missense_variant",
        "impact":      "MODERATE",
    },
]


def draw_gene_track(ax, case):
    """Publication-style gene track: thin intron line with directional chevrons,
    narrow exon rectangles, small triangle marker for variant."""
    exons  = case["exons"]
    var    = case["var_pos"]
    strand = case["strand"]
    accent = IMPACT_ACCENT[case["impact"]]

    gene_start = min(e[0] for e in exons)
    gene_end   = max(e[1] for e in exons)
    pad        = (gene_end - gene_start) * 0.02

    # Intron line — thin, medium gray
    ax.plot([gene_start, gene_end], [0, 0],
            color=CLR_INTRON, lw=0.7, solid_capstyle="butt", zorder=1)

    # Directional chevrons (subtle, sparse)
    n_chev = 4
    for i in range(1, n_chev):
        x = gene_start + (gene_end - gene_start) * i / n_chev
        mk = "$>$" if strand == "+" else "$<$"
        ax.text(x, 0, mk, fontsize=6, color=CLR_INTRON,
                ha="center", va="center", zorder=1)

    # Exons — thin dark rectangles
    exon_h = 0.32
    for i, (s, e) in enumerate(exons, start=1):
        ax.add_patch(Rectangle(
            (s, -exon_h / 2), max(e - s, (gene_end - gene_start) * 0.002),
            exon_h,
            facecolor=CLR_EXON, edgecolor="none", zorder=3,
        ))
        # Only number the exons at the ends + those hosting/adjacent-to variant
        # to avoid clutter (paper convention)
        show_num = (i == 1 or i == len(exons))
        # Also show the exon closest to the variant
        # (heuristic: within 3-exon window of variant containing/adjacent exon)
        if not show_num:
            for j, (js, je) in enumerate(exons, start=1):
                if js - 500 <= var <= je + 500 and abs(i - j) <= 1:
                    show_num = True; break
        if show_num:
            ax.text((s + e) / 2, -0.42, str(i), fontsize=6.5,
                    color=CLR_MUTED, ha="center", va="top", zorder=4)

    # Variant marker — thin vertical line + small triangle
    ax.plot([var, var], [-exon_h * 0.55, exon_h * 1.05],
            color=accent, lw=1.1, zorder=5)
    ax.plot(var, exon_h * 1.15, marker="v", markersize=5.5,
            markerfacecolor=accent, markeredgecolor=accent, zorder=6)

    ax.set_xlim(gene_start - pad, gene_end + pad)
    ax.set_ylim(-0.75, 0.75)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(False)


def draw_panel(fig, case, gs_slot):
    """One case panel: letter, gene name, coord line, track, annotation."""
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    sub = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_slot,
                                  height_ratios=[1.0, 1.1, 1.4],
                                  hspace=0.10)

    # ── Row 1: panel letter + gene name + coord ──────────────────────────────
    ax0 = fig.add_subplot(sub[0])
    ax0.set_xlim(0, 10); ax0.set_ylim(0, 1); ax0.axis("off")
    ax0.text(0.00, 0.85, case["letter"],
             fontsize=13, fontweight="bold", color=CLR_TEXT,
             ha="left", va="top")
    ax0.text(0.55, 0.85, case["gene"],
             fontsize=12, style="italic", color=CLR_TEXT,
             ha="left", va="top")
    gene_s = min(e[0] for e in case["exons"])
    gene_e = max(e[1] for e in case["exons"])
    ax0.text(0.55, 0.35,
             f"{case['chrom']}:{gene_s:,}–{gene_e:,}  "
             f"({(gene_e - gene_s) / 1000:.1f} kb, {case['strand']} strand)",
             fontsize=8, color=CLR_MUTED, ha="left", va="top")

    # ── Row 2: gene track ────────────────────────────────────────────────────
    ax_track = fig.add_subplot(sub[1])
    draw_gene_track(ax_track, case)

    # ── Row 3: annotation block (4 fields, left-aligned "label  value") ─────
    ax2 = fig.add_subplot(sub[2])
    ax2.set_xlim(0, 10); ax2.set_ylim(0, 1); ax2.axis("off")

    label_x = 0.05
    value_x = 1.35
    ys = [0.88, 0.66, 0.42, 0.15]

    # Line 1 — Variant coordinate
    ax2.text(label_x, ys[0], "Variant",  fontsize=8, color=CLR_MUTED,
             ha="left", va="center")
    ax2.text(value_x, ys[0],
             f"{case['chrom']}:{case['var_pos']:,}   ({case['var_context']})",
             fontsize=8.5, color=CLR_TEXT, ha="left", va="center",
             family="monospace")

    # Line 2 — VEP consequence
    ax2.text(label_x, ys[1], "VEP",      fontsize=8, color=CLR_MUTED,
             ha="left", va="center")
    ax2.text(value_x, ys[1], case["consequence"],
             fontsize=8.5, color=CLR_TEXT, ha="left", va="center",
             family="monospace")

    # Line 3 — Impact (colored subtly with accent, italic)
    ax2.text(label_x, ys[2], "Impact",   fontsize=8, color=CLR_MUTED,
             ha="left", va="center")
    ax2.text(value_x, ys[2], case["impact"].lower(),
             fontsize=8.5, color=IMPACT_ACCENT[case["impact"]],
             style="italic", ha="left", va="center")

    # Line 4 — SpliceAI score
    ax2.text(label_x, ys[3], "SpliceAI", fontsize=8, color=CLR_MUTED,
             ha="left", va="center")
    ax2.text(value_x, ys[3], f"{case['spliceai']:.2f}",
             fontsize=8.5, color=CLR_TEXT, ha="left", va="center",
             family="monospace")


# ── Build figure ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13.5, 3.6))
gs = GridSpec(1, 3, figure=fig, wspace=0.14,
              left=0.03, right=0.98, top=0.94, bottom=0.06)

for i, case in enumerate(CASES):
    draw_panel(fig, case, gs[0, i])

out_png = OUT_DIR / "vep_schematic_3cases.png"
out_pdf = OUT_DIR / "vep_schematic_3cases.pdf"
fig.savefig(str(out_png), dpi=220, bbox_inches="tight", facecolor="white")
fig.savefig(str(out_pdf),                bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved -> {out_png}")
print(f"Saved -> {out_pdf}")
