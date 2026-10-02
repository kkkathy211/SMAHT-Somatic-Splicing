# -*- coding: utf-8 -*-
# ============================================================================
# Phase 3 - Prediction-first candidate selection
#
# Per-donor sashimi: one row per donor, carrier highlighted, every arc labelled with its read count.
#
# Original location in the analysis project:
#     0727/code/01_tierA_per_donor_sashimi.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
01_tierA_per_donor_sashimi.py
==============================
Per-donor sashimi for Tier A candidates (BAMs already local on Broad):
  1. UNC13A (SMHT039)  chr19:17,639,526  splice_acceptor  VAF 0.055
  2. MEI1   (SMHT039)  chr22:41,763,173  splice_acceptor  VAF 0.060
  3. COPB2  (SMHT016)  chr3:139,359,179  splice_acceptor  VAF 0.049

All three carriers already have 3AD BAM in
`/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/`.  No download needed.

Method: samtools + pysam (no LeafCutter / no pyqtl dependency).
  - MANE Select exons from v47 GTF
  - Per-donor junction reads via pysam CIGAR-N + XS-tag strand filter
  - Per-donor sashimi rows stacked, carrier highlighted (red row background)
  - Junction arcs labeled with actual read counts

Run on Broad:
    use Samtools
    cd /medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/0727update
    python3 code/01_tierA_per_donor_sashimi.py
"""

import shutil
import subprocess
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Arc
import pysam

if not shutil.which("samtools"):
    raise SystemExit("ERROR: samtools not in PATH. Run 'use Samtools' first.")

# ── Paths ────────────────────────────────────────────────────────────────────
BAM_DIR  = Path("/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin")
OUT_DIR  = Path("/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/0727update/output")
GTF_PATH = Path("/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/0720update/gencode.v47.annotation.gtf")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Nature-style palette ─────────────────────────────────────────────────────
COL_CARRIER = "#B91C1C"     # deep red for carrier
COL_NONCARR = "#94A3B8"     # slate for non-carriers
COL_EXON    = "#1F2937"     # near-black exons
COL_INTRON  = "#94A3B8"
COL_TEXT    = "#111827"
COL_MUTED   = "#6B7280"
BG_CARRIER  = "#FEF2F2"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42,
    "axes.linewidth": 0.6,
})

# ── Tier-A candidates ────────────────────────────────────────────────────────
CANDIDATES = [
    {"gene": "UNC13A", "donor": "SMHT039", "spliceai": 1.00, "vaf": 0.055,
     "var_pos": 17639526, "vep": "splice_acceptor (HIGH)"},
    {"gene": "MEI1",   "donor": "SMHT039", "spliceai": 0.99, "vaf": 0.060,
     "var_pos": 41763173, "vep": "splice_acceptor (HIGH)"},
    {"gene": "COPB2",  "donor": "SMHT016", "spliceai": 0.98, "vaf": 0.049,
     "var_pos": 139359179, "vep": "splice_acceptor (HIGH)"},
]

# ── Filter tuning ────────────────────────────────────────────────────────────
MIN_JCT_READS   = 5    # skip junctions with fewer than this in any donor
MIN_JCT_LENGTH  = 25   # skip tiny (likely artifact) junctions
ZOOM_BUFFER     = 500  # bp padding around visible exons


# ── Parse GTF for MANE Select exons ──────────────────────────────────────────
def parse_mane_exons(gtf_path, gene_names):
    genes = {}
    mane_tid_by_gene = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript": continue
            attr = parts[8]
            if 'gene_name "' not in attr: continue
            gname = attr.split('gene_name "')[1].split('"')[0]
            if gname not in gene_names: continue
            if 'tag "MANE_Select"' not in attr: continue
            tid = attr.split('transcript_id "')[1].split('"')[0]
            mane_tid_by_gene[gname] = tid
            genes[gname] = {"chr": parts[0], "strand": parts[6], "tid": tid,
                            "gene_start": int(parts[3]), "gene_end": int(parts[4]),
                            "exons": []}
    mane_tids = set(mane_tid_by_gene.values())
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon": continue
            attr = parts[8]
            tid = attr.split('transcript_id "')[1].split('"')[0]
            if tid not in mane_tids: continue
            gname = attr.split('gene_name "')[1].split('"')[0]
            if gname in genes:
                genes[gname]["exons"].append((int(parts[3]), int(parts[4])))
    for g in genes:
        genes[g]["exons"].sort()
    return genes


# ── BAM helpers ──────────────────────────────────────────────────────────────
def pick_samples_for_tissue(bam_dir, tissue_tag="3AD"):
    all_bams = sorted(f for f in bam_dir.glob("*.bam") if not str(f).endswith(".bai"))
    picked = []
    for bam in all_bams:
        if f"-{tissue_tag}-" not in bam.name: continue
        if "-A101-" not in bam.name: continue
        parts = bam.name.split("-")
        picked.append({"name": bam.stem, "donor": parts[0], "bam": bam})
    return picked


def junctions_stranded(bam_path, chrom, start, end, gene_strand):
    """Junctions via pysam CIGAR N with XS-tag strand filter."""
    jcts = defaultdict(int)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if (read.is_unmapped or read.is_secondary or
                read.is_supplementary or read.is_duplicate or read.is_qcfail):
                continue
            if read.mapping_quality < 1: continue
            xs = read.get_tag("XS") if read.has_tag("XS") else None
            if xs is None or xs != gene_strand: continue
            rpos = read.reference_start
            for op, length in (read.cigartuples or []):
                if op == 0: rpos += length
                elif op == 3:
                    if start <= rpos and (rpos+length) <= end:
                        jcts[(rpos, rpos+length)] += 1
                    rpos += length
                elif op in (2, 7, 8): rpos += length
    return jcts


def compute_zoom(exons, var_pos, target_kb=15):
    """Zoom around variant: include ≥2 exons on each side, cap ~target_kb wide."""
    all_ex = sorted(exons)
    # Find exon indices flanking the variant
    left_idx  = None
    right_idx = None
    for i, (s, e) in enumerate(all_ex):
        if e < var_pos:
            left_idx = i
        elif s > var_pos and right_idx is None:
            right_idx = i
    # Grab 2 exons on each side (or fewer if at edge)
    lo_idx = max(0, (left_idx or 0) - 1)
    hi_idx = min(len(all_ex) - 1, (right_idx or len(all_ex) - 1) + 1)
    lo = all_ex[lo_idx][0] - ZOOM_BUFFER
    hi = all_ex[hi_idx][1] + ZOOM_BUFFER
    # Cap width at target_kb
    if hi - lo > target_kb * 1000:
        lo = var_pos - target_kb * 500
        hi = var_pos + target_kb * 500
    return lo, hi


# ── Parse GTF once ───────────────────────────────────────────────────────────
print("Parsing v47 GTF for MANE Select exons ...")
GENE_STRUCTS = parse_mane_exons(GTF_PATH, [c["gene"] for c in CANDIDATES])
for g in [c["gene"] for c in CANDIDATES]:
    if g in GENE_STRUCTS:
        s = GENE_STRUCTS[g]
        print(f"  {g}: {s['chr']}:{s['gene_start']}-{s['gene_end']} "
              f"({s['strand']} strand, {len(s['exons'])} exons)")
    else:
        print(f"  {g}: NOT FOUND in v47 GTF !!")


# ── Load BAMs ────────────────────────────────────────────────────────────────
samples = pick_samples_for_tissue(BAM_DIR, "3AD")
print(f"\nBAMs available for 3AD (A101): {len(samples)}")
for s in samples:
    print(f"  [{s['donor']}]  {s['bam'].name}")


# ── Per-candidate figure ─────────────────────────────────────────────────────
for cand in CANDIDATES:
    gname   = cand["gene"]
    donor   = cand["donor"]
    var_pos = cand["var_pos"]

    print(f"\n{'='*70}")
    print(f"{gname}  carrier={donor}  SpliceAI={cand['spliceai']:.2f}  VAF={cand['vaf']:.3f}")
    print("="*70)

    if gname not in GENE_STRUCTS:
        print(f"  !! skipping (gene not in v47 GTF)")
        continue
    gs = GENE_STRUCTS[gname]
    chrom = gs["chr"]

    # Zoom
    lo, hi = compute_zoom(gs["exons"], var_pos, target_kb=15)
    print(f"  Zoom: {chrom}:{lo:,}-{hi:,}  ({(hi-lo)/1000:.1f} kb)")

    # Check carrier
    car_samps = [s for s in samples if s["donor"] == donor]
    if not car_samps:
        print(f"  !! carrier {donor} BAM missing, skipping")
        continue

    # Extract junctions per sample
    per_donor_jcts = {}
    for s in samples:
        print(f"  {s['donor']} ...", end=" ", flush=True)
        jc = junctions_stranded(s["bam"], chrom, lo, hi, gs["strand"])
        per_donor_jcts[s["name"]] = jc
        print(f"n_jcts={len(jc)}, total_reads={sum(jc.values())}")

    # Collect junctions to display: any junction with ≥MIN_JCT_READS in ≥1 donor
    all_jcts = set()
    for name, jc in per_donor_jcts.items():
        for k, v in jc.items():
            if v >= MIN_JCT_READS and (k[1] - k[0]) >= MIN_JCT_LENGTH:
                all_jcts.add(k)
    all_jcts = sorted(all_jcts)
    print(f"  Displaying {len(all_jcts)} junctions (≥{MIN_JCT_READS} reads in ≥1 donor)")

    # ── Draw figure ──────────────────────────────────────────────────────────
    N = len(samples)
    row_h = 0.38
    fig_h = row_h * N + 1.8
    fig_w = 15
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs_fig = fig.add_gridspec(N + 1, 1, height_ratios=[1] * N + [0.9],
                              hspace=0.15, left=0.10, right=0.98,
                              top=0.92, bottom=0.06)

    # Max read count for arc height scaling
    max_reads_global = max(
        [per_donor_jcts[s["name"]].get(k, 0) for s in samples for k in all_jcts]
        + [1]
    )

    for i, s in enumerate(samples):
        ax = fig.add_subplot(gs_fig[i, 0])
        is_carrier = (s["donor"] == donor)
        color = COL_CARRIER if is_carrier else COL_NONCARR

        if is_carrier:
            ax.set_facecolor(BG_CARRIER)

        # Exon boxes (compact, near baseline)
        box_h = 0.28
        for es, ee in gs["exons"]:
            if ee < lo or es > hi: continue
            es_c, ee_c = max(es, lo), min(ee, hi)
            ax.add_patch(Rectangle((es_c, -box_h / 2), ee_c - es_c, box_h,
                                   facecolor=COL_EXON, edgecolor="none", zorder=3))
        # Intron baseline
        vis_exons = [(es, ee) for es, ee in gs["exons"] if not (ee < lo or es > hi)]
        if len(vis_exons) >= 2:
            ax.plot([vis_exons[0][1], vis_exons[-1][0]], [0, 0],
                    color=COL_INTRON, lw=0.5, zorder=1)

        # Junction arcs
        for jk in all_jcts:
            n_reads = per_donor_jcts[s["name"]].get(jk, 0)
            if n_reads == 0: continue
            j0, j1 = jk
            w = j1 - j0
            h_arc = 0.28 + 1.20 * (min(n_reads, max_reads_global) / max_reads_global) ** 0.55
            lw = 0.5 + 2.8 * (min(n_reads, max_reads_global) / max_reads_global) ** 0.45
            arc = Arc(((j0 + j1) / 2, 0), w, h_arc * 2,
                      angle=0, theta1=0, theta2=180,
                      color=color, lw=lw, alpha=0.85, zorder=5)
            ax.add_patch(arc)
            # Label with read count
            ax.text((j0 + j1) / 2, h_arc + 0.03, f"{n_reads:,}",
                    fontsize=7, fontweight="bold", color=color,
                    ha="center", va="bottom", zorder=10,
                    bbox=dict(boxstyle="round,pad=0.08", facecolor="white",
                              edgecolor=color, lw=0.3, alpha=0.95))

        # Sample label on the left
        ax.text(lo - (hi - lo) * 0.015, 0, s["donor"],
                fontsize=10, fontweight="bold" if is_carrier else "normal",
                color=color if is_carrier else COL_TEXT,
                ha="right", va="center", zorder=10)
        if is_carrier:
            ax.text(lo - (hi - lo) * 0.075, 0, "★",
                    fontsize=13, color=COL_CARRIER, ha="right", va="center")

        # Variant vertical line
        ax.axvline(var_pos, color="#111827", ls=":", lw=0.6, alpha=0.5, zorder=2)

        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.4, 1.9)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ("top", "right", "bottom", "left"):
            ax.spines[sp].set_visible(False)

    # ── Bottom gene track ────────────────────────────────────────────────────
    ax_g = fig.add_subplot(gs_fig[N, 0])
    box_h = 0.55
    # Exon labels: number by MANE exon order (strand-aware)
    strand = gs["strand"]
    all_exons_sorted = sorted(gs["exons"])
    n_ex_total = len(all_exons_sorted)
    for i, (es, ee) in enumerate(all_exons_sorted, start=1):
        if ee < lo or es > hi: continue
        es_c, ee_c = max(es, lo), min(ee, hi)
        exon_num = i if strand == "+" else n_ex_total - i + 1
        ax_g.add_patch(Rectangle((es_c, -box_h / 2), ee_c - es_c, box_h,
                                 facecolor=COL_EXON, edgecolor="none", zorder=3))
        ax_g.text((es_c + ee_c) / 2, box_h / 2 + 0.12, f"e{exon_num}",
                  fontsize=9, fontweight="bold", ha="center", va="bottom",
                  color=COL_EXON)
    # Intron line
    vis_ex = [(es, ee) for es, ee in all_exons_sorted if not (ee < lo or es > hi)]
    if len(vis_ex) >= 2:
        ax_g.plot([vis_ex[0][1], vis_ex[-1][0]], [0, 0],
                  color=COL_MUTED, lw=0.9)

    # Variant marker with label
    ax_g.axvline(var_pos, color="#B91C1C", ls="--", lw=1.2, alpha=0.7)
    ax_g.text(var_pos, -0.9,
              f"variant\n{chrom}:{var_pos:,}",
              fontsize=8.5, color="#B91C1C", fontweight="bold",
              ha="center", va="top",
              bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                        edgecolor="#B91C1C", lw=0.6))

    ax_g.set_xlim(lo, hi)
    ax_g.set_ylim(-1.5, 1.2)
    ax_g.set_xticks([]); ax_g.set_yticks([])
    for sp in ("top", "right", "bottom", "left"):
        ax_g.spines[sp].set_visible(False)

    # ── Title & legend ───────────────────────────────────────────────────────
    fig.suptitle(
        f"{gname}  |  per-donor junction usage  (3AD Skin Calf, n={N} samples)  |  "
        f"carrier {donor} ★  |  SpliceAI = {cand['spliceai']:.2f}  |  "
        f"VEP HIGH {cand['vep']}  |  VAF = {cand['vaf']:.3f}",
        fontsize=11, fontweight="bold", x=0.05, ha="left", y=0.98,
        color=COL_TEXT)

    legend_ax = fig.add_axes([0.83, 0.94, 0.15, 0.04])
    legend_ax.axis("off")
    legend_ax.text(0.0, 0.5,
                   "arc label = # supporting reads\n"
                   "arc thickness ∝ reads",
                   transform=legend_ax.transAxes, fontsize=8, color=COL_MUTED,
                   style="italic", ha="left", va="center")

    tag = f"{gname}_{donor}_3AD_per_donor_sashimi"
    out_png = OUT_DIR / (tag + ".png")
    out_pdf = OUT_DIR / (tag + ".pdf")
    fig.savefig(str(out_png), dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(str(out_pdf), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved -> {out_png}")

print(f"\n{'='*70}")
print(f"Done. Output in {OUT_DIR}")
print("="*70)
