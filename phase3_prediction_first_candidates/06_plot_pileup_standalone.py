# -*- coding: utf-8 -*-
# ============================================================================
# Phase 3 - Prediction-first candidate selection
#
# Same figure as 05 with no pyqtl dependency: samtools depth + pysam with XS-tag strand filter.
#
# Original location in the analysis project:
#     0720/code/19_tier1_3AD_custom_pileup.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
19_tier1_3AD_custom_pileup.py
==============================
Fallback pileup script (no pyqtl dependency) for the 11 Tier-1 3AD candidates.
Uses only samtools + pysam + matplotlib — all pre-installed on Broad.

Method (mirrors Laura's method as closely as possible without pyqtl):
  1. Parse gencode.v47 GTF locally for MANE Select exons per gene
  2. samtools depth for coverage curves
  3. pysam CIGAR-N for junction reads, WITH strand filter (XS tag matching
     gene.strand) — same rigor as Laura's regtools call
  4. samtools idxstats for library size normalization (per BAM)
  5. matplotlib: coverage curves + junction arcs + isoform stack

Layout follows Laura's convention: coverage on top, isoforms on bottom,
Carrier red / Non-carrier blue, arcs sized by junction read count.

Run on Broad:
    use Samtools
    cd /medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/0720update
    python3 code/19_tier1_3AD_custom_pileup.py
"""

import shutil
import subprocess
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Rectangle
import pysam

if not shutil.which("samtools"):
    raise SystemExit("ERROR: samtools not in PATH. Run 'use Samtools' first.")

# ── Paths ────────────────────────────────────────────────────────────────────
BAM_DIR   = Path("/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin")
OUT_DIR   = Path("/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/0720update/output")
GTF_PATH  = Path("/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/0720update/gencode.v47.annotation.gtf")
OUT_DIR.mkdir(parents=True, exist_ok=True)

LIBSIZE_CACHE = {}

COL_CARRIER = plt.cm.Reds(0.6)
COL_NONCARR = plt.cm.Blues(0.6)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size":   11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})

# ── Tier-1 3AD candidates ────────────────────────────────────────────────────
CANDIDATES = [
    {"gene":"MEI1",    "donor":"SMHT039", "spliceai":0.99, "vaf":0.060, "var_pos":41763173,  "vep":"splice_acceptor (HIGH)"},
    {"gene":"UNC13A",  "donor":"SMHT039", "spliceai":1.00, "vaf":0.055, "var_pos":17639526,  "vep":"splice_acceptor (HIGH)"},
    {"gene":"COPB2",   "donor":"SMHT016", "spliceai":0.98, "vaf":0.049, "var_pos":139359179, "vep":"splice_acceptor (HIGH)"},
    {"gene":"PLEKHA3", "donor":"SMHT015", "spliceai":0.99, "vaf":0.041, "var_pos":178503759, "vep":"splice_acceptor (HIGH)"},
    {"gene":"FAT1",    "donor":"SMHT039", "spliceai":0.90, "vaf":0.036, "var_pos":186609182, "vep":"splice_donor (HIGH)"},
    {"gene":"KAT2A",   "donor":"SMHT015", "spliceai":0.99, "vaf":0.036, "var_pos":42117677,  "vep":"splice_donor (HIGH)"},
    {"gene":"ZDHHC2",  "donor":"SMHT040", "spliceai":0.99, "vaf":0.033, "var_pos":17198379,  "vep":"splice_acceptor (HIGH)"},
    {"gene":"TMC7",    "donor":"SMHT016", "spliceai":0.99, "vaf":0.032, "var_pos":19016599,  "vep":"splice_donor (HIGH)"},
    {"gene":"NOTCH1",  "donor":"SMHT018", "spliceai":1.00, "vaf":0.032, "var_pos":136513020, "vep":"splice_donor (HIGH)"},
    {"gene":"TP53",    "donor":"SMHT040", "spliceai":1.00, "vaf":0.031, "var_pos":7675237,   "vep":"splice_acceptor (HIGH)"},
    {"gene":"SPSB4",   "donor":"SMHT039", "spliceai":0.99, "vaf":0.031, "var_pos":141051994, "vep":"splice_donor (HIGH)"},
]


# ── 1. Parse v47 GTF for MANE Select exon structure of each candidate ───────
def parse_mane_exons(gtf_path, gene_names):
    """Return {gene_name: {'chr':..., 'strand':..., 'exons':[(s,e),...],
                            'gene_start':..., 'gene_end':..., 'tid':...}}"""
    genes = {}
    mane_tid_by_gene = {}
    # First pass: find MANE Select transcript IDs for each gene
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
            genes[gname] = {
                "chr":      parts[0],
                "strand":   parts[6],
                "tid":      tid,
                "gene_start": int(parts[3]),
                "gene_end":   int(parts[4]),
                "exons":    [],
            }
    # Second pass: pull exons of those MANE transcripts
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

print("Parsing v47 GTF for MANE Select exons ...")
gene_names_needed = list({c["gene"] for c in CANDIDATES})
GENE_STRUCTS = parse_mane_exons(GTF_PATH, gene_names_needed)
print(f"  Got structure for {len(GENE_STRUCTS)}/{len(gene_names_needed)} genes")
for g in gene_names_needed:
    if g in GENE_STRUCTS:
        s = GENE_STRUCTS[g]
        print(f"    {g}: {s['chr']}:{s['gene_start']}-{s['gene_end']} ({s['strand']}), "
              f"{len(s['exons'])} exons")


# ── 2. Helpers (same as previous laura-style scripts) ────────────────────────
def pick_samples_for_tissue(bam_dir, tissue_tag="3AD"):
    """Return list of {name, donor, bam} matching tissue tag and A101 library."""
    all_bams = sorted(f for f in bam_dir.glob("*.bam") if not str(f).endswith(".bai"))
    picked = []
    for bam in all_bams:
        if f"-{tissue_tag}-" not in bam.name: continue
        if "-A101-" not in bam.name: continue
        parts = bam.name.split("-")
        picked.append({"name": bam.stem, "donor": parts[0], "bam": bam})
    return picked


def libsize(bam_path):
    key = bam_path.name
    if key in LIBSIZE_CACHE: return LIBSIZE_CACHE[key]
    r = subprocess.run(["samtools","idxstats",str(bam_path)],
                       capture_output=True, text=True, check=True)
    total = 0
    for l in r.stdout.strip().split("\n"):
        p = l.split("\t")
        if len(p) >= 3:
            try: total += int(p[2])
            except ValueError: pass
    LIBSIZE_CACHE[key] = max(total, 1)
    return LIBSIZE_CACHE[key]


def samtools_depth(bam_path, chrom, start, end):
    n = end - start
    cov = np.zeros(n, dtype=np.float64)
    region = f"{chrom}:{start+1}-{end}"
    try:
        r = subprocess.run(["samtools","depth","-a","-r",region,str(bam_path)],
                           capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"    samtools depth failed for {bam_path.name}: {e}")
        return cov
    for l in r.stdout.strip().split("\n"):
        if not l: continue
        p = l.split("\t")
        pos, depth = int(p[1]), int(p[2])
        idx = pos - 1 - start
        if 0 <= idx < n: cov[idx] = depth
    return cov


def junctions_stranded(bam_path, chrom, start, end, gene_strand):
    """Junctions with strand filter — only keep reads whose XS tag matches
    gene_strand (matches Laura's method)."""
    jcts = defaultdict(int)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if (read.is_unmapped or read.is_secondary or
                read.is_supplementary or read.is_duplicate or read.is_qcfail):
                continue
            if read.mapping_quality < 1: continue
            # Strand filter via XS tag
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


def draw_isoform_stack(ax, exons, lo, hi):
    y = 0.5; box_h = 0.55
    visible = [(s,e) for s,e in exons if e >= lo and s <= hi]
    if visible:
        line_lo = max(min(s for s,_ in visible), lo)
        line_hi = min(max(e for _,e in visible), hi)
        ax.plot([line_lo, line_hi], [y, y], color="black", lw=0.8, zorder=2)
        for s, e in visible:
            ax.add_patch(Rectangle((max(s,lo), y-box_h/2),
                                   min(e,hi)-max(s,lo), box_h,
                                   facecolor="black", zorder=3))
    ax.set_xlim(lo, hi); ax.set_ylim(-0.3, 1.3)
    ax.set_yticks([]); ax.set_xticks([])
    ax.set_ylabel("Isoform", fontsize=10, labelpad=12, rotation=90,
                  va="center", ha="center")
    for sp in ("top","right","bottom","left"):
        ax.spines[sp].set_visible(False)


# ── 3. Loop over candidates ──────────────────────────────────────────────────
samples = pick_samples_for_tissue(BAM_DIR, "3AD")
print(f"\nBAMs available for 3AD (A101): {len(samples)}")
for s in samples:
    print(f"  [{s['donor']}]  {s['bam'].name}")

for i, cand in enumerate(CANDIDATES, start=1):
    gname   = cand["gene"]
    donor   = cand["donor"]
    print(f"\n{'='*70}")
    print(f"[{i}/{len(CANDIDATES)}]  {gname}  carrier={donor}  "
          f"(SpliceAI={cand['spliceai']:.2f}, VAF={cand['vaf']:.3f})")
    print("="*70)

    if gname not in GENE_STRUCTS:
        print(f"  !! {gname} not found in GTF, skipping")
        continue
    gs = GENE_STRUCTS[gname]
    chrom = gs["chr"]
    lo, hi = gs["gene_start"], gs["gene_end"]

    # Check carrier BAM
    car_samps = [s for s in samples if s["donor"] == donor]
    if not car_samps:
        print(f"  !! Carrier {donor}-3AD BAM not in folder. Download it first:")
        print(f"     gcloud storage cp \\")
        print(f"       gs://smaht-p25-bam-files/{donor}-3AD-*.bam \\")
        print(f"       gs://smaht-p25-bam-files/{donor}-3AD-*.bam.bai \\")
        print(f"       {BAM_DIR}/")
        continue
    nc_samps = [s for s in samples if s["donor"] != donor]
    if not nc_samps:
        print(f"  !! No non-carrier BAMs. Skipping.")
        continue

    # Per-sample coverage + junctions
    cov_rpm, jcts_rpm = {}, {}
    for s in samples:
        print(f"  {s['name']} ...", end=" ", flush=True)
        lib = libsize(s["bam"])
        cov = samtools_depth(s["bam"], chrom, lo, hi)
        jc  = junctions_stranded(s["bam"], chrom, lo, hi, gs["strand"])
        cov_rpm[s["name"]]  = cov / lib * 1e6
        jcts_rpm[s["name"]] = {k: v/lib*1e6 for k,v in jc.items()}
        print(f"max_RPM={cov_rpm[s['name']].max():.2f}, "
              f"n_jcts={len(jc)}")

    car_cov = np.stack([cov_rpm[s["name"]] for s in car_samps]).mean(axis=0)
    nc_cov  = np.stack([cov_rpm[s["name"]] for s in nc_samps]).mean(axis=0)

    all_keys = set()
    for name in cov_rpm: all_keys.update(jcts_rpm[name].keys())
    car_jct = {k: float(np.mean([jcts_rpm[s["name"]].get(k,0.0) for s in car_samps]))
               for k in all_keys}
    nc_jct  = {k: float(np.mean([jcts_rpm[s["name"]].get(k,0.0) for s in nc_samps]))
               for k in all_keys}

    # Plot
    n_iso = 1
    iso_h = 0.6
    fig, (ax, ax_iso) = plt.subplots(
        2, 1, figsize=(13, 3.5),
        gridspec_kw={"height_ratios": [3.0, iso_h], "hspace": 0.05},
        sharex=False,
    )
    x = np.arange(lo, hi)
    ax.plot(x, car_cov, lw=1.4, color=COL_CARRIER, label=f"Carrier {donor}")
    ax.plot(x, nc_cov,  lw=1.4, color=COL_NONCARR,
            label=f"Non-carrier mean (n={len(nc_samps)})")

    y_top = max(float(car_cov.max()), float(nc_cov.max()), 1.0)
    arc_h_max = y_top * 0.55

    def draw_arcs(jdict, color, alpha=0.85):
        items = [(k,v) for k,v in jdict.items() if v>0 and (k[1]-k[0])>25]
        if not items: return
        max_v = max(v for _,v in items)
        cutoff = max_v * 0.02
        items = [(k,v) for k,v in items if v >= cutoff]
        for (j0,j1), v in items:
            w  = j1-j0
            h  = arc_h_max * (v/max_v) * 2
            lw = 0.6 + 1.8 * (v/max_v)
            arc = Arc(((j0+j1)/2, 0), w, h, angle=0, theta1=0, theta2=180,
                      color=color, lw=lw, alpha=alpha, zorder=5)
            ax.add_patch(arc)
    draw_arcs(nc_jct, COL_NONCARR)
    draw_arcs(car_jct, COL_CARRIER)

    ax.set_xlim(lo, hi); ax.set_ylim(bottom=0)
    ax.set_ylabel("Mean RPM", fontsize=11)
    ax.set_xticks([])
    ax.legend(fontsize=9.5, frameon=True, loc="upper right",
              edgecolor="lightgray", facecolor="white")

    # Variant marker
    ymax_now = ax.get_ylim()[1]
    ax.axvline(cand["var_pos"], color="black", ls="--", lw=1.1, alpha=0.7, zorder=4)
    ax.text(cand["var_pos"], ymax_now*0.97, "variant",
            ha="center", va="top", fontsize=9, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.20", facecolor="white",
                      edgecolor="black", lw=0.5), zorder=11)

    ax.set_title(
        f"{gname}  ({donor}, Skin Calf)  |  SpliceAI={cand['spliceai']:.2f}  |  "
        f"{cand['vep']}  |  VAF={cand['vaf']:.3f}  |  "
        f"variant {chrom}:{cand['var_pos']:,}",
        fontsize=10, fontweight="bold", loc="left", color="#1F2937", pad=10)

    draw_isoform_stack(ax_iso, gs["exons"], lo, hi)
    ax_iso.axvline(cand["var_pos"], color="black", ls="--", lw=0.8, alpha=0.6, zorder=4)

    fig.tight_layout()
    tag = f"{gname}_{donor}_3AD_custom"
    out_png = OUT_DIR / (tag + ".png")
    out_pdf = OUT_DIR / (tag + ".pdf")
    fig.savefig(str(out_png), dpi=160, bbox_inches="tight")
    fig.savefig(str(out_pdf), bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved -> {out_png}")

print(f"\n{'='*70}")
print(f"All done. Output in {OUT_DIR}")
print("="*70)
