# -*- coding: utf-8 -*-
# ============================================================================
# Phase 3 - Prediction-first candidate selection
#
# Pileup via pyqtl (Broad internal). Template; edit the CASE dict per gene.
#
# Original location in the analysis project:
#     0716/code/16_laura_method_generic_pileup.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
16_laura_method_generic_pileup.py
==================================
Generic Laura-style pileup — for any new candidate gene from the VEP splice list.
Uses Laura's EXACT method (as in SMaHT_P25_splicing_pileup_plots.ipynb).

Key differences from earlier custom scripts:
  - qtl.annotation.Annotation to parse the FULL gene structure from GTF
  - qtl.pileup.samtools_depth (library wrapper) — matches Laura's normalization
  - qtl.pileup.regtools_extract_junctions with strand filtering
  - Pre-computed libsize file (libsize_skin_calf.txt)
  - qtl.pileup.plot for the visualization (handles isoform layout, intron
    compression, arc rendering — battle-tested)

Configure the CASE dict at the top for the gene you want to plot.

Run on Broad:
    cd /medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/0716update
    python3 16_laura_method_generic_pileup.py
"""

import copy
import os
import subprocess
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Laura's qtl library (should already be installed on Broad)
import qtl.annotation
import qtl.pileup

plt.rcParams.update({
    "font.family": "Helvetica",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

# ── Configure here for each new candidate ────────────────────────────────────
# EDIT this dict per gene. Example uses TP53 in SMHT040 (from the top-25 list).
CASE = {
    "gene":            "TP53",
    "carrier_sample":  "SMHT040-3AD-...",   # full sample id — set after checking fc_df
    "tissue_folder":   "3AD_Skin_Calf",
    # SpliceAI + VEP metadata for the title
    "spliceai":        1.00,
    "vep":             "splice_acceptor_variant (HIGH)",
    "var_pos":         7675237,
    # ymax + manual carrier scale — tune after first render if needed
    "ymax":            5,
    "carrier_scale":   1.0,
    # optional isoform whitelist (transcript IDs) — leave empty to auto-pick
    # (MANE Select + top expressed protein_coding). Fill in after first inspect.
    "isoform_ids":     [],
}

# Paths on Broad
GTF_PATH        = "/home/ldomenec/references/gencode.v47.annotation.gtf"
SAMPLE_TSV      = "sample_tissue.tsv"       # {sample -> tissue folder}
LIBSIZE_TSV     = f"libsize_{CASE['tissue_folder'].lower()}.txt"   # {sample -> libsize}
BAM_INDEX_DIR   = "/tmp/bam_indexes"
OUT_DIR         = "."


# ── 1. Load annotation ───────────────────────────────────────────────────────
print(f"Loading GTF ...")
annot = qtl.annotation.Annotation(GTF_PATH)

gene = annot.get_gene(CASE["gene"])
print(f"\nGene {gene.name}  ({gene.id})  {gene.chr}:{gene.start_pos}-{gene.end_pos}  "
      f"strand={gene.strand}")
print(f"Available transcripts ({len(gene.transcripts)}):")
for t in gene.transcripts:
    print(f"  {t.id}  {t.name}")

# Deep copy so we can subset transcripts
g = copy.deepcopy(gene)

# Isoform selection — either whitelist or default filter
if CASE["isoform_ids"]:
    kept = [t for t in g.transcripts if t.id in CASE["isoform_ids"]]
    # Reverse so MANE Select is at top of stack (Laura convention)
    g.set_transcripts(kept[::-1])
    print(f"\nUsing {len(kept)} whitelisted isoforms")
else:
    # Default: exclude junk biotypes (Laura's default in gene.plot())
    kept = [t for t in g.transcripts
            if t.type not in ("retained_intron",
                              "protein_coding_CDS_not_defined",
                              "nonsense_mediated_decay")]
    g.set_transcripts(kept[::-1])
    print(f"\nUsing {len(kept)} isoforms after biotype filter "
          f"(exclude retained_intron/NMD/etc)")


# ── 2. Load BAM list + tissue filter ─────────────────────────────────────────
print(f"\nListing BAMs from gs://smaht-p25-bam-files/ ...")
result = subprocess.run(["gsutil", "ls", "gs://smaht-p25-bam-files/"],
                        capture_output=True, text=True)
import re
bam_files = [l for l in result.stdout.strip().split("\n") if l.endswith(".bam")]
fc_df = pd.DataFrame({"bam_file": bam_files})
fc_df["sample_id"] = fc_df["bam_file"].apply(
    lambda x: re.sub(r"-SMAFI\w+.*", "", os.path.basename(x))
)
fc_df = fc_df.set_index("sample_id")

tissue_df = pd.read_csv(SAMPLE_TSV, sep="\t", index_col="sample")
fc_df = fc_df.join(tissue_df)
fc_df = fc_df[fc_df["tissue"] == CASE["tissue_folder"]]
print(f"BAMs in {CASE['tissue_folder']}: {len(fc_df)}")
for sid in fc_df.index:
    marker = "  ⭐ CARRIER" if sid == CASE["carrier_sample"] else ""
    print(f"  {sid}{marker}")


# ── 3. Extract coverage + junctions (Laura's wrappers) ───────────────────────
os.makedirs(BAM_INDEX_DIR, exist_ok=True)
region_str = f"{g.chr}:{g.start_pos}-{g.end_pos}"
print(f"\nRegion: {region_str}")

qtl.refresh_gcs_token()
bam_s = fc_df["bam_file"]

print("Running samtools depth ...")
pileups_df = qtl.pileup.samtools_depth(region_str, bam_s,
                                       bam_index_dir=BAM_INDEX_DIR, num_threads=8)
print("Running regtools junctions extract ...")
junctions_df = qtl.pileup.regtools_extract_junctions(region_str, bam_s,
                                                     bam_index_dir=BAM_INDEX_DIR,
                                                     strand="XS", num_threads=8)

# Strand filter — only junctions matching gene's strand
junctions_df = junctions_df[
    junctions_df.index.map(lambda x: x.rsplit(":", 1)[1]) == gene.strand
].copy()
junctions_df.index = junctions_df.index.map(lambda x: x.rsplit(":", 1)[0])


# ── 4. Normalize by libsize ──────────────────────────────────────────────────
libsize_df = pd.read_csv(LIBSIZE_TSV, sep="\t", index_col="sample")
libsize_df["libsize"] = libsize_df["libsize"].astype(str).str.replace(",", "").astype(int)
fc_df = fc_df.join(libsize_df)
libsize_s = fc_df["libsize"]

pileups_rpm = qtl.pileup.norm_pileups(pileups_df, libsize_s,
                                       covariates_df=None,
                                       id_map=lambda x: x.split("-")[0])
junctions_rpm = qtl.pileup.norm_pileups(junctions_df, libsize_s,
                                         covariates_df=None,
                                         id_map=lambda x: x.split("-")[0])


# ── 5. Carrier vs Non-carrier groups ─────────────────────────────────────────
fc_df["group"] = "Non-carrier"
fc_df.loc[CASE["carrier_sample"], "group"] = "Carrier"
cohort_s = fc_df["group"]
cohort_s.index = cohort_s.index.map(lambda x: x.split("-")[0])

groups = ["Carrier", "Non-carrier"]
group_df = pd.concat([pileups_rpm[cohort_s[cohort_s == i].index].mean(axis=1).rename(i)
                      for i in groups], axis=1)
junc_group_df = pd.concat([junctions_rpm[cohort_s[cohort_s == i].index].mean(axis=1).rename(i)
                            for i in groups], axis=1)
junc_group_df.index = junc_group_df.index.map(lambda x: x.split(":")[1])

# Optional manual carrier scaling
if CASE["carrier_scale"] != 1.0:
    group_df["Carrier"] *= CASE["carrier_scale"]
    print(f"Applied manual carrier scale: * {CASE['carrier_scale']}")


# ── 6. Plot ──────────────────────────────────────────────────────────────────
colors = {"Carrier":     plt.cm.Reds(0.6),
          "Non-carrier": plt.cm.Blues(0.6)}

max_intron = 100
g.set_plot_coords(max_intron=max_intron)

ax = qtl.pileup.plot(group_df, g, outline=True,
                     junctions_df=junc_group_df[
                         junc_group_df.max(1) > junc_group_df.max().max() * 0.01
                     ],
                     colors=colors,
                     lw=1.25, db=0.25, dr=1.5,
                     ymax=CASE["ymax"],
                     aw=8, ah=2,
                     max_intron=max_intron)

ax[0].set_title(
    f"{gene.name}  ({CASE['carrier_sample'].split('-')[0]}, {CASE['tissue_folder']})"
    f"  |  SpliceAI={CASE['spliceai']:.2f}  |  {CASE['vep']}"
    f"  |  variant chr{gene.chr[3:]}:{CASE['var_pos']:,}",
    fontsize=11)
ax[1].set_ylabel("Isoforms", fontsize=11, labelpad=25)

out_pdf = os.path.join(OUT_DIR, f"{gene.name}_{CASE['carrier_sample'].split('-')[0]}.pileup.pdf")
out_png = out_pdf.replace(".pdf", ".png")
plt.savefig(out_pdf, bbox_inches="tight")
plt.savefig(out_png, dpi=180, bbox_inches="tight")
print(f"\nSaved -> {out_pdf}")
print(f"Saved -> {out_png}")
