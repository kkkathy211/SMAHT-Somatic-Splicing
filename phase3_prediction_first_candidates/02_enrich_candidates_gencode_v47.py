# -*- coding: utf-8 -*-
# ============================================================================
# Phase 3 - Prediction-first candidate selection
#
# Enrich candidates with GENCODE v47: MANE Select transcript, exon context, pileup region string.
#
# Original location in the analysis project:
#     0716/code/17_enrich_candidates_v47.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
17_enrich_candidates_v47.py
============================
Enrich the VEP splice candidates (SpliceAI ≥ 0.2 + Consequence contains 'splice')
with gencode.v47 gene annotation so each row is a fully-prepared spec for
running Laura's Laura-style pileup (16_laura_method_generic_pileup.py).

Adds per gene:
  - v47 gene_id, chr, start, end, strand
  - v47 MANE Select transcript id (for isoform whitelist)
  - Number of transcripts (protein-coding, non-NMD)
  - Gene length (kb)

Sorts by priority (has RNA-seq, not yet scanned, SpliceAI, IMPACT) and outputs
a single clean CSV / TSV.
"""

from pathlib import Path
import pandas as pd

GTF = "/Users/spectremac/Desktop/Differential_splicing/0720/gencode.v47.annotation.gtf"
INPUT_CSV = "/Users/spectremac/Desktop/Differential_splicing/0716/output/vep_splice_candidates.csv"
OUT_CSV   = "/Users/spectremac/Desktop/Differential_splicing/0716/output/vep_splice_candidates_v47_enriched.csv"
OUT_TSV   = OUT_CSV.replace(".csv", ".tsv")

# ── 1. Parse GTF for gene + MANE Select transcript info ─────────────────────
print("Parsing gencode.v47 GTF (gene + MANE Select transcript rows) ...")
gene_info  = {}    # gene_name -> dict
trans_by_g = {}    # gene_name -> list of (tid, biotype, mane_select_bool)

with open(GTF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        feature = parts[2]
        chrom, start, end, strand, attr = parts[0], int(parts[3]), int(parts[4]), parts[6], parts[8]

        # Extract gene_name from attributes
        if 'gene_name "' not in attr:
            continue
        gname = attr.split('gene_name "')[1].split('"')[0]

        if feature == "gene":
            biotype = attr.split('gene_type "')[1].split('"')[0] if 'gene_type "' in attr else ""
            gene_info[gname] = {
                "gene_id":   attr.split('gene_id "')[1].split('"')[0],
                "chr":       chrom,
                "start":     start,
                "end":       end,
                "strand":    strand,
                "biotype":   biotype,
                "length_kb": round((end - start) / 1000, 2),
            }
        elif feature == "transcript":
            tid = attr.split('transcript_id "')[1].split('"')[0]
            tbio = attr.split('transcript_type "')[1].split('"')[0] if 'transcript_type "' in attr else ""
            is_mane = 'tag "MANE_Select"' in attr
            trans_by_g.setdefault(gname, []).append((tid, tbio, is_mane))

print(f"Parsed {len(gene_info)} genes, {sum(len(v) for v in trans_by_g.values())} transcripts")

# ── 2. Build per-gene enrichment fields ──────────────────────────────────────
enrich = {}
for gname, gi in gene_info.items():
    ts = trans_by_g.get(gname, [])
    mane = [t for (t, b, m) in ts if m]
    pc_non_nmd = [t for (t, b, m) in ts
                  if b == "protein_coding"]
    enrich[gname] = {
        "gene_id_v47":     gi["gene_id"],
        "chr_v47":         gi["chr"],
        "gene_start_v47":  gi["start"],
        "gene_end_v47":    gi["end"],
        "strand_v47":      gi["strand"],
        "gene_biotype_v47":gi["biotype"],
        "gene_length_kb":  gi["length_kb"],
        "mane_select_tid": mane[0] if mane else "",
        "n_transcripts":   len(ts),
        "n_pc_transcripts":len(pc_non_nmd),
    }

# ── 3. Load candidate list and merge ─────────────────────────────────────────
print(f"\nLoading candidate list: {INPUT_CSV}")
cand = pd.read_csv(INPUT_CSV)
print(f"  {len(cand)} rows")

# Left-join by Gene name
cand["enriched"] = cand["Gene"].map(lambda g: enrich.get(g))
missing = cand[cand["enriched"].isna()]
if len(missing):
    print(f"⚠️  {len(missing)} rows missing v47 annotation (name mismatch?):")
    print(f"    genes: {missing['Gene'].unique()[:20].tolist()}")

# Expand enriched dict into columns
for col in ["gene_id_v47", "chr_v47", "gene_start_v47", "gene_end_v47",
            "strand_v47", "gene_biotype_v47", "gene_length_kb",
            "mane_select_tid", "n_transcripts", "n_pc_transcripts"]:
    cand[col] = cand["enriched"].apply(lambda d: d[col] if isinstance(d, dict) else "")
cand = cand.drop(columns=["enriched"])

# Sanity check: chr from VEP vs chr from v47 should match
mismatch = cand[(cand["chr_v47"] != "") & (cand["chr_v47"] != cand["chr"])]
if len(mismatch):
    print(f"⚠️  {len(mismatch)} chromosome mismatches between VEP and v47:")
    print(mismatch[["Gene","chr","chr_v47"]].head())

# ── 4. Compute variant → nearest intron/exon boundary (from MANE Select) ────
# For each row: find which intron the variant falls in (based on MANE transcript
# exons in v47). This gives us a v47-authoritative context.
print("\nParsing MANE Select exons for each candidate gene ...")
mane_exons_by_tid = {}
mane_tids_needed = set(cand["mane_select_tid"]) - {""}
if mane_tids_needed:
    with open(GTF) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attr = parts[8]
            tid = attr.split('transcript_id "')[1].split('"')[0]
            if tid not in mane_tids_needed:
                continue
            mane_exons_by_tid.setdefault(tid, []).append((int(parts[3]), int(parts[4])))
    # sort exon lists
    for tid in mane_exons_by_tid:
        mane_exons_by_tid[tid].sort()

def var_context(row):
    """Return string like 'intron 3/N' or 'exon K/N' based on MANE exons."""
    tid = row["mane_select_tid"]
    if not tid or tid not in mane_exons_by_tid:
        return ""
    exons = mane_exons_by_tid[tid]
    n = len(exons)
    strand = row["strand_v47"]
    pos = row["pos"]
    # Determine plus-strand exon index containing pos, or intron index
    for i, (s, e) in enumerate(exons, start=1):
        if s <= pos <= e:
            # inside an exon
            idx = i if strand == "+" else (n - i + 1)
            return f"exon {idx}/{n}"
    # Otherwise find which intron
    for i in range(len(exons) - 1):
        left_end = exons[i][1]
        right_start = exons[i + 1][0]
        if left_end < pos < right_start:
            idx = (i + 1) if strand == "+" else (n - i - 1)
            return f"intron {idx}/{n-1}"
    return "outside_transcript"

cand["v47_var_context"] = cand.apply(var_context, axis=1)

# ── 5. Compose the "run-ready" region string for each candidate ─────────────
cand["region_v47"] = cand.apply(
    lambda r: f"{r['chr_v47']}:{r['gene_start_v47']}-{r['gene_end_v47']}"
    if r["chr_v47"] else "", axis=1
)

# ── 6. Reorder columns for readability + sort ───────────────────────────────
priority_cols = [
    "Gene", "SpliceAI", "Impact", "Consequence",
    "Donor", "TissueID", "Tissue",
    "chr", "pos", "ref", "alt",
    "VAF",
    "v47_var_context",
    "region_v47",
    "mane_select_tid",
    "gene_length_kb", "n_transcripts", "n_pc_transcripts", "gene_biotype_v47",
    "strand_v47", "chr_v47", "gene_start_v47", "gene_end_v47",
    "has_RNAseq", "in_scan", "Already_case", "priority",
    "SpliceAI_band",
]
cand = cand[[c for c in priority_cols if c in cand.columns]]
cand = cand.sort_values(["priority", "SpliceAI"], ascending=[False, False])

cand.to_csv(OUT_CSV, index=False)
cand.to_csv(OUT_TSV, index=False, sep="\t")
print(f"\nSaved -> {OUT_CSV}")
print(f"Saved -> {OUT_TSV}")

# ── 7. Print top-20 NEW candidates as a compact preview ──────────────────────
print("\n" + "=" * 100)
print("TOP 20 NEW candidates (has RNA-seq + not yet scanned + high priority)")
print("=" * 100)
new_top = cand[(cand["has_RNAseq"]=="✓") & (cand["in_scan"]=="")].head(20)
show_cols = ["Gene", "SpliceAI", "Impact", "Donor", "TissueID",
             "v47_var_context", "region_v47", "mane_select_tid",
             "gene_length_kb"]
print(new_top[show_cols].to_string(index=False))
