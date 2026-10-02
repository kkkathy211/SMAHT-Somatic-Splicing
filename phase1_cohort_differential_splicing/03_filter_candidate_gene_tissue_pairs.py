#!/usr/bin/env python3
# ============================================================================
# Phase 1 - Cohort differential splicing (group test)
#
# Apply carrier-count filters to produce Strategy A / B / C candidate lists.
#
# Original location in the analysis project:
#     0520_analysis/code/03_candidate_filter.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
Phase 0.3 — Filter candidates for Strategies A, B, C.

Inputs:
    output/tables/00_carriers_per_tissue_gene.csv      (Phase 0.2)
    output/tables/00_annotated_mutations.csv           (Phase 0.1)

Outputs:
    output/tables/00_candidates_A.csv   (tissue,gene)  pairs for within-tissue
    output/tables/00_candidates_B.csv   (gene)         cross-tissue genes
    output/tables/00_candidates_C.csv   (variant)      priority SF / rare variants
    output/logs/03_candidate_filter.log

Filtering rules:
    A: n_carriers >= 4  AND  n_non_carriers >= 4
        + gene_type == "protein_coding" (drop lincRNA / pseudogene noise)
    B: gene satisfies A in >= 2 tissues
    C: variant in SF-gene panel (any carrier, no n threshold)
       OR very high SR_VAF (>= 0.20)  AND  protein_coding gene
"""

import sys
from pathlib import Path
import pandas as pd

ROOT      = Path(__file__).resolve().parent.parent
PER_TG    = ROOT / "output" / "tables" / "00_carriers_per_tissue_gene.csv"
MUTS      = ROOT / "output" / "tables" / "00_annotated_mutations.csv"
OUT_A     = ROOT / "output" / "tables" / "00_candidates_A.csv"
OUT_B     = ROOT / "output" / "tables" / "00_candidates_B.csv"
OUT_C     = ROOT / "output" / "tables" / "00_candidates_C.csv"
LOG       = ROOT / "output" / "logs" / "03_candidate_filter.log"

SF_GENES = {
    "SRSF2", "SF3B1", "U2AF1", "U2AF2", "RBM10", "RBM39",
    "PRPF8", "SF1", "ZRSR2", "DDX3X", "SRSF7", "SRSF3",
    "SF3A1", "SF3A2", "FUBP1",
}

MIN_N_CARRIERS    = 4
MIN_N_NONCARRIERS = 4
MIN_TISSUES_FOR_B = 2
HIGH_VAF_THRESH   = 0.20


def main():
    log_fh = open(LOG, "w")
    def log(m):
        print(m); log_fh.write(m + "\n"); log_fh.flush()

    log("== Phase 0.3 — Candidate filtering (A/B/C) ==\n")
    per_tg = pd.read_csv(PER_TG)
    muts = pd.read_csv(MUTS)
    log(f"Loaded {len(per_tg):,} (tissue, gene) rows")
    log(f"Loaded {len(muts):,} HC mutation rows")

    # =====================================================================
    # CANDIDATE LIST A — within-tissue, ≥4 carriers + ≥4 non-carriers
    # =====================================================================
    log("\n--------------------------------------")
    log("[A] Strategy A: within-tissue carrier-vs-non-carrier")
    log("--------------------------------------")
    log(f"Filters:")
    log(f"   n_carriers     >= {MIN_N_CARRIERS}")
    log(f"   n_non_carriers >= {MIN_N_NONCARRIERS}")
    log(f"   gene_type      == protein_coding")

    A = per_tg[
        (per_tg["n_carriers"] >= MIN_N_CARRIERS) &
        (per_tg["n_non_carriers"] >= MIN_N_NONCARRIERS) &
        (per_tg["gene_type"] == "protein_coding")
    ].copy().sort_values(["n_carriers", "median_VAF_carriers"],
                         ascending=[False, False])

    log(f"\nCandidate A list size: {len(A):,} (tissue, gene) pairs")
    log(f"Top 20 by n_carriers:")
    cols = ["tissue", "gene_name", "n_carriers", "n_non_carriers",
            "median_VAF_carriers", "max_VAF_carriers", "is_SF_gene"]
    for _, r in A.head(20).iterrows():
        log(f"   {r['tissue']:28s} {r['gene_name']:15s} "
            f"carriers={r['n_carriers']:2d}  non={r['n_non_carriers']:2d}  "
            f"median_VAF={r['median_VAF_carriers']:.3f}  "
            f"SF={'Y' if r['is_SF_gene'] else 'N'}")

    log(f"\nBreakdown by tissue:")
    for tissue, sub in A.groupby("tissue"):
        log(f"   {tissue:28s}  {len(sub):4d} candidates")

    A.to_csv(OUT_A, index=False)
    log(f"\nWrote {OUT_A}")

    # =====================================================================
    # CANDIDATE LIST B — same gene, >=2 tissues meeting A's filter
    # =====================================================================
    log("\n--------------------------------------")
    log("[B] Strategy B: cross-tissue (same gene across tissues)")
    log("--------------------------------------")
    log(f"Filter: gene appears in candidates_A in >= {MIN_TISSUES_FOR_B} tissues")

    gene_tissue_count = A.groupby("gene_name")["tissue"].nunique()
    B_genes = gene_tissue_count[gene_tissue_count >= MIN_TISSUES_FOR_B].index

    B = A[A["gene_name"].isin(B_genes)].copy()
    # Add the per-gene tissue count
    B["n_eligible_tissues"] = B["gene_name"].map(gene_tissue_count)
    B = B.sort_values(["n_eligible_tissues", "gene_name", "n_carriers"],
                      ascending=[False, True, False])

    log(f"\nCandidate B list:")
    log(f"   Unique genes with >=2 tissue eligibility: "
        f"{B['gene_name'].nunique():,}")
    log(f"   Total (gene, tissue) pair rows:           {len(B):,}")

    log(f"\nGenes that show up in the most tissues:")
    for g, n in B[["gene_name", "n_eligible_tissues"]].drop_duplicates() \
                  .sort_values("n_eligible_tissues", ascending=False) \
                  .head(20).values:
        carriers_per_t = (
            B[B["gene_name"] == g][["tissue", "n_carriers"]]
            .to_dict("records")
        )
        carrier_str = ', '.join(
            f"{r['tissue']}:{r['n_carriers']}" for r in carriers_per_t[:5]
        )
        log(f"   {g:18s}  {n} tissues  carriers: {carrier_str}")

    B.to_csv(OUT_B, index=False)
    log(f"\nWrote {OUT_B}")

    # =====================================================================
    # CANDIDATE LIST C — priority variants (SF panel OR high VAF)
    # =====================================================================
    log("\n--------------------------------------")
    log("[C] Strategy C: priority variants (SF genes + high-VAF)")
    log("--------------------------------------")
    log(f"Filter (OR logic — keep variant if ANY of these):")
    log(f"   (1) gene in SF panel ({len(SF_GENES)} genes)")
    log(f"   (2) VAF >= {HIGH_VAF_THRESH:.2f} AND protein-coding gene")

    # Variant-level filter (work on individual mutations, not aggregate)
    is_sf       = muts["gene_name"].isin(SF_GENES)
    is_high_vaf = (muts["vaf"] >= HIGH_VAF_THRESH) & \
                  (muts["gene_type"] == "protein_coding")
    C = muts[is_sf | is_high_vaf].copy()

    # Annotate priority reason
    C["is_SF_gene"]     = C["gene_name"].isin(SF_GENES)
    C["is_high_vaf"]    = (C["vaf"] >= HIGH_VAF_THRESH) & \
                          (C["gene_type"] == "protein_coding")
    C["priority_class"] = "SF_only"
    C.loc[C["is_high_vaf"] & ~C["is_SF_gene"], "priority_class"] = "high_VAF_only"
    C.loc[C["is_SF_gene"]  &  C["is_high_vaf"], "priority_class"] = "SF_AND_high_VAF"
    C = C.sort_values(["is_SF_gene", "vaf"], ascending=[False, False])

    log(f"\nCandidate C list size: {len(C):,} variants")
    log(f"   SF-gene variants (any VAF):  {C['is_SF_gene'].sum():,}")
    log(f"   High-VAF protein-coding:     {C['is_high_vaf'].sum():,}")
    log(f"   Class breakdown:")
    for cls, n in C["priority_class"].value_counts().items():
        log(f"      {cls:25s}  {n:,}")

    log(f"\nSF-gene variants in detail:")
    for _, r in C[C["is_SF_gene"]].iterrows():
        log(f"   {r['donor']:12s}  {r['tissue']:28s}  "
            f"{r['gene_name']:8s}  {r['chrom']}:{r['pos']:>11d}  "
            f"VAF={r['vaf']:.3f}")

    log(f"\nTop 20 high-VAF non-SF variants:")
    high_only = C[C["priority_class"] == "high_VAF_only"].head(20)
    for _, r in high_only.iterrows():
        log(f"   {r['donor']:12s}  {r['tissue']:28s}  "
            f"{r['gene_name']:18s}  VAF={r['vaf']:.3f}")

    C.to_csv(OUT_C, index=False)
    log(f"\nWrote {OUT_C}")

    log("\n=== ALL DONE ===")
    log_fh.close()


if __name__ == "__main__":
    main()
