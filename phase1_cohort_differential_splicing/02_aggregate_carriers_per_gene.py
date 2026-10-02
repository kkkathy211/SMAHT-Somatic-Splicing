#!/usr/bin/env python3
# ============================================================================
# Phase 1 - Cohort differential splicing (group test)
#
# Aggregate mutations to per-(tissue, gene) carrier statistics.
#
# Original location in the analysis project:
#     0520_analysis/code/02_carrier_aggregate.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
Phase 0.2 — Aggregate to per-(tissue, gene) carrier statistics.

Inputs:
    output/tables/00_annotated_mutations.csv       (Phase 0.1 output)
    ../v2_batch_corrected/NEW/data/tierA_sf_carriers.csv  (SF-gene panel)

Outputs:
    output/tables/00_carriers_per_tissue_gene.csv
    output/logs/02_carrier_aggregate.log

For each (tissue, gene) we compute:
    n_carriers          = number of donors with >=1 HC mutation in this gene
                          in this tissue
    n_total_donors      = number of donors with ANY data in this tissue
    n_non_carriers      = n_total_donors - n_carriers
    median_VAF_carriers = median of carrier donors' max VAF in this gene
    max_VAF_carriers    = max  of carrier donors' max VAF in this gene
    n_mutations_total   = total HC mutations in this gene in this tissue
    is_SF_gene          = boolean — does this gene match the SF-gene panel?
"""

import sys
from pathlib import Path
import pandas as pd

ROOT      = Path(__file__).resolve().parent.parent
IN_TABLE  = ROOT / "output" / "tables" / "00_annotated_mutations.csv"
SF_TABLE  = ROOT.parent / "v2_batch_corrected" / "NEW" / "data" / "tierA_sf_carriers.csv"
OUT_TABLE = ROOT / "output" / "tables" / "00_carriers_per_tissue_gene.csv"
LOG       = ROOT / "output" / "logs" / "02_carrier_aggregate.log"

# Curated SF-gene panel — same 15 genes used in v1
SF_GENES = {
    "SRSF2", "SF3B1", "U2AF1", "U2AF2", "RBM10", "RBM39",
    "PRPF8", "SF1", "ZRSR2", "DDX3X", "SRSF7", "SRSF3",
    "SF3A1", "SF3A2", "FUBP1",
}


def main():
    log_fh = open(LOG, "w")
    def log(m):
        print(m); log_fh.write(m + "\n"); log_fh.flush()

    log("== Phase 0.2 — Per-(tissue, gene) carrier aggregation ==")
    log(f"Input: {IN_TABLE}")
    log(f"SF gene list: {len(SF_GENES)} curated genes")
    log("")

    # Load annotated mutations
    mut = pd.read_csv(IN_TABLE)
    log(f"Loaded {len(mut):,} mutations  ({mut['donor'].nunique()} donors, "
        f"{mut['tissue'].nunique()} tissues)")

    # Drop INTERGENIC mutations — they can't be aggregated by gene
    mut_in_genes = mut[mut["gene_name"] != "INTERGENIC"].copy()
    log(f"After dropping INTERGENIC: {len(mut_in_genes):,} mutations")
    log(f"   In {mut_in_genes['gene_name'].nunique():,} unique genes")
    log("")

    # -------------------------------------------------------------------------
    # Per (tissue, gene, donor): donor's max VAF in this gene
    #   — collapse multiple mutations from same donor in same gene to max VAF
    # -------------------------------------------------------------------------
    log("[1/3] Collapsing multiple mutations per (donor, tissue, gene) ...")
    per_donor_gene = (
        mut_in_genes
        .groupby(["donor", "tissue", "gene_name", "gene_type"], as_index=False)
        .agg(donor_max_vaf=("vaf", "max"),
             donor_n_muts=("vaf", "count"))
    )
    log(f"   Unique (donor, tissue, gene) rows: {len(per_donor_gene):,}")

    # -------------------------------------------------------------------------
    # Per (tissue, gene): aggregate across donors
    # -------------------------------------------------------------------------
    log("\n[2/3] Aggregating to (tissue, gene) ...")
    per_tg = (
        per_donor_gene
        .groupby(["tissue", "gene_name", "gene_type"], as_index=False)
        .agg(n_carriers=("donor", "nunique"),
             median_VAF_carriers=("donor_max_vaf", "median"),
             max_VAF_carriers=("donor_max_vaf", "max"),
             min_VAF_carriers=("donor_max_vaf", "min"),
             n_mutations_total=("donor_n_muts", "sum"))
    )
    log(f"   Unique (tissue, gene) rows: {len(per_tg):,}")

    # -------------------------------------------------------------------------
    # Per-tissue donor count
    # -------------------------------------------------------------------------
    log("\n[3/3] Computing per-tissue donor cohort size ...")
    donors_per_tissue = (
        mut[["donor", "tissue"]]
        .drop_duplicates()
        .groupby("tissue")["donor"]
        .nunique()
        .rename("n_total_donors")
        .reset_index()
    )
    log("   Donors per tissue:")
    for _, r in donors_per_tissue.iterrows():
        log(f"      {r['tissue']:32s}  {r['n_total_donors']:3d}")
    log("")

    # Merge total donor counts
    per_tg = per_tg.merge(donors_per_tissue, on="tissue", how="left")
    per_tg["n_non_carriers"] = per_tg["n_total_donors"] - per_tg["n_carriers"]
    per_tg["is_SF_gene"] = per_tg["gene_name"].isin(SF_GENES)

    # Reorder columns
    final = per_tg[[
        "tissue", "gene_name", "gene_type",
        "n_carriers", "n_total_donors", "n_non_carriers",
        "median_VAF_carriers", "max_VAF_carriers", "min_VAF_carriers",
        "n_mutations_total", "is_SF_gene",
    ]]
    final = final.sort_values(["tissue", "n_carriers", "n_mutations_total"],
                                ascending=[True, False, False])

    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUT_TABLE, index=False)
    log(f"   Wrote {OUT_TABLE} ({len(final):,} rows)")

    # -------------------------------------------------------------------------
    # Quick summary stats
    # -------------------------------------------------------------------------
    log("\n=== Quick summary ===")
    log(f"Total unique (tissue, gene) pairs: {len(final):,}")

    log("\nPairs with >= 4 carriers (Strategy A eligibility, by tissue):")
    a_elig = final[final["n_carriers"] >= 4].copy()
    a_elig["both_groups_ok"] = a_elig["n_non_carriers"] >= 4
    log(f"   Total (tissue, gene) with n_carriers>=4:       {len(a_elig):,}")
    log(f"   With n_non_carriers also >=4 (both groups):    "
        f"{a_elig['both_groups_ok'].sum():,}")
    log("")
    log("   By tissue:")
    for tissue, sub in a_elig[a_elig["both_groups_ok"]].groupby("tissue"):
        log(f"      {tissue:32s}  {len(sub):4d} eligible (tissue,gene) pairs")

    log("\nSF gene presence per tissue (any carrier):")
    sf_in_tissue = (final[(final["is_SF_gene"]) & (final["n_carriers"] >= 1)]
                     .groupby("tissue")["gene_name"]
                     .apply(list))
    for tissue, genes in sf_in_tissue.items():
        log(f"   {tissue:32s}  SF carriers: {', '.join(sorted(set(genes)))}")

    log("\nGenes mutated in >= 2 tissues with >=4 carriers each "
        "(Strategy B candidates preview):")
    b_pre = a_elig[a_elig["both_groups_ok"]]
    cross = (b_pre.groupby("gene_name")["tissue"]
             .nunique().rename("n_tissues_eligible"))
    cross = cross[cross >= 2].sort_values(ascending=False)
    log(f"   Total such genes: {len(cross):,}")
    log(f"   Top 15:")
    for g, n in cross.head(15).items():
        log(f"      {g:20s}  {n} tissues")

    log("\nDONE.")
    log_fh.close()


if __name__ == "__main__":
    main()
