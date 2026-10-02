#!/usr/bin/env python3
# ============================================================================
# Phase 1 - Cohort differential splicing (group test)
#
# Join somatic mutations to SpliceAI scores; build the LOOSE / Q95 carrier tiers.
#
# Original location in the analysis project:
#     0520_analysis/code/02b_annotate_splice_pred.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
Phase 0.2b — Annotate every somatic mutation with SpliceAI prediction, then
aggregate to per-(donor, tissue, gene) splice-affecting-carrier flags.

Inputs:
    output/tables/00_annotated_mutations.csv     (Phase 0.1)
    spliceAI_VAF_quantile_fraction.txt           (Yilin V1 list)

Outputs:
    output/tables/02b_mutations_with_spliceAI.csv
        per-mutation row, plus columns:
            spliceAI_score    (0..1, NA if not in Yilin's table)
            spliceAI_type     (Q95 / Q99 / other / NA)
    output/tables/02b_carrier_per_donor_gene_v2.csv
        per (donor, tissue, gene) row, with:
            donor_max_vaf
            donor_n_muts
            donor_max_spliceAI            (max across mutations in gene)
            donor_n_high_splice_variants  (count with score >= 0.20)
            donor_n_Q95plus_variants      (count with type in {Q95, Q99})
            is_carrier_ANY                (1 = donor has any mut in gene)
            is_carrier_LOOSE              (1 = donor_max_spliceAI >= 0.20)
            is_carrier_Q95                (1 = donor_n_Q95plus_variants >= 1)
    output/logs/02b_annotate_splice_pred.log

The SpliceAI score is per-variant (function of sequence), so the score lookup
is keyed on (chrom, pos, ref, alt). Multiple Yilin rows for the same variant
(different donor/tissue) carry the same score — we collapse them with .first().

Carrier hierarchy (NESTED):
    Q95   ⊂   LOOSE   ⊂   ANY

For lm A/B v2: non-carriers = donors with NO somatic mutation in gene.
Donors that are ANY but not LOOSE are SKIPPED in LOOSE / Q95 runs (ambiguous).
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
MUT_IN  = ROOT / "output" / "tables" / "00_annotated_mutations.csv"
AI_IN   = ROOT / "spliceAI_VAF_quantile_fraction.txt"
MUT_OUT = ROOT / "output" / "tables" / "02b_mutations_with_spliceAI.csv"
PDG_OUT = ROOT / "output" / "tables" / "02b_carrier_per_donor_gene_v2.csv"
LOG     = ROOT / "output" / "logs" / "02b_annotate_splice_pred.log"

# Score thresholds (configurable)
LOOSE_THRESH = 0.20
Q95_TYPES = {"Q95", "Q99"}


def main():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    log_fh = open(LOG, "w")
    def L(m):
        print(m); log_fh.write(m + "\n"); log_fh.flush()

    L("== Phase 0.2b — annotate mutations with SpliceAI ==\n")

    # -------------------------------------------------------------------
    # Load mutation table (canonical)
    # -------------------------------------------------------------------
    mut = pd.read_csv(MUT_IN)
    L(f"Loaded mutations:   {len(mut):>8,} rows  "
      f"({mut['donor'].nunique()} donors, "
      f"{mut['tissue'].nunique()} tissues, "
      f"{(mut['gene_name'] != '-1').sum():,} in-gene)")

    # -------------------------------------------------------------------
    # Load SpliceAI table (variant-level prediction)
    # -------------------------------------------------------------------
    ai = pd.read_csv(AI_IN, sep="\t", low_memory=False)
    L(f"Loaded spliceAI:    {len(ai):>8,} rows")
    L(f"   columns: {list(ai.columns)}")

    # Per-variant lookup keyed by chrom/pos/ref/alt (multiple Yilin rows
    # for same variant just repeat the same prediction)
    ai_key_cols = ["chr", "pos", "ref", "alt"]
    ai_uniq = (
        ai[ai_key_cols + ["score", "spliceAI_type"]]
        .drop_duplicates(subset=ai_key_cols, keep="first")
        .rename(columns={"chr": "chrom",
                          "score": "spliceAI_score",
                          "spliceAI_type": "spliceAI_class"})
    )
    L(f"   unique variants:  {len(ai_uniq):,}")

    # Score distribution
    sc = ai_uniq["spliceAI_score"]
    L(f"   score distribution:")
    L(f"     max:        {sc.max():.3f}")
    L(f"     >= 0.5:     {(sc >= 0.5).sum():,}")
    L(f"     >= 0.2:     {(sc >= 0.20).sum():,}")
    L(f"     >= 0.1:     {(sc >= 0.1).sum():,}")
    L(f"     == 0:       {(sc == 0).sum():,}")
    L(f"   class breakdown: {dict(ai_uniq['spliceAI_class'].value_counts())}")
    L("")

    # -------------------------------------------------------------------
    # Join: mutation × spliceAI by (chrom, pos, ref, alt)
    # -------------------------------------------------------------------
    L("[1/3] Joining mutation × spliceAI ...")
    mut2 = mut.merge(ai_uniq,
                      on=["chrom", "pos", "ref", "alt"],
                      how="left")
    n_matched = mut2["spliceAI_score"].notna().sum()
    L(f"   {n_matched:,} / {len(mut2):,} ({100*n_matched/len(mut2):.1f}%) "
      f"mutations matched a spliceAI prediction")

    # Fill NA score with 0 (no prediction → assume no splice effect)
    mut2["spliceAI_score"] = mut2["spliceAI_score"].fillna(0.0)
    mut2["spliceAI_class"] = mut2["spliceAI_class"].fillna("none")

    mut2.to_csv(MUT_OUT, index=False)
    L(f"   Wrote {MUT_OUT}")

    # -------------------------------------------------------------------
    # Drop INTERGENIC (gene_name == "-1") for aggregation
    # -------------------------------------------------------------------
    in_gene = mut2[mut2["gene_name"] != "-1"].copy()
    L(f"   In-gene mutations: {len(in_gene):,}")

    # -------------------------------------------------------------------
    # Aggregate per (donor, tissue, gene)
    # -------------------------------------------------------------------
    L("\n[2/3] Aggregating to (donor, tissue, gene) ...")
    pdg = (
        in_gene
        .groupby(["donor", "tissue", "gene_name", "gene_type"],
                  as_index=False)
        .agg(donor_max_vaf=("vaf", "max"),
             donor_n_muts=("vaf", "count"),
             donor_max_spliceAI=("spliceAI_score", "max"),
             donor_n_high_splice_variants=("spliceAI_score",
                                            lambda s: int((s >= LOOSE_THRESH).sum())),
             donor_n_Q95plus_variants=("spliceAI_class",
                                        lambda s: int(s.isin(Q95_TYPES).sum())),
             variants_summary=("chrom",
                                lambda _: ""))  # placeholder, filled below
    )

    # Re-compute a compact variants_summary string (for breakdown later)
    L("   Building per-(donor, tissue, gene) variant summary strings ...")
    summary_map = {}
    for (donor, tis, gene, _gt), sub in in_gene.groupby(
            ["donor", "tissue", "gene_name", "gene_type"]):
        parts = []
        for _, r in sub.sort_values("spliceAI_score", ascending=False).iterrows():
            parts.append(f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} "
                          f"(VAF={r['vaf']:.2f}, AI={r['spliceAI_score']:.2f}"
                          f"{', '+r['spliceAI_class'] if r['spliceAI_class'] in Q95_TYPES else ''})")
        summary_map[(donor, tis, gene)] = " | ".join(parts)
    pdg["variants_summary"] = pdg.apply(
        lambda r: summary_map.get((r["donor"], r["tissue"], r["gene_name"]),
                                    ""),
        axis=1)

    # Carrier flags
    pdg["is_carrier_ANY"]   = 1
    pdg["is_carrier_LOOSE"] = (pdg["donor_max_spliceAI"]
                                >= LOOSE_THRESH).astype(int)
    pdg["is_carrier_Q95"]   = (pdg["donor_n_Q95plus_variants"]
                                >= 1).astype(int)

    L(f"   {len(pdg):,} (donor, tissue, gene) rows")
    L(f"   carriers ANY:    {(pdg['is_carrier_ANY']==1).sum():,}")
    L(f"   carriers LOOSE:  {(pdg['is_carrier_LOOSE']==1).sum():,}  "
      f"(SpliceAI ≥ {LOOSE_THRESH})")
    L(f"   carriers Q95:    {(pdg['is_carrier_Q95']==1).sum():,}  "
      f"(spliceAI_class in {sorted(Q95_TYPES)})")

    pdg.to_csv(PDG_OUT, index=False)
    L(f"   Wrote {PDG_OUT}")

    # -------------------------------------------------------------------
    # Quick per-tissue sanity print
    # -------------------------------------------------------------------
    L("\n[3/3] Per-tissue carrier counts (LOOSE / Q95):")
    by_tis = (pdg.groupby("tissue")
                  .agg(n_ANY=("is_carrier_ANY", "sum"),
                       n_LOOSE=("is_carrier_LOOSE", "sum"),
                       n_Q95=("is_carrier_Q95", "sum")))
    for tis, r in by_tis.iterrows():
        L(f"   {tis:32s}  ANY={int(r['n_ANY']):>5d}  "
          f"LOOSE={int(r['n_LOOSE']):>4d}  Q95={int(r['n_Q95']):>3d}")

    L("\nDONE.")
    log_fh.close()


if __name__ == "__main__":
    main()
