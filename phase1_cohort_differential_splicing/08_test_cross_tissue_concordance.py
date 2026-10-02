#!/usr/bin/env python3
# ============================================================================
# Phase 1 - Cohort differential splicing (group test)
#
# Strategy B: does the same gene move the same direction in >=2 tissues?
#
# Original location in the analysis project:
#     0520_analysis/code/10_strategy_B_concordance.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
Strategy B — Cross-tissue concordance.

Two levels of concordance:

  Level 1 — GENE-LEVEL (primary, always works):
     For each cross-tissue gene G, look at the *top* junction per tissue
     (from A_summary).  Ask:
        - Do tissues with q<0.10 have consistent β direction?
        - Do all tested tissues (even non-sig) show same direction?
        - How many tissues show signal?

  Level 2 — JUNCTION-LEVEL (when same junction is sig in ≥2 tissues):
     Match stripped junction IDs (chr:start:end:strand) across tissues.
     For shared junctions, classify as concordant/discordant.

Reads:
    output/tables/A_summary.csv
    output/tables/A_top_hits.csv
    output/tables/00_candidates_B.csv

Writes:
    output/tables/B_summary.csv           (per gene, one row)
    output/tables/B_per_tissue.csv        (long, gene × tissue × top stats)
    output/tables/B_junction_concordant.csv  (junction-level matches if any)
    output/logs/10_strategy_B.log
"""

import re
import sys
from pathlib import Path
import pandas as pd
import numpy as np

ROOT  = Path(__file__).resolve().parent.parent
A_SUM = ROOT / "output" / "tables" / "A_summary.csv"
A_HIT = ROOT / "output" / "tables" / "A_top_hits.csv"

OUT_GENE  = ROOT / "output" / "tables" / "B_summary.csv"
OUT_TIS   = ROOT / "output" / "tables" / "B_per_tissue.csv"
OUT_JCONC = ROOT / "output" / "tables" / "B_junction_concordant.csv"
LOG       = ROOT / "output" / "logs" / "10_strategy_B.log"

JUNC_STRIP_RE = re.compile(r"^(.+):clu_\d+_([+-])$")


def strip_cluster(jid):
    if pd.isna(jid): return None
    m = JUNC_STRIP_RE.match(jid)
    if not m: return jid
    return f"{m.group(1)}:{m.group(2)}"


def sign_str(x):
    if pd.isna(x): return "?"
    if x > 0:      return "+"
    if x < 0:      return "−"
    return "0"


def classify_gene(per_tissue):
    """per_tissue is a DataFrame of tissues for one gene, with A3_top_q,
    A3_top_beta, A1_top_q, A1_top_beta etc."""
    # signs of A.3 top β where defined
    signs = [np.sign(b) for b in per_tissue["A3_top_beta"]
              if not pd.isna(b) and b != 0]
    sig_mask = per_tissue["A3_top_q"] < 0.10
    n_sig = int(sig_mask.sum())
    sig_signs = [np.sign(b) for b, s in zip(per_tissue["A3_top_beta"], sig_mask)
                  if s and not pd.isna(b)]

    if n_sig == 0:
        return "NO_SIGNAL", 0, 0
    if n_sig == 1:
        return "TISSUE_SPECIFIC", 0, 0
    # ≥2 tissues sig
    n_concord_sig = sum(s == sig_signs[0] for s in sig_signs)
    if n_concord_sig == n_sig:
        return "CONCORDANT (same sign in all sig tissues)", n_concord_sig, 0
    if n_concord_sig == 1:
        return "DISCORDANT", 0, n_sig - 1
    return "MIXED", n_concord_sig, n_sig - n_concord_sig


def main():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    log_fh = open(LOG, "w")
    def L(m):
        print(m); log_fh.write(m + "\n"); log_fh.flush()

    L("== Strategy B — Cross-tissue concordance ==\n")

    A_summary = pd.read_csv(A_SUM)
    A_hits    = pd.read_csv(A_HIT)
    L(f"A_summary:   {len(A_summary)} rows")
    L(f"A_top_hits:  {len(A_hits)} rows\n")

    # Genes tested in ≥2 tissues
    gene_tissue_count = (A_summary.groupby("gene")["tissue"].nunique()
                                    .rename("n_tissues_tested"))
    multi_tissue_genes = gene_tissue_count[gene_tissue_count >= 2].index.tolist()
    L(f"Genes tested in ≥2 tissues: {len(multi_tissue_genes)}")
    L(f"  {multi_tissue_genes}\n")

    # ------------------------------------------------------------------
    # B.1 — per-tissue table (one row per gene × tissue)
    # ------------------------------------------------------------------
    L("[1/3] Building per-tissue gene table ...")
    keep_cols = [
        "tissue", "gene", "n_carrier_cores", "n_non_carr_cores",
        "n_carrier_donors", "n_total_donors", "n_junctions",
        "A1_top_jid", "A1_top_q", "A1_top_beta",
        "A2_top_jid", "A2_top_q", "A2_top_beta",
        "A3_top_jid", "A3_top_q", "A3_top_beta",
        "A1_n_q10", "A2_n_q10", "A3_n_q10",
        "A1_n_q05", "A2_n_q05", "A3_n_q05",
    ]
    per_tis = A_summary[A_summary["gene"].isin(multi_tissue_genes)][
        [c for c in keep_cols if c in A_summary.columns]
    ].copy()
    per_tis["A3_sig_q10"] = per_tis["A3_top_q"] < 0.10
    per_tis["A3_top_sign"] = per_tis["A3_top_beta"].apply(sign_str)
    per_tis = per_tis.sort_values(["gene", "A3_top_q"])
    per_tis.to_csv(OUT_TIS, index=False)
    L(f"   Wrote {OUT_TIS}  ({len(per_tis)} rows)")

    # ------------------------------------------------------------------
    # B.2 — per-gene concordance summary
    # ------------------------------------------------------------------
    L("\n[2/3] Per-gene concordance classification ...")
    summary_rows = []
    for gene in multi_tissue_genes:
        gsub = A_summary[A_summary["gene"] == gene]
        label, n_conc_sig, n_disc_sig = classify_gene(gsub)
        n_tested = gsub["tissue"].nunique()
        n_sig_q10 = int((gsub["A3_top_q"] < 0.10).sum())
        n_sig_q05 = int((gsub["A3_top_q"] < 0.05).sum())

        tissues_tested = sorted(gsub["tissue"].unique().tolist())
        sig_tissues = sorted(gsub.loc[gsub["A3_top_q"] < 0.10, "tissue"].tolist())
        per_tissue_str = "; ".join(
            f"{r['tissue']}: A3_q={r['A3_top_q']:.2g} β={r['A3_top_beta']:+.2f}"
            for _, r in gsub.sort_values("A3_top_q").iterrows()
            if not pd.isna(r["A3_top_q"])
        )
        signs_str = "/".join(
            sign_str(r["A3_top_beta"]) for _, r in gsub.iterrows()
        )

        summary_rows.append({
            "gene": gene,
            "n_tissues_tested": n_tested,
            "n_tissues_sig_q10": n_sig_q10,
            "n_tissues_sig_q05": n_sig_q05,
            "tissues_tested": ",".join(tissues_tested),
            "sig_tissues": ",".join(sig_tissues),
            "A3_top_beta_signs": signs_str,
            "n_concordant_sig": n_conc_sig,
            "n_discordant_sig": n_disc_sig,
            "concordance_label": label,
            "per_tissue_top_stats": per_tissue_str,
        })

    summary_df = pd.DataFrame(summary_rows).sort_values(
        ["n_concordant_sig", "n_tissues_sig_q10"], ascending=[False, False])
    summary_df.to_csv(OUT_GENE, index=False)
    L(f"   Wrote {OUT_GENE}  ({len(summary_df)} rows)\n")

    # ------------------------------------------------------------------
    # B.3 — Junction-level concordance (when same junction shared)
    # ------------------------------------------------------------------
    L("[3/3] Junction-level concordance (stripped junction IDs) ...")
    A_hits["jid_strip"] = A_hits["junction_id"].apply(strip_cluster)
    multi_hits = A_hits[A_hits["gene"].isin(multi_tissue_genes)]
    L(f"   Top-hits restricted to multi-tissue genes: {len(multi_hits)} rows")

    junction_rows = []
    for (gene, jid_s), sub in multi_hits.groupby(["gene", "jid_strip"]):
        if sub["tissue"].nunique() < 2:
            continue
        # This jid is sig (q<0.10) in ≥2 tissues
        tissues = sub["tissue"].tolist()
        betas = sub["A3_beta"].tolist()
        qs    = sub["A3_q"].tolist()
        signs = [np.sign(b) for b in betas if not pd.isna(b) and b != 0]
        if len(signs) < 2:
            label = "AMBIGUOUS"
        elif all(s == signs[0] for s in signs):
            label = "CONCORDANT"
        else:
            label = "DISCORDANT"
        junction_rows.append({
            "gene": gene,
            "junction_stripped": jid_s,
            "n_tissues": sub["tissue"].nunique(),
            "tissues": ",".join(tissues),
            "A3_betas": ",".join(f"{b:+.3f}" if not pd.isna(b) else "NA"
                                  for b in betas),
            "A3_qs": ",".join(f"{q:.2g}" if not pd.isna(q) else "NA"
                               for q in qs),
            "concordance": label,
        })

    if junction_rows:
        junc_df = pd.DataFrame(junction_rows).sort_values("n_tissues",
                                                            ascending=False)
    else:
        junc_df = pd.DataFrame(columns=["gene", "junction_stripped",
                                          "n_tissues", "tissues",
                                          "A3_betas", "A3_qs", "concordance"])
    junc_df.to_csv(OUT_JCONC, index=False)
    L(f"   Wrote {OUT_JCONC}  ({len(junc_df)} junction-level rows)")

    if junc_df.empty:
        L("\n   ⚠ NO single junction reached q<0.10 in ≥2 tissues.")
        L("     Conclusion: cross-tissue replication at the EXACT junction "
          "level is absent in our data. Gene-level concordance is the more "
          "informative analysis (see B_summary.csv).")
    else:
        for _, r in junc_df.iterrows():
            L(f"     {r['gene']:10s}  {r['junction_stripped']:40s}  "
              f"{r['tissues']}  βs={r['A3_betas']}  → {r['concordance']}")

    # ------------------------------------------------------------------
    L("\n=== Strategy B summary ===\n")
    L("Per-gene concordance (gene-level):\n")
    L(f"{'gene':<14s} {'tested':>6s} {'sig_q10':>7s} {'sig_q05':>7s} "
      f"{'signs':>10s}  conclusion")
    L("-" * 80)
    for _, r in summary_df.iterrows():
        L(f"{r['gene']:<14s} {r['n_tissues_tested']:>6d} "
          f"{r['n_tissues_sig_q10']:>7d} {r['n_tissues_sig_q05']:>7d} "
          f"{r['A3_top_beta_signs']:>10s}  {r['concordance_label']}")

    L("\nDONE.")
    log_fh.close()


if __name__ == "__main__":
    main()
