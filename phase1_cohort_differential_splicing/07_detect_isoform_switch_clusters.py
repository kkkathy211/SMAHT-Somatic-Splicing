#!/usr/bin/env python3
# ============================================================================
# Phase 1 - Cohort differential splicing (group test)
#
# Roll junction hits up to clusters; detect mirror pairs (isoform-switch signature).
#
# Original location in the analysis project:
#     0520_analysis/code/06_aggregate_to_clusters.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
Strategy A Step 6 — Aggregate per-junction hits to cluster level
                    + detect mirror pairs.

Reads:
    output/tables/A_top_hits.csv         per-junction long format
    output/tables/A_cluster_to_gene.csv  cluster ↔ gene mapping

Produces:
    output/tables/A_cluster_hits.csv     per-cluster summary with:
        - n_junctions_sig (q<0.10 in any test)
        - has_mirror_pair (cluster has ≥2 junctions with opposite-sign β
                            of similar magnitude)
        - confidence_stars   (0–3)
    output/logs/06_aggregate_to_clusters.log

Mirror pair criteria:
    Within a cluster, find any pair (j1, j2) such that:
       - both q < 0.10  (in test A.3, our primary)
       - β(j1) and β(j2) have opposite signs
       - |β(j1)| within 50% of |β(j2)|  (i.e. magnitude ratio 0.5–2.0)
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
HITS = ROOT / "output" / "tables" / "A_top_hits.csv"
C2G  = ROOT / "output" / "tables" / "A_cluster_to_gene.csv"
SUMR = ROOT / "output" / "tables" / "A_summary.csv"
OUT  = ROOT / "output" / "tables" / "A_cluster_hits.csv"
LOG  = ROOT / "output" / "logs" / "06_aggregate_to_clusters.log"

Q_THRESH       = 0.10
MIRROR_TOL     = 0.50    # |β| ratio must be within 0.5 to 2.0


def log_open(path):
    fh = open(path, "w")
    def L(m):
        print(m); fh.write(m + "\n"); fh.flush()
    return L, fh


def detect_mirror_pair(sub_df: pd.DataFrame) -> bool:
    """Sub_df = one cluster's rows. Returns True iff mirror pair exists."""
    rows = sub_df[(sub_df["A3_q"] < Q_THRESH) & sub_df["A3_beta"].notna()]
    if len(rows) < 2:
        return False
    betas = rows["A3_beta"].to_numpy()
    abs_betas = np.abs(betas)
    for i in range(len(betas)):
        for j in range(i + 1, len(betas)):
            if betas[i] * betas[j] < 0:                # opposite sign
                if abs_betas[i] == 0 or abs_betas[j] == 0:
                    continue
                ratio = abs_betas[i] / abs_betas[j]
                if MIRROR_TOL <= ratio <= 1 / MIRROR_TOL:
                    return True
    return False


def main():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    L, fh = log_open(LOG)
    L("== Strategy A Step 6 — Aggregate to cluster level ==\n")

    if not HITS.exists() or pd.read_csv(HITS).empty:
        L("No A_top_hits to aggregate.  Writing empty cluster_hits.")
        pd.DataFrame(columns=[
            "tissue","gene","cluster_id","n_junctions","n_sig_q10",
            "n_sig_q05","has_mirror_pair","top_q","top_beta",
            "confidence_stars"]).to_csv(OUT, index=False)
        fh.close()
        return

    hits = pd.read_csv(HITS)
    L(f"Loaded {len(hits):,} per-junction hit rows")

    summ = pd.read_csv(SUMR)
    L(f"Loaded {len(summ):,} (tissue, gene) summary rows")
    L("")

    # Add cluster_id column already exists in hits
    # Group per (tissue, gene, cluster) and compute stats
    rows = []
    for (t, g, clu), sub in hits.groupby(["tissue", "gene", "cluster_id"]):
        n_jct      = len(sub)
        n_sig_q10  = int((sub["A3_q"] < 0.10).sum() +
                          (sub["A1_q"] < 0.10).sum() +
                          (sub["A2_q"] < 0.10).sum())
        n_sig_q10_distinct = int(((sub[["A1_q","A2_q","A3_q"]].fillna(1) < 0.10)
                                    .any(axis=1)).sum())
        n_sig_q05_distinct = int(((sub[["A1_q","A2_q","A3_q"]].fillna(1) < 0.05)
                                    .any(axis=1)).sum())
        mirror     = detect_mirror_pair(sub)
        top_idx    = sub["A3_q"].idxmin()
        if pd.isna(top_idx):
            top_idx = sub["A1_q"].idxmin()
        if pd.isna(top_idx):
            top_idx = sub.index[0]
        top_row = sub.loc[top_idx]

        # Confidence stars
        if mirror:
            stars = 3
        elif n_sig_q05_distinct >= 2:
            stars = 3
        elif n_sig_q05_distinct >= 1:
            stars = 2
        elif n_sig_q10_distinct >= 1:
            stars = 1
        else:
            stars = 0

        rows.append({
            "tissue": t,
            "gene": g,
            "cluster_id": clu,
            "n_junctions_in_cluster_with_hit": n_jct,
            "n_sig_q10_distinct": n_sig_q10_distinct,
            "n_sig_q05_distinct": n_sig_q05_distinct,
            "has_mirror_pair": mirror,
            "top_junction": top_row.get("junction_id"),
            "top_A3_q": top_row.get("A3_q"),
            "top_A3_beta": top_row.get("A3_beta"),
            "top_A1_q": top_row.get("A1_q"),
            "top_A1_beta": top_row.get("A1_beta"),
            "confidence_stars": stars,
        })

    cluster_hits = pd.DataFrame(rows).sort_values(
        ["confidence_stars", "top_A3_q"],
        ascending=[False, True])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cluster_hits.to_csv(OUT, index=False)
    L(f"Wrote {OUT}  ({len(cluster_hits):,} clusters with ≥1 hit)")
    L("")

    # Summary statistics
    L("=== Cluster-level confidence breakdown ===")
    for stars in [3, 2, 1, 0]:
        n = (cluster_hits["confidence_stars"] == stars).sum()
        L(f"   {stars} stars:  {n:3d}")
    L("")

    L("Mirror-pair clusters:")
    for _, r in cluster_hits[cluster_hits["has_mirror_pair"]].iterrows():
        L(f"   {r['tissue']:25s} {r['gene']:12s} {r['cluster_id']:15s}  "
          f"top q={r['top_A3_q']:.2e}  β={r['top_A3_beta']:+.3f}")
    L("")

    L("Top 15 (tissue, gene, cluster) by confidence:")
    for _, r in cluster_hits.head(15).iterrows():
        L(f"   ⭐{r['confidence_stars']}  {r['tissue']:25s} {r['gene']:12s}  "
          f"{r['cluster_id']:15s}  q={r['top_A3_q']:.2e}  "
          f"sig_q05={r['n_sig_q05_distinct']}  mirror={r['has_mirror_pair']}")

    L("\nDONE.")
    fh.close()


if __name__ == "__main__":
    main()
