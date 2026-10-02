# -*- coding: utf-8 -*-
# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# Re-rank the 71 triplets by a visibility heuristic; three candidate-selection strategies.
#
# Original location in the analysis project:
#     0630/code/01_rescan_for_visible_signal.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
01_rescan_for_visible_signal.py
================================
Re-rank the 71 testable triplets to find candidates whose splicing difference
is most likely to be VISUALLY VISIBLE in carrier vs non-carrier pileup figures.

KEY FINDING from data exploration:
  - Max VAF in entire dataset = 0.058 (PLBD2 is already the highest!)
  - 0 candidates with VAF >= 0.10
  - Median VAF = 0.020
  - => SMaHT P25 cohort is fundamentally low-VAF
  - => Need to rank by COMBINED metrics, not just VAF

Three strategies:
  Option A — Conservative: SpliceAI >= 0.2, VAF >= 0.05 (3 hits)
  Option B — Aggressive:   SpliceAI >= 0.1, |z| >= 3 (8 hits)
  Option C — Visibility:   rank ALL 71 by visibility score, show top 10

Visibility score:
  visibility = |z| * VAF * |delta_psi| * log10(cluster_reads + 1)

Higher visibility => more likely the carrier-aberrant junction will appear
as a clearly enriched arc in a Cummings-style pileup figure.

Usage:
    python3 01_rescan_for_visible_signal.py
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT      = Path("/Users/spectremac/Desktop/Differential_splicing")
SCAN_CSV  = ROOT / "0619" / "output" / "tables" / "02_outlier_scan_all.csv"
OUT_DIR   = ROOT / "0630" / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("Re-ranking outlier triplets for VISIBLE splicing signal")
print("=" * 70)

df = pd.read_csv(SCAN_CSV)
df["abs_z"]         = df["z"].abs()
df["abs_delta_psi"] = df["psi_shift"].abs()

# Visibility score for everyone
df["visibility"] = (
    df["abs_z"] *
    df["SR_VAF"] *
    df["abs_delta_psi"] *
    np.log10(df["cluster_reads_median"] + 1)
)

print("\nLoaded " + str(len(df)) + " testable triplets")
print("\nVAF range: " + str(round(df['SR_VAF'].min(), 4)) +
      " - " + str(round(df['SR_VAF'].max(), 4)) +
      " (median " + str(round(df['SR_VAF'].median(), 4)) + ")")
print("Note: max VAF = " + str(round(df['SR_VAF'].max(), 4)) +
      " — there are NO variants with VAF >= 0.1 in this cohort")

DISPLAY_COLS = ["gene", "donor", "tissue", "SpliceAI", "SR_VAF",
                "abs_z", "abs_delta_psi", "cluster_reads_median",
                "visibility"]

# ─────────────────────────────────────────────────────────────────────────────
# Option A — Conservative
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("OPTION A — Conservative (SpliceAI >= 0.2, VAF >= 0.05, |z| >= 2)")
print("=" * 70)
mask_A = (
    (df["SpliceAI"] >= 0.2)  &
    (df["SR_VAF"]   >= 0.05) &
    (df["abs_z"]    >= 2)
)
A = df[mask_A].sort_values("visibility", ascending=False)
print("Pass: " + str(len(A)) + " triplets")
print(A[DISPLAY_COLS].to_string(index=False))
A.to_csv(OUT_DIR / "optionA_conservative.csv", index=False)

# ─────────────────────────────────────────────────────────────────────────────
# Option B — Aggressive
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("OPTION B — Aggressive (SpliceAI >= 0.1, |z| >= 3) — no VAF floor")
print("=" * 70)
mask_B = (
    (df["SpliceAI"] >= 0.1) &
    (df["abs_z"]    >= 3)
)
B = df[mask_B].sort_values("visibility", ascending=False)
print("Pass: " + str(len(B)) + " triplets")
print(B[DISPLAY_COLS].to_string(index=False))
B.to_csv(OUT_DIR / "optionB_aggressive.csv", index=False)

# ─────────────────────────────────────────────────────────────────────────────
# Option C — Top 10 by visibility (NO filter — pure ranking)
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("OPTION C — Top 10 by VISIBILITY SCORE (no pre-filter)")
print("=" * 70)
C = df.sort_values("visibility", ascending=False).head(10)
print(C[DISPLAY_COLS].to_string(index=False))
C.to_csv(OUT_DIR / "optionC_top10_visibility.csv", index=False)

# ─────────────────────────────────────────────────────────────────────────────
# Combined recommendation
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("RECOMMENDED CANDIDATES TO VISUALIZE")
print("=" * 70)

# Union of A, B, C top hits, dedup by (gene, donor, tissue)
combined = pd.concat([A, B, C.head(5)])
combined = combined.drop_duplicates(subset=["gene", "donor", "tissue"])
combined = combined.sort_values("visibility", ascending=False)
print(combined[DISPLAY_COLS + ["chr", "pos"]].to_string(index=False))
combined.to_csv(OUT_DIR / "recommended_candidates.csv", index=False)

# ─────────────────────────────────────────────────────────────────────────────
# Comparison: visibility scores of the original 3 hero cases
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("COMPARISON — original 3 hero cases vs new candidates")
print("=" * 70)
heroes_mask = df["gene"].isin(["PLBD2", "ISCU", "C11orf54"])
heroes = df[heroes_mask].sort_values("visibility", ascending=False)
print("\nOriginal 3 heroes (visibility scores):")
print(heroes[DISPLAY_COLS].to_string(index=False))

new_candidates = df[~heroes_mask].sort_values("visibility", ascending=False).head(5)
print("\nNew top candidates by visibility (NOT in original 3):")
print(new_candidates[DISPLAY_COLS].to_string(index=False))

print("\n" + "=" * 70)
print("Output files written to: " + str(OUT_DIR))
print("=" * 70)
