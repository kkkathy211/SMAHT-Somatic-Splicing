# -*- coding: utf-8 -*-
# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# Six-way filter that removes low-count and distance artifacts. Exactly one candidate survives.
#
# Original location in the analysis project:
#     0630/code/04_strict_visibility_filter.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
04_strict_visibility_filter.py
===============================
Apply a STRICTER visibility filter to weed out candidates whose "large Δψ"
comes from expression artifacts rather than true splicing changes.

Filter:
  1. SpliceAI          >= 0.2      (predicted splice-disrupting)
  2. VAF               >= 0.02     (some clonal signal)
  3. |z|               >= 3        (statistical outlier)
  4. cluster_reads_median >= 200   (sufficient cohort expression — avoids
                                     low-count artifacts like CAPN11 with 18 reads)
  5. |Δψ|              >= 0.03     (real ratio shift)
  6. variant inside/near outlier junction span
     (proxy for "variant inside cluster span" — ensures direct causal link)

Compare results to original 3 hero cases and to relaxed filter results.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT     = Path("/Users/spectremac/Desktop/Differential_splicing")
SCAN_CSV = ROOT / "0619" / "output" / "tables" / "02_outlier_scan_all.csv"
OUT_DIR  = ROOT / "0630" / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SCAN_CSV)
df["abs_z"]         = df["z"].abs()
df["abs_delta_psi"] = df["psi_shift"].abs()

# Distance of variant to outlier junction (0 if inside the junction span)
df["dist_to_junction"] = np.where(
    (df["pos"] >= df["jct_start"]) & (df["pos"] <= df["jct_end"]),
    0,
    np.minimum((df["pos"] - df["jct_start"]).abs(),
               (df["pos"] - df["jct_end"]).abs()),
)

# Visibility (kept for reference)
df["visibility"] = (
    df["abs_z"] *
    df["SR_VAF"] *
    df["abs_delta_psi"] *
    np.log10(df["cluster_reads_median"] + 1)
)

DISPLAY = ["gene", "donor", "tissue", "SpliceAI", "SR_VAF",
           "abs_z", "abs_delta_psi", "cluster_reads_median",
           "dist_to_junction", "visibility"]

print("=" * 70)
print("STRICT visibility filter — remove expression artifacts")
print("=" * 70)

# Apply filter step by step to see attrition
print("\nAttrition per criterion:")
n = len(df); print("  Start:                     " + str(n))
n = (df["SpliceAI"] >= 0.2).sum(); print("  SpliceAI >= 0.2:           " + str(n))
n = ((df["SpliceAI"] >= 0.2) & (df["SR_VAF"] >= 0.02)).sum(); print("  + VAF >= 0.02:             " + str(n))
n = ((df["SpliceAI"] >= 0.2) & (df["SR_VAF"] >= 0.02) &
     (df["abs_z"] >= 3)).sum(); print("  + |z| >= 3:                " + str(n))
n = ((df["SpliceAI"] >= 0.2) & (df["SR_VAF"] >= 0.02) &
     (df["abs_z"] >= 3) & (df["cluster_reads_median"] >= 200)).sum()
print("  + cluster_reads >= 200:    " + str(n))
n = ((df["SpliceAI"] >= 0.2) & (df["SR_VAF"] >= 0.02) &
     (df["abs_z"] >= 3) & (df["cluster_reads_median"] >= 200) &
     (df["abs_delta_psi"] >= 0.03)).sum()
print("  + |Δψ| >= 0.03:            " + str(n))

mask = (
    (df["SpliceAI"]              >= 0.2)  &
    (df["SR_VAF"]                >= 0.02) &
    (df["abs_z"]                 >= 3)    &
    (df["cluster_reads_median"]  >= 200)  &
    (df["abs_delta_psi"]         >= 0.03) &
    (df["dist_to_junction"]      <= 500)      # variant inside/very near junction
)
n = mask.sum()
print("  + variant near/inside jct: " + str(n) + "  <== FINAL")

# ─────────────────────────────────────────────────────────────────────────────
# Show what passes
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("CANDIDATES PASSING STRICT FILTER")
print("=" * 70)
strict = df[mask].sort_values("visibility", ascending=False)
if len(strict) == 0:
    print("None. Filter is too strict.")
else:
    print(strict[DISPLAY].to_string(index=False))
strict.to_csv(OUT_DIR / "strict_visibility_filtered.csv", index=False)

# ─────────────────────────────────────────────────────────────────────────────
# What was DROPPED that had passed the relaxed filter?
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("EXPLAINING WHY NEW CANDIDATES WERE DROPPED")
print("=" * 70)
new_ones = df[df["gene"].isin(["MYO15B", "CAPN11", "MCM7", "CADPS2",
                                "SLC43A3", "SGIP1", "SHROOM2",
                                "SNX11", "KDM4C", "PRKAG1", "PFKP",
                                "P4HB", "WDR1"])].copy()
new_ones["passes_strict"] = mask[new_ones.index].values
new_ones["fail_reason"]   = ""
new_ones.loc[new_ones["SpliceAI"]              < 0.2,  "fail_reason"] += "SpliceAI<0.2 "
new_ones.loc[new_ones["SR_VAF"]                < 0.02, "fail_reason"] += "VAF<0.02 "
new_ones.loc[new_ones["abs_z"]                 < 3,    "fail_reason"] += "|z|<3 "
new_ones.loc[new_ones["cluster_reads_median"]  < 200,  "fail_reason"] += "reads<200 "
new_ones.loc[new_ones["abs_delta_psi"]         < 0.03, "fail_reason"] += "Δψ<0.03 "
new_ones.loc[new_ones["dist_to_junction"]      > 500,  "fail_reason"] += "dist>500 "

print("\nWhy each new candidate is DROPPED:")
print(new_ones[["gene", "donor", "SpliceAI", "SR_VAF", "abs_z",
                "abs_delta_psi", "cluster_reads_median",
                "dist_to_junction", "passes_strict", "fail_reason"]]
      .to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# Compare to original 3 hero cases
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ORIGINAL 3 HEROES — do they pass the strict filter?")
print("=" * 70)
heroes = df[df["gene"].isin(["PLBD2", "ISCU", "C11orf54"])].copy()
heroes["passes_strict"] = mask[heroes.index].values
heroes["fail_reason"]   = ""
heroes.loc[heroes["SpliceAI"]              < 0.2,  "fail_reason"] += "SpliceAI<0.2 "
heroes.loc[heroes["SR_VAF"]                < 0.02, "fail_reason"] += "VAF<0.02 "
heroes.loc[heroes["abs_z"]                 < 3,    "fail_reason"] += "|z|<3 "
heroes.loc[heroes["cluster_reads_median"]  < 200,  "fail_reason"] += "reads<200 "
heroes.loc[heroes["abs_delta_psi"]         < 0.03, "fail_reason"] += "Δψ<0.03 "
heroes.loc[heroes["dist_to_junction"]      > 500,  "fail_reason"] += "dist>500 "
print(heroes[["gene", "donor", "SpliceAI", "SR_VAF", "abs_z",
              "abs_delta_psi", "cluster_reads_median",
              "dist_to_junction", "passes_strict", "fail_reason"]]
      .to_string(index=False))

print("\nDone.")
