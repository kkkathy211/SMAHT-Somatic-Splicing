# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# Inventory donors/junctions/clusters per tissue; count testable (donor, tissue, variant) triplets.
#
# Original location in the analysis project:
#     0619/code/01_inspect_data.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
01_inspect_data.py
==================

Research question
-----------------
Before any expensive analysis: do all 10 perind files load cleanly, how many
donors does each tissue have, and — most importantly — how many of the V2
SpliceAI (donor, tissue) pairs can we actually test (i.e., the donor has RNA
in the matched tissue's perind)?

Output
------
- output/tables/01_donor_tissue_inventory.csv
- console: per-tissue summary + V2-vs-perind match summary
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import pandas as pd
from _utils import TISSUES, TISSUE_PREFIX, load_perind, load_v2

OUT = ROOT / "output" / "tables" / "01_donor_tissue_inventory.csv"

print("=" * 70)
print("STEP 01 — perind & V2 inventory")
print("=" * 70)

# --- 1. Per-tissue inventory ---
rows = []
perind_donors = {}  # tissue → set of donors
for tissue in TISSUES:
    p = load_perind(tissue)
    donors = sorted(set(p["donors"]))
    n_samples = len(p["samples"])
    n_donors  = len(donors)
    n_jct     = len(p["junctions"])
    n_clusters = len({j["clu"] for j in p["junctions"]})
    perind_donors[tissue] = set(donors)
    rows.append({
        "tissue": tissue, "n_donors": n_donors, "n_samples": n_samples,
        "cores_per_donor": round(n_samples / n_donors, 2),
        "n_junctions": n_jct, "n_clusters": n_clusters,
        "donors": ",".join(donors),
    })
    print(f"  {tissue:30s}  {n_donors:3d} donors  {n_samples:3d} samples  "
          f"{n_clusters:6d} clusters  {n_jct:7d} junctions")

inv_df = pd.DataFrame(rows)

# --- 2. V2 vs perind cross-reference ---
print("\n" + "=" * 70)
print("V2 SpliceAI ≥ 0.2 → match to perind tissues")
print("=" * 70)

v2 = load_v2(min_spliceai=0.2)
print(f"  Total V2 variants with SpliceAI ≥ 0.2: {len(v2)}")
print(f"  Unique donors in V2: {v2['SampleDonors'].nunique()}")
print(f"  Unique TissueIDs in V2: {v2['TissueID'].nunique()}")

# Matchable: V2 row has both (a) tissue with perind, (b) donor in that perind
v2["tissue_full"] = v2["TissueID"].map(TISSUE_PREFIX)
v2["has_perind"]  = v2["tissue_full"].notna()
v2["donor_in_perind"] = False
for i, row in v2.iterrows():
    if row["has_perind"]:
        v2.at[i, "donor_in_perind"] = row["SampleDonors"] in perind_donors[row["tissue_full"]]

testable = v2[v2["donor_in_perind"]]
print(f"\n  Testable (V2 row + donor has RNA in matched tissue): {len(testable)}")
print(f"  Unique testable (donor, tissue, variant) triplets: "
      f"{testable.groupby(['SampleDonors','TissueID','id']).ngroups}")

print("\n  Testable count by SpliceAI bin:")
for lo, hi, label in [(0.2, 0.5, "LOOSE 0.2-0.5"),
                       (0.5, 0.8, "STRICT 0.5-0.8"),
                       (0.8, 1.01, "high-conf ≥0.8")]:
    n = ((testable["score"] >= lo) & (testable["score"] < hi)).sum()
    print(f"    {label:20s}: {n}")

print("\n  Testable count by VAF tier:")
for lo, hi, label in [(0,    0.05, "VAF<0.05"),
                       (0.05, 0.10, "0.05≤VAF<0.10"),
                       (0.10, 1.01, "VAF≥0.10")]:
    n = ((testable["SR_VAF"] >= lo) & (testable["SR_VAF"] < hi)).sum()
    print(f"    {label:20s}: {n}")

print("\n  SMHT029 specific:")
s29 = testable[testable["SampleDonors"] == "SMHT029"]
print(f"    Testable SMHT029 variant-tissue pairs: {len(s29)}")
for _, r in s29.sort_values(["score", "SR_VAF"], ascending=False).iterrows():
    print(f"    {r['gene']:12s}  {r['TissueID']:5s}  "
          f"SpliceAI={r['score']:.2f}  VAF={r['SR_VAF']:.3f}  {r['id']}")

# --- 3. Save inventory ---
inv_df.to_csv(OUT, index=False)
print(f"\n  → wrote {OUT}")
