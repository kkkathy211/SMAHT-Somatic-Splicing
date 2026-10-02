# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# CORE ANALYSIS. MAD robust z per junction + the |dPSI| < VAF physical-plausibility filter.
#
# Original location in the analysis project:
#     0619/code/02_outlier_scan.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
02_outlier_scan.py
==================

Research question
-----------------
For every V2 SpliceAI variant (score ≥ 0.2) whose carrier has matched RNA-seq
in the variant's tissue: is the carrier's intron-usage ratio (psi) in the
LeafCutter cluster(s) overlapping that variant an outlier vs the cohort?

Approach
--------
1. For each testable triplet (donor, tissue, variant):
   a. Find LeafCutter clusters overlapping the variant ± window
   b. For each cluster, require min cohort coverage (median cluster reads ≥ 10)
   c. For each junction in the cluster, compute per-donor psi (mean across cores)
   d. Compute MAD-based |z| for the carrier on each junction
   e. Record the junction with max |z| as the cluster's signal
2. Output: ranked table of (variant, donor, tissue) with best cluster's |z|

VAF-aware reporting
-------------------
Carrier signal is physically capped at ~VAF magnitude in psi (e.g., VAF=0.05
→ max ~5% psi shift even if SpliceAI prediction is perfect). We report:
  - z (MAD-based)
  - psi_shift = psi_carrier - psi_others_median  (signed)
  - shift_ratio = |psi_shift| / VAF  (signal vs physical max; > 1 = suspicious
    of false positive or cohort variance issues)

Output
------
- output/tables/02_outlier_scan_all.csv  (one row per testable triplet)
"""
import sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import numpy as np
import pandas as pd
from _utils import (TISSUES, TISSUE_PREFIX, load_perind, load_v2,
                     per_donor_psi, find_clusters_overlapping,
                     junctions_of_cluster, robust_zscore)

OUT = ROOT / "output" / "tables" / "02_outlier_scan_all.csv"
WINDOW          = 2000   # bp around variant to search for cluster
MIN_CLUSTER_RD  = 10     # require cohort median cluster reads ≥ this

print("=" * 70)
print("STEP 02 — per-donor outlier scan (V2 ≥ 0.2)")
print("=" * 70)

# ---- 1. Load V2 + filter to testable ----
v2 = load_v2(min_spliceai=0.2)
v2["tissue_full"] = v2["TissueID"].map(TISSUE_PREFIX)
v2 = v2[v2["tissue_full"].notna()].reset_index(drop=True)
print(f"  V2 rows w/ matched-tissue perind available: {len(v2)}")

# ---- 2. Cache perind loads (one per tissue) ----
print("\n  Loading perind files...")
caches = {}
for t in TISSUES:
    t0 = time.time()
    p = load_perind(t)
    psi_d, donors_d = per_donor_psi(p)
    caches[t] = {"perind": p, "psi_donor": psi_d, "donors": donors_d}
    print(f"    {t:30s}  loaded in {time.time()-t0:.1f}s  "
          f"({len(donors_d)} donors, {p['psi'].shape[0]} junctions)")

# ---- 3. Scan ----
print("\n  Scanning testable triplets...")
rows = []
skipped = {"donor_not_in_perind": 0, "no_cluster_overlap": 0,
            "low_coverage": 0, "ok": 0}

for idx, r in v2.iterrows():
    tissue = r["tissue_full"]
    donor  = r["SampleDonors"]
    c = caches[tissue]
    if donor not in c["donors"]:
        skipped["donor_not_in_perind"] += 1; continue
    d_i = c["donors"].index(donor)

    clusters = find_clusters_overlapping(c["perind"], r["chr"], int(r["pos"]),
                                            window=WINDOW)
    if not clusters:
        skipped["no_cluster_overlap"] += 1; continue

    # Pick the best (highest |z|) cluster
    best = None
    for clu in clusters:
        jcts = junctions_of_cluster(c["perind"], clu)
        j_idx = [j["idx"] for j in jcts]
        # cohort coverage check (median cluster total reads across all samples)
        cluster_totals_per_sample = c["perind"]["b"][j_idx[0], :]
        if np.median(cluster_totals_per_sample) < MIN_CLUSTER_RD: continue

        # per-junction z on this donor
        psi_d = c["psi_donor"][j_idx, :]  # n_jct × n_donors
        carrier = psi_d[:, d_i]
        others_mask = np.array([i != d_i for i in range(psi_d.shape[1])])
        for k, j in enumerate(jcts):
            others = psi_d[k, others_mask]
            z = robust_zscore(carrier[k], others)
            if np.isnan(z): continue
            shift = carrier[k] - np.nanmedian(others)
            ratio = abs(shift) / r["SR_VAF"] if r["SR_VAF"] > 0 else np.nan
            cand = {
                "cluster": clu, "jct_start": j["start"], "jct_end": j["end"],
                "jct_len": j["end"] - j["start"],
                "psi_carrier": float(carrier[k]),
                "psi_others_median": float(np.nanmedian(others)),
                "psi_shift": float(shift), "z": float(z),
                "shift_over_vaf": float(ratio),
                "n_donors_cohort": int(others_mask.sum()),
                "cluster_reads_median": float(np.median(cluster_totals_per_sample)),
            }
            if best is None or abs(cand["z"]) > abs(best["z"]):
                best = cand
    if best is None:
        skipped["low_coverage"] += 1; continue

    rows.append({
        "gene": r["gene"], "donor": donor, "tissue": tissue,
        "variant_id": r["id"], "chr": r["chr"], "pos": int(r["pos"]),
        "SpliceAI": float(r["score"]), "SR_VAF": float(r["SR_VAF"]),
        **best,
    })
    skipped["ok"] += 1

print(f"\n  Done. Outcomes:")
for k, v in skipped.items():
    print(f"    {k:25s}: {v}")

# ---- 4. Save ranked output ----
out_df = pd.DataFrame(rows).sort_values("z", key=lambda s: s.abs(), ascending=False)
out_df.to_csv(OUT, index=False)
print(f"\n  → wrote {OUT}  ({len(out_df)} rows)")

# ---- 5. Print TOP 20 + SMHT029 highlights ----
print("\n  TOP 20 by |z|:")
cols_show = ["gene", "donor", "tissue", "SpliceAI", "SR_VAF",
              "psi_carrier", "psi_others_median", "z", "shift_over_vaf"]
print(out_df[cols_show].head(20).to_string(index=False))

print("\n  SMHT029 results (all triplets):")
s29 = out_df[out_df["donor"] == "SMHT029"].sort_values(
    "z", key=lambda s: s.abs(), ascending=False)
print(s29[cols_show].to_string(index=False))
