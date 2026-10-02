# -*- coding: utf-8 -*-
# ============================================================================
# Phase 3 - Prediction-first candidate selection
#
# NEGATIVE CONTROL. Outlier scan of 4 large genes with no SpliceAI filter. All 273 hits have SpliceAI 0.00.
#
# Original location in the analysis project:
#     0716/code/08_unsupervised_biggene_scan.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
08_unsupervised_biggene_scan.py
================================
Unsupervised per-junction outlier scan of 4 big genes, WITHOUT SpliceAI filter.

Motivation: SpliceAI may miss splice-disrupting variants (esp. deep intronic,
regulatory-element, or splicing-factor knock-down effects). For each big gene:
  For each junction in the gene:
    Compute psi per donor
    Compute robust z (MAD)
    Flag donors with |z| >= 2
  For each outlier donor:
    Check if this donor has ANY somatic variant in this gene
    If yes -> HIT (candidate causal variant)

Genes:
  CNTNAP2  chr7:146,116,002-148,420,998   2.30 Mb
  CSMD1    chr8:2,935,353-4,994,972       2.06 Mb
  RBFOX1   chr16:5,239,802-7,713,340      2.47 Mb  (splicing factor)
  DMD      chrX:31,097,677-33,339,609     2.24 Mb
"""

import gzip
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

ROOT = Path("/Users/spectremac/Desktop/Differential_splicing")
LEAF_DIR = ROOT / "leafcutter_results"
V2_FILE  = ROOT / "0619" / "data" / "V2_spliceAI_VAF_quantile_fraction.txt"
OUT_DIR  = ROOT / "0716" / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Gene regions (GRCh38, from gencode.v44)
GENES = {
    "CNTNAP2": ("chr7",  146116002, 148420998),
    "CSMD1":   ("chr8",    2935353,   4994972),
    "RBFOX1":  ("chr16",   5239802,   7713340),
    "DMD":     ("chrX",   31097677,  33339609),
}

# Tissues we have RNA-seq perind for
TISSUES = [
    "3AD_Skin_Calf", "3AF_Skin_Abdomen", "3AH_Muscle",
    "3AK_Brain_Frontal_Lobe", "3AL_Brain_Temporal_Lobe",
    "3AM_Brain_Cerebellum", "3C_Esophagus", "3G_Colon_Desc",
    "3Q_Lung", "3S_Heart", "3U_Testis_L",
]
# Map to TissueID used in V2 file
TISSUE_ID = {t: t.split("_")[0] for t in TISSUES}


# ─────────────────────────────────────────────────────────────────────────────
# 1. Load V2 somatic variants for the 4 genes: which (donor, tissue) has any?
# ─────────────────────────────────────────────────────────────────────────────
print("=== Loading V2 somatic variants for 4 big genes ===")
v2 = pd.read_csv(V2_FILE, sep="\t")
v2_big = v2[v2["gene"].isin(GENES.keys())].copy()
print(f"Total somatic variants in these 4 genes: {len(v2_big)}")

# Build (gene, donor, tissue) -> list of variants
carrier_map = defaultdict(list)   # (gene, donor, tissueID) -> [(pos, score, ref, alt), ...]
for _, r in v2_big.iterrows():
    key = (r["gene"], r["SampleDonors"], r["TissueID"])
    carrier_map[key].append((int(r["pos"]), float(r["score"]),
                             r["ref"], r["alt"], float(r["SR_VAF"])))
print(f"Unique (gene, donor, tissue) carrier keys: {len(carrier_map)}")


# ─────────────────────────────────────────────────────────────────────────────
# 2. Robust z-score (MAD-based, matching prior scan methodology)
# ─────────────────────────────────────────────────────────────────────────────
def robust_z(carrier_val, cohort_vals):
    """MAD-based robust z-score; returns 0 if MAD == 0 or too few values."""
    cohort_vals = np.asarray([v for v in cohort_vals if not np.isnan(v)])
    if len(cohort_vals) < 5:
        return np.nan
    med = np.median(cohort_vals)
    mad = np.median(np.abs(cohort_vals - med))
    if mad == 0:
        return np.nan
    return 0.6745 * (carrier_val - med) / mad


def psi_from_kn(s):
    try:
        k, n = s.split("/")
        k, n = float(k), float(n)
        return k / n if n > 0 else np.nan
    except Exception:
        return np.nan


# ─────────────────────────────────────────────────────────────────────────────
# 3. For each tissue, scan each big gene
# ─────────────────────────────────────────────────────────────────────────────
Z_CUT = 2.0
MIN_CLUSTER_READS = 20    # avoid low-count junctions

all_hits = []
per_tissue_stats = []

for tissue in TISSUES:
    tid = TISSUE_ID[tissue]
    perind = LEAF_DIR / tissue / f"{tissue}_perind.counts.gz"
    if not perind.exists():
        continue

    print(f"\n=== {tissue} (TissueID={tid}) ===")

    with gzip.open(str(perind), "rt") as f:
        header = f.readline().strip().split()
        sample_cols = header[1:]
        donors = [s.split("-")[0] for s in sample_cols]

        tissue_stats = {"tissue": tissue, "genes": {}}
        for line in f:
            parts = line.strip().split()
            jid = parts[0]           # e.g. chr7:146200000:146201000:clu_1234_+
            fields = jid.split(":")
            if len(fields) < 4:
                continue
            jchr, jstart, jend, _clu = fields[0], int(fields[1]), int(fields[2]), fields[3]

            # Which gene does this junction fall in?
            for gene, (gchr, glo, ghi) in GENES.items():
                if jchr != gchr:
                    continue
                if not (glo <= jstart and jend <= ghi):
                    continue

                # Parse K/N per donor
                kn = parts[1:]
                psi = np.array([psi_from_kn(v) for v in kn])
                # cohort read count sanity
                n_reads = np.array([
                    float(v.split("/")[1]) if v not in (".", "NA") else np.nan
                    for v in kn
                ])
                if np.nanmedian(n_reads) < MIN_CLUSTER_READS:
                    continue

                # For each donor, compute z vs the rest
                for i, d in enumerate(donors):
                    if np.isnan(psi[i]):
                        continue
                    rest = np.delete(psi, i)
                    z = robust_z(psi[i], rest)
                    if np.isnan(z) or abs(z) < Z_CUT:
                        continue
                    # Cross-ref: does this (gene, donor, tissue) have variants?
                    variants = carrier_map.get((gene, d, tid), [])
                    if not variants:
                        continue
                    # HIT
                    all_hits.append({
                        "gene":        gene,
                        "tissue":      tissue,
                        "donor":       d,
                        "junction":    f"{jchr}:{jstart}-{jend}",
                        "psi_carrier": round(float(psi[i]), 4),
                        "psi_median":  round(float(np.nanmedian(rest)), 4),
                        "z":           round(float(z), 3),
                        "n_reads_med": int(np.nanmedian(n_reads)),
                        "n_variants":  len(variants),
                        "max_score":   max(v[1] for v in variants),
                        "max_VAF":     max(v[4] for v in variants),
                        "variant_positions": ",".join(str(v[0]) for v in variants[:5]),
                    })
                tissue_stats["genes"].setdefault(gene, 0)
                tissue_stats["genes"][gene] += 1
    per_tissue_stats.append(tissue_stats)


# ─────────────────────────────────────────────────────────────────────────────
# 4. Output
# ─────────────────────────────────────────────────────────────────────────────
if all_hits:
    hits_df = pd.DataFrame(all_hits).sort_values("z", key=abs, ascending=False)
    out_csv = OUT_DIR / "biggene_unsupervised_hits.csv"
    hits_df.to_csv(out_csv, index=False)
    print(f"\n{'='*70}")
    print(f"TOTAL HITS (|z|>={Z_CUT}, carrier has variant in gene): {len(hits_df)}")
    print(f"{'='*70}")
    print(hits_df.to_string(index=False))
    print(f"\nSaved: {out_csv}")

    # Summary per gene
    print(f"\n{'='*70}")
    print("SUMMARY BY GENE")
    print(f"{'='*70}")
    summary = hits_df.groupby("gene").agg(
        n_hits=("z", "size"),
        max_abs_z=("z", lambda x: max(abs(x))),
        n_donors=("donor", "nunique"),
        n_tissues=("tissue", "nunique"),
    ).sort_values("n_hits", ascending=False)
    print(summary.to_string())
else:
    print("\n=== NO HITS ===")
    print("No junction in any of the 4 big genes had a carrier |z|>=2 outlier.")

print("\nDone.")
