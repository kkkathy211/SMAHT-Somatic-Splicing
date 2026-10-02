# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# SHARED ENGINE. perind parsing, per-donor psi, cluster overlap, MAD robust z-score.
#
# Original location in the analysis project:
#     0619/code/_utils.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
Shared utilities for V2 SpliceAI × LeafCutter per-donor outlier analysis.

Conventions
-----------
- A "perind" file is LeafCutter's per-individual junction count file.
  Row 1: header with sample names (full filenames).
  Row 2+: first column = "chrom:start:end:clu_<id>_<strand>",
          remaining columns = "a/b" where a = junction reads, b = cluster total reads.
- A "donor" is identified by SMHT### prefix; multiple cores per donor are
  averaged for the donor-level psi.
"""

from __future__ import annotations
import gzip, re
from pathlib import Path
import numpy as np
import pandas as pd

DATA_DIR   = Path("/Users/spectremac/Desktop/Differential_splicing/0619/data")
PERIND_DIR = DATA_DIR / "perind"
V2_FILE    = DATA_DIR / "V2_spliceAI_VAF_quantile_fraction.txt.gz"

# All tissues with available perind (Laura's filtered_junc_files output)
TISSUES = [
    "3AD_Skin_Calf", "3AF_Skin_Abdomen", "3AH_Muscle",
    "3AK_Brain_Frontal_Lobe", "3AM_Brain_Cerebellum", "3C_Esophagus",
    "3G_Colon_Desc", "3Q_Lung", "3S_Heart", "3U_Testis_L",
]
# Map TissueID prefix → full tissue name
TISSUE_PREFIX = {t.split("_")[0]: t for t in TISSUES}

DONOR_RE = re.compile(r"(SMHT\d+)")


def donor_of(sample_name: str) -> str:
    """Extract SMHT### from a full perind sample column name."""
    m = DONOR_RE.match(sample_name)
    return m.group(1) if m else sample_name


def load_perind(tissue: str) -> dict:
    """Load a perind.counts.gz. Returns dict with:
        - junctions: list of dicts {chrom, start, end, clu, idx}
        - samples: list of full sample column names
        - donors: list of donor SMHT### per sample
        - counts: junction × sample matrix of (a, b) tuples
        - psi: junction × sample matrix of float (a/b)
    """
    path = PERIND_DIR / f"{tissue}_perind.counts.gz"
    with gzip.open(path, "rt") as fh:
        header = fh.readline().strip().split()
        samples = header[1:]
        junctions, rows_a, rows_b = [], [], []
        for ln in fh:
            parts = ln.strip().split()
            chrom, start, end, clu = parts[0].split(":")
            junctions.append({
                "chrom": chrom, "start": int(start), "end": int(end),
                "clu": clu, "idx": len(junctions),
            })
            a_row, b_row = [], []
            for v in parts[1:]:
                a, b = v.split("/")
                a_row.append(int(a))
                b_row.append(int(b))
            rows_a.append(a_row)
            rows_b.append(b_row)
    a_mat = np.array(rows_a, dtype=np.int32)
    b_mat = np.array(rows_b, dtype=np.int32)
    with np.errstate(divide="ignore", invalid="ignore"):
        psi = np.where(b_mat > 0, a_mat / b_mat, np.nan)
    donors = [donor_of(s) for s in samples]
    return {
        "tissue": tissue, "junctions": junctions, "samples": samples,
        "donors": donors, "a": a_mat, "b": b_mat, "psi": psi,
    }


def per_donor_psi(perind: dict) -> tuple[np.ndarray, list[str]]:
    """Collapse per-sample (multiple cores per donor) → per-donor psi (mean across cores).
    Returns (psi: junction × donor, donor_list)."""
    donors = sorted(set(perind["donors"]))
    n_j = perind["psi"].shape[0]
    out = np.full((n_j, len(donors)), np.nan)
    for d_i, d in enumerate(donors):
        sample_idx = [i for i, s_d in enumerate(perind["donors"]) if s_d == d]
        out[:, d_i] = np.nanmean(perind["psi"][:, sample_idx], axis=1)
    return out, donors


def find_clusters_overlapping(perind: dict, chrom: str, pos: int,
                                window: int = 2000) -> list[str]:
    """Return list of cluster IDs whose junction span includes (pos ± window).
    Note: matches any junction in the cluster that brackets the variant — the
    variant doesn't have to be at a junction endpoint, it just has to fall
    within the intron defined by that junction (or close to it)."""
    hits = set()
    for j in perind["junctions"]:
        if j["chrom"] != chrom: continue
        if j["start"] <= pos + window and j["end"] >= pos - window:
            hits.add(j["clu"])
    return sorted(hits)


def junctions_of_cluster(perind: dict, clu: str) -> list[dict]:
    """All junctions belonging to a cluster."""
    return [j for j in perind["junctions"] if j["clu"] == clu]


def robust_zscore(donor_val: float, others: np.ndarray) -> float:
    """MAD-based robust z; falls back to mean/SD if MAD is 0/nan."""
    others = others[~np.isnan(others)]
    if len(others) < 2: return np.nan
    med = np.nanmedian(others)
    mad = np.nanmedian(np.abs(others - med)) * 1.4826
    if mad == 0 or np.isnan(mad):
        sd = np.nanstd(others, ddof=1)
        if sd == 0 or np.isnan(sd): return np.nan
        return (donor_val - med) / sd
    return (donor_val - med) / mad


def load_v2(min_spliceai: float = 0.0) -> pd.DataFrame:
    """Load V2 SpliceAI table, optionally filtered by min score."""
    df = pd.read_csv(V2_FILE, sep="\t")
    if min_spliceai > 0:
        df = df[df["score"] >= min_spliceai].reset_index(drop=True)
    return df
