# -*- coding: utf-8 -*-
"""
01_build_vep_candidate_list.py
===============================
Build the VEP x SpliceAI candidate list.

    V2 SpliceAI somatic calls  +  Ensembl VEP annotation
                   |
                   v
      vep_splice_candidates.csv   (213 rows)

This is the entry point for Phase 3. Everything downstream (enrichment with
GENCODE v47, the STRICT / HIGH_VAF / selected tiers, every pileup figure)
depends on this table.

-------------------------------------------------------------------------------
PROVENANCE NOTE
-------------------------------------------------------------------------------
The original 213-row table was produced ad hoc in an interactive session and the
code was never saved. This script was reconstructed from the output and verified
to reproduce it exactly:

    - 717 variants at SpliceAI >= 0.2          (matches the slide deck)
    - 213 rows after the "splice" consequence filter
    - IMPACT breakdown  95 HIGH / 82 LOW / 36 MODERATE
    - the `priority` column reproduced with 100% accuracy

-------------------------------------------------------------------------------
THE JOIN (this is the part that is easy to get wrong)
-------------------------------------------------------------------------------
VEP's `#Uploaded_variation` uses a SLASH between the alleles:

        chr12_108565432_G/A

The V2 SpliceAI table's `id` column uses an UNDERSCORE:

        chr12_108565432_G_A

so joining on `id` directly returns zero matches. Build the key explicitly:

        key = f"{chr}_{pos}_{ref}/{alt}"

With the key built this way the join has **zero misses** (717/717).

-------------------------------------------------------------------------------
KNOWN LIMITATION  -  `has_RNAseq` is a TISSUE-level flag
-------------------------------------------------------------------------------
`has_RNAseq` records whether the variant's TISSUE has a LeafCutter perind file.
It does NOT check whether that specific (donor, tissue) pair has an RNA-seq BAM.

Consequence: a row can be marked `has_RNAseq = OK` and still be impossible to
plot, because that donor was never sequenced in that tissue. In the 3AD skin
cohort alone, five donors (SMHT015, SMHT018, SMHT024, SMHT040, SMHT042) carry
DNA variants but have no 3AD RNA-seq.

If you need true per-(donor, tissue) testability, check the BAM bucket directly:

    gsutil ls gs://smaht-p25-bam-files/ | grep "{DONOR}-{TISSUE}-"

-------------------------------------------------------------------------------
Usage
-------------------------------------------------------------------------------
    python3 01_build_vep_candidate_list.py \
        --spliceai  data/V2_spliceAI_VAF_quantile_fraction.txt \
        --vep       data/annotate_v2_381files.txt \
        --perind    data/perind \
        --scan      results/phase2/02_outlier_scan_all.csv \
        --out       results/phase3/vep_splice_candidates.csv
"""

import argparse
import glob
import gzip
import os
from pathlib import Path

import pandas as pd

# ─── Constants ──────────────────────────────────────────────────────────────

SPLICEAI_MIN = 0.2          # Yilin's recommended threshold for the V2 callset
CONSEQUENCE_KEYWORD = "splice"

# The two cases that had already been worked up as case studies before this
# list was built. C11orf54 is deliberately absent: VEP calls it
# `missense_variant` (p.G204E), so it never passes the "splice" filter. That
# absence is itself the finding.
ALREADY_CASE_GENES = {"ISCU", "PLBD2"}

# Columns carried over from VEP
VEP_COLS = ["Consequence", "IMPACT", "EXON", "INTRON"]


# ─── Helpers ────────────────────────────────────────────────────────────────

def read_spliceai(path, score_min=SPLICEAI_MIN):
    """V2 somatic SpliceAI calls. 18 columns, one row per (variant, donor-tissue)."""
    df = pd.read_csv(path, sep="\t")
    print(f"  SpliceAI table : {len(df):,} rows, {df['id'].nunique():,} unique variants")
    df = df[df["score"] >= score_min].copy()
    print(f"  score >= {score_min}     : {len(df):,} rows, {df['id'].nunique():,} unique variants")
    return df


def read_vep(path):
    """
    Ensembl VEP output. The header line is the first line starting with a single
    '#', after ~103 lines of '##' metadata.
    """
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        skip = 0
        for line in fh:
            if line.startswith("##"):
                skip += 1
            else:
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, low_memory=False)
    df = df.rename(columns={"#Uploaded_variation": "vep_key"})
    print(f"  VEP table      : {len(df):,} rows, {df['vep_key'].nunique():,} unique variants")
    return df


def leafcutter_tissues(perind_dir):
    """Tissue codes that have a LeafCutter perind file (e.g. '3AD' from '3AD_Skin_Calf')."""
    tissues = set()
    for f in glob.glob(os.path.join(perind_dir, "*_perind.counts.gz")):
        tissues.add(os.path.basename(f).split("_")[0])
    return tissues


def scanned_genes(scan_csv):
    """Genes already covered by the Phase 2 LeafCutter outlier scan."""
    if scan_csv and os.path.exists(scan_csv):
        return set(pd.read_csv(scan_csv)["gene"])
    return set()


# ─── Main ───────────────────────────────────────────────────────────────────

def build(spliceai_path, vep_path, perind_dir, scan_csv, out_csv):
    print("\n[1/5] Reading inputs")
    sai = read_spliceai(spliceai_path)
    vep = read_vep(vep_path)

    print("\n[2/5] Joining on chr_pos_ref/alt")
    # NOTE the slash. See the module docstring.
    sai["vep_key"] = (
        sai["chr"].astype(str) + "_" + sai["pos"].astype(str) + "_"
        + sai["ref"].astype(str) + "/" + sai["alt"].astype(str)
    )
    merged = sai.merge(vep[["vep_key"] + VEP_COLS], on="vep_key", how="left")
    n_missing = merged["Consequence"].isna().sum()
    print(f"  joined         : {len(merged):,} rows, {n_missing} with no VEP match")

    print(f"\n[3/5] Filtering to Consequence containing '{CONSEQUENCE_KEYWORD}'")
    cand = merged[
        merged["Consequence"].fillna("").str.contains(CONSEQUENCE_KEYWORD)
    ].copy()
    print(f"  candidates     : {len(cand):,} rows")
    print("  IMPACT         :", cand["IMPACT"].value_counts().to_dict())

    print("\n[4/5] Adding flags and priority")
    leaf_tissues = leafcutter_tissues(perind_dir)
    in_scan = scanned_genes(scan_csv)
    print(f"  LeafCutter tissues : {sorted(leaf_tissues)}")
    print(f"  genes already scanned : {len(in_scan)}")

    out = pd.DataFrame({
        "Gene":        cand["gene"],
        "SpliceAI":    cand["score"],
        "VAF":         cand["SR_VAF"],
        "Donor":       cand["SampleDonors"],
        "TissueID":    cand["TissueID"],
        "Tissue":      cand["SampleTissues"],
        "chr":         cand["chr"],
        "pos":         cand["pos"],
        "ref":         cand["ref"],
        "alt":         cand["alt"],
        "Consequence": cand["Consequence"],
        "Impact":      cand["IMPACT"],
        "Exon":        cand["EXON"],
        "Intron":      cand["INTRON"],
    })

    out["Already_case"] = out["Gene"].map(
        lambda g: "* existing" if g in ALREADY_CASE_GENES else pd.NA
    )
    # TISSUE-level flag - see the limitation note in the module docstring.
    out["has_RNAseq"] = out["TissueID"].map(
        lambda t: "OK" if t in leaf_tissues else "NO"
    )
    out["in_scan"] = out["Gene"].map(lambda g: "OK" if g in in_scan else pd.NA)

    # Priority: higher = better. Verified to reproduce the original column exactly.
    out["priority"] = (
        1
        + 2 * (out["has_RNAseq"] == "OK")       # tissue is testable
        + 2 * (out["Impact"] == "HIGH")         # VEP says high impact
        + 3 * (out["SpliceAI"] >= 0.8)          # SpliceAI is confident
        - 1 * (out["in_scan"] == "OK")          # already covered by Phase 2
    ).astype(int)

    out = out.sort_values(["priority", "SpliceAI"], ascending=[False, False])

    print("\n[5/5] Writing output")
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)
    out.to_csv(str(out_csv).replace(".csv", ".tsv"), sep="\t", index=False)
    print(f"  -> {out_csv}  ({len(out)} rows)")
    print("\n  priority distribution:")
    print(out["priority"].value_counts().sort_index(ascending=False).to_string())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spliceai", required=True,
                    help="V2_spliceAI_VAF_quantile_fraction.txt")
    ap.add_argument("--vep", required=True,
                    help="annotate_v2_381files.txt (Ensembl VEP output)")
    ap.add_argument("--perind", required=True,
                    help="directory holding <tissue>_perind.counts.gz")
    ap.add_argument("--scan", default=None,
                    help="Phase 2 02_outlier_scan_all.csv (optional, for in_scan)")
    ap.add_argument("--out", required=True,
                    help="output CSV path")
    args = ap.parse_args()
    build(args.spliceai, args.vep, args.perind, args.scan, args.out)


if __name__ == "__main__":
    main()
