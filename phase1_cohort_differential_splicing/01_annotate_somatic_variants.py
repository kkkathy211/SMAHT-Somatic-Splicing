#!/usr/bin/env python3
# ============================================================================
# Phase 1 - Cohort differential splicing (group test)
#
# Parse 351 somatic VCFs (HighConf only), overlap with GENCODE v44 gene bodies.
#
# Original location in the analysis project:
#     0520_analysis/code/01_extract_and_annotate.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
Phase 0.1 — Extract HC mutations from all VCFs + annotate with gene_symbol.

Inputs:
    ../../somatic_snvIndel_p25/vcf_p25/*.vcf.gz   (351 VCFs)
    data/gencode.v44.basic.gtf.gz                  (GENCODE gene annotation)

Outputs:
    output/tables/00_annotated_mutations.csv      (one row per HC mutation,
                                                    with donor / tissue / gene)
    output/logs/01_extract_annotate.log

Notes:
    - HC = FILTER==HighConf  (the official HighConf tier per VCF header)
    - Gene assignment = first overlapping protein-coding gene from GENCODE basic
    - "Intergenic" variants are kept with gene = None
"""

import gzip, os, re, sys, time
from pathlib import Path
import pandas as pd
import pyranges as pr

# -----------------------------------------------------------------------------
# Paths
ROOT       = Path(__file__).resolve().parent.parent  # 0520_analysis/
VCF_DIR    = ROOT.parent / "somatic_snvIndel_p25" / "vcf_p25"
GTF        = ROOT / "data" / "gencode.v44.basic.gtf.gz"
OUT_TABLE  = ROOT / "output" / "tables" / "00_annotated_mutations.csv"
LOG        = ROOT / "output" / "logs" / "01_extract_annotate.log"

# Sample-ID parsing pattern  (SMHTXXX-TIS-CORE-AgeSex-...)
# Example: SMHT020-3AH-MAMC-M81-XX-dac-SMAFIKQATZJV
SAMPLE_RE = re.compile(r"^(SMHT\d+)-(\w+?)-")

# Map tissue codes to full names (matching our covariate matrix convention)
TISSUE_MAP = {
    "3AD": "3AD_Skin_Calf",
    "3AF": "3AF_Skin_Abdomen",
    "3AH": "3AH_Muscle",
    "3AK": "3AK_Brain_Frontal_Lobe",
    "3AL": "3AL_Brain_Temporal_Lobe",
    "3AM": "3AM_Brain_Cerebellum",
    "3C":  "3C_Esophagus",
    "3G":  "3G_Colon_Desc",
    "3Q":  "3Q_Lung",
    "3S":  "3S_Heart",
    "3U":  "3U_Testis_L",
}

KEEP_FILTERS = {"HighConf"}   # 只保留 HC tier
INFO_FIELDS  = {"SR_VAF"}     # 只取需要的 INFO 字段


# -----------------------------------------------------------------------------
def parse_info(info_str: str) -> dict:
    """Parse VCF INFO column → dict of key/value."""
    out = {}
    for kv in info_str.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            if k in INFO_FIELDS:
                out[k] = v
        else:
            out[kv] = True   # flag
    return out


def extract_vcf(vcf_path: Path):
    """Yield (chrom, pos, ref, alt, vaf, filter_str) for variants in a VCF."""
    fname = vcf_path.name
    m = SAMPLE_RE.match(fname)
    if not m:
        return
    donor = m.group(1)
    tcode = m.group(2)
    tissue = TISSUE_MAP.get(tcode, None)
    if tissue is None:
        return  # 不在我们 11 个 tissue 列表里

    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt, _qual, filt, info = parts[:8]
            if filt not in KEEP_FILTERS:
                continue
            info_d = parse_info(info)
            vaf = float(info_d.get("SR_VAF", "nan"))
            yield {
                "donor": donor,
                "tissue": tissue,
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "vaf": vaf,
                "filter": filt,
            }


def main():
    log_fh = open(LOG, "w")
    def log(msg):
        print(msg)
        log_fh.write(msg + "\n")
        log_fh.flush()

    t0 = time.time()
    log("== Phase 0.1 — Extract HC mutations + annotate with gene_symbol ==")
    log(f"GTF: {GTF}")
    log(f"VCF dir: {VCF_DIR}")
    log(f"Tissue map: {len(TISSUE_MAP)} tissues")

    # -------------------------------------------------------------------------
    # 1) Read GTF: keep only protein_coding gene rows, build PyRanges
    # -------------------------------------------------------------------------
    log("\n[1/3] Loading GENCODE GTF ...")
    rows = []
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1   # GTF 1-based inclusive → 0-based half-open
            end   = int(parts[4])
            # parse attributes
            attrs = parts[8]
            gname = re.search(r'gene_name "([^"]+)"', attrs)
            gtype = re.search(r'gene_type "([^"]+)"', attrs)
            gid   = re.search(r'gene_id "([^"]+)"', attrs)
            hgnc  = re.search(r'hgnc_id "([^"]+)"', attrs)
            is_readthrough = 'tag "readthrough_gene"' in attrs
            if not gname:
                continue
            gene_name = gname.group(1)
            gene_type = gtype.group(1) if gtype else ""
            gene_id   = gid.group(1) if gid else ""
            has_hgnc  = hgnc is not None
            # has_real_name = gene_name doesn't look like "ENSGxxxxxxx"
            has_real_name = not gene_name.startswith("ENSG")
            rows.append((chrom, start, end, gene_name, gene_type, gene_id,
                          has_hgnc, has_real_name, is_readthrough))

    gtf_df = pd.DataFrame(rows,
                           columns=["Chromosome", "Start", "End",
                                    "gene_name", "gene_type", "gene_id",
                                    "has_hgnc", "has_real_name", "is_readthrough"])
    log(f"   Total gene records: {len(gtf_df):,}")
    # Keep all for now (so we can also identify lincRNA / SF genes downstream)
    gtf_pr = pr.PyRanges(gtf_df)
    log(f"   Built PyRanges with {len(gtf_pr):,} intervals")

    # -------------------------------------------------------------------------
    # 2) Loop over all VCFs, accumulate HC mutations
    # -------------------------------------------------------------------------
    log("\n[2/3] Parsing VCFs ...")
    vcf_files = sorted(VCF_DIR.glob("*.vcf.gz"))
    log(f"   Found {len(vcf_files)} VCFs")

    all_mut_rows = []
    n_skipped = 0
    for i, vcf in enumerate(vcf_files):
        for rec in extract_vcf(vcf):
            all_mut_rows.append(rec)
        if (i + 1) % 50 == 0:
            log(f"   parsed {i+1}/{len(vcf_files)} VCFs, "
                f"so far {len(all_mut_rows):,} HC mutations")
    log(f"   Done. Total HC mutations across all VCFs: {len(all_mut_rows):,}")

    mut_df = pd.DataFrame(all_mut_rows)
    log(f"   Unique (donor,tissue) pairs:  "
        f"{mut_df[['donor','tissue']].drop_duplicates().shape[0]}")
    log(f"   Unique donors:               {mut_df['donor'].nunique()}")
    log(f"   Unique tissues:              {mut_df['tissue'].nunique()}")
    log(f"   Tissue counts:")
    for t, c in mut_df['tissue'].value_counts().items():
        log(f"      {t:32s}  {c:6,}")

    # -------------------------------------------------------------------------
    # 3) Overlap mutations with gene intervals (PyRanges)
    # -------------------------------------------------------------------------
    log("\n[3/3] Annotating mutations with gene_symbol ...")
    mut_pr = pr.PyRanges(
        mut_df.rename(columns={"chrom": "Chromosome", "pos": "Start"})
              .assign(End=lambda d: d["Start"] + 1)
    )

    # Overlap: each mutation gets attached to any overlapping gene
    overlap = mut_pr.join(gtf_pr, how="left")
    over_df = overlap.df

    # Reduce: if a mutation hits multiple genes, prefer protein_coding,
    # else first encountered.
    log(f"   Raw overlap rows (before dedup): {len(over_df):,}")
    if "Start_b" in over_df.columns:
        over_df = over_df.drop(columns=["Start_b", "End_b"], errors="ignore")

    # Re-cast bool cols (pyranges may pass them as strings)
    for c in ("has_hgnc", "has_real_name", "is_readthrough"):
        over_df[c] = over_df[c].map(lambda v: bool(v) if isinstance(v, bool)
                                              else str(v).lower() == "true")
    # Composite ranking (lower = better priority):
    #   1. has HGNC ID  (real curated gene)
    #   2. has real gene_name (not "ENSGxxx")
    #   3. NOT a readthrough_gene
    #   4. protein_coding
    over_df["rank"] = (
        (~over_df["has_hgnc"]).astype(int) * 1000   # missing HGNC = +1000 (worst)
        + (~over_df["has_real_name"]).astype(int) * 100  # ENSG-style name = +100
        + over_df["is_readthrough"].astype(int) * 10    # readthrough = +10
        + (over_df["gene_type"] != "protein_coding").astype(int)
    )
    over_df = (over_df.sort_values(["donor", "tissue", "Chromosome",
                                      "Start", "rank"])
                       .drop_duplicates(subset=["donor", "tissue",
                                                  "Chromosome", "Start",
                                                  "ref", "alt"],
                                          keep="first"))
    over_df = over_df.drop(columns=["rank", "has_hgnc", "has_real_name",
                                      "is_readthrough"], errors="ignore")

    # Tidy column names back to lowercase
    final = over_df.rename(columns={"Chromosome": "chrom", "Start": "pos"}) \
                    .drop(columns=["End"], errors="ignore")
    final["gene_name"] = final["gene_name"].fillna("INTERGENIC")
    final["gene_type"] = final["gene_type"].fillna("")
    final["gene_id"]   = final["gene_id"].fillna("")

    # Final column order
    cols = ["donor", "tissue", "chrom", "pos", "ref", "alt", "vaf",
             "filter", "gene_name", "gene_type", "gene_id"]
    final = final[cols].sort_values(["donor", "tissue", "chrom", "pos"])

    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUT_TABLE, index=False)
    log(f"\n   Wrote {OUT_TABLE} ({len(final):,} rows)")
    log(f"\n   Gene summary:")
    log(f"      INTERGENIC:               "
        f"{(final['gene_name']=='INTERGENIC').sum():,}")
    log(f"      protein_coding annotated: "
        f"{(final['gene_type']=='protein_coding').sum():,}")
    log(f"      Unique annotated genes:   "
        f"{final[final['gene_name']!='INTERGENIC']['gene_name'].nunique():,}")

    log(f"\nTotal runtime: {(time.time()-t0)/60:.1f} min")
    log("DONE.")
    log_fh.close()


if __name__ == "__main__":
    main()
