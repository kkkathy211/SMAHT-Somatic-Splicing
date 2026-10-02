#!/usr/bin/env python3
# ============================================================================
# Phase 1 - Cohort differential splicing (group test)
#
# Map each LeafCutter cluster to a gene by majority vote over its junctions.
#
# Original location in the analysis project:
#     0520_analysis/code/04_cluster_to_gene_mapping.py
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
"""
Strategy A Pre-step — Map LeafCutter clusters to GENCODE gene_symbols.

For each tissue's perind.counts.gz, every junction in the file is a row.
Junction IDs look like:  chr10:86683311:86683509:clu_5278_-

We:
  1. Parse the header row to get all junction IDs
  2. Parse junction → (chrom, start, end, strand, cluster_id)
  3. Overlap with GENCODE protein-coding genes
  4. Per cluster: take majority overlapping gene (if any).
  5. Cluster with no protein-coding overlap → marked INTERGENIC.

Inputs:
   ../v2_batch_corrected/NEW/data/leafcutter_results/<tissue>/<tissue>_perind.counts.gz
   data/gencode.v44.basic.gtf.gz

Output:
   output/tables/A_cluster_to_gene.csv
   output/logs/04_cluster_to_gene.log
"""

import gzip, re, sys, time
from pathlib import Path
from collections import Counter
import pandas as pd
import pyranges as pr

ROOT  = Path(__file__).resolve().parent.parent
LEAF  = ROOT.parent / "v2_batch_corrected" / "NEW" / "data" / "leafcutter_results"
GTF   = ROOT / "data" / "gencode.v44.basic.gtf.gz"
OUT   = ROOT / "output" / "tables" / "A_cluster_to_gene.csv"
LOG   = ROOT / "output" / "logs" / "04_cluster_to_gene.log"

TISSUES = [
    "3AD_Skin_Calf", "3AF_Skin_Abdomen", "3AH_Muscle",
    "3AK_Brain_Frontal_Lobe", "3AL_Brain_Temporal_Lobe",
    "3AM_Brain_Cerebellum", "3C_Esophagus", "3G_Colon_Desc",
    "3Q_Lung", "3S_Heart", "3U_Testis_L",
]

# Junction ID regex:  chrXX:start:end:clu_NNN_[+-]
JUNC_RE = re.compile(r"^(chr[^:]+):(\d+):(\d+):(clu_\d+_[+-])$")


def parse_gtf_gene_intervals(gtf_path: Path):
    rows = []
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            attrs = parts[8]
            gname = re.search(r'gene_name "([^"]+)"', attrs)
            gtype = re.search(r'gene_type "([^"]+)"', attrs)
            gid = re.search(r'gene_id "([^"]+)"', attrs)
            hgnc = 'hgnc_id "' in attrs
            is_rt = 'tag "readthrough_gene"' in attrs
            if not gname:
                continue
            gene_name = gname.group(1)
            gene_type = gtype.group(1) if gtype else ""
            gene_id = gid.group(1) if gid else ""
            has_real = not gene_name.startswith("ENSG")
            rows.append((chrom, start, end, gene_name, gene_type,
                         gene_id, hgnc, has_real, is_rt))
    return pd.DataFrame(rows, columns=[
        "Chromosome", "Start", "End", "gene_name", "gene_type",
        "gene_id", "has_hgnc", "has_real_name", "is_readthrough"])


def extract_junctions(perind_path: Path) -> pd.DataFrame:
    """Read just the header → get all junction IDs."""
    with gzip.open(perind_path, "rt") as fh:
        # First line in perind has format:  chrom SAMPLE1 SAMPLE2 ...
        # then each data row starts with junction_id
        next(fh)  # skip header
        rows = []
        for line in fh:
            jid = line.split(" ", 1)[0]
            m = JUNC_RE.match(jid)
            if not m:
                continue
            chrom, s, e, clu = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
            rows.append((jid, chrom, s, e, clu))
    return pd.DataFrame(rows,
                        columns=["jid", "Chromosome", "Start", "End",
                                 "cluster_id"])


def main():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    log_fh = open(LOG, "w")
    def L(m):
        print(m); log_fh.write(m + "\n"); log_fh.flush()

    t0 = time.time()
    L("== Strategy A Pre-step — LeafCutter cluster → gene mapping ==\n")

    L("[1/3] Load GENCODE gene intervals ...")
    gtf_df = parse_gtf_gene_intervals(GTF)
    L(f"   Total genes in GTF: {len(gtf_df):,}")
    L(f"   protein_coding only:  "
      f"{(gtf_df['gene_type']=='protein_coding').sum():,}")

    # Composite priority: prefer real HGNC-curated protein-coding genes
    gtf_df["rank"] = (
        (~gtf_df["has_hgnc"]).astype(int) * 1000 +
        (~gtf_df["has_real_name"]).astype(int) * 100 +
        gtf_df["is_readthrough"].astype(int) * 10 +
        (gtf_df["gene_type"] != "protein_coding").astype(int)
    )
    gtf_pr = pr.PyRanges(gtf_df[["Chromosome", "Start", "End", "gene_name",
                                  "gene_type", "gene_id", "rank"]])

    all_rows = []
    for tname in TISSUES:
        L(f"\n[2/3] Process {tname} ...")
        perind = LEAF / tname / f"{tname}_perind.counts.gz"
        if not perind.exists():
            L(f"   ! missing: {perind}"); continue

        junc_df = extract_junctions(perind)
        L(f"   Junctions: {len(junc_df):,}")
        L(f"   Clusters:  {junc_df['cluster_id'].nunique():,}")

        # Overlap junctions with gene intervals
        junc_pr = pr.PyRanges(junc_df)
        ov = junc_pr.join(gtf_pr, how="left").df
        # Sort to prefer top-ranked gene per junction
        ov = ov.sort_values(["jid", "rank"])
        # Recast bool cols not needed (just rank int)
        ov["rank"] = ov["rank"].fillna(99999).astype(int)
        # Per junction: pick best (lowest rank) gene
        ov_best = ov.drop_duplicates(subset=["jid"], keep="first")
        ov_best = ov_best.assign(gene_name=ov_best["gene_name"].fillna("INTERGENIC"),
                                  gene_type=ov_best["gene_type"].fillna(""))

        # Per cluster: pick majority gene across its junctions
        clu_gene = []
        for clu, sub in ov_best.groupby("cluster_id"):
            counts = Counter(sub["gene_name"])
            # Drop INTERGENIC from counter if a real gene also tagged
            if "INTERGENIC" in counts and len(counts) > 1:
                del counts["INTERGENIC"]
            majority = counts.most_common(1)[0][0]
            gene_row = sub[sub["gene_name"] == majority].iloc[0]
            clu_gene.append({
                "tissue": tname,
                "cluster_id": clu,
                "n_junctions": len(sub),
                "gene_name": majority,
                "gene_type": gene_row["gene_type"],
                "gene_id": gene_row["gene_id"] if pd.notna(gene_row["gene_id"]) else "",
            })
        L(f"   Clusters mapped: {len(clu_gene):,}")
        all_rows.extend(clu_gene)

    L(f"\n[3/3] Write merged table ...")
    final = pd.DataFrame(all_rows)
    final = final.sort_values(["tissue", "gene_name", "cluster_id"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUT, index=False)
    L(f"   Wrote {OUT}")
    L(f"   Total rows: {len(final):,}")

    # Summary per tissue
    L("\nClusters per tissue:")
    for tissue, sub in final.groupby("tissue"):
        n_pc = (sub["gene_type"] == "protein_coding").sum()
        n_inter = (sub["gene_name"] == "INTERGENIC").sum()
        L(f"   {tissue:30s}  total={len(sub):5d}  "
          f"protein_coding={n_pc:5d}  intergenic={n_inter:5d}")

    L(f"\nTotal runtime: {(time.time()-t0)/60:.1f} min")
    L("DONE.")
    log_fh.close()


if __name__ == "__main__":
    main()
