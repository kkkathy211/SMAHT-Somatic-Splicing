# Phase 1 — Cohort differential splicing (group test)

## Research question

> **For gene *G* in tissue *T*, do donors carrying a somatic mutation in *G*
> splice *G* differently from donors who don't?**

Two complementary arms:

- **Strategy A** — within-tissue. Per-junction regression of intron usage on
  carrier status.
- **Strategy B** — cross-tissue. If a gene is mutated in several tissues, does
  the splicing change go the same direction each time?

## Test type

**Group test.** Carriers (n ≥ 4) vs non-carriers (n ≥ 4) within one tissue.
This is a conventional case/control design — and it is exactly why the phase
eventually stalls (see *Why this phase ended*).

---

## Pipeline

```
351 somatic VCFs                 GENCODE v44 GTF
        │                               │
        └──────────┬────────────────────┘
                   ▼
    01_annotate_somatic_variants.py          [pyranges interval overlap]
                   ▼
        00_annotated_mutations.csv  (34,907 rows)
                   ▼
    02_aggregate_carriers_per_gene.py
                   ▼
    00_carriers_per_tissue_gene.csv  (15,644 rows)
                   ▼
    03_filter_candidate_gene_tissue_pairs.py    [>=4 carriers, >=4 non-carriers,
                   ▼                             protein_coding]
        00_candidates_A.csv  (182)   00_candidates_B.csv  (125)

  LeafCutter perind ──► 04_map_leafcutter_clusters_to_genes.py
                                    ▼
                        A_cluster_to_gene.csv  (516,762 rows)

  SpliceAI calls   ──► 05_join_spliceai_annotations.py
                                    ▼
                        02b_mutations_with_spliceAI.csv

        candidates + clusters + covariates + PSI
                   ▼
    06_test_carrier_vs_noncarrier_lm.R           <- the core engine
                   ▼
        A_summary.csv (62)   A_top_hits.csv (24)
                   ▼
    07_detect_isoform_switch_clusters.py         [mirror-pair detection]
                   ▼
        A_cluster_hits.csv (15)

    08_test_cross_tissue_concordance.py
                   ▼
        B_summary.csv (10)   B_junction_concordant.csv (0 rows)

    09_plot_sashimi_leafcutter_ds.R              [figure for the DMD vignette]
```

---

## File I/O

### Inputs

| File | Format | What it holds |
|---|---|---|
| `vcf_p25/*.vcf.gz` (351) | VCFv4.2, **sites-only** | Somatic SNVs. No FORMAT column — everything is in `FILTER` (`HighConf` / `LowConf` / `LikelyArtifact`) and `INFO`. Key fields: `SR_VAF` (this tissue's VAF), `TISSUE_SR_VAFS` (pipe-delimited VAF in the donor's *other* tissues), `CALLERS`, `PB_PHASING`. |
| `gencode.v44.basic.gtf.gz` | GTF | Gene body coordinates + HGNC IDs |
| `<tissue>_perind.counts.gz` | LeafCutter | Space-delimited. Row key `chr:start:end:clu_N_strand`, cells `numer/denom`. ψ = numer/denom. |
| `<tissue>.leafcutter.PCs.txt` | TSV | Splicing PC1–PC10 × samples |
| `covariate_matrix.csv` | CSV, 167 rows | `sample_id, donor, tissue, core, sex, age, protocol, center` parsed from sample IDs |
| `spliceAI_VAF_quantile_fraction.txt` | TSV | V1 SpliceAI somatic calls (only 33.5% coverage of HC mutations) |

### Outputs

| File | Rows | Contents |
|---|---|---|
| `00_annotated_mutations.csv` | 34,907 | `donor, tissue, chrom, pos, ref, alt, vaf, filter, gene_name, gene_type, gene_id` |
| `00_carriers_per_tissue_gene.csv` | 15,644 | `tissue, gene_name, n_carriers, n_non_carriers, median_VAF_carriers, max_spliceai, is_SF_gene` |
| `00_candidates_A.csv` | 182 | gene × tissue pairs passing the carrier-count filter |
| `A_cluster_to_gene.csv` | 516,762 | `tissue, cluster_id, n_junctions, gene_name, gene_type` |
| `A_summary.csv` | 62 | 30 cols — per candidate, the A.1/A.2/A.3 top hit, q, β, and the covariate set actually used |
| `A_top_hits.csv` | 24 | junction-level significant hits |
| `A_cluster_hits.csv` | 15 | cluster-level roll-up + `has_mirror_pair` + confidence stars |
| `B_summary.csv` | 10 | per gene: tissues tested, β signs, concordance label |

---

## Model

Three regressions per candidate gene × tissue, all on logit-transformed intron
usage:

```
A.1   logit(PSI) ~ carrier_flag + center + sex + age + PC1..PCk     # binary
A.2   logit(PSI) ~ VAF          + center + sex + age + PC1..PCk     # carriers only
A.3   logit(PSI) ~ VAF          + center + sex + age + PC1..PCk     # all donors, non-carriers = 0
```

A fourth test was added after a reviewer asked whether the signal is driven by
predicted-splice variants specifically:

```
A.4   logit(PSI) ~ donor_max_SpliceAI + center + sex + age + PC1..PCk
```

**Adaptive covariate selection.** Small tissues (n = 9–12) saturate the design
matrix. The engine drops terms in a fixed order — `PC5 → PC4 → PC3 → PC2 → PC1 →
age → sex → center` — until `df_resid ≥ 3`. Without this, four tissues returned
all-`NA` p-values that were silently reported as "zero hits".

**Multiple testing.** Benjamini–Hochberg **within each (tissue, gene, test)** —
typically ~20 junctions per gene. This is far more permissive than a genome-wide
correction; hit counts are not comparable to analyses that correct across all
junctions.

**Cluster-level test** (`09_plot_sashimi_leafcutter_ds.R` consumes its output).
LeafCutter's Dirichlet-Multinomial `leafcutter_ds.R`, run on the Broad cluster
as a qsub array. It would not install locally — Stan/rstan compile failure on
macOS — which is why per-junction `lm` is the primary engine throughout.

**Isoform-switch detection.** Within a LeafCutter cluster, intron usage sums to
1, so a genuine isoform switch produces a *mirror pair*: two junctions with
opposite-sign β of similar magnitude. `07_detect_isoform_switch_clusters.py`
flags pairs where both q < 0.10, β signs differ, and |β| ratio is in 0.5–2.0.

---

## Tools

| Tool | Where |
|---|---|
| `pyranges` | interval overlap of variants against gene bodies, and junctions against genes |
| R `lm` | per-junction regression (all three tests) |
| `leafcutter_ds.R` | cluster-level Dirichlet-Multinomial (Broad cluster) |
| `data.table` | PSI matrix handling |
| `ggplot2` + `patchwork` | sashimi figure |

---

## Carrier definitions

Three were run in parallel to test whether the signal is driven by
splice-predicted variants:

| Definition | Criterion |
|---|---|
| **ANY** | donor has any HighConf somatic mutation in the gene |
| **LOOSE** | donor has a mutation with SpliceAI ≥ 0.2 |
| **STRICT / Q95** | donor has a mutation with `spliceAI_type ∈ {Q95, Q99}` |

Donors with a sub-threshold mutation are **excluded** from the LOOSE/STRICT runs
rather than used as controls — they are neither clean controls nor convincing
carriers.

A fourth definition was added later:

| **ANY_noGerm** | ANY, minus any donor carrying a high-confidence *germline* SpliceAI variant in the same gene |

---

## Results

```
34,907  HighConf somatic mutations
19,938  inside an annotated gene        (13,852 protein-coding)
 4,805  unique mutated genes            (3,319 protein-coding)
15,644  (tissue, gene) pairs with >=1 carrier
   182  Candidate A pairs
    62  testable after RNA matching
 1,458  junction-level lm tests
    10  pairs at FDR < 10%
     6  pairs at FDR < 5%
     5  clusters with mirror-pair architecture
     0  genes concordant across tissues
```

**Top hits (A.3):**

| Gene | Tissue | q | β |
|---|---|---|---|
| DMD | 3C Esophagus | 0.0013 | −1.64 |
| RBFOX1 | 3AH Muscle | 0.0016 | −0.33 |
| DUSP22 | 3AK Brain frontal | 0.0051 | — |
| PTPRD | 3Q Lung | 0.013 | +4.46 |
| NRG1 | 3Q Lung | 0.05 | +7.21 |

**lm vs Dirichlet-Multinomial agreement** (64 candidates re-tested at cluster
level):

| Category | n | Genes |
|---|---|---|
| Both non-significant | 92 | — |
| lm only | 9 | RBFOX1/Muscle, DUSP22 ×3, PTPRD/Lung, … |
| DS only | 2 | MACROD2/Lung, KMT2C/Muscle |
| **Both significant** | **1** | **DMD / 3C Esophagus** |

**The DMD vignette.** Five carriers vs nine non-carriers. One carrier also
carried a germline DMD variant with SpliceAI > 0.5. Excluding them:

| Test | q | effect |
|---|---|---|
| lm, 5 carriers | 0.0013 | β = −1.64 |
| DM, 5 carriers | 0.016 | Δψ = −0.11 |
| **DM, 4 carriers (no germline)** | **0.0042** | **Δψ = −0.24** |

q improves ~4×, effect size doubles. The germline-suspect carrier was *diluting*
the signal, not driving it. Mean ψ on the long-range junction: 0.13 in
non-carriers → 0.37 in carriers.

**Strategy B is a clean negative.** Zero junctions reach q < 0.10 at the same
chromosomal position in ≥2 tissues. DMD flips sign between tissues
(esophagus −1.64, skin +2.17, lung +1.42) — direct empirical support for
analyzing tissues separately rather than pooling.

---

## Why this phase ended

Two independent walls:

1. **Carrier counts.** Most high-SpliceAI somatic variants have **exactly one
   carrier per gene × tissue**. One carrier against ten non-carriers has no power
   under any group test — not `lm`, not Dirichlet-Multinomial.

2. **SpliceAI coverage.** The V1 callset annotated only **33.5%** of HighConf
   mutations. Of 182 candidates, 173 had *zero* SpliceAI variation among their
   carriers. Test A.4 returned **zero hits** (best q = 0.12), and the
   SpliceAI-stratified carrier definitions produced **empty candidate lists**.
   Worse: all 10 Strategy-A hits had carriers whose maximum SpliceAI was **0.00**,
   and none of those 17 mutations sat within 50 bp of an annotated splice site.

Phase 2 abandons the group test entirely.

---

## Paths to edit

These scripts carry hardcoded paths from the original environment:

| Script | Hardcoded path |
|---|---|
| `01_annotate_somatic_variants.py` | `../somatic_snvIndel_p25/vcf_p25/`, `data/gencode.v44.basic.gtf.gz` |
| `04_map_leafcutter_clusters_to_genes.py` | `../v2_batch_corrected/NEW/data/leafcutter_results/` |
| `06_test_carrier_vs_noncarrier_lm.R` | `../v2_batch_corrected/NEW/output/tables/covariate_matrix.csv` |
| `09_plot_sashimi_leafcutter_ds.R` | `3C_Esophagus_perind.counts.gz`, `DMD.ANY_noGerm.txt` |

---

## Caveats

- **The PSI values are quantile-normalized, not covariate-residualized.** The
  `.leafcutter.bed.gz` files come out of `cluster_prepare_fastqtl.py`, which only
  applies `qqnorm`. Batch effects are still in the phenotype; that is what the
  explicit `center` term in every model is for.
- **`pyranges` fills unmatched rows with `-1`,** not `NaN`, so the
  `fillna("INTERGENIC")` calls are no-ops. 14,969 mutations carry
  `gene_name == "-1"`. Harmless downstream because of the `protein_coding`
  filter, but do not reuse those tables without handling it.
- **Splicing PC columns are keyed by `donor-tissue` and are non-unique** when a
  donor has two cores. The `merge` + `match` pattern gives both cores the *first*
  PC vector.
- **DUSP22 is probably an artifact** — present in all 11 tissues with 6–9 carriers
  each at VAF ≈ 0.25. Suspected callability or pseudogene-region issue. Needs IGV
  verification at chr6:304,600–345,900.
