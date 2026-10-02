# Somatic Splicing Analysis 

Code for testing whether **somatic mutations alter RNA splicing** in SMaHT P25 cohort.

---

## The question

Somatic mutations accumulate in normal tissue throughout life. Some of them land
on or near splice sites, and computational predictors (SpliceAI, VEP) flag them
as splice-disrupting. **Does that prediction show up in the matched RNA?**

The cohort:

| | |
|---|---|
| Donors | 23 (somatic calls), 25 (germline) |
| Tissues with DNA | 21 |
| Tissues with matched RNA-seq | 11 (10 after re-clustering QC) |
| RNA samples ("cores") | 167 |
| Sequencing centers | 4 (broad, nygc, bcm, washu) |
| Somatic VCFs | 351 |
| **Median somatic VAF** | **0.02** |

That last number is the one that matters. A variant present in 2% of cells can
shift bulk intron-usage by at most ~2% — at or below the noise floor of
short-read bulk RNA-seq. Every analysis in this repo eventually runs into it.

---

## Three framings

The project changed its statistical framing twice. Each phase is a self-contained
pipeline in its own directory.

| Phase | Framing | Test type | Statistical unit | Key output |
|---|---|---|---|---|
| **1** | Carrier vs non-carrier, per gene × tissue | **Group** | junction (lm), cluster (DM) | 10 hits @ FDR<10% |
| **2** | One carrier vs the cohort distribution | **Individual** | donor | 8 outliers → 1 after strict filter |
| **3** | Prediction-first; RNA only confirms | Individual | — (visual) | 213 candidates → 77 STRICT |

```
Phase 1  group test  ──┐  most variants have only ONE carrier
                       │  per gene-tissue → no statistical power
                       ▼
Phase 2  individual outlier  ──┐  max VAF in the whole cohort is 0.058;
                               │  only 1 candidate survives artifact filtering
                               ▼
Phase 3  prediction-first  ──→  rank by SpliceAI + VEP agreement,
                                use RNA as confirmation, not discovery
```

---

## Repository layout

```
.
├── README.md
├── requirements.txt
├── .gitignore
│
├── docs/
│   ├── methods.md            statistical detail for all three phases
│   ├── data_dictionary.md    every input file format, column by column
│   └── results.md            findings, case studies, known artifacts
│
├── phase1_cohort_differential_splicing/     group test (carrier vs non-carrier)
│   ├── README.md
│   ├── 01_annotate_somatic_variants.py
│   ├── 02_aggregate_carriers_per_gene.py
│   ├── 03_filter_candidate_gene_tissue_pairs.py
│   ├── 04_map_leafcutter_clusters_to_genes.py
│   ├── 05_join_spliceai_annotations.py
│   ├── 06_test_carrier_vs_noncarrier_lm.R
│   ├── 07_detect_isoform_switch_clusters.py
│   ├── 08_test_cross_tissue_concordance.py
│   └── 09_plot_sashimi_leafcutter_ds.R
│
├── phase2_per_donor_outlier_scan/           individual test (1 vs cohort)
│   ├── README.md
│   ├── utils.py                             shared engine: perind parsing, psi, MAD z
│   ├── 01_inventory_testable_triplets.py
│   ├── 02_scan_psi_outliers.py              <- the core analysis
│   ├── 03_rank_candidates_by_visibility.py
│   ├── 04_apply_strict_artifact_filter.py
│   ├── 05_plot_case_study_arcs.py
│   ├── 06_plot_vaf_spliceai_landscape.py
│   ├── 07_plot_psi_stripplot.py
│   ├── 08_plot_pileup_cummings_style.py
│   └── 09_fraser2_cross_validation.R
│
└── phase3_prediction_first_candidates/      prediction-first candidate selection
    ├── README.md
    ├── 01_build_vep_candidate_list.py       <- entry point (reconstructed)
    ├── 02_enrich_candidates_gencode_v47.py
    ├── 03_unsupervised_biggene_scan.py
    ├── 04_plot_candidates_by_tissue.py
    ├── 05_plot_pileup_pyqtl.py
    ├── 06_plot_pileup_standalone.py
    ├── 07_plot_per_donor_sashimi.py
    ├── 08_plot_vep_schematic.py
    └── qc/
        ├── check_housekeeping_coverage.sh
        └── check_intron_exon_ratio.sh
```

---

## Inputs (not in this repo)

All inputs are controlled-access SMaHT data and are **not** committed here.
See [`docs/data_dictionary.md`](docs/data_dictionary.md) for the exact format of
each file.

| File | Source | Used by |
|---|---|---|
| `vcf_p25/*.vcf.gz` (351) | SMaHT DAC somatic SNV release | Phase 1 |
| `<tissue>_perind.counts.gz` | LeafCutter, run by Laura| Phases 1, 2, 3 |
| `<tissue>.leafcutter.PCs.txt` | same | Phase 1 |
| `V2_spliceAI_VAF_quantile_fraction.txt` | V2 somatic SpliceAI | Phases 2, 3 |
| `annotate_v2_381files.txt` | Ensembl VEP v113 | Phase 3 |
| `gencode.v47.annotation.gtf` | GENCODE | Phase 3 |
| `gs://smaht-p25-bam-files/*.bam` | SMaHT RNA-seq alignments | Phases 2, 3 (pileups) |

---

## Dependencies

**Python ≥3.9**

```
pandas  numpy  scipy  matplotlib  seaborn  pysam  pyranges
```

**R ≥4.3**

```
data.table  dplyr  ggplot2  patchwork  matrixStats  broom  lme4  lmerTest
leafcutter  FRASER
```

**Command-line**

```
samtools  regtools  gsutil
```

**Optional** — `pyqtl` (Broad internal) for `phase3/05_plot_pileup_pyqtl.py`.
A dependency-free equivalent is provided as `06_plot_pileup_standalone.py`.

See [`requirements.txt`](requirements.txt).

---

## Quick start

Each phase runs independently. Scripts are numbered in execution order.

```bash
# Phase 1 — group test
cd phase1_cohort_differential_splicing
python 01_annotate_somatic_variants.py
python 02_aggregate_carriers_per_gene.py
python 03_filter_candidate_gene_tissue_pairs.py
python 04_map_leafcutter_clusters_to_genes.py
Rscript 06_test_carrier_vs_noncarrier_lm.R

# Phase 2 — per-donor outlier scan
cd ../phase2_per_donor_outlier_scan
python 01_inventory_testable_triplets.py
python 02_scan_psi_outliers.py
python 04_apply_strict_artifact_filter.py

# Phase 3 — prediction-first
cd ../phase3_prediction_first_candidates
python 01_build_vep_candidate_list.py --spliceai ... --vep ... --perind ... --out ...
python 02_enrich_candidates_gencode_v47.py
```

Most scripts still carry hardcoded paths from the original analysis environment
(either a local Mac or the Broad cluster). Each phase README lists the paths you
need to edit.

---

## Results in one table

| Phase | Funnel | Outcome |
|---|---|---|
| **1** | 34,907 HC mutations → 182 candidate gene×tissue pairs → 62 testable → 1,458 junction tests | **10 pairs at FDR<10%**, 6 at <5%, 5 isoform-switch clusters, **0 replicating across tissues** |
| **2** | 180,162 variants → 717 at SpliceAI≥0.2 → 71 testable triplets | **8 pass 3 filters** → **1 survives strict artifact filtering (PLBD2)** |
| **3** | 717 → 213 VEP-splice candidates | **77 STRICT** (SpliceAI≥0.8 + VEP HIGH) → 17 selected → **3 Tier A with BAMs available** |

The headline finding from Phase 1 is **DMD in esophagus**: significant under both
per-junction `lm` (q = 0.0013) and LeafCutter's Dirichlet-Multinomial cluster
test (q = 0.016), and the effect *doubles* (Δψ −0.11 → −0.24) once a carrier with
a confounding germline splice variant is excluded.

Full detail, including the negative results and the artifacts we identified, is in
[`docs/results.md`](docs/results.md).

---

## Known issues

These are documented rather than silently fixed, because the published figures
were generated with the code as-is.

1. **The Phase 3 entry point was missing.** `vep_splice_candidates.csv` had no
   producing script — it was built in an interactive session that was never saved.
   `01_build_vep_candidate_list.py` is a reconstruction, verified to reproduce the
   original 213 rows and the `priority` column exactly.

2. **Gene-symbol version skew.** VEP was run against cache `105_GRCh38` while the
   GTF is GENCODE v47. 20 of 213 candidates (12 symbols) get no v47 annotation —
   e.g. `FAM49A` was renamed `CYRIA`, which silently dropped it from one batch.

3. **`has_RNAseq` is a tissue-level flag,** not per-(donor, tissue). A row can be
   marked testable and still have no BAM for that donor in that tissue. Check the
   bucket directly for real testability.

4. **LeafCutter cluster IDs are not stable across runs.** Match junctions on
   `chr:start:end`, never on `clu_N`. `phase2/07_plot_psi_stripplot.py` does this
   deliberately.

5. **Intergenic sentinel.** `pyranges.join(how="left")` fills unmatched rows with
   `-1`, not `NaN`, so the `fillna("INTERGENIC")` calls in Phase 1 are no-ops.
   14,969 mutations carry `gene_name == "-1"`. Downstream filters on
   `gene_type == "protein_coding"` make this harmless, but it is not cosmetic if
   you reuse those tables.

6. **FDR scope is not comparable between phases.** Phase 1 corrects within each
   (tissue, gene, test) — typically ~20 junctions. An earlier cohort-wide version
   corrected across ~30,000 junctions per tissue. Hit counts from the two are not
   on the same scale.

---

## Where this is going

Short-read bulk RNA-seq has a hard ceiling at SMaHT-level VAF. Two directions
follow from that:

- **Long-read (PacBio Kinnex)** — full-length isoforms resolve structures that
  junction-level inference can only hint at, and reads are long enough to phase a
  variant onto the transcript carrying the aberrant splice.
- **Unsupervised outlier detection (FRASER2 / DROP)** — per-sample normalized, so
  it sidesteps the library-size and batch artifacts that dominate the RPM
  comparisons here, and it does not depend on a SpliceAI candidate list.
  A cross-validation harness is already in
  `phase2/09_fraser2_cross_validation.R`.
