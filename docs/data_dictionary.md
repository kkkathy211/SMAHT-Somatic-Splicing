# Data dictionary

Every input file this code reads, with its exact format. None of these are
committed to the repo — they are controlled-access SMaHT data.

---

## 1. Somatic variant calls

### `vcf_p25/*.vcf.gz` — 351 files

One VCF per (donor × tissue). Naming:

```
SMHT001-3A-MAMC-M42-XX-dac-SMAFIFT86X7Z-SmahtSNV_1.0.0_GRCh38.filtered.snv.vcf.gz
   │     │                     │
 donor  tissue            file accession
```

Every script parses donor/tissue by splitting on `-` and taking fields 1 and 2.

**These VCFs are sites-only** — there is no sample or FORMAT column. All
information lives in `FILTER` and `INFO`.

| Field | Values / meaning |
|---|---|
| `FILTER` | `HighConf` = CrossTech **or** (CrossTissue **and** CrossCaller)<br>`LowConf` = only one of CrossTissue / CrossCaller<br>`LikelyArtifact` = no cross-evidence tag |
| `SR_VAF` | Short-read VAF **in this VCF's tissue**. This is "the VAF" everywhere downstream. |
| **`TISSUE_SR_VAFS`** | Pipe-delimited map of VAF in the donor's *other* tissues:<br>`3G:0.066667\|3E:0.030303\|3A:0.203540`<br>**The most load-bearing field in the project** — it is how cross-tissue presence is reconstructed without re-reading every VCF. |
| `CALLERS` | RUFUS / Strelka2 / TNhaplotyper2 / longcallD |
| `TISSUE_PB_VAF`, `POOLED_PB_VAF`, `POOLED_ONT_VAF` | long-read VAFs |
| `SR_ADF`, `SR_ADR`, `PB_*`, `ONT_*` | strand-split depths |
| `SB_PVAL` | strand bias |
| `GERMLINE_PVAL*` | germline exclusion |
| `PB_PHASING` | e.g. `MOSAIC_PHASED` |

Variant counts vary widely per VCF: SMHT001-3A has 167 total (72 HighConf);
SMHT005-3AD has 2,785 total but only 70 HighConf.

**Coverage:** 23 donors × 21 tissue codes
(`3A, 3AD, 3AF, 3AH, 3AK, 3AL, 3AM, 3AN, 3AO, 3B, 3C, 3E, 3G, 3I, 3K, 3M, 3O, 3Q, 3S, 3U, 3Y`).
Only 11 of these have matched RNA.

---

## 2. SpliceAI predictions

### `V2_spliceAI_VAF_quantile_fraction.txt` — 180,162 rows, 18 columns

V2 somatic SpliceAI callset. One row per (variant, donor-tissue).

| Column | Meaning |
|---|---|
| `gene` | gene symbol |
| `score` | SpliceAI max delta score, 0–1 |
| `id` | `chr_pos_ref_alt` — **underscore separator** (see the join note below) |
| `chr`, `pos`, `ref`, `alt` | coordinates |
| `type` | VCF FILTER (`PASS` etc.) |
| `donor_tissue` | e.g. `SMHT029-3A` |
| `VAF_each_core` | per-core VAF, semicolon-delimited |
| `SR_VAF` | short-read VAF (the one used) |
| `CrossTissue`, `inTissues` | cross-tissue evidence |
| `SampleDonors` | donor ID |
| `TissueID` | short tissue code (`3AD`) |
| `SampleTissues` | full tissue name (`3AD - Skin, Calf`) |
| `vaf_bin`, `spliceAI_type` | quantile bins (`Q95`, `Q99`) |

**717 rows have `score >= 0.2`** (620 unique variants). That is the working set
for Phases 2 and 3.

### `smaht_germline_SNV_spliceai_type.txt.gz` — 37M rows, 161 MB

Germline SpliceAI, all 25 donors. Columns:
`id, chr, pos, ref, alt, FileAccession, SampleDonors, spliceAI_type, score, gene`

Used once, in Phase 1, to define the `ANY_noGerm` carrier tier for the DMD
vignette.

---

## 3. VEP annotation

### `annotate_v2_381files.txt` — 998,650 rows, 79 columns, 289 MB

Ensembl **VEP v113.0**, run 2026-06-15, cache **`105_GRCh38`**.

- 103 lines of `##` metadata, then the header at line 104
- 998,546 unique variants in 998,649 data rows → effectively **one row per
  variant** (VEP was run with `--pick`-style canonical selection)
- 671 of the 717 SpliceAI ≥ 0.2 variants carry a `MANE_SELECT` transcript

Key columns:

| Column | Meaning |
|---|---|
| `#Uploaded_variation` | `chr12_108565432_G/A` — **slash between alleles** |
| `Consequence` | comma-separated SO terms, e.g. `splice_donor_variant` |
| `IMPACT` | `HIGH` / `MODERATE` / `LOW` / `MODIFIER` |
| `SYMBOL`, `HGNC_ID`, `BIOTYPE` | gene identity |
| `CANONICAL`, `MANE`, `MANE_SELECT` | transcript selection |
| `EXON`, `INTRON` | e.g. `8/14` |
| `SIFT`, `PolyPhen` | missense predictors |
| ~30 `gnomAD*_AF`, `MAX_AF` | population frequency |
| `CLIN_SIG`, `SOMATIC`, `PHENO` | clinical annotation |

### ⚠️ The join

The SpliceAI `id` column uses an **underscore** (`chr12_108565432_G_A`); VEP uses
a **slash** (`chr12_108565432_G/A`). Joining on `id` returns **zero** matches.
Build the key explicitly:

```python
key = f"{chr}_{pos}_{ref}/{alt}"
```

With the key built this way the join has **zero misses** (717/717).

### ⚠️ Version skew

VEP cache is `105_GRCh38` but the GTF used downstream is **v47**. 20 of 213
candidates (12 symbols) therefore get no v47 annotation. See the Phase 3 README.

---

## 4. LeafCutter output

Produced by Laura Domenech. Cromwell command:

```
python3 cluster_prepare_fastqtl.py <junc_list> gencode.v47.GRCh38.exons.txt.gz \
    gencode.v47.genes.gtf <TISSUE> sample_participant.tsv \
    --min_clu_reads 30 --max_intron_len 500000 --num_pcs 10
```

Junctions from `regtools`; alignments STAR 2.7.10b against GENCODE v47.

### `<tissue>_perind.counts.gz` — the main phenotype

Space-delimited. First column is the junction key, remaining columns are samples.

```
chrom                          SMHT001-3AD-001A3-M42-A101-broad.regtools_junc.txt.gz  ...
chr1:16765:16854:clu_1_-       17/106                                                 ...
```

- Row key: `chr:start:end:clu_N_strand`
- Cells: `numerator/denominator` = junction reads / cluster total reads
- **ψ = numerator / denominator**, and within a cluster ψ sums to 1

**Two versions exist and they are not interchangeable:**

| Version | Location | 3AD clusters |
|---|---|---|
| v1 (raw) | `leafcutter_results/` | 163,106 |
| **v2 (re-clustered after mappability + annotated-junction QC)** | used from June onward | **68,121** (42% retention) |

Phases 2 and 3 use v2. `3AL_Brain_Temporal_Lobe` did not come back from the
re-cluster, which is why those phases cover 10 tissues rather than 11.

### `<tissue>.leafcutter.bed.gz`

bgzipped BED, tabix-indexed. Columns 1–4 are `#chr, start, end, ID`; columns 5+
are one per sample.

- `ID` = `chr1:14829:14970:clu_2387_-:ENSG00000310526.1`
- Values are **quantile-normalized (qqnorm)** intron excision ratios, roughly
  N(0,1)

> **⚠️ These are not covariate-residualized.** `cluster_prepare_fastqtl.py` only
> applies `qqnorm`. Sequencing-center batch effects are still in the phenotype.
> That is why every Phase 1 model carries an explicit `center` term.

### `<tissue>.leafcutter.PCs.txt`

TSV, rows PC1–PC10, columns = samples.

> **⚠️ Column headers are truncated to `donor-tissue`** (`SMHT001-3AH`), so they
> are **non-unique** when a donor has two cores. The standard
> `merge(by="donor_tissue")` + `match(colnames)` pattern silently gives both
> replicate cores the *first* PC vector.

### `<tissue>.leafcutter.phenotype_groups.txt`

Two-column TSV: junction ID → ENSG. Produced for grouped QTL testing. Can be used
as an authoritative cluster→gene map instead of building one with `pyranges`.

### Cohort sizes

| Tissue | samples | introns |
|---|---|---|
| 3AH_Muscle | 32 | 131,164 |
| 3AF_Skin_Abdomen | 17 | 148,992 |
| 3AM_Brain_Cerebellum | 15 | 122,914 |
| 3S_Heart | 15 | 107,279 |
| 3AK_Brain_Frontal_Lobe | 14 | 139,504 |
| 3C_Esophagus | 14 | 105,200 |
| 3Q_Lung | 14 | 137,241 |
| 3AL_Brain_Temporal_Lobe | 12 | 115,608 |
| 3G_Colon_Desc | 12 | 100,310 |
| 3AD_Skin_Calf | 11 | 104,260 |
| 3U_Testis_L | 11 | 188,376 |

`3A_Whole_Blood` is **absent** — which blocks HMHA1, RIMS1 and DDX10.

---

## 5. Covariates

### `covariate_matrix.csv` — 167 rows

Built by parsing sample IDs, not from a metadata file:

```
SMHT020-3AH-001C5-M81-A101-bcm
   │     │    │    │ │   │    │
 donor  tissue core sex age protocol center
```

| Column | Example |
|---|---|
| `sample_id` | `SMHT020-3AH-001C5-M81-A101-bcm` |
| `donor` | `SMHT020` |
| `tissue_code` | `3AH` |
| `core` | `001C5` |
| `sex` | `M` |
| `age` | `81` |
| `protocol` | `A101` |
| `center` | `bcm` |
| `donor_tissue` | `SMHT020-3AH` |
| `has_multi_center_replicate` | TRUE/FALSE |

Cohort composition: 4 centers (broad 49, nygc 50, bcm 45, washu 23), 53 F /
114 M, ages 42–89, and **34 donor-tissue pairs sequenced at two centers** — the
empirical handle on batch effects.

Measured batch effect: within-donor-different-center ρ = **0.954** vs
between-donor-same-tissue ρ = **0.946**. Δ = 0.008. Tissue identity dominates;
center batch is real but small at the bulk PSI level.

---

## 6. Reference annotation

| File | Used by | Note |
|---|---|---|
| `gencode.v44.basic.gtf.gz` (28 MB) | Phase 1 | gene bodies for variant annotation |
| `gencode.v47.annotation.gtf` (1.8 GB) | Phase 3 | MANE Select transcripts + exon structure |

The v44/v47 split is historical — Phase 1 predates the switch. Laura's LeafCutter
run already used v47, so Phase 1's cluster→gene map and the perind files were
built against different annotation versions.

---

## 7. RNA-seq alignments

`gs://smaht-p25-bam-files/*.bam`

```
SMHT039-3AD-006A2-M74-A101-nygc-SMAFIWDGZQQV-star_2.7.10b_GRCh38_gencode_v47.aligned.sorted.bam
```

Sample ID is everything before `-SMAFI`. Index files (`.bai`) are alongside.

Also mirrored on the Broad cluster at
`/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/`.

Used for coverage (`samtools depth`), junction extraction (`pysam` CIGAR-N with
`XS` tag strand filtering), and library size (`samtools idxstats`).

---

## 8. Files present in the original project but never read

Worth noting so they are not mistaken for dependencies:

| File | Size | Status |
|---|---|---|
| `germline_phased_snvIndel/release_*/*.vcf.gz` | **13 GB** | 25 donors of phased germline WGS (PacBio+ONT, WhatsHap 2.8, ~1 switch error / 30 Mb). **No script in the project reads these.** The germline work used a separate SpliceAI-annotated table instead. |
| `SMaHT_P25_RNA_WM_samples.tsv` | 214 rows | Full RNA sample sheet with `FinalQCStatus` and `QCNotes`. Metadata is parsed from sample ID strings instead. |
| `leafcutter_results/10_tissue_groups.zip` | — | **Not splicing data.** An LDSC partitioned-heritability run from an unrelated project. Misfiled. |
