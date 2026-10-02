# Phase 2 — Per-donor outlier scan (individual test)

## Research question

> **When SpliceAI flags a somatic variant as splice-disrupting, can we detect the
> effect in that same donor's RNA-seq?**
>
> Reframed: *in this one person, is their junction usage abnormal compared to
> everyone else in the same tissue?*

## Test type

**Individual test.** One carrier against the distribution of all other donors in
the same tissue. This replaces the Phase 1 group test, because most
high-SpliceAI variants have only a single carrier per gene × tissue — which
leaves a conventional case/control design with no power.

---

## Pipeline

```
V2 SpliceAI calls (180,162 variants)      LeafCutter perind (10 tissues)
            │                                       │
            │  score >= 0.2  →  717 variants        │
            └───────────────┬───────────────────────┘
                            │  join on (donor, tissue, position)
                            ▼
            01_inventory_testable_triplets.py
                            ▼
                01_donor_tissue_inventory.csv

            02_scan_psi_outliers.py                 <- the core analysis
                            ▼
            02_outlier_scan_all.csv  (71 testable triplets)
                            │
            ┌───────────────┴───────────────┐
            ▼                               ▼
  03_rank_candidates_by_visibility.py   04_apply_strict_artifact_filter.py
            ▼                               ▼
  recommended_candidates.csv (15)      strict_visibility_filtered.csv (1)

  Figures:  05_plot_case_study_arcs.py        per-case arc + psi scatter
            06_plot_vaf_spliceai_landscape.py SpliceAI x VAF landscape
            07_plot_psi_stripplot.py          per-donor psi strip (best for low VAF)
            08_plot_pileup_cummings_style.py  Cummings et al. Fig.2 layout

  Orthogonal check:
            09_fraser2_cross_validation.R     FRASER2 / DROP, then overlap
```

---

## File I/O

### Inputs

| File | Format | What it holds |
|---|---|---|
| `V2_spliceAI_VAF_quantile_fraction.txt` | TSV, 180,162 rows | 18 columns: `gene, score, id, chr, pos, ref, alt, type, donor_tissue, VAF_each_core, SR_VAF, CrossTissue, inTissues, SampleDonors, TissueID, SampleTissues, vaf_bin, spliceAI_type`. **717 rows have `score >= 0.2`.** |
| `perind/<tissue>_perind.counts.gz` × 10 | LeafCutter | Laura's **re-clustered** output, after mappability + annotated-junction QC (45.6% junction retention). Row key `chr:start:end:clu_N_strand`, cells `numer/denom`. |
| `*.bam` (`gs://smaht-p25-bam-files/`) | BAM | Only for the pileup figures (scripts 05–08) |

Tissue coverage: 9–21 donors each, 59,820–100,818 junctions, 19,766–31,498
clusters. `3A_Whole_Blood` is **absent** from Laura's output, which blocks
HMHA1, RIMS1 and DDX10.

### Outputs

| File | Rows | Contents |
|---|---|---|
| `01_donor_tissue_inventory.csv` | 10 | `tissue, n_donors, n_samples, cores_per_donor, n_junctions, n_clusters, donors` |
| **`02_outlier_scan_all.csv`** | **71** | `gene, donor, tissue, variant_id, chr, pos, SpliceAI, SR_VAF, cluster, jct_start, jct_end, jct_len, psi_carrier, psi_others_median, psi_shift, z, shift_over_vaf, n_donors_cohort, cluster_reads_median` |
| `recommended_candidates.csv` | 15 | the scan table + `abs_z, abs_delta_psi, visibility` |
| `strict_visibility_filtered.csv` | **1** | same + `dist_to_junction`, after 6-way filtering |

---

## Model

This phase uses **no regression**. The statistic is a robust outlier score.

### ψ (intron usage ratio)

```
psi = junction reads / cluster total reads
```

Note this is LeafCutter's intron-usage ratio, **not** classical cassette-exon
PSI. Within a cluster, ψ sums to 1.

Multiple cores per donor are collapsed by `nanmean` before scoring, so a donor
with two sequencing runs does not get double weight.

### MAD-based robust z

```python
z = (psi_carrier - median(psi_others)) / (1.4826 * MAD(psi_others))
```

Falls back to mean/SD when MAD is 0 or NaN. Returns NaN with fewer than 5 cohort
values. Implemented once in `utils.py:robust_zscore()` and used unchanged for
the rest of the project.

### The three filters

A triplet is called a hit only if all three hold:

| Filter | Rationale |
|---|---|
| `\|z\| >= 2` | statistical outlier |
| **`\|Δψ\| < VAF`** | **physical plausibility** — a variant in 5% of cells cannot shift bulk ψ by more than ~5%. Values above 1 indicate a false positive or a cohort-variance problem. |
| variant within 2 kb of the cluster | causal plausibility. The variant need not be at a junction endpoint; it may sit inside the intron. |

The second filter is the methodological contribution of this phase. It converts
VAF from a nuisance parameter into a hard upper bound on the observable effect.

### Visibility score (script 03)

Introduced once it became clear the figures were not showing much:

```
visibility = |z| x VAF x |Δψ| x log10(cluster_reads_median + 1)
```

### Strict artifact filter (script 04)

Six conditions, all required:

```
SpliceAI            >= 0.2
VAF                 >= 0.02
|z|                 >= 3
cluster_reads_median >= 200      # kills low-count artifacts
|Δψ|                >= 0.03
dist_to_junction    <= 500       # variant inside/adjacent to the outlier junction
```

The `cluster_reads >= 200` term matters: CAPN11's spectacular Δψ = 0.565 comes
from a cluster with **18 reads**.

---

## Tools

| Tool | Where |
|---|---|
| `pandas` / `numpy` | perind parsing, ψ, MAD z |
| `pysam` | junction extraction via CIGAR `N` operations |
| `samtools depth` | coverage. Replaced `pysam.pileup`, which smears introns instead of dropping them to zero. |
| `matplotlib.patches.Arc` | junction arcs |
| `FRASER` (Bioconductor) | orthogonal outlier detection, script 09 |

---

## Results

```
180,162  V2 somatic variants
    717  SpliceAI >= 0.2
     71  testable triplets (carrier has matched RNA, cluster has >=10 median reads)
      8  pass all three filters
      3  also have the variant inside/adjacent to the cluster
      1  survives the strict artifact filter
```

**The 8 that pass the three filters:**

| Gene | Donor | Tissue | SpliceAI | z | Δψ | VAF |
|---|---|---|---|---|---|---|
| PFKP | SMHT040 | 3AF | — | +5.25 | — | — |
| **ISCU** | SMHT039 | 3AD | 0.99 | **+4.31** | +0.011 | 0.019 |
| **C11orf54** | SMHT039 | 3AD | 0.38 | −3.84 | −0.044 | 0.056 |
| **PLBD2** | SMHT016 | 3AD | **1.00** | −3.55 | −0.038 | 0.058 |
| PON2 | SMHT039 | 3AD | — | +2.93 | — | — |
| MFNG | SMHT016 | 3AD | — | +2.74 | — | — |
| ROCK2 | SMHT027 | 3AM | — | −2.50 | — | — |
| MARK2 | SMHT039 | 3AD | — | +2.07 | — | — |

All three headline cases are in 3AD skin — the tissue with the highest somatic
mutation burden.

**After the strict filter: exactly one survivor, PLBD2.** ISCU fails on
|Δψ| = 0.011 < 0.03. C11orf54 fails. Every new high-visibility candidate fails on
read depth or variant-to-junction distance.

---

## The number that defines the project

Stated at the top of `03_rank_candidates_by_visibility.py`:

> **Maximum VAF in the entire dataset is 0.058 — PLBD2 is already the highest.**
> **Zero candidates have VAF ≥ 0.10. Median VAF is 0.020.**

The cohort is fundamentally low-VAF. Ranking candidates by VAF alone is useless,
and no amount of figure polish will make a 2% effect visible against 98%
wild-type transcript.

---

## An instructive negative: RSRC2

`RSRC2 chr12:122,519,030 C>T` — VEP `splice_acceptor_variant`, HIGH impact
(c.208-1G>A, the canonical AG at −1), **SpliceAI = 1.00**, VAF 0.014.

LeafCutter gives **|z| = 0.77, |Δψ| = 0.02** — no outlier signal at all.

Two non-exclusive explanations:

1. **NMD.** The aberrant transcript carries a premature stop and is degraded
   before it can be sequenced, so the observed fraction is *lower* than the
   cellular fraction.
2. **VAF.** At 1.4%, the aberrant isoform is simply invisible against the intact
   one.

A perfect prediction from both tools, and zero RNA evidence. This case is the
clearest statement of the short-read detection ceiling.

---

## Figure development

`08_plot_pileup_cummings_style.py` is the last of seven iterations. The
progression, and why each step was needed:

| Change | Reason |
|---|---|
| region-local RPM → **global library-size RPM** | the original normalization made cross-donor comparison meaningless |
| grey spaghetti → **blue mean + min–max band** | 10 overlapping non-carrier traces are unreadable |
| add **isoform track** at the bottom | you cannot interpret an arc without knowing where the exons are |
| stack **carrier above control** (Cummings Fig. 2) | overlaid traces hide which one is higher |
| **label every arc with its read count** | a 15-read cryptic arc is invisible next to a 1,000-read canonical one, but the number is not |
| `pysam.pileup` → **`samtools depth`** | pysam smears coverage across introns instead of dropping to zero |

`07_plot_psi_stripplot.py` is the honest alternative: one dot per donor, carrier
as a star. Its docstring says it directly — *"compared to sashimi/pileup which
fails to show low-VAF signal, this strip plot dramatically shows carrier vs
cohort separation."* For a 2% effect, a strip plot works where a pileup does not.

That script also matches junctions on `chr:start:end` only, **deliberately
ignoring the cluster ID**, because LeafCutter renumbers clusters between runs.
The PLBD2 cluster is `clu_5571_+` in one run and a different ID in another.

---

## Paths to edit

| Script | Hardcoded path |
|---|---|
| `utils.py`, `01`–`04` | `ROOT = /Users/spectremac/Desktop/Differential_splicing` |
| `05`–`08` | `/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/*.bam` (Broad) |
| `05`–`08` | `LIBSIZES` dict — hardcoded from `libsize_skincalf.txt` |
| `09_fraser2_cross_validation.R` | `/medpop/esp2/yxliu/SMAHT/fraser2_output/` |

---

## Caveats

- **ψ here is LeafCutter's intron-usage ratio, not cassette-exon PSI.** Stated
  explicitly because the two are routinely confused, and the DMD vignette in
  Phase 1 depends on the distinction.
- **LeafCutter cluster IDs are not stable across runs.** Always match on
  `chr:start:end`.
- **All three headline cases are in one tissue** (3AD skin). That is a mutation
  burden effect, not a biological statement about skin splicing.
- **The per-donor sashimi reveals that the ISCU cryptic junction is not
  carrier-exclusive**: SMHT039 (carrier) has 15 reads, but SMHT007 has 14 and
  SMHT022 has 12. The carrier is the highest but not an isolated outlier in raw
  counts — it is the ψ normalization that produces |z| = 4.31. See
  `phase3/07_plot_per_donor_sashimi.py`.
