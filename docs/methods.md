# Methods

Statistical detail for all three phases. For file formats see
[`data_dictionary.md`](data_dictionary.md); for findings see
[`results.md`](results.md).

---

## The quantity being modelled

Everything in this project models **ψ (psi), LeafCutter's intron usage ratio**:

```
psi = reads supporting this junction / total reads in the junction's cluster
```

A LeafCutter *cluster* is a set of junctions sharing splice sites — one
alternative splicing unit. Within a cluster, ψ sums to 1.

> **This is not classical cassette-exon PSI.** The two are routinely confused.
> Classical PSI is inclusion / (inclusion + skipping) for one exon; LeafCutter's
> ψ is one junction's share of its whole cluster. The distinction matters for
> interpretation — a cluster-level Dirichlet-Multinomial test can detect
> multi-exon events that an exon-level PSI cannot.

When a donor has multiple sequencing cores, ψ is collapsed by `nanmean` (Phase 2)
or `median` (Phase 1) before any statistic is computed, so a donor with two runs
does not get double weight.

---

## Phase 1 — Group test

### Design

Carriers vs non-carriers within one tissue, for one gene.

- **Carrier**: donor has ≥1 HighConf somatic mutation in gene *G*
- **Non-carrier**: donor has no somatic mutation in *G* in that tissue
- **Inclusion**: ≥4 carriers and ≥4 non-carriers (relaxed to ≥3 RNA cores per
  group at test time, because donors with mutations often lack matched RNA)

### Models

```
A.1   logit(PSI) ~ carrier_flag + center + sex + age + PC1..PCk
A.2   logit(PSI) ~ VAF          + center + sex + age + PC1..PCk     # carriers only
A.3   logit(PSI) ~ VAF          + center + sex + age + PC1..PCk     # all donors, NC = 0
A.4   logit(PSI) ~ donor_max_SpliceAI + center + sex + age + PC1..PCk
```

A.3 is the primary test. Setting non-carriers to VAF = 0 rather than dropping
them keeps the full cohort and treats clonal fraction as a dose.

### Adaptive covariate selection

Small tissues (n = 9–12) saturate the design matrix: 12 parameters at n = 11
leaves `df_resid = 0` and every p-value `NA`. The engine drops terms in a fixed
order until `df_resid ≥ 3`:

```
PC5 → PC4 → PC3 → PC2 → PC1 → age → sex → center
```

Without this, four tissues (3AD, 3AK, 3AL, 3G) silently returned "zero hits"
when in fact no test had run. With it, 3AD recovered 39 hits and 3G recovered 25.

### Why continuous VAF matters

A median split on VAF found almost nothing (9 hits across all tissues). The same
data with VAF as a continuous predictor found 13–28 hits per tissue. Dichotomizing
a continuous dose destroyed the signal.

### Cluster-level test

LeafCutter's `leafcutter_ds.R` — a **Dirichlet-Multinomial** likelihood-ratio
test on the full count vector of a cluster. This is the statistically correct
test for compositional junction data and was Laura's explicit recommendation,
because it detects multi-exon events that per-junction PSI regression misses.

It requires **two clean groups with ≥4 samples each**. It would not install
locally (Stan/rstan compile failure on macOS, three install strategies attempted)
and was run on the Broad cluster as a qsub array — ~12 minutes per candidate.

### Isoform-switch detection

Within a cluster ψ sums to 1, so a genuine isoform switch leaves a **mirror
pair**: two junctions with opposite-sign β of comparable magnitude. Flagged when:

- both junctions have q < 0.10
- β signs differ
- |β| ratio is within 0.5–2.0

Example from 3S Heart cluster `clu_5278_-`: β = −4.43 (q = 2.5e−13) and
β = +4.43 (q = 2.9e−13). That structure is hard to produce by noise and is used
as evidence the signal is real splicing rather than a coverage artifact.

### Multiple testing

Benjamini–Hochberg **within each (tissue, gene, test)** — typically ~20 junctions
per gene.

> This is far more permissive than correcting genome-wide. An earlier
> cohort-wide version of this analysis corrected across ~30,000 junctions per
> tissue. **Hit counts from the two are not on the same scale** and should not be
> compared.

### Carrier definitions

Run in parallel to ask whether the signal is driven by splice-predicted variants:

| Definition | Criterion |
|---|---|
| **ANY** | any HighConf somatic mutation in the gene |
| **LOOSE** | a mutation with SpliceAI ≥ 0.2 |
| **STRICT / Q95** | a mutation with `spliceAI_type ∈ {Q95, Q99}` |
| **ANY_noGerm** | ANY, minus donors carrying a high-confidence *germline* SpliceAI variant in the same gene |

Donors with a sub-threshold mutation are **excluded** from LOOSE/STRICT runs
rather than used as controls — they are neither clean controls nor convincing
carriers.

The `ANY_noGerm` tier exists to answer a specific reviewer concern: *is this
actually a germline splicing QTL leaking through?* For DMD it answered cleanly —
excluding the one germline-suspect carrier **improved** q by ~4× and **doubled**
the effect size, so that carrier was diluting the signal, not creating it.

---

## Phase 2 — Individual test

### Design

One carrier against the distribution of all other donors in the same tissue.
No group comparison, because most high-SpliceAI variants have exactly one carrier
per gene × tissue.

### Statistic

A MAD-based robust z-score:

```python
z = (psi_carrier - median(psi_others)) / (1.4826 * MAD(psi_others))
```

The 1.4826 factor makes MAD a consistent estimator of σ for normally distributed
data. Falls back to mean/SD when MAD is 0 or NaN; returns NaN with fewer than 5
cohort values.

MAD rather than SD because the cohort is small (9–21 donors) and a single aberrant
donor would otherwise inflate the denominator and hide itself.

### The three filters

A triplet is a hit only if all three hold:

| Filter | Threshold | Rationale |
|---|---|---|
| Statistical | `\|z\| >= 2` | the carrier is an outlier |
| **Physical** | **`\|Δψ\| < VAF`** | **a variant present in a fraction *f* of cells cannot shift bulk ψ by more than ~*f*.** Values above 1 indicate a false positive or a cohort-variance problem. |
| Causal | variant within 2 kb of the cluster | plausible mechanism. The variant need not sit at a junction endpoint; it may be inside the intron. |

The physical filter is the methodological contribution of this phase. It turns
VAF from a nuisance parameter into a **hard upper bound on the observable
effect**, and it is what makes a 1-vs-cohort comparison interpretable at all.

### Cluster coverage requirement

`median(cluster reads across cohort) >= 10` — below that, ψ is dominated by
sampling noise.

### Visibility score

Introduced once it was clear the figures were not showing much:

```
visibility = |z| x VAF x |Δψ| x log10(cluster_reads_median + 1)
```

Not a formal statistic — a heuristic for ranking which candidates are most likely
to be *visible* in a coverage figure. Combines statistical strength, clonal
fraction, effect size and read depth.

### Strict artifact filter

Six conditions, all required:

```
SpliceAI             >= 0.2
VAF                  >= 0.02
|z|                  >= 3
cluster_reads_median >= 200
|Δψ|                 >= 0.03
dist_to_junction     <= 500
```

The `cluster_reads >= 200` term is the one that does the work: CAPN11's apparently
spectacular Δψ = 0.565 comes from a cluster with **18 reads**.

### Orthogonal validation — FRASER2

`09_fraser2_cross_validation.R` is a deliberately faithful reimplementation of the
DROP `aberrantSplicing` module, using FRASER2 with DROP's default parameters:

```
implementation           = "PCA"
minExpressionInOneSample = 20
quantileForFiltering     = 0.75
padjCutoff               = 0.1
deltaPsiCutoff           = 0.1
metric                   = jaccard
```

Encoding dimension `q` via `optimHyperParams`, falling back to q = 3 for the small
cohort.

FRASER2 is attractive here precisely because it is **per-sample normalized and
phenotype-agnostic** — it does not need a SpliceAI candidate list and it sidesteps
the library-size artifacts that dominate RPM comparisons. A custom section merges
FRASER2 outliers with the LeafCutter scan on (gene, donor) and asks directly
whether the three hero cases are recovered.

---

## Phase 3 — Prediction-first

**No statistical model.** This phase is annotation, filtering and visualization.
RNA evidence is read visually.

The shift is deliberate. Phases 1 and 2 let the RNA pick candidates — significance
in a regression, or an outlier ψ. That produced artifacts: an unsupervised scan of
four large genes returned 273 hits whose **maximum SpliceAI was 0.00**, with one
donor contributing 126 of them. Phase 3 inverts it: **predictions pick the
candidates, RNA confirms.**

### Candidate filter

```
SpliceAI >= 0.2  AND  Consequence contains "splice"
```

The substring test is deliberately loose, so it admits weak classes
(`splice_polypyrimidine_tract_variant`, `splice_donor_5th_base_variant`,
`splice_donor_region_variant`, `missense_variant,splice_region_variant`). That is
why 82 of 213 rows are LOW impact — and why the downstream STRICT tier adds
`IMPACT == HIGH`.

### Tiers

| Tier | Filter | n | Axis |
|---|---|---|---|
| enriched | SpliceAI ≥ 0.2 + VEP splice + GENCODE v47 annotation | 213 | annotation |
| **STRICT** | **SpliceAI ≥ 0.8 AND IMPACT == HIGH** | **77** | confidence (two tools agree) |
| HIGH_VAF | VAF ≥ 0.05, SpliceAI threshold dropped | 19 | detectability |
| selected | HIGH + SpliceAI ≥ 0.5 + VAF ≥ 0.041, minus 1 manual drop | 17 (16 distinct) | VAF-maximizing |
| Tier A | STRICT + carrier BAM already available | 3 | practical |

**STRICT and enriched are different axes, not nested tiers of one criterion.**

### Priority score

```
priority = 1
         + 2 * has_RNAseq        # tissue has a LeafCutter run
         + 2 * (IMPACT == HIGH)
         + 3 * (SpliceAI >= 0.8)
         - 1 * in_scan           # already covered by Phase 2
```

Range 1–8. Reverse-engineered from the original table and verified to reproduce
it exactly. The 49 rows at priority 8 are the practical working set.

### Junction extraction for figures

`pysam`, CIGAR `N` operations, with an **explicit `XS` tag strand check** against
the gene strand — reads lacking a matching `XS` are discarded. This matches the
rigor of Laura's `regtools junctions extract -s XS` call.

Filters: `MIN_JCT_READS = 5` (in ≥1 donor), `MIN_JCT_LENGTH = 25`, mapq ≥ 1,
excluding duplicate / secondary / supplementary / QC-fail reads.

### Coverage for figures

`samtools depth -a` via subprocess, **not** `pysam.pileup`. The latter smears
coverage across introns instead of dropping it to zero, which makes exon
boundaries unreadable.

### Normalization

Global library-size RPM:

```
RPM = coverage / library_size * 1e6
```

with library size from `samtools idxstats` **computed per BAM at runtime**.
Hardcoding library sizes is wrong across sequencing centers.

An earlier version used region-local normalization (`cov / cov.sum() * 1e6`),
which makes cross-donor comparison meaningless. That was a real bug, fixed in
the figure-development sequence.

---

## Why each phase ended

| Phase | Wall |
|---|---|
| **1** | Most high-SpliceAI variants have **one carrier per gene × tissue**. 1 vs 10 has no power under any group test. Separately, V1 SpliceAI covered only 33.5% of HighConf mutations, so the splice-stratified carrier definitions produced empty candidate lists and test A.4 returned zero hits. |
| **2** | **Maximum VAF in the entire cohort is 0.058; median 0.020.** After strict artifact filtering, exactly one candidate survives. No ranking heuristic can manufacture signal that is below the noise floor. |
| **3** | Ongoing. Tier A figures were generated on the cluster but the outputs were never pulled back. |

---

## The constraint that governs everything

> A somatic variant present in 2% of cells can shift bulk intron usage by at most
> ~2%. Short-read bulk RNA-seq cannot resolve that against 98% wild-type
> transcript — and if the aberrant transcript carries a premature stop codon,
> nonsense-mediated decay reduces the observed fraction below even the cellular
> fraction.

Two directions follow:

- **Long-read (PacBio Kinnex)** — full-length isoforms resolve structures that
  junction-level inference can only infer, and reads are long enough to phase a
  variant onto the transcript carrying the aberrant splice.
- **Unsupervised outlier detection (FRASER2 / DROP)** — per-sample normalized,
  phenotype-agnostic, does not depend on a SpliceAI candidate list. The
  cross-validation harness is already in `phase2/09_fraser2_cross_validation.R`.
