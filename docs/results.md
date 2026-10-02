# Results

Findings across all three phases, including the negative results and the
artifacts identified along the way.

---

## Phase 1 — Group test

### Funnel

```
34,907  HighConf somatic mutations
19,938  inside an annotated gene        (13,852 protein-coding)
 4,805  unique mutated genes            (3,319 protein-coding)
15,644  (tissue, gene) pairs with >=1 carrier
   182  Candidate A pairs  (>=4 carriers, >=4 non-carriers, protein-coding)
    62  testable after RNA matching and >=3 cores per group
 1,458  junction-level lm tests
    10  pairs at FDR < 10%
     6  pairs at FDR < 5%
     5  clusters with mirror-pair architecture
     0  genes concordant across tissues
```

### Top hits (test A.3)

| Gene | Tissue | q | β | Note |
|---|---|---|---|---|
| **DMD** | 3C Esophagus | 0.0013 | −1.64 | the only hit significant under both lm and DM |
| **RBFOX1** | 3AH Muscle | 0.0016 | −0.33 | itself an RNA-binding splicing factor |
| DUSP22 | 3AK Brain frontal | 0.0051 | — | probable artifact, see below |
| **PTPRD** | 3Q Lung | 0.013 | +4.46 | near-binary (ψ 0.97 ↔ 0.03) |
| DUSP22 | 3U Testis | 0.03 | — | |
| NRG1 | 3Q Lung | 0.05 | +7.21 | |

**Effect sizes are modest.** RBFOX1 in muscle is a ~3% absolute PSI difference
(β = −0.33 per unit VAF), though consistent in direction across all 32 samples.

### lm vs Dirichlet-Multinomial

64 candidates were re-tested at the cluster level with `leafcutter_ds`:

| Category | n | Genes |
|---|---|---|
| Both non-significant | 92 | — |
| lm only | 9 | RBFOX1/Muscle, DUSP22 ×3, PTPRD/Lung, DMD/Skin, NRG1/Lung, … |
| DS only | 2 | MACROD2/Lung (q = 0.061), KMT2C/Muscle (q = 0.075) |
| **Both significant** | **1** | **DMD / 3C Esophagus** |

RBFOX1 — the headline of an earlier presentation — is **lm-only**: DS q = 0.350,
Δψ = 0.0064. It was honestly downgraded.

### The DMD vignette

Five carriers vs nine non-carriers in esophagus, cluster `clu_20135_-`.

One of the five carriers also had a **germline** DMD variant with SpliceAI > 0.5.
Re-running with that donor excluded:

| Test | n carriers | q | effect |
|---|---|---|---|
| lm (A.3) | 5 | 0.0013 | β = −1.64 |
| DM, ANY | 5 | 0.016 | Δψ = −0.11 |
| **DM, ANY_noGerm** | **4** | **0.0042** | **Δψ = −0.24** |

**q improves ~4×, effect size doubles.** The germline-suspect carrier was
*diluting* the signal, not driving it — which directly answers the concern that
this might be a germline splicing QTL leaking through.

Mean ψ on the long-range junction: **0.13 in non-carriers → 0.37 in carriers** —
a 24-point shift corresponding to an isoform switch that skips two exons.

### Strategy B — a clean negative

| Classification | n | Genes |
|---|---|---|
| CONCORDANT | **0** | — |
| DISCORDANT | 2 | DMD, DUSP22 |
| TISSUE_SPECIFIC | 3 | NRG1, PTPRD, RBFOX1 |
| NO_SIGNAL | 5 | ERBB4, GPRIN2, KMT2C, VPS13B, ZFHX3 |

> **Zero junctions reach q < 0.10 at the same chromosomal position in ≥2 tissues.**

DMD flips sign between tissues: esophagus β = −1.64 (q = 1.3e−3), skin β = +2.17
(q = 0.07), lung β = +1.42 (n.s.). This is direct empirical support for analyzing
tissues separately rather than pooling them.

### The SpliceAI negative result

This is what drove the pivot to Phase 2.

- V1 SpliceAI annotated only **33.5%** of HighConf mutations (11,696 / 34,907)
- Of 182 candidates, **173 had zero SpliceAI variation** among their carriers
- Test A.4 (SpliceAI as predictor) returned **zero hits**; best q = 0.12, and only
  9 of 62 candidates had any SpliceAI variation at all
- The SpliceAI-stratified carrier definitions produced **empty candidate lists**
  (LOOSE ≥3 carriers: 0; LOOSE ≥2: 0; Q95 ≥2: 0)
- **All 10 Strategy-A hits had carriers whose maximum SpliceAI was 0.00** — 17
  distinct mutations across RBFOX1 (5 donors), DUSP22 (6), DMD (6), **none within
  50 bp of an annotated splice site**

Three honest interpretations were offered at the time: (a) V1 coverage gap,
(b) SpliceAI is blind to exonic-splicing-enhancer, regulatory and trans effects,
(c) germline sQTL leakage. The germline screening in the DMD vignette addressed
(c) for that one gene.

---

## Phase 2 — Individual outlier scan

### Funnel

```
180,162  V2 somatic variants
    717  SpliceAI >= 0.2
     71  testable triplets
      8  pass all three filters
      3  also have the variant inside/adjacent to the cluster
      1  survives the strict artifact filter
```

### The 8 that pass

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
mutation burden, not a biological statement about skin.

### After the strict filter: one survivor

**PLBD2 / SMHT016 / 3AD** is the only candidate passing all six conditions.

- ISCU fails on |Δψ| = 0.011 < 0.03
- C11orf54 fails
- Every new high-visibility candidate fails on read depth or variant-to-junction
  distance

### The defining number

> **Maximum VAF in the entire dataset is 0.058 — PLBD2 is already the highest.**
> **Zero candidates have VAF ≥ 0.10. Median VAF is 0.020.**

The cohort is fundamentally low-VAF. No ranking heuristic can manufacture signal
that is below the noise floor.

### Candidates that looked good and were not

Top visibility scores beat the original three heroes, but each failed on
inspection:

| Gene | visibility | Why it fails |
|---|---|---|
| MYO15B | 0.0636 | — |
| **CAPN11** | 0.0334 | **Δψ = 0.565 from a cluster with 18 reads.** This is the case that motivated the `cluster_reads >= 200` filter. |
| MCM7 | 0.0318 | SpliceAI 0.94 but VEP LOW |
| SLC43A3 | 0.0315 | — |
| SGIP1 | — | \|z\| = 7.40, but low depth |

---

## Phase 3 — Prediction-first

### Candidate list

```
717  variants at SpliceAI >= 0.2
213  with VEP Consequence containing "splice"
       95 HIGH / 82 LOW / 36 MODERATE
       24 donors; tissue-skewed: 3AD 106/213, 3A 21, 3AF 16, 3I 13
       has_RNAseq: 153 yes / 60 no
 77  STRICT  (SpliceAI >= 0.8 AND IMPACT == HIGH)
 19  HIGH_VAF (VAF >= 0.05)
 17  selected (16 distinct variants)
  3  Tier A (carrier BAM available)
```

### The three case studies, cross-validated against VEP

| Case | SpliceAI | VEP | Verdict |
|---|---|---|---|
| **ISCU** | 0.99 | `splice_donor_variant`, **HIGH** | **Both tools agree — the most credible case.** |
| **PLBD2** | **1.00** | **LOW** (polypyrimidine tract, deep intronic, not canonical ±2) | **Informative discordance.** This is what SpliceAI adds over VEP's hard rules. |
| **C11orf54** | 0.38 | **`missense_variant` (p.G204E)** | **Refuted.** Not a splicing variant. Downgraded to "mechanism unclear". |

C11orf54 is **absent from the 213-row list by construction** — its VEP consequence
contains no `"splice"` substring. That absence is itself the finding.

### ISCU in detail

Variant at chr12:108,565,432, the exon 3 splice donor +1. Zooming to the
downstream region and labelling every arc with its read count:

| Junction | carrier | cohort | ratio |
|---|---|---|---|
| canonical e4→e5 | 1156 | 763 | 1.5× (baseline expression) |
| **cryptic donor → e5** | **15** | **4** | **3.75×** |

The cryptic excess over baseline is ~2.5×, and that is the variant-specific
signal behind |z| = 4.31.

> **⚠️ Important caveat.** The per-donor sashimi shows the cryptic junction is
> **not carrier-exclusive**:
>
> | Donor | cryptic reads |
> |---|---|
> | **SMHT039 (carrier)** | **15** |
> | SMHT007 | 14 |
> | SMHT022 | 12 |
> | SMHT027 | 7 |
>
> The carrier is highest, but not an isolated outlier in raw counts. The |z| = 4.31
> comes from ψ normalization. Any writeup should state this.

### Donor 005 — both leads negative

The PI flagged this donor as a likely smoking gun. Two variants investigated:

- **RSRC2** — splice acceptor −1 (c.208-1G>A), VEP **HIGH**, SpliceAI **1.00**,
  VAF 0.014. LeafCutter |z| = **0.77**, |Δψ| = 0.02. **No RNA signal despite a
  perfect prediction from both tools.** Two non-exclusive explanations: the
  aberrant transcript carries a premature stop and is degraded by NMD before
  sequencing, or VAF 1.4% is simply invisible.
- **DMD** — an unsupervised-scan hit with SpliceAI ≈ 0. Carrier coverage is
  systematically low across the whole gene — a low-expression-tissue artifact,
  since DMD is barely expressed in skin.

> **Take-home:** even a perfect prediction (VEP HIGH + SpliceAI 1.00) can be
> short-read invisible at ~2% VAF with bulk RNA and NMD. This is a hard ceiling on
> short-read validation alone.

### SLC10A6 tandem repeat — expression, not splicing

Following a long-read haplotype-specific TR finding in intron 1, a cryptic
acceptor at chr4:86,837,573 looked promising: cohort ψ 0–0.032, with SMHT029
ranked first (0.032) and SMHT039 second (0.018).

The intron/exon coverage ratio check closed it:

| Donor | intron 1 / exon 1 |
|---|---|
| SMHT029 | **0.12** |
| SMHT039 | **0.11** |
| cohort median | **0.15** |

**Both candidates are below the cohort**, which rejects intron retention. The real
explanation is that SMHT029 expresses SLC10A6 at roughly twice cohort levels in
skin, so all coverage scales up together.

> **Not a splicing phenomenon — an expression phenomenon.** Allele-specific
> expression is the right framework; the repeat acts on transcription, not
> splicing.

### Unsupervised scan — the negative control

Dropping the SpliceAI filter entirely and scanning four very large genes
(CNTNAP2, CSMD1, RBFOX1, DMD) for ψ outliers:

```
273 hits:  DMD 230, RBFOX1 30, CSMD1 10, CNTNAP2 3
Maximum SpliceAI across all 273:  0.00
One donor (SMHT029) accounts for 126 of 273
```

That concentration is the tell. These are expression and QC artifacts, not
splicing events — and it is the direct justification for the prediction-first
pivot.

---

## Case study status

| Case | Status | Evidence |
|---|---|---|
| **ISCU** | ✅ credible | SpliceAI 0.99 + VEP HIGH agree; cryptic donor chr12:108,567,757; caveat above |
| **PLBD2** | ⚠️ informative discordance | SpliceAI 1.00 vs VEP LOW; only candidate surviving the strict filter |
| **DMD / esophagus** | ✅ strongest group-level result | lm + DM both significant; germline exclusion doubles the effect |
| **C11orf54** | ❌ refuted | VEP `missense_variant` p.G204E |
| **RSRC2** | ❌ perfect prediction, zero signal | VEP HIGH + SpliceAI 1.00, \|z\| = 0.77 |
| **MCM7** | 🟡 contrast case | SpliceAI 0.94, VEP LOW, \|z\| = 2.78, \|Δψ\| = 0.213 |
| **RBFOX1 / muscle** | ⚠️ downgraded | lm q = 0.0016 but DM q = 0.350, Δψ = 0.0064 |
| **SLC10A6** | ❌ expression, not splicing | intron/exon ratio below cohort |
| **DMD / Donor 005** | ❌ low-expression artifact | systematically low carrier coverage |
| **DUSP22** | ⚠️ probable artifact | present in all 11 tissues, 6–9 carriers each, VAF ≈ 0.25; needs IGV check at chr6:304,600–345,900 |

---

## Technical findings worth keeping

These are measurements, not analysis results, but they shaped every decision.

**Sequencing-center batch effect is real but small at bulk PSI level.** Using 34
same-donor-different-center replicate pairs:

| Comparison | Pearson ρ | n |
|---|---|---|
| Within donor, different center | **0.954** | 34 |
| Between donors, same tissue | 0.946 | 318 |

Δ = 0.008, median RMSE Δψ ≈ 0.10. **Tissue identity dominates.** Donor identity
adds ~1%. Center batch is under 1% — contrary to what was feared.

But covariate adjustment still moved ~6,000 junctions across p = 0.05 in a pilot
test, so `center` stays in every model.

**Donor random effects absorb some apparent signal.** Adding `(1|donor)` wiped the
burden signal in 3AH Muscle and 3S Heart entirely — those were within-donor
pseudo-replication. A cross-tissue mixed model returned **zero hits**, meaning any
burden effect is tissue-specific rather than a global program.

**Satterthwaite degrees of freedom matter at this n.** A naive Wald z-test
reported 788 significant junctions at n = 9. The same model with `lmerTest`
Satterthwaite df reported **5**.

---

## Known artifacts and gotchas

| Issue | Impact |
|---|---|
| **LeafCutter cluster IDs are not stable across runs** | PLBD2 is `clu_5571_+` in one run and a different ID in another. Always match junctions on `chr:start:end`. |
| **VEP cache 105 vs GENCODE v47** | 20 of 213 candidates lose their v47 annotation. `FAM49A` (renamed `CYRIA`) was silently dropped from a batch; `SSPO` and `MSH5-SAPCD1` cannot be plotted as written. |
| **`has_RNAseq` is tissue-level** | A row can be marked testable with no BAM for that donor. In 3AD alone, 5 donors have DNA variants but no RNA-seq. |
| **`pyranges` fills `-1`, not `NaN`** | 14,969 mutations carry `gene_name == "-1"`. Harmless downstream due to the `protein_coding` filter, but do not reuse those tables blind. |
| **Splicing PC columns are non-unique** | Keyed by `donor-tissue`; donors with two cores get the first PC vector for both. |
| **`selected_17` has duplicates** | 16 distinct variants. `MSH5` / `MSH5-SAPCD1` are the same variant; `CYP2W1 ×3` and `BRCC3 ×2` are single variants across tissues. |
| **The PSI matrices are qqnorm'd, not residualized** | Despite comments calling them "residualized PSI". Batch is still in the phenotype. |

---

## Summary

| | Phase 1 | Phase 2 | Phase 3 |
|---|---|---|---|
| **Test** | group | individual | individual, prediction-first |
| **Model** | `logit(PSI) ~ carrier/VAF + covariates`; DM at cluster level | MAD robust z + 3 filters | none (visual) |
| **Funnel** | 182 → 62 testable | 717 → 71 testable | 717 → 213 → 77 |
| **Hits** | 10 @ FDR<10%, 1 double-significant | 8 → 1 after strict filter | 77 STRICT → 3 Tier A |
| **Headline** | DMD / esophagus | PLBD2 | ISCU (both tools agree) |
| **Wall** | 1 carrier per gene×tissue | max VAF 0.058 | outputs not yet retrieved |

> **One sentence:** median somatic VAF in this cohort is 0.02, so a variant can
> shift bulk intron usage by at most ~2% — at or below the noise floor of
> short-read bulk RNA-seq, which is why every framing above eventually hits a
> ceiling and why the project is moving toward long-read and unsupervised
> outlier detection.
