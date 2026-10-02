# Phase 3 — Prediction-first candidate selection

## Research question

> **Cross-validate SpliceAI against an independent, rule-based annotator (VEP),
> build a systematic candidate list, and use RNA only to confirm — not to
> discover.**

The strategic shift is recorded as a one-line note at the end of the Phase 2
presentation script:

> *"prediction score + VEP rather than z score"*

Phases 1 and 2 let the RNA pick the candidates (significance in a regression, or
an outlier ψ). That produced artifacts: the unsupervised big-gene scan returned
273 hits whose **maximum SpliceAI was 0.00**, with one donor (SMHT029)
contributing 126 of them. Phase 3 inverts the logic — predictions pick the
candidates, RNA is only asked to confirm.

## Test type

**Individual** (carrying on from Phase 2), but with the candidate source changed
from RNA-driven to prediction-driven. **No statistical model is applied at this
stage** — this phase is annotation, filtering and visualization. RNA evidence is
read visually from pileups and per-donor sashimi plots.

---

## Pipeline

```
V2 SpliceAI (180,162)            Ensembl VEP v113 (998,650 variants)
         │                                    │
         │  score >= 0.2  →  717              │
         └────────────┬───────────────────────┘
                      │  join on "chr_pos_ref/alt"   <- note the SLASH
                      ▼
      01_build_vep_candidate_list.py
                      ▼
      vep_splice_candidates.csv  (213)
                      ▼
      02_enrich_candidates_gencode_v47.py     + gencode.v47.annotation.gtf
                      ▼
      vep_splice_candidates_v47_enriched.csv  (213 x 27)
                      │
      ┌───────────────┼───────────────┬──────────────────┐
      ▼               ▼               ▼                  ▼
  STRICT (77)    HIGH_VAF (19)   selected (17)      Tier A (3)
  SA>=0.8 &      VAF>=0.05       HIGH & SA>=0.5     STRICT &
  IMPACT HIGH                    & VAF>=0.041       BAM available

  Figures:  04_plot_candidates_by_tissue.py    candidate counts, testability gap
            05_plot_pileup_pyqtl.py            Laura's toolchain
            06_plot_pileup_standalone.py       dependency-free equivalent
            07_plot_per_donor_sashimi.py       stacked per-donor rows
            08_plot_vep_schematic.py           publication gene schematic

  Sanity:   03_unsupervised_biggene_scan.py    the negative control
            qc/check_housekeeping_coverage.sh  RNA quality vs real biology
            qc/check_intron_exon_ratio.sh      intron retention vs high expression
```

---

## File I/O

### Inputs

| File | Format | What it holds |
|---|---|---|
| `V2_spliceAI_VAF_quantile_fraction.txt` | TSV, 180,162 rows | SpliceAI score + VAF per somatic variant |
| `annotate_v2_381files.txt` | VEP tab output, 998,650 rows, **79 columns** | Ensembl VEP **v113.0**, cache `105_GRCh38`, run 2026-06-15. 103 `##` metadata lines before the header. Key columns: `#Uploaded_variation`, `Consequence`, `IMPACT`, `SYMBOL`, `MANE_SELECT`, `EXON`, `INTRON`, `BIOTYPE`, plus ~30 population-frequency columns. |
| `gencode.v47.annotation.gtf` | GTF, 1.8 GB | GENCODE v47 / Ensembl 113, GRCh38. Used for MANE Select transcripts and exon structure. |
| `perind/<tissue>_perind.counts.gz` | LeafCutter | only for the `has_RNAseq` flag and the unsupervised scan |
| `gs://smaht-p25-bam-files/*.bam` | BAM | all pileup and sashimi figures |

### The join

VEP's `#Uploaded_variation` uses a **slash** between alleles:

```
chr12_108565432_G/A
```

The SpliceAI table's `id` column uses an **underscore**:

```
chr12_108565432_G_A
```

Joining on `id` returns zero matches. Build the key explicitly:

```python
key = f"{chr}_{pos}_{ref}/{alt}"
```

With the key built this way, **all 717 variants match** — zero misses.

### Outputs

| File | Rows | Filter | Notes |
|---|---|---|---|
| `vep_splice_candidates.csv` | **213** | SpliceAI ≥ 0.2 **and** `Consequence` contains `"splice"` | 95 HIGH / 82 LOW / 36 MODERATE. 24 donors. Tissue-skewed: 3AD 106/213. |
| `vep_splice_candidates_v47_enriched.csv` | 213 × 27 | + GENCODE v47 annotation | adds `v47_var_context` (strand-aware `intron K/N`), `region_v47`, `mane_select_tid`, `gene_length_kb`, transcript counts |
| `vep_splice_candidates_STRICT.tsv` | **77** | `SpliceAI >= 0.8` **AND** `IMPACT == HIGH` | the two-tool-agreement list. SpliceAI 0.81–1.00, VAF 0.0105–0.0621, 18 donors |
| `vep_splice_HIGH_VAF.tsv` | **19** | `VAF >= 0.05` (SpliceAI threshold dropped) | 11 HIGH / 5 LOW / 3 MODERATE. **10 of 19 have no RNA-seq.** |
| `selected_17_candidates.tsv` | **17** | `IMPACT == HIGH` **and** `SpliceAI >= 0.5` **and** `VAF >= ~0.041`, minus one manual drop | VAF-maximizing, not testability-maximizing — only 9/17 have RNA-seq |
| `biggene_unsupervised_hits.csv` | 273 | the negative control, see below | |

**`STRICT` and `enriched` are different axes, not nested tiers.** `enriched` is
the annotation axis (everything VEP calls splice-related, with gene structure
attached). `STRICT` is the confidence axis (both tools agree strongly).

### The `priority` column

Reverse-engineered from the original table and verified to reproduce it exactly:

```
priority = 1
         + 2 * (has_RNAseq)        # tissue has a LeafCutter run
         + 2 * (IMPACT == HIGH)    # VEP high impact
         + 3 * (SpliceAI >= 0.8)   # SpliceAI confident
         - 1 * (in_scan)           # already covered by Phase 2
```

Range 1–8. The 49 rows at priority 8 are the practical working set.

---

## Tools

| Tool | Where |
|---|---|
| Ensembl **VEP v113** | upstream annotation (not run here — output is an input) |
| `pyqtl` (Broad internal) | `05_plot_pileup_pyqtl.py` — `qtl.pileup.samtools_depth`, `regtools_extract_junctions(strand="XS")`, `norm_pileups`, `plot(max_intron=100)` |
| `samtools depth` + `idxstats` | `06_plot_pileup_standalone.py` — the dependency-free path |
| `pysam` + **XS-tag strand filter** | junction extraction. Reads without an `XS` tag matching the gene strand are discarded, matching the rigor of Laura's `regtools` call. |
| `python-pptx` | deck assembly (not included here) |

`05` and `06` produce the same figure from the same data. Use `05` on Broad where
`pyqtl` is available; use `06` anywhere else.

---

## Results

### The three case studies, cross-validated

| Case | SpliceAI | VEP | Verdict |
|---|---|---|---|
| **ISCU** | 0.99 | `splice_donor_variant`, **HIGH** | **Both tools agree — the most credible case.** Cryptic donor at chr12:108,567,757. |
| **PLBD2** | **1.00** | **LOW** — variant is in the polypyrimidine tract, deep intronic, not at the canonical ±2 | **Informative discordance.** This is exactly what SpliceAI adds over VEP's hard rules. |
| **C11orf54** | 0.38 | **`missense_variant` (p.G204E)** | **Refuted.** Not a splicing variant. Downgraded to "mechanism unclear". |

C11orf54 is **absent from the 213-row candidate list by construction** — its VEP
consequence contains no `"splice"` substring, so the filter drops it. That
absence *is* the finding.

### Two leads investigated and closed

**Donor 005** (flagged by the PI as a likely smoking gun) — both negative:

- `RSRC2`: splice acceptor −1, VEP HIGH + SpliceAI 1.00. LeafCutter |z| = 0.77,
  |Δψ| = 0.02. **No RNA signal despite a perfect prediction.**
- `DMD`: an unsupervised-scan hit with SpliceAI ≈ 0. Carrier coverage is
  systematically low across the gene — a low-expression-tissue artifact, since
  DMD is barely expressed in skin.

**SLC10A6 tandem-repeat expansion** (following a long-read haplotype finding).
A cryptic acceptor at chr4:86,837,573 shows cohort ψ 0–0.032, with SMHT029
ranked first (0.032) and SMHT039 second (0.018). Looks promising — until
`qc/check_intron_exon_ratio.sh`:

| Donor | intron 1 / exon 1 coverage |
|---|---|
| SMHT029 | **0.12** |
| SMHT039 | **0.11** |
| cohort median | **0.15** |

Both candidates are **below** the cohort. That rejects intron retention. The real
explanation is that SMHT029 expresses SLC10A6 at roughly twice cohort levels in
skin, so all coverage scales up together.

> **Not a splicing phenomenon — an expression phenomenon.** Allele-specific
> expression is the right framework; the repeat acts on transcription, not
> splicing.

### The unsupervised scan as a negative control

`03_unsupervised_biggene_scan.py` drops the SpliceAI filter entirely and scans
four very large genes (CNTNAP2, CSMD1, RBFOX1, DMD) for ψ outliers, then checks
whether the outlier donor carries any somatic variant in that gene.

**273 hits. Maximum SpliceAI across all of them: 0.00. One donor (SMHT029)
accounts for 126 of 273.**

That concentration is the tell — these are expression and QC artifacts, not
splicing events. `qc/check_housekeeping_coverage.sh` exists to adjudicate exactly
this: compare the suspect donor's coverage in GAPDH / ACTB / GUSB against the
gene of interest. Low everywhere → RNA quality problem, discard. Normal in
housekeeping genes but low in the target → real biology.

This result is the direct justification for the prediction-first pivot.

---

## Tier A — where this phase currently stands

Three candidates whose carrier BAMs are already available, no download needed:

| Gene | Donor | Position | VEP | SpliceAI | VAF |
|---|---|---|---|---|---|
| UNC13A | SMHT039 | chr19:17,639,526 | splice_acceptor, HIGH | 1.00 | 0.055 |
| MEI1 | SMHT039 | chr22:41,763,173 | splice_acceptor, HIGH | 0.99 | 0.060 |
| COPB2 | SMHT016 | chr3:139,359,179 | splice_acceptor, HIGH | 0.98 | 0.049 |

This is the intersection of highest-confidence prediction and actual data
availability — the practical answer to the testability gap that
`04_plot_candidates_by_tissue.py` makes visible (hatched bars = no LeafCutter
perind for that tissue).

`07_plot_per_donor_sashimi.py` renders these. **The outputs were never pulled
back from the cluster**, so the conclusion for Tier A is not in this repo.

---

## Paths to edit

| Script | Hardcoded path |
|---|---|
| `01_build_vep_candidate_list.py` | none — fully parameterized via CLI args |
| `02_enrich_candidates_gencode_v47.py` | `GTF`, `INPUT_CSV`, `OUT_CSV` at the top of the file |
| `03_unsupervised_biggene_scan.py` | `ROOT`, `LEAF_DIR`, `V2_FILE` |
| `04_plot_candidates_by_tissue.py` | input CSV path; `SPLICEAI_MIN = 0.5` |
| `05_plot_pileup_pyqtl.py` | `/home/ldomenec/references/gencode.v47.annotation.gtf`, `libsize_<tissue>.txt`, the `CASE` dict |
| `06`, `07` | `/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin/` |
| `qc/*.sh` | `BAM_DIR` |

---

## Known issues

### 1. The entry point was missing

`vep_splice_candidates.csv` had **no producing script** anywhere in the original
project — it was built interactively and never saved. Everything from this phase
onward depends on it.

`01_build_vep_candidate_list.py` is a reconstruction. It was verified to
reproduce the original exactly: 717 variants at SpliceAI ≥ 0.2, 213 rows after
the consequence filter, the 95/82/36 IMPACT breakdown, and the `priority` column
with 100% accuracy.

### 2. Gene-symbol version skew — an active bug

VEP was run against cache **`105_GRCh38`**; the GTF is **v47**. **20 of 213 rows
(12 symbols) get no v47 annotation**:

```
C6orf183   CCBL1    FAM129C   FAM49A    HMHA1     LRRC16A
MSH5-SAPCD1  MUM1   RP11-507M3.1  RP11-545J16.1  SLCO1B7  SSPO
```

Their `mane_select_tid`, `region_v47` and `gene_length_kb` are empty.
Consequences:

- `FAM49A` was **silently dropped** from a 16-candidate batch — it was renamed
  `CYRIA` in v47, so `parse_mane_exons` found nothing.
- **`SSPO` and `MSH5-SAPCD1` — 2 of the 17 selected candidates — cannot be
  plotted as written**, because there is no MANE transcript to draw.

A symbol-alias map would fix this.

### 3. `has_RNAseq` is a tissue-level flag

It records whether the **tissue** has a LeafCutter perind file, not whether that
specific **(donor, tissue)** pair has a BAM. A row can be marked testable and
still be impossible to plot. In the 3AD skin cohort alone, five donors
(SMHT015, SMHT018, SMHT024, SMHT040, SMHT042) carry DNA variants with no 3AD
RNA-seq.

For real testability, query the bucket:

```bash
gsutil ls gs://smaht-p25-bam-files/ | grep "SMHT039-3AD-"
```

### 4. `selected_17` contains duplicate variants

It is really **16 distinct variants**:

- `MSH5` and `MSH5-SAPCD1` are the **same variant** (chr6:31,745,321, SMHT040,
  3A, VAF 0.0635) counted twice under two gene symbols
- `CYP2W1 ×3` is one variant (chr7:987,246, SMHT022) across three brain regions
- `BRCC3 ×2` is one variant across 3A and 3O

### 5. Per-donor sashimi exposes a caveat the slides gloss over

For ISCU, the cryptic junction is **not carrier-exclusive**:

| Donor | cryptic reads |
|---|---|
| **SMHT039 (carrier)** | **15** |
| SMHT007 | 14 |
| SMHT022 | 12 |
| SMHT027 | 7 |

The carrier is highest but not an isolated outlier in raw counts. The |z| = 4.31
comes from ψ normalization, not from a raw excess. Worth stating explicitly in
any writeup.
