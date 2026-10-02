# ============================================================================
# Phase 2 - Per-donor outlier scan (individual test)
#
# FRASER2 / DROP aberrant-splicing module, then cross-validate against the LeafCutter scan.
#
# Original location in the analysis project:
#     0630/code/05_fraser2_standalone.R
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
# ============================================================================
# 05_fraser2_standalone.R
# ----------------------------------------------------------------------------
# Standalone FRASER2 splicing outlier detection for SMaHT P25 cohort
# 3AD_Skin_Calf (12 BAMs, 10 unique donors)
#
# STRICT ADHERENCE TO DROP PIPELINE
# ----------------------------------
# - Uses only FRASER Bioconductor package functions (no modifications)
# - Uses DROP default parameters from
#   https://github.com/gagneurlab/drop/blob/master/drop/template/config.yaml
# - Follows DROP aberrantSplicing module R script order:
#     01_*_countRNA*  -> Sections 3-4  (counting)
#     02_psi_calculation.R -> Section 5
#     03_filter_expression.R -> Section 5
#     04_2_fit_encdim_jaccard.R + 05_fit_finalPCA_jaccard.R -> Section 6
#     06_stats_calculation_jaccard.R + 07_1_results.R -> Section 7
#
# CUSTOM ADDITION
# ---------------
# Section 9: cross-validation with LeafCutter outlier_scan_all.csv
# (explicitly marked with "# CUSTOM:" comments)
#
# Run on Broad:
#   Rscript 05_fraser2_standalone.R
# ============================================================================


# ─── Section 1: Setup ──────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(FRASER)
  library(data.table)
  library(GenomicFeatures)
  library(magrittr)
})

# Paths
BAM_DIR    <- "/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin"
GTF_FILE   <- "/medpop/esp2/yxliu/SMAHT/gencode.v47.annotation.gtf"
WORK_DIR   <- "/medpop/esp2/yxliu/SMAHT/fraser2_output"
LEAFCUTTER_CSV <- "/medpop/esp2/yxliu/SMAHT/02_outlier_scan_all.csv"   # optional
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

# Parameters — ALL from DROP config.yaml default
# https://github.com/gagneurlab/drop/blob/master/drop/template/config.yaml
FRASER_VERSION           <- "FRASER2"
IMPLEMENTATION           <- "PCA"        # DROP default (autoencoder via PCA)
MIN_EXPRESSION_ONE_SAMPLE <- 20          # DROP default
MIN_DELTA_PSI            <- 0.0          # DROP default (filter step)
PADJ_CUTOFF              <- 0.1          # DROP default
DELTA_PSI_CUTOFF         <- 0.1          # DROP default
QUANTILE_FOR_FILTERING   <- 0.75         # DROP default (FRASER2-specific)
MAX_TESTED_DIM_PROPORTION <- 6           # DROP default (autoencoder tuning)

cat("=========================================================\n")
cat("FRASER2 standalone — SMaHT P25 3AD_Skin_Calf\n")
cat("=========================================================\n")
cat("BAM dir:  ", BAM_DIR,  "\n")
cat("GTF file: ", GTF_FILE, "\n")
cat("Work dir: ", WORK_DIR, "\n")
cat("---------------------------------------------------------\n")


# ─── Section 2: Sample annotation table ────────────────────────────────────
# Required by FRASER: sampleID, bamFile
# Required by DROP:   pairedEnd, strand, group, etc.

bam_files <- list.files(BAM_DIR,
                        pattern = "3AD-.*\\.bam$",
                        full.names = TRUE)

sample_table <- data.table(
  sampleID  = tools::file_path_sans_ext(basename(bam_files)),
  bamFile   = bam_files,
  pairedEnd = TRUE,       # STAR-aligned RNA-seq (paired-end)
  strand    = "no",       # unstranded / default DROP setting
  group     = "3AD_Skin_Calf"
)

# Extract donor ID from sampleID
sample_table[, donor := sub("^(SMHT[0-9]+).*", "\\1", sampleID)]

cat("Sample table (", nrow(sample_table), " samples):\n", sep = "")
print(sample_table[, .(sampleID, donor)])

fwrite(sample_table,
       file.path(WORK_DIR, "sample_annotation.tsv"),
       sep = "\t")


# ─── Section 3-4: Counting (DROP 01_*.R equivalents) ───────────────────────
# FRASER internally handles split-read and non-split-read counting when
# countRNAData() is called. This matches DROP's 01_1..01_4 sequence.

cat("\n[Section 3-4] Creating FraserDataSet and counting reads...\n")
cat("(This is the slowest step; takes ~30-90 min for 12 BAMs)\n")

fds <- FraserDataSet(
  colData    = sample_table,
  workingDir = WORK_DIR,
  name       = "SMaHT_P25_3AD_Skin_Calf"
)

# Counting — equivalent to DROP 01_2_countRNA_splitReads.R and 01_3_countRNA_nonSplitReads.R
fds <- countRNAData(fds,
                    minAnchor = 5,        # FRASER default
                    recount   = FALSE)    # DROP default

saveRDS(fds, file.path(WORK_DIR, "fds_after_counting.rds"))
cat("Saved fds after counting.\n")


# ─── Section 5: PSI + Filter (DROP 02, 03) ─────────────────────────────────
cat("\n[Section 5] Computing PSI values and filtering...\n")

# Corresponds to DROP 02_psi_calculation.R
fds <- calculatePSIValues(fds)

# Corresponds to DROP 03_filter_expression.R
fds <- filterExpressionAndVariability(
  fds,
  minExpressionInOneSample = MIN_EXPRESSION_ONE_SAMPLE,   # 20
  minDeltaPsi              = MIN_DELTA_PSI,               # 0.0
  filter                   = TRUE
)

# Keep only junctions passing filter
fds <- fds[mcols(fds, type = "j")[, "passed"], ]

cat("Junctions passing filter:", nrow(fds), "\n")

saveRDS(fds, file.path(WORK_DIR, "fds_after_filter.rds"))


# ─── Section 6: FRASER2 model fit (DROP 04, 05, 06) ────────────────────────
# For FRASER2, the type is "jaccard" (Intron Jaccard Index — single metric)

cat("\n[Section 6] Fitting FRASER2 autoencoder + statistical model...\n")

# 6.1 — Optimize autoencoder dimension q
# Corresponds to DROP 04_2_fit_encdim_jaccard.R
# For small cohorts (n=12), auto-optim may be unstable; try/catch with fallback
best_q <- tryCatch({
  fds <- optimHyperParams(
    fds,
    type                = "jaccard",
    implementation      = IMPLEMENTATION,
    q_param             = seq(2, min(10, ceiling(nrow(sample_table) / MAX_TESTED_DIM_PROPORTION * 2)), by = 1),
    minDeltaPsi         = MIN_DELTA_PSI,
    plot                = FALSE
  )
  bestQ(fds, type = "jaccard")
}, error = function(e) {
  cat("optimHyperParams failed:", conditionMessage(e), "\n")
  cat("Falling back to q = 3 (default for small cohorts)\n")
  3L
})
cat("Best autoencoder dimension q =", best_q, "\n")

# 6.2 — Final FRASER2 fit (autoencoder + beta-binomial p-values)
# Corresponds to DROP 05_fit_finalPCA_jaccard.R and 06_stats_calculation_jaccard.R
fds <- FRASER(
  fds,
  q              = c(jaccard = best_q),
  type           = "jaccard",
  implementation = IMPLEMENTATION
)

saveRDS(fds, file.path(WORK_DIR, "fds_final.rds"))
cat("Saved fitted FRASER2 model.\n")


# ─── Section 7: Annotation + Results (DROP 07) ─────────────────────────────
cat("\n[Section 7] Annotating with gene names and extracting outliers...\n")

# Build TxDb from user's GENCODE v47 GTF
txdb <- makeTxDbFromGFF(GTF_FILE, format = "gtf")

# Annotate junctions with gene names
fds <- annotateRanges(fds, txdb = txdb)

# Extract outliers — equivalent to DROP 07_1_results.R
res <- results(
  fds,
  padjCutoff     = PADJ_CUTOFF,
  deltaPsiCutoff = DELTA_PSI_CUTOFF
)

res_dt <- as.data.table(res)
cat("Outliers passing cutoffs (padj <", PADJ_CUTOFF, ", |ΔJaccard| >",
    DELTA_PSI_CUTOFF, "):", nrow(res_dt), "\n")

fwrite(res_dt, file.path(WORK_DIR, "fraser2_outliers.tsv"), sep = "\t")


# ─── Section 8: FRASER built-in plots ──────────────────────────────────────
cat("\n[Section 8] Generating FRASER built-in plots...\n")

pdf(file.path(WORK_DIR, "fraser2_plots.pdf"), width = 10, height = 7)

# Aberrant events per sample (bar chart)
plotAberrantPerSample(fds,
                      padjCutoff     = PADJ_CUTOFF,
                      deltaPsiCutoff = DELTA_PSI_CUTOFF)

# Q-Q plot — global model calibration
plotQQ(fds, type = "jaccard")

# Volcano plot per sample (all 12 samples)
for (sid in sample_table$sampleID) {
  tryCatch({
    print(plotVolcano(fds, sampleID = sid, type = "jaccard",
                      deltaPsiCutoff = DELTA_PSI_CUTOFF,
                      padjCutoff = PADJ_CUTOFF))
  }, error = function(e) {
    cat("Volcano failed for", sid, ":", conditionMessage(e), "\n")
  })
}

# Sample-sample correlation heatmap (before vs after autoencoder correction)
plotCountCorHeatmap(fds, type = "jaccard", normalized = FALSE)
plotCountCorHeatmap(fds, type = "jaccard", normalized = TRUE)

# Per-outlier expression scatter (top 6 outliers by padj)
if (nrow(res_dt) > 0) {
  top_out <- head(res_dt[order(padjust)], 6)
  for (i in seq_len(nrow(top_out))) {
    tryCatch({
      print(plotExpression(fds,
                           result = res[i],
                           type   = "jaccard"))
    }, error = function(e) {
      cat("Expression plot failed for outlier", i, "\n")
    })
  }
}

dev.off()
cat("Saved plots to fraser2_plots.pdf\n")


# ============================================================================
# CUSTOM ADDITION — NOT part of DROP/FRASER
# ============================================================================
# Cross-validate FRASER2 outliers against our LeafCutter per-donor outlier scan
# ============================================================================

cat("\n[Section 9] CUSTOM: cross-validation with LeafCutter results\n")

if (file.exists(LEAFCUTTER_CSV)) {
  lc <- fread(LEAFCUTTER_CSV)
  # CUSTOM: LeafCutter has sampleID as donor (SMHTxxx); FRASER2 has full sampleID
  # match on (gene, donor)
  res_dt[, donor := sub("^(SMHT[0-9]+).*", "\\1", sampleID)]

  common <- merge(
    res_dt[, .(gene = hgncSymbol, donor, seqnames, start, end,
               fraser_padj = padjust, fraser_delta = deltaPsi,
               fraser_zScore = zScore)],
    lc[, .(gene, donor, tissue, variant_id, SpliceAI, SR_VAF,
           lc_z = abs(z), lc_delta = abs(psi_shift))],
    by = c("gene", "donor")
  )

  cat("LeafCutter hits:                  ", nrow(lc), "\n")
  cat("FRASER2 outliers:                 ", nrow(res_dt), "\n")
  cat("Overlapping (same gene + donor):  ", nrow(common), "\n")

  fwrite(common,
         file.path(WORK_DIR, "fraser2_vs_leafcutter_overlap.tsv"),
         sep = "\t")

  # Specifically check the 3 hero cases
  heroes <- c("PLBD2", "ISCU", "C11orf54")
  hero_check <- res_dt[hgncSymbol %in% heroes]
  cat("\nOriginal 3 hero cases detected by FRASER2?\n")
  if (nrow(hero_check) > 0) {
    print(hero_check[, .(hgncSymbol, sampleID, seqnames, start, end,
                          padjust, deltaPsi, zScore)])
  } else {
    cat("None of PLBD2/ISCU/C11orf54 passed FRASER2 cutoffs.\n")
    cat("Note: this may still be biologically consistent with low VAF.\n")
  }
} else {
  cat("LeafCutter CSV not found at", LEAFCUTTER_CSV, "\n")
  cat("Skipping cross-validation. To enable:\n")
  cat("  scp <local>/02_outlier_scan_all.csv to", LEAFCUTTER_CSV, "\n")
}

cat("\n=========================================================\n")
cat("DONE. Output directory:", WORK_DIR, "\n")
cat("Key files:\n")
cat("  - fraser2_outliers.tsv          : main results\n")
cat("  - fraser2_plots.pdf             : QC + outlier plots\n")
cat("  - fraser2_vs_leafcutter_overlap.tsv : cross-validation\n")
cat("  - fds_final.rds                 : saved FRASER2 object\n")
cat("=========================================================\n")
