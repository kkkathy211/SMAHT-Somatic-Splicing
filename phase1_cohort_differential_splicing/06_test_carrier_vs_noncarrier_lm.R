#!/usr/bin/env Rscript
# ============================================================================
# Phase 1 - Cohort differential splicing (group test)
#
# CORE ENGINE. Per-junction lm: tests A.1 (binary), A.2 (VAF in carriers), A.3 (VAF all donors).
#
# Original location in the analysis project:
#     0520_analysis/code/05_strategy_A_lm.R
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
# =============================================================================
# Strategy A — Per-junction lm tests A.1 / A.2 / A.3
# -----------------------------------------------------------------------------
# For each (tissue, gene) candidate in candidates_A.csv:
#   1. Load PSI matrix from LeafCutter perind.counts.gz
#   2. Filter junctions to clusters mapped to this gene (via A_cluster_to_gene)
#   3. Define groups:
#        carriers       = donors with HC mutation in this gene
#        non_carriers   = donors without
#        VAF_per_donor  = max VAF of HC mutations in this gene (0 if non-carrier)
#   4. Run 3 lm tests per junction:
#        Test A.1 — PSI ~ carrier_flag + center + sex + age + PCs[adaptive]
#        Test A.2 — PSI ~ VAF + ...   (carriers only)
#        Test A.3 — PSI ~ VAF + ...   (all donors, NC=0)
#   5. BH-FDR adjust within each (tissue, gene, test)
#   6. Aggregate stats
#
# Outputs:
#   output/tables/A_summary.csv
#   output/tables/A_top_hits.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(stringr)
})

# ---- paths -----------------------------------------------------------------
# When run via Rscript from code/ dir
ROOT       <- ".."
LEAF        <- file.path(ROOT, "..", "v2_batch_corrected", "NEW",
                          "data", "leafcutter_results")
COV_PATH    <- file.path(ROOT, "..", "v2_batch_corrected", "NEW",
                          "output", "tables", "covariate_matrix.csv")
CAND_PATH   <- file.path(ROOT, "output", "tables", "00_candidates_A.csv")
MUT_PATH    <- file.path(ROOT, "output", "tables", "00_annotated_mutations.csv")
C2G_PATH    <- file.path(ROOT, "output", "tables", "A_cluster_to_gene.csv")
OUT_SUM     <- file.path(ROOT, "output", "tables", "A_summary.csv")
OUT_HITS    <- file.path(ROOT, "output", "tables", "A_top_hits.csv")
LOG_PATH    <- file.path(ROOT, "output", "logs", "05_strategy_A_lm.log")

dir.create(dirname(LOG_PATH), recursive = TRUE, showWarnings = FALSE)
log_con <- file(LOG_PATH, open = "wt")
sink(log_con, split = TRUE)

cat("== Strategy A — Per-junction lm tests A.1 / A.2 / A.3 ==\n\n")
t_global <- Sys.time()

# ---- helpers ---------------------------------------------------------------
read_psi_matrix <- function(tissue_full) {
  cp <- file.path(LEAF, tissue_full, paste0(tissue_full, "_perind.counts.gz"))
  raw <- fread(cmd = sprintf("gunzip -c %s", cp), header = TRUE, sep = " ")
  sample_cols <- setdiff(colnames(raw), "chrom")
  clean <- sub("\\.regtools_junc\\.txt\\.gz$", "", sample_cols)
  setnames(raw, sample_cols, clean)

  parse_xy <- function(v) {
    s <- strsplit(v, "/", fixed = TRUE)
    list(num = as.integer(vapply(s, `[`, character(1), 1)),
         den = as.integer(vapply(s, `[`, character(1), 2)))
  }

  M  <- matrix(NA_real_, nrow(raw), length(clean),
               dimnames = list(raw$chrom, clean))
  Dn <- matrix(NA_integer_, nrow(raw), length(clean),
               dimnames = list(raw$chrom, clean))
  for (j in seq_along(clean)) {
    p <- parse_xy(raw[[clean[j]]])
    ok <- p$den > 0
    M[ok, j] <- p$num[ok] / p$den[ok]
    Dn[, j]  <- p$den
  }
  list(psi = M, denom = Dn)
}

read_splicing_pcs <- function(tissue_full) {
  pp <- file.path(LEAF, tissue_full,
                   paste0(tissue_full, ".leafcutter.PCs.txt"))
  if (!file.exists(pp)) return(NULL)
  pc <- fread(pp, header = TRUE)
  rn <- pc$ID
  pc[, ID := NULL]
  out <- as.data.table(t(as.matrix(pc)))
  setnames(out, rn)
  out[, donor_tissue := colnames(pc)]
  out
}

choose_covar_set <- function(n_samples, n_centers, n_sex_levels,
                              available_pcs = paste0("PC", 1:5),
                              fixed_term, target_df_resid = 3L) {
  fixed_terms_full <- c(fixed_term,
                         if (n_centers >= 2) "center" else NULL,
                         if (n_sex_levels >= 2) "sex" else NULL,
                         "age")
  drop_order <- c(rev(available_pcs), "age", "sex", "center")
  current <- c(fixed_terms_full, available_pcs)
  n_params <- function(terms) {
    p <- 1L
    if (fixed_term %in% terms)  p <- p + 1L
    if ("center" %in% terms)    p <- p + (n_centers - 1L)
    if ("sex" %in% terms)       p <- p + (n_sex_levels - 1L)
    if ("age" %in% terms)       p <- p + 1L
    p <- p + sum(terms %in% paste0("PC", 1:5))
    p
  }
  while (n_params(current) > n_samples - target_df_resid) {
    drop <- intersect(drop_order, current)
    if (length(drop) == 0) break
    current <- setdiff(current, drop[1])
  }
  current
}

safe_extract <- function(m, term) {
  if (is.null(m)) return(c(beta = NA_real_, p = NA_real_))
  co <- suppressWarnings(summary(m)$coefficients)
  if (!is.matrix(co) || !(term %in% rownames(co)))
    return(c(beta = NA_real_, p = NA_real_))
  c(beta = co[term, "Estimate"], p = co[term, "Pr(>|t|)"])
}

# ---- Load global tables ----------------------------------------------------
cat("[1/4] Loading global tables ...\n")
covm <- fread(COV_PATH)
cand <- fread(CAND_PATH)
mut  <- fread(MUT_PATH)
c2g  <- fread(C2G_PATH)
cat(sprintf("   Covariate matrix:  %d rows\n", nrow(covm)))
cat(sprintf("   Candidates A:      %d rows\n", nrow(cand)))
cat(sprintf("   Annotated muts:    %d rows\n", nrow(mut)))
cat(sprintf("   Cluster→gene:      %d rows\n", nrow(c2g)))
cat("\n")

# Build a fast per-(donor, tissue, gene) max VAF lookup
mut_per_dtg <- mut[, .(donor_max_vaf = max(vaf, na.rm = TRUE)),
                   by = .(donor, tissue, gene_name)]

# ---- Main loop over (tissue, gene) candidates ------------------------------
all_results_summary <- list()
all_top_hits        <- list()
psi_cache           <- list()    # cache PSI matrix per tissue
pc_cache            <- list()    # cache PC table per tissue

tissues_processed <- unique(cand$tissue)
cat(sprintf("[2/4] Tissues to process: %s\n", paste(tissues_processed, collapse = ", ")))
cat("[3/4] Running per-(tissue, gene) tests ...\n\n")

n_done <- 0
for (i in seq_len(nrow(cand))) {
  tn <- cand$tissue[i]
  gn <- cand$gene_name[i]
  n_done <- n_done + 1

  # Cache PSI matrix
  if (is.null(psi_cache[[tn]])) {
    psi_cache[[tn]] <- read_psi_matrix(tn)
    pc_cache[[tn]]  <- read_splicing_pcs(tn)
  }
  M_full  <- psi_cache[[tn]]$psi
  Dn_full <- psi_cache[[tn]]$denom
  pcs     <- pc_cache[[tn]]

  # Get clusters mapped to this gene in this tissue
  this_clusters <- c2g[tissue == tn & gene_name == gn, cluster_id]
  if (length(this_clusters) == 0) next

  # All junctions in those clusters
  all_jids <- rownames(M_full)
  jid_clu  <- sub("^.*:(clu_\\d+_[+-])$", "\\1", all_jids)
  in_gene  <- all_jids[jid_clu %in% this_clusters]
  if (length(in_gene) < 2) next

  M  <- M_full[in_gene, , drop = FALSE]
  Dn <- Dn_full[in_gene, , drop = FALSE]
  N  <- ncol(M)

  # Filter junctions: coverage + variance
  ok <- rowSums(Dn >= 5, na.rm = TRUE) >= ceiling(0.80 * N) &
        rowVars(M, na.rm = TRUE) > 1e-6 &
        rowSums(is.finite(M)) >= ceiling(0.80 * N)
  M <- M[ok, , drop = FALSE]
  if (nrow(M) < 2) next
  n_junc_tested <- nrow(M)

  # Per-tissue covariate
  this_cov <- covm[tissue_full == tn]
  this_cov <- this_cov[match(colnames(M), sample_id)]
  this_cov[, center := factor(center)]
  this_cov[, sex    := factor(sex)]
  this_cov[, donor_tissue := paste(donor, tissue_code, sep = "-")]

  # Merge splicing PCs
  pc_cols_avail <- character()
  if (!is.null(pcs)) {
    pc_cols_avail <- intersect(paste0("PC", 1:5), colnames(pcs))
    this_cov <- merge(this_cov, pcs[, c("donor_tissue", pc_cols_avail), with = FALSE],
                       by = "donor_tissue", all.x = TRUE, sort = FALSE)
    this_cov <- this_cov[match(colnames(M), sample_id)]
  }

  # Per-donor carrier / VAF flags for this gene in this tissue
  donor_vaf_tbl <- mut_per_dtg[tissue == tn & gene_name == gn,
                                .(donor, donor_max_vaf)]
  this_cov <- merge(this_cov, donor_vaf_tbl, by = "donor",
                     all.x = TRUE, sort = FALSE)
  this_cov <- this_cov[match(colnames(M), sample_id)]
  this_cov[, vaf := fifelse(is.na(donor_max_vaf), 0, donor_max_vaf)]
  this_cov[, carrier_flag := as.integer(vaf > 0)]

  n_carriers     <- sum(this_cov$carrier_flag == 1)
  n_non_carriers <- sum(this_cov$carrier_flag == 0)
  # Relax to >=3 cores per group — many donors with HC mutation don't have
  # matched RNA-seq cores in the same tissue, so the effective n in the PSI
  # matrix is smaller than the candidate file's donor-level count.
  if (n_carriers < 3 || n_non_carriers < 3) next
  n_carrier_donors <- length(unique(this_cov[carrier_flag == 1, donor]))
  n_total_donors   <- length(unique(this_cov$donor))

  n_centers <- nlevels(droplevels(this_cov$center))
  n_sex     <- nlevels(droplevels(this_cov$sex))

  # ===========================================================
  # Test A.1 — Binary (carrier vs non-carrier), all donors
  # ===========================================================
  chosen_A1 <- choose_covar_set(N, n_centers, n_sex, pc_cols_avail,
                                  fixed_term = "carrier_flag")
  fmla_A1   <- as.formula(paste("psi ~", paste(chosen_A1, collapse = " + ")))
  res_A1 <- vector("list", nrow(M))
  for (j in seq_len(nrow(M))) {
    y <- as.numeric(M[j, ])
    d <- copy(this_cov); d[, psi := y]
    d <- d[is.finite(psi)]
    if (nrow(d) < 8 || length(unique(d$psi)) < 2) {
      res_A1[[j]] <- c(beta = NA_real_, p = NA_real_); next
    }
    fit <- tryCatch(suppressWarnings(lm(fmla_A1, data = d)),
                    error = function(e) NULL)
    res_A1[[j]] <- safe_extract(fit, "carrier_flag")
  }
  beta_A1 <- vapply(res_A1, function(r) r["beta"], numeric(1))
  p_A1    <- vapply(res_A1, function(r) r["p"], numeric(1))
  q_A1    <- p.adjust(p_A1, method = "BH")

  # ===========================================================
  # Test A.2 — Continuous VAF within carriers
  # ===========================================================
  carr_idx <- which(this_cov$carrier_flag == 1)
  if (length(carr_idx) >= 4) {
    M2       <- M[, carr_idx, drop = FALSE]
    cov2     <- this_cov[carr_idx]
    nc_centers <- nlevels(droplevels(cov2$center))
    nc_sex     <- nlevels(droplevels(cov2$sex))
    chosen_A2 <- choose_covar_set(length(carr_idx), nc_centers, nc_sex,
                                    pc_cols_avail, fixed_term = "vaf")
    fmla_A2 <- as.formula(paste("psi ~", paste(chosen_A2, collapse = " + ")))
    res_A2 <- vector("list", nrow(M2))
    for (j in seq_len(nrow(M2))) {
      y <- as.numeric(M2[j, ])
      d <- copy(cov2); d[, psi := y]
      d <- d[is.finite(psi)]
      if (nrow(d) < 8 || sd(d$vaf, na.rm = TRUE) == 0 ||
          length(unique(d$psi)) < 2) {
        res_A2[[j]] <- c(beta = NA_real_, p = NA_real_); next
      }
      fit <- tryCatch(suppressWarnings(lm(fmla_A2, data = d)),
                      error = function(e) NULL)
      res_A2[[j]] <- safe_extract(fit, "vaf")
    }
    beta_A2 <- vapply(res_A2, function(r) r["beta"], numeric(1))
    p_A2    <- vapply(res_A2, function(r) r["p"], numeric(1))
    q_A2    <- p.adjust(p_A2, method = "BH")
  } else {
    beta_A2 <- p_A2 <- q_A2 <- rep(NA_real_, nrow(M))
    chosen_A2 <- character()
  }

  # ===========================================================
  # Test A.3 — Continuous VAF all donors (non-carriers = 0)
  # ===========================================================
  chosen_A3 <- choose_covar_set(N, n_centers, n_sex, pc_cols_avail,
                                  fixed_term = "vaf")
  fmla_A3   <- as.formula(paste("psi ~", paste(chosen_A3, collapse = " + ")))
  res_A3 <- vector("list", nrow(M))
  for (j in seq_len(nrow(M))) {
    y <- as.numeric(M[j, ])
    d <- copy(this_cov); d[, psi := y]
    d <- d[is.finite(psi)]
    if (nrow(d) < 8 || length(unique(d$psi)) < 2 ||
        sd(d$vaf, na.rm = TRUE) == 0) {
      res_A3[[j]] <- c(beta = NA_real_, p = NA_real_); next
    }
    fit <- tryCatch(suppressWarnings(lm(fmla_A3, data = d)),
                    error = function(e) NULL)
    res_A3[[j]] <- safe_extract(fit, "vaf")
  }
  beta_A3 <- vapply(res_A3, function(r) r["beta"], numeric(1))
  p_A3    <- vapply(res_A3, function(r) r["p"], numeric(1))
  q_A3    <- p.adjust(p_A3, method = "BH")

  # ---------------------------------------------------------
  # Summarize per (tissue, gene)
  # ---------------------------------------------------------
  top_A1_idx <- which.min(p_A1); top_A2_idx <- which.min(p_A2); top_A3_idx <- which.min(p_A3)
  jids <- rownames(M)
  cluster_in_gene <- unique(sub("^.*:(clu_\\d+_[+-])$", "\\1", jids))

  summ_row <- data.table(
    tissue            = tn,
    gene              = gn,
    n_carrier_cores   = n_carriers,
    n_non_carr_cores  = n_non_carriers,
    n_carrier_donors  = n_carrier_donors,
    n_total_donors    = n_total_donors,
    n_junctions       = nrow(M),
    n_clusters        = length(cluster_in_gene),
    cluster_ids       = paste(cluster_in_gene, collapse = ";"),

    A1_n_q05          = sum(q_A1 < 0.05, na.rm = TRUE),
    A1_n_q10          = sum(q_A1 < 0.10, na.rm = TRUE),
    A1_top_jid        = if (length(top_A1_idx)) jids[top_A1_idx] else NA_character_,
    A1_top_p          = if (length(top_A1_idx)) p_A1[top_A1_idx]  else NA_real_,
    A1_top_q          = if (length(top_A1_idx)) q_A1[top_A1_idx]  else NA_real_,
    A1_top_beta       = if (length(top_A1_idx)) beta_A1[top_A1_idx] else NA_real_,
    A1_covars         = paste(setdiff(chosen_A1, "carrier_flag"), collapse = "+"),

    A2_n_q05          = sum(q_A2 < 0.05, na.rm = TRUE),
    A2_n_q10          = sum(q_A2 < 0.10, na.rm = TRUE),
    A2_top_jid        = if (length(top_A2_idx)) jids[top_A2_idx] else NA_character_,
    A2_top_p          = if (length(top_A2_idx)) p_A2[top_A2_idx]  else NA_real_,
    A2_top_q          = if (length(top_A2_idx)) q_A2[top_A2_idx]  else NA_real_,
    A2_top_beta       = if (length(top_A2_idx)) beta_A2[top_A2_idx] else NA_real_,
    A2_covars         = paste(setdiff(chosen_A2, "vaf"), collapse = "+"),

    A3_n_q05          = sum(q_A3 < 0.05, na.rm = TRUE),
    A3_n_q10          = sum(q_A3 < 0.10, na.rm = TRUE),
    A3_top_jid        = if (length(top_A3_idx)) jids[top_A3_idx] else NA_character_,
    A3_top_p          = if (length(top_A3_idx)) p_A3[top_A3_idx]  else NA_real_,
    A3_top_q          = if (length(top_A3_idx)) q_A3[top_A3_idx]  else NA_real_,
    A3_top_beta       = if (length(top_A3_idx)) beta_A3[top_A3_idx] else NA_real_,
    A3_covars         = paste(setdiff(chosen_A3, "vaf"), collapse = "+")
  )
  all_results_summary[[length(all_results_summary) + 1]] <- summ_row

  # Top-hits long-format (only q<0.10 in any test)
  cluster_ids_per_junc <- sub("^.*:(clu_\\d+_[+-])$", "\\1", jids)
  any_sig <- (q_A1 < 0.10) | (q_A2 < 0.10) | (q_A3 < 0.10)
  any_sig[is.na(any_sig)] <- FALSE
  if (sum(any_sig) > 0) {
    hits <- data.table(
      tissue = tn, gene = gn,
      junction_id = jids[any_sig],
      cluster_id  = cluster_ids_per_junc[any_sig],
      A1_beta = beta_A1[any_sig], A1_q = q_A1[any_sig],
      A2_beta = beta_A2[any_sig], A2_q = q_A2[any_sig],
      A3_beta = beta_A3[any_sig], A3_q = q_A3[any_sig]
    )
    all_top_hits[[length(all_top_hits) + 1]] <- hits
  }

  if (n_done %% 20 == 0 || n_done == nrow(cand)) {
    cat(sprintf("   [%3d/%3d]  %s/%s  carriers=%d  junc=%d  A1_q05=%d  A2_q05=%d  A3_q05=%d\n",
                n_done, nrow(cand), tn, gn,
                n_carriers, nrow(M),
                summ_row$A1_n_q05, summ_row$A2_n_q05, summ_row$A3_n_q05))
  }
}

# ---- Save outputs ----------------------------------------------------------
cat("\n[4/4] Saving outputs ...\n")
A_summary <- rbindlist(all_results_summary, fill = TRUE)
A_summary <- A_summary[order(-pmax(A1_n_q10, A2_n_q10, A3_n_q10, na.rm = TRUE),
                              A3_top_q, A2_top_q, A1_top_q)]
fwrite(A_summary, OUT_SUM)
cat(sprintf("   Wrote %s (%d rows)\n", OUT_SUM, nrow(A_summary)))

A_top_hits <- rbindlist(all_top_hits, fill = TRUE)
A_top_hits <- A_top_hits[order(tissue, gene, A3_q, A1_q)]
fwrite(A_top_hits, OUT_HITS)
cat(sprintf("   Wrote %s (%d rows)\n", OUT_HITS, nrow(A_top_hits)))

# ---- Quick summary ---------------------------------------------------------
cat("\n=== SUMMARY ===\n")
cat(sprintf("Total (tissue,gene) tested: %d\n", nrow(A_summary)))
cat(sprintf("Pairs with >=1 q<0.10 in any test:  %d\n",
            sum((A_summary$A1_n_q10 > 0) | (A_summary$A2_n_q10 > 0) |
                  (A_summary$A3_n_q10 > 0), na.rm = TRUE)))
cat(sprintf("Pairs with >=1 q<0.05 in any test:  %d\n",
            sum((A_summary$A1_n_q05 > 0) | (A_summary$A2_n_q05 > 0) |
                  (A_summary$A3_n_q05 > 0), na.rm = TRUE)))

cat("\nTop 15 (tissue, gene) by A.3 significance:\n")
top15 <- A_summary[!is.na(A3_top_q)][order(A3_top_q)][1:15]
for (i in seq_len(nrow(top15))) {
  r <- top15[i]
  cat(sprintf("   %-30s %-12s  n_carr=%2d  n_junc=%3d  A3 q05=%2d  A3 top_q=%.2g  A3 top_beta=%.3f\n",
              r$tissue, r$gene, r$n_carrier_cores, r$n_junctions,
              r$A3_n_q05, r$A3_top_q, r$A3_top_beta))
}

cat(sprintf("\nTotal runtime: %.1f min\n",
            as.numeric(Sys.time() - t_global, units = "mins")))
cat("DONE.\n")

sink()
close(log_con)
