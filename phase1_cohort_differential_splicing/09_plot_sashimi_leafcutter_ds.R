#!/usr/bin/env Rscript
# ============================================================================
# Phase 1 - Cohort differential splicing (group test)
#
# Sashimi for the DMD / esophagus vignette, carrier vs non-carrier, germline-excluded.
#
# Original location in the analysis project:
#     0611/sashimi_DMD_noGerm_v4.R
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
# Adapted from 09_sashimi_plots.R — generates the same A_sashimi style
# for the NEW DM-sig DMD cluster in 3C_Esophagus (noGerm carrier definition).
#
# Inputs (in this folder):
#   3C_Esophagus_perind.counts.gz   — Laura's re-clustered counts
#   DMD.ANY_noGerm.txt              — 4 carriers vs 9 non_carriers
#   ../0520_analysis/data/gencode.v44.basic.gtf.gz  — exon annotation
#
# Output:
#   A_sashimi_3C_Esophagus_DMD_clu_20135_v4.png

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(reshape2)
})

ROOT       <- "/Users/spectremac/Desktop/Differential_splicing/0611"
COUNTS     <- file.path(ROOT, "3C_Esophagus_perind.counts.gz")
GROUPS     <- file.path(ROOT, "DMD.ANY_noGerm.txt")
GTF        <- "/Users/spectremac/Desktop/Differential_splicing/0520_analysis/data/gencode.v44.basic.gtf.gz"
OUT_DIR    <- ROOT

TISSUE     <- "3C_Esophagus"
GENE       <- "DMD"
TARGET_CLU <- "clu_20135_-"     # NEW DM-sig cluster (q=0.0042)

# ====== plot styling (matches v3 deck) ======
COL_CARR <- "#C53030"
COL_NC   <- "#1F77B4"
COL_EXON <- "#1A2A4F"
COL_INTR <- "#9CA3AF"

theme_sashimi <- theme_void() +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10),
        legend.position = "top",
        axis.title.y  = element_text(angle = 90, size = 10),
        axis.title.x  = element_text(size = 10),
        axis.text     = element_text(size = 9),
        axis.text.y   = element_text(size = 9),
        axis.line.x   = element_line(colour = "grey50", linewidth = 0.3),
        plot.margin   = margin(8, 14, 8, 14))

# ====== helpers ======
arc_points <- function(x1, x2, height, n = 50) {
  cx <- (x1 + x2) / 2
  rx <- abs(x2 - x1) / 2
  t  <- seq(pi, 0, length.out = n)
  data.frame(x = cx + rx * cos(t),
             y = height * sin(pi - t))
}

read_psi_matrix <- function(counts_path) {
  # leafcutter perind.counts.gz format: junction_id sample1 sample2 ...
  # Each cell is "numer/denom"
  raw <- fread(cmd = sprintf("gunzip -c %s", counts_path), header = TRUE, sep = " ")
  sample_cols <- setdiff(colnames(raw), "chrom")
  clean <- sub("\\.regtools_junc\\.txt\\.gz$", "", sample_cols)
  setnames(raw, sample_cols, clean)

  parse_xy <- function(v) {
    s <- strsplit(v, "/", fixed = TRUE)
    list(num = as.integer(vapply(s, `[`, character(1), 1)),
         den = as.integer(vapply(s, `[`, character(1), 2)))
  }
  M <- matrix(NA_real_, nrow(raw), length(clean),
              dimnames = list(raw$chrom, clean))
  for (j in seq_along(clean)) {
    p <- parse_xy(raw[[clean[j]]])
    ok <- p$den > 0
    M[ok, j] <- p$num[ok] / p$den[ok]
  }
  M
}

get_gene_exons <- function(gene_symbol) {
  # Use system grep to filter GTF fast
  cmd <- sprintf(
    "gunzip -c '%s' | awk '$3==\"exon\"' | grep -F 'gene_name \"%s\"'",
    GTF, gene_symbol)
  raw <- tryCatch(
    fread(cmd = cmd, header = FALSE, sep = "\t"),
    error = function(e) NULL)
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  setnames(raw, c("seqnames", "source", "feature", "start", "end",
                   "score", "strand", "frame", "attrs"))
  raw[, .(seqnames = as.character(seqnames),
           start = as.integer(start),
           end   = as.integer(end))]
}

# ====== load data ======
cat("Loading PSI matrix...\n")
M_full <- read_psi_matrix(COUNTS)
cat(sprintf("  %d junctions × %d samples\n", nrow(M_full), ncol(M_full)))

cat("Loading groups...\n")
groups <- fread(GROUPS, header = FALSE, col.names = c("sample", "group"))
cat(sprintf("  groups: %s\n", paste(table(groups$group), names(table(groups$group)),
                                       sep = " ", collapse = ", ")))

# Strip suffix from sample names in groups.txt to match M_full column names
groups[, sample_clean := sub("\\.regtools_junc\\.txt\\.gz$", "", sample)]

# Match groups to PSI matrix columns
sample_to_group <- setNames(groups$group, groups$sample_clean)

cat("Loading GENCODE exons for DMD...\n")
exon_df <- get_gene_exons(GENE)
if (!is.null(exon_df)) {
  cat(sprintf("  %d exon rows (will collapse to union)\n", nrow(exon_df)))
} else {
  cat("  no exon annotation\n")
}

# ====== filter to target cluster ======
jids <- rownames(M_full)
jid_clu <- sub("^.*:(clu_\\d+_[+-])$", "\\1", jids)
in_clu <- jids[jid_clu == TARGET_CLU]
cat(sprintf("Target cluster %s: %d junctions\n", TARGET_CLU, length(in_clu)))

if (length(in_clu) == 0) {
  cat("\nERROR: cluster not found. Available clusters near 20135:\n")
  near <- unique(jid_clu[grepl("clu_201[2-4]", jid_clu)])
  print(head(near, 20))
  stop("cluster not found")
}

# Parse coords
parsed <- regmatches(in_clu, regexec("^([^:]+):(\\d+):(\\d+):.*$", in_clu))
jn_dt <- rbindlist(lapply(seq_along(parsed), function(i) {
  m <- parsed[[i]]
  if (length(m) < 4) return(NULL)
  data.table(jid = in_clu[i],
             chrom = m[2],
             start = as.integer(m[3]),
             end   = as.integer(m[4]))
}))
cat(sprintf("Cluster span: %s:%d-%d\n",
            jn_dt$chrom[1], min(jn_dt$start), max(jn_dt$end)))

# Mean PSI per junction × group
psi_long <- as.data.table(reshape2::melt(
  M_full[in_clu, , drop = FALSE],
  varnames = c("jid", "sample_id"), value.name = "psi"))
psi_long[, sample_id := as.character(sample_id)]
psi_long[, jid       := as.character(jid)]
psi_long[, group := sample_to_group[sample_id]]
psi_long <- psi_long[is.finite(psi) & !is.na(group)]

mean_psi <- psi_long[, .(mean_psi = mean(psi, na.rm = TRUE)),
                       by = .(jid, group)]
mean_psi <- merge(mean_psi, jn_dt[, .(jid, start, end)], by = "jid")
cat(sprintf("\nMean PSI per junction × group:\n"))
print(mean_psi)

# Build arc traces
arcs <- rbindlist(lapply(seq_len(nrow(mean_psi)), function(i) {
  r <- mean_psi[i]
  a <- arc_points(r$start, r$end, r$mean_psi)
  a$jid   <- r$jid
  a$group <- r$group
  a$mean_psi <- r$mean_psi
  a
}))

xmin <- min(jn_dt$start) - 200
xmax <- max(jn_dt$end)   + 200

# ====== Panel A: junction arcs ======
p_arc <- ggplot() +
  geom_path(data = arcs,
            aes(x = x, y = y, group = interaction(jid, group),
                colour = group, alpha = mean_psi,
                linewidth = mean_psi)) +
  scale_colour_manual(values = c(carrier = COL_CARR, non_carrier = COL_NC),
                       name = "Group") +
  scale_alpha_continuous(range = c(0.45, 1), guide = "none") +
  scale_linewidth_continuous(range = c(0.4, 2.2),
                               name = "mean PSI",
                               breaks = c(0.1, 0.3, 0.6, 0.9),
                               limits = c(0, max(arcs$mean_psi))) +
  coord_cartesian(xlim = c(xmin, xmax), ylim = c(0, max(arcs$mean_psi)*1.15),
                   expand = FALSE) +
  labs(title = sprintf("%s in %s — cluster %s (DM noGerm: q=0.0042)",
                        GENE, TISSUE, TARGET_CLU),
       subtitle = sprintf("4 carriers vs 9 non_carriers — sQTL-suspect carrier excluded | %s",
                          "arc height & width = mean PSI"),
       x = NULL, y = "junction PSI") +
  theme_sashimi

# ====== Panel B: exon structure ======
if (!is.null(exon_df) && nrow(exon_df) > 0) {
  exon_df <- exon_df[start <= xmax & end >= xmin]
  if (nrow(exon_df) > 0) {
    exon_iv <- exon_df[order(start)][, .(start, end)]
    merged <- list()
    cur_s <- exon_iv$start[1]; cur_e <- exon_iv$end[1]
    if (nrow(exon_iv) > 1) {
      for (k in 2:nrow(exon_iv)) {
        if (exon_iv$start[k] <= cur_e + 1) {
          cur_e <- max(cur_e, exon_iv$end[k])
        } else {
          merged[[length(merged) + 1]] <- c(cur_s, cur_e)
          cur_s <- exon_iv$start[k]; cur_e <- exon_iv$end[k]
        }
      }
    }
    merged[[length(merged) + 1]] <- c(cur_s, cur_e)
    merged_df <- as.data.frame(do.call(rbind, merged))
    colnames(merged_df) <- c("start", "end")
    junction_ends <- unique(c(jn_dt$start, jn_dt$end))
    merged_df$at_junction <-
      sapply(seq_len(nrow(merged_df)), function(k) {
        any(abs(junction_ends - merged_df$start[k]) < 50) ||
          any(abs(junction_ends - merged_df$end[k]) < 50)
      })

    p_gene <- ggplot() +
      geom_segment(aes(x = xmin, xend = xmax, y = 1, yend = 1),
                   colour = COL_INTR, linewidth = 0.5) +
      geom_rect(data = merged_df,
                 aes(xmin = start, xmax = end,
                     ymin = 0.65, ymax = 1.35,
                     fill = at_junction, colour = at_junction),
                 linewidth = 0.4) +
      scale_fill_manual(values = c(`TRUE` = "#C53030", `FALSE` = COL_EXON),
                         guide = "none") +
      scale_colour_manual(values = c(`TRUE` = "#7A1818", `FALSE` = COL_EXON),
                            guide = "none") +
      coord_cartesian(xlim = c(xmin, xmax), ylim = c(0.3, 1.7), expand = FALSE) +
      labs(x = sprintf("Genomic position on %s  (red exons = junction endpoints)",
                       jn_dt$chrom[1]),
           y = "Gene") +
      theme_sashimi +
      theme(axis.text.y = element_blank())
  } else {
    p_gene <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5,
                label = "(no exons in window)", size = 4)
  }
} else {
  p_gene <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5,
              label = sprintf("(no GENCODE annotation for %s)", GENE), size = 4)
}

# ====== Panel C: per-sample PSI box plot ======
psi_long[, jid_short := sub("^[^:]+:(\\d+):(\\d+):.*$", "\\1-\\2", jid)]
p_box <- ggplot(psi_long, aes(jid_short, psi, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               position = position_dodge(width = 0.7)) +
  geom_jitter(aes(colour = group),
              position = position_jitterdodge(jitter.width = 0.18,
                                                dodge.width = 0.7),
              size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = c(carrier = "#FBCFCF", non_carrier = "#CCE0F5"),
                     guide = "none") +
  scale_colour_manual(values = c(carrier = COL_CARR, non_carrier = COL_NC),
                       guide = "none") +
  labs(x = "junction (start–end)", y = "PSI per sample") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# ====== combine + save ======
combined <- (p_arc / p_gene / p_box) +
  plot_layout(heights = c(2.0, 0.9, 1.3))

out_png <- file.path(OUT_DIR,
                     sprintf("A_sashimi_%s_%s_%s_v4.png",
                             TISSUE, GENE, sub("[+-]$", "", TARGET_CLU)))
ggsave(out_png, combined, width = 12, height = 8, dpi = 150)
cat(sprintf("\nWrote: %s\n", out_png))
