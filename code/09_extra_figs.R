# ------------------------------------------------------------
# 09_extra_figs.R
# Additional figures specifically used by the paper:
#   - vine pair-copula contour matrix (example stratum)
#   - dependence comparison: empirical Spearman vs vine-implied
#   - CV-vs-in-sample CRPS_rel scatter (overfitting check)
#   - per-sex calibration density (PIT)
# ------------------------------------------------------------

.this_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(fa)) dirname(fa[1])
  else if (!is.null(sf <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL))) dirname(sf)
  else getwd()
}
source(file.path(.this_dir, "00_setup.R"))
source(file.path(CODE_DIR, "helpers.R"))

fits <- readRDS(file.path(RESULTS_DIR, "04_marginals_fits.rds"))
mt   <- read.csv(file.path(RESULTS_DIR, "04_marginals_table.csv"),
                 stringsAsFactors = FALSE)
cv   <- read.csv(file.path(RESULTS_DIR, "07_cv_marginals.csv"),
                 stringsAsFactors = FALSE)

# 1) CV vs in-sample CRPS_rel (overfitting diagnostic) ---------
mer <- merge(mt[, c("sex", "sport", "test_type", "metric",
                     "CRPS_rel", "KS_D")],
             cv[, c("sex", "sport", "test_type", "metric",
                     "cv_CRPS_rel", "cv_KS_D")],
             by = c("sex", "sport", "test_type", "metric"))
n_test_levels <- length(unique(mer$test_type))
shape_pal <- (15:25)[seq_len(n_test_levels)]
p1 <- ggplot(mer, aes(pmin(CRPS_rel, 3), pmin(cv_CRPS_rel, 3),
                      color = sex, shape = test_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
  scale_shape_manual(values = shape_pal, name = "Test") +
  labs(title = "Out-of-sample vs in-sample CRPS / IQR",
       subtitle = "Points near the diagonal indicate negligible overfitting",
       x = "In-sample CRPS / IQR", y = "5-fold CV CRPS / IQR")
ggsave(file.path(PAPER_FIGS, "fig_cv_vs_insample.pdf"), p1,
       width = 9, height = 5.5)

# 2) Vine dependence: example stratum --------------------------
vine_dirs <- list.dirs(file.path(RESULTS_DIR, "strata"), recursive = FALSE)
big <- vine_dirs[which.max(vapply(vine_dirs, function(d) {
  fp <- file.path(d, "vine.rds")
  if (!file.exists(fp)) return(0L)
  obj <- readRDS(fp); nrow(obj$U)
}, integer(1)))]
if (length(big)) {
  obj <- readRDS(file.path(big, "vine.rds"))
  pdf(file.path(PAPER_FIGS, "fig_vine_pair_contours.pdf"),
      width = 8, height = 8)
  tryCatch({
    par(mar = c(2, 2, 2, 2))
    VineCopula::contour(obj$vine,
                        margins = "norm",
                        cex = 0.6)
    title(sprintf("Pair-copula contours: %s",
                  basename(big)))
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Vine contour plot unavailable:", e$message))
  })
  dev.off()

  emp  <- cor(obj$U, method = "spearman")
  rho  <- 2 * sin(pi * VineCopula::TauMatrix(obj$U) / 6)  # Spearman from Kendall (approx)
  emp_long <- reshape2::melt(emp)
  p_dep <- ggplot(emp_long, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 2.7) +
    scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#cb181d",
                          midpoint = 0, limits = c(-1, 1),
                          name = "Spearman ρ") +
    labs(title = paste("Empirical Spearman dependence:", basename(big)),
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))
  ggsave(file.path(PAPER_FIGS, "fig_dependence_heatmap.pdf"), p_dep,
         width = 8, height = 6)
}

# 3) PIT-density per sex --------------------------------------
pit_rows <- list()
for (sid in names(fits)) {
  s <- fits[[sid]]
  for (mn in names(s$metrics)) {
    pit <- s$metrics[[mn]]$pit
    pit_rows[[length(pit_rows) + 1]] <- data.frame(
      sex = s$sex, pit = pit
    )
  }
}
pit_df <- do.call(rbind, pit_rows)
p_pit <- ggplot(pit_df, aes(pit, fill = sex, color = sex)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
  scale_fill_manual(values  = c(Female = PAL$female, Male = PAL$male)) +
  labs(title = "Probability-integral-transform residuals across all marginals",
       subtitle = "Uniform reference shown as the dashed line at y = 1",
       x = "PIT", y = "Density")
ggsave(file.path(PAPER_FIGS, "fig_pit_density.pdf"), p_pit,
       width = 8, height = 4.5)

cat("[09] extra figures written.\n")
