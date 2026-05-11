# ------------------------------------------------------------
# 08_paper_assets.R
# Build the final tables and figures consumed by paper/main.tex.
# Outputs:
#   paper/figures/fig_*.pdf
#   paper/tables/tab_*.tex
#   paper/numbers.tex          \newcommand{\NumXxx}{...} for inline numbers
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

mt <- read.csv(file.path(RESULTS_DIR, "04_marginals_table.csv"),
               stringsAsFactors = FALSE)
js <- if (file.exists(file.path(RESULTS_DIR, "05_joint_summary.csv")))
  read.csv(file.path(RESULTS_DIR, "05_joint_summary.csv"),
           stringsAsFactors = FALSE) else NULL
cv <- if (file.exists(file.path(RESULTS_DIR, "07_cv_marginals.csv")))
  read.csv(file.path(RESULTS_DIR, "07_cv_marginals.csv"),
           stringsAsFactors = FALSE) else NULL

# ---- Figure: per-stratum mean CRPS_rel heatmap (sex facet) ----
mean_tbl <- mt |>
  dplyr::group_by(sex, sport, test_type) |>
  dplyr::summarise(mean_CRPS_rel = mean(CRPS_rel, na.rm = TRUE),
                   mean_KS_D     = mean(KS_D, na.rm = TRUE),
                   mean_AIC      = mean(AIC, na.rm = TRUE),
                   n             = max(n),
                   .groups = "drop")

p_heat <- ggplot(mean_tbl, aes(test_type, sport,
                               fill = pmin(mean_CRPS_rel, 5))) +
  geom_tile(color = "white") +
  geom_text(aes(label = formatC(mean_CRPS_rel, format = "g", digits = 2)),
            size = 3) +
  facet_wrap(~ sex) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey90",
                       name = "Mean CRPS / IQR") +
  labs(title = "Marginal fit quality across cohorts",
       subtitle = "Each cell averages 4 metrics; lower is better.",
       x = "Test type", y = "Sport") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"))
ggsave(file.path(PAPER_FIGS, "fig_crps_heatmap.pdf"), p_heat,
       width = 10, height = 5.5)
ggsave(file.path(RESULTS_DIR, "08_crps_heatmap.png"), p_heat,
       width = 10, height = 5.5, dpi = 130)

# ---- Figure: KS_D vs CRPS_rel (calibration vs sharpness) ----
mt_clip <- mt
mt_clip$CRPS_rel_c <- pmin(mt$CRPS_rel, 3)
n_test_levels <- length(unique(mt_clip$test_type))
shape_pal <- (15:25)[seq_len(n_test_levels)]
p_cs <- ggplot(mt_clip, aes(KS_D, CRPS_rel_c, color = sex, shape = test_type)) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
  scale_shape_manual(values = shape_pal, name = "Test") +
  labs(title = "Calibration vs sharpness",
       subtitle = "Lower KS_D = better calibration; lower CRPS = sharper.",
       x = "KS_D on PIT", y = "CRPS / IQR (clipped at 3)")
ggsave(file.path(PAPER_FIGS, "fig_calibration_sharpness.pdf"),
       p_cs, width = 9, height = 5.5)
ggsave(file.path(RESULTS_DIR, "08_calibration_sharpness.png"),
       p_cs, width = 9, height = 5.5, dpi = 130)

# ---- Figure: female vs male CRPS_rel paired comparison ----
paired <- mt |>
  dplyr::group_by(sport, test_type, metric) |>
  dplyr::summarise(
    F = ifelse(any(sex == "Female"), CRPS_rel[sex == "Female"][1], NA_real_),
    M = ifelse(any(sex == "Male"),   CRPS_rel[sex == "Male"][1],   NA_real_),
    .groups = "drop"
  ) |>
  dplyr::filter(!is.na(F), !is.na(M))
if (nrow(paired) > 0) {
  p_pair <- ggplot(paired, aes(M, F, color = sport, shape = test_type)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
    geom_point(size = 2.4, alpha = 0.85) +
    scale_shape_manual(values = shape_pal, name = "Test") +
    coord_equal() +
    labs(title = "Female vs male CRPS / IQR (matched on metric)",
         subtitle = "Above the dashed line: model fits worse for women.",
         x = "Male  CRPS / IQR", y = "Female  CRPS / IQR")
  ggsave(file.path(PAPER_FIGS, "fig_paired_sex.pdf"),
         p_pair, width = 7, height = 5.5)
  ggsave(file.path(RESULTS_DIR, "08_paired_sex.png"),
         p_pair, width = 7, height = 5.5, dpi = 130)
}

# ---- Figure: example multivariable centile panel ----
fits_path <- file.path(RESULTS_DIR, "04_marginals_fits.rds")
example_pdf <- file.path(PAPER_FIGS, "fig_example_centiles.pdf")
if (file.exists(fits_path)) {
  fits <- readRDS(fits_path)
  # Largest-stratum heuristic (e.g. Football CMJ Male)
  pick_score <- vapply(fits, function(x) {
    if (length(x$metrics) == 0) return(0L)
    yv <- x$metrics[[1]]$y
    if (is.null(yv)) 0L else as.integer(length(yv) * length(x$metrics))
  }, integer(1))
  pick <- names(fits)[which.max(pick_score)]
  s <- fits[[pick]]
  if (length(s$metrics)) {
    src <- file.path(RESULTS_DIR, "strata", pick, "centiles.pdf")
    if (file.exists(src)) file.copy(src, example_pdf, overwrite = TRUE)
  }
}

# ---- Tables (LaTeX) ----
top_strata <- mean_tbl |>
  dplyr::arrange(mean_CRPS_rel) |>
  dplyr::slice_head(n = 12)
worst_strata <- mean_tbl |>
  dplyr::arrange(dplyr::desc(mean_CRPS_rel)) |>
  dplyr::slice_head(n = 8)

write_latex_table <- function(df, path, caption, label) {
  cn <- colnames(df)
  cn <- gsub("_", "\\\\_", cn)
  align <- paste(c("l", rep("r", ncol(df) - 1)), collapse = "")
  con <- file(path, "w"); on.exit(close(con))
  cat("\\begin{table}[t]\n\\centering\n\\small\n", file = con)
  cat(sprintf("\\caption{%s}\n\\label{%s}\n", caption, label), file = con)
  cat(sprintf("\\begin{tabular}{%s}\n\\toprule\n", align), file = con)
  cat(paste(cn, collapse = " & "), " \\\\\n\\midrule\n", file = con)
  for (i in seq_len(nrow(df))) {
    row <- vapply(df[i, ], function(v) {
      if (is.numeric(v)) formatC(v, format = "g", digits = 3)
      else as.character(v)
    }, character(1))
    cat(paste(row, collapse = " & "), " \\\\\n", file = con)
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n", file = con)
}

write_latex_table(
  top_strata, file.path(PAPER_TABS, "tab_top_strata.tex"),
  caption = "Cohorts with the best mean fit quality (lower CRPS/IQR is better). $n$ is the largest stratum-level sample size used by any of the four selected metrics.",
  label = "tab:top"
)
write_latex_table(
  worst_strata, file.path(PAPER_TABS, "tab_worst_strata.tex"),
  caption = "Cohorts where marginal fit is most challenged: small $n$ and heavy-tailed metrics dominate.",
  label = "tab:worst"
)

# ---- Inline numbers for the paper ----
fits <- if (file.exists(fits_path)) readRDS(fits_path) else list()
n_athletes <- length(unique(unlist(lapply(fits, function(s)
  unlist(lapply(s$metrics, function(m) length(m$y)))))))
n_strata   <- length(fits)
n_metrics  <- nrow(mt)
n_combos   <- nrow(mean_tbl)
sports     <- length(unique(mt$sport))
tests      <- length(unique(mt$test_type))
mean_crps  <- round(mean(mt$CRPS_rel, na.rm = TRUE), 3)
median_ks  <- round(median(mt$KS_D, na.rm = TRUE), 3)
best_str   <- top_strata$sport[1]
best_test  <- top_strata$test_type[1]
best_sex   <- top_strata$sex[1]
best_val   <- round(top_strata$mean_CRPS_rel[1], 3)
worst_str  <- worst_strata$sport[1]
worst_test <- worst_strata$test_type[1]
worst_sex  <- worst_strata$sex[1]
worst_val  <- round(worst_strata$mean_CRPS_rel[1], 3)

inline_path <- file.path(PROJECT_DIR, "paper", "numbers.tex")
con <- file(inline_path, "w"); on.exit(close(con), add = TRUE)
def <- function(name, value) cat(sprintf("\\newcommand{\\%s}{%s}\n", name, value), file = con)
def("NumStrata",      n_strata)
def("NumCombos",      n_combos)
def("NumMetricRows",  n_metrics)
def("NumSports",      sports)
def("NumTests",       tests)
def("MeanCRPSrel",    sprintf("%.3f", mean_crps))
def("MedianKS",       sprintf("%.3f", median_ks))
def("BestStratum",    sprintf("%s %s (%s)", best_str, best_test, best_sex))
def("BestCRPSrel",    sprintf("%.3f", best_val))
def("WorstStratum",   sprintf("%s %s (%s)", worst_str, worst_test, worst_sex))
def("WorstCRPSrel",   sprintf("%.3f", worst_val))

if (!is.null(js)) {
  js_n <- nrow(js)
  med_aic <- round(median(js$vine_AIC, na.rm = TRUE), 1)
  def("NumVines", js_n)
  def("MedVineAIC", sprintf("%.1f", med_aic))
}

log_step("Paper assets written.")
