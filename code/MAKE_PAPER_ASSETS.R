# =============================================================================
#                       MAKE_PAPER_ASSETS.R  (v2)
# Builds every figure and table used in paper/main.tex from the RDS/CSV files
# produced by RUN_PAPER_FINAL.R.
#
# Design choices for v2 (per author request):
#   - Figures saved as 300-dpi JPEG (not PDF). Cleaner, less "AI generated".
#   - Muted blue/grey palette only.  No rainbow/viridis.
#   - Larger axis text and tile labels so numbers are readable in print.
#   - Single panels by default (no busy facet/subfig cocktails in body figures).
# =============================================================================

PROJECT_DIR <- normalizePath("~/Downloads/VALD_Multivariate_Centiles_STAIX2026",
                              mustWork = TRUE)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(gridExtra); library(scales); library(fmsb)
  library(scoringRules); library(VineCopula)
  library(gamlss); library(gamlss.dist); library(reshape2)
})

# ---- shared theme ----------------------------------------------------------
theme_paper <- function() {
  theme_minimal(base_size = 12, base_family = "") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          axis.title       = element_text(colour = "grey20"),
          axis.text        = element_text(colour = "grey25"),
          strip.text       = element_text(colour = "grey15", face = "bold"),
          plot.title       = element_text(face = "bold", colour = "grey15"),
          legend.position  = "bottom",
          legend.title     = element_text(colour = "grey25"),
          legend.text      = element_text(colour = "grey25"))
}
theme_set(theme_paper())

BLUE_DARK   <- "#1F4F7F"
BLUE_MED    <- "#3F7AAD"
BLUE_LIGHT  <- "#9DC3E6"
GREY_DARK   <- "#444444"
GREY_LIGHT  <- "#CFCFCF"
ORANGE      <- "#C6612A"
SEX_PAL     <- c(Female = "#A93226", Male = BLUE_DARK)

JPG <- function(path, plot, w = 6.5, h = 4.0, dpi = 300) {
  ggsave(path, plot = plot, width = w, height = h, dpi = dpi,
         units = "in", device = "jpeg")
}

dir.create("paper/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("paper/tables",  showWarnings = FALSE, recursive = TRUE)

# ---- Wipe stale assets (R4) -------------------------------------------------
# Asset-generation must be authoritative. Old figure/table files left over
# from prior pipeline iterations are deleted so that pdflatex never silently
# pulls a stale artefact. The audit script (99_check_paper_assets.R) is the
# enforcement layer; this is the cleanup layer.
cat("Wiping stale paper assets...\n")
.stale_dirs <- list(
  "paper/figures" = c("\\.jpg$", "\\.jpeg$", "\\.png$", "\\.pdf$"),
  "paper/tables"  = c("\\.tex$")
)
.n_wiped <- 0L
for (.d in names(.stale_dirs)) {
  for (.pat in .stale_dirs[[.d]]) {
    .ff <- list.files(.d, pattern = .pat, full.names = TRUE)
    if (length(.ff)) {
      file.remove(.ff)
      .n_wiped <- .n_wiped + length(.ff)
    }
  }
}
cat(sprintf("  -> wiped %d files from paper/figures and paper/tables\n",
             .n_wiped))

# ---- per-artefact write logger ---------------------------------------------
wrote <- function(path) {
  fi <- file.info(path)
  sz <- if (is.finite(fi$size)) {
    if (fi$size > 1024*1024)
      sprintf("%.1f MB", fi$size/1024/1024)
    else if (fi$size > 1024)
      sprintf("%.1f KB", fi$size/1024)
    else sprintf("%d B", as.integer(fi$size))
  } else "?"
  cat(sprintf("  -> wrote %s (%s)\n", path, sz))
}

# ---- artefacts -------------------------------------------------------------
data_all <- read.csv("data/analysis_dataset.csv",
                      stringsAsFactors = FALSE, check.names = FALSE)
g <- toupper(trimws(as.character(data_all$gender)))
data_all$gender <- ifelse(g %in% c("M","MALE"), "Male",
                   ifelse(g %in% c("F","FEMALE"), "Female", NA))
data_all$sport <- trimws(gsub("\\s+(Overall(\\s+Data)?)$", "",
                              data_all$sport, ignore.case = TRUE))
data_all <- data_all[!is.na(data_all$gender) & !is.na(data_all$sport) &
                       !is.na(data_all$age_at_test) &
                       data_all$age_at_test >= 8 & data_all$age_at_test <= 60, ]
data_all$bm <- suppressWarnings(as.numeric(data_all$weight))

cohort_fits   <- readRDS("results/final/cohort_fits.rds")
margin_tab    <- read.csv("results/final/marginal_competition.csv",
                            stringsAsFactors = FALSE)
joint_tab     <- read.csv("results/final/joint_competition.csv",
                            stringsAsFactors = FALSE)
cohort_tab    <- read.csv("results/final/cohort_table.csv",
                            stringsAsFactors = FALSE)
rel_tab       <- read.csv("results/final/reliability_table.csv",
                            stringsAsFactors = FALSE)
rel_summary   <- read.csv("results/final/reliability_summary.csv",
                            stringsAsFactors = FALSE)

# =============================================================================
# 1) EDA — counts heatmap (muted, dark text on light fill)
# =============================================================================
cat("EDA figures...\n")
eligible <- cohort_tab %>% transmute(sex, sport, test, n)

p_counts <- ggplot(eligible, aes(sport, test, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n), size = 3.0, colour = "grey15", fontface = "bold") +
  scale_fill_gradient(low = "#EEF3F8", high = BLUE_DARK,
                       trans = "sqrt", name = expression("n   ")) +
  facet_wrap(~ sex) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid  = element_blank(),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(1.2, "cm"))
JPG("paper/figures/fig_eda_counts.jpg", p_counts, w = 7.0, h = 3.0)

# 1b) Age distribution by athlete
ages <- data_all %>% group_by(profile_id) %>% slice(1) %>% ungroup()
p_age <- ggplot(ages, aes(age_at_test, fill = gender)) +
  geom_histogram(bins = 28, position = "identity",
                  alpha = 0.65, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = SEX_PAL, name = NULL) +
  labs(x = "Age (years)", y = "Athletes") +
  theme(legend.position = c(0.85, 0.85))
JPG("paper/figures/fig_eda_age.jpg", p_age, w = 4.5, h = 2.8)

# 1c) Cumulative tests
cum <- data_all %>%
  mutate(test_date = as.Date(test_date)) %>%
  filter(!is.na(test_date)) %>%
  group_by(test_date) %>% summarise(n = n_distinct(test_id)) %>%
  arrange(test_date) %>% mutate(cum = cumsum(n))
p_cum <- ggplot(cum, aes(test_date, cum)) +
  geom_step(colour = BLUE_DARK, linewidth = 0.9) +
  labs(x = "Date", y = "Cumulative tests")
JPG("paper/figures/fig_eda_cumulative.jpg", p_cum, w = 4.5, h = 2.8)

# 1d) Missingness map (panel availability)
panel_slugs <- list(
  CMJ=c("JumpHeight","PeakPower_BM","mRSI","EccMeanForce"),
  DJ =c("RSI_FTCT","ContactTime","JumpHeight","PeakLandF_BW"),
  SJ =c("JumpHeight","PeakPower_BM","ConcMF_BM","ForceAtPP"),
  SLJ=c("ConcImp_R","ConcImp_L","ConcMF_R","ConcMF_L"),
  SLDJ=c("ConcImp_R","ConcImp_L","MeanLandPow_R","MeanLandPow_L"),
  SLLAH=c("DropLand_R","DropLand_L","PeakDLF_R"),
  SLHAR=c("PeakTakeoffF_R","PeakTakeoffF_L","PeakLandF_R","PeakLandF_L"))
miss_rows <- list()
for (cf in cohort_fits) {
  panel <- panel_slugs[[cf$test]]; if (is.null(panel)) next
  for (m in panel) {
    miss_rows[[length(miss_rows)+1]] <- data.frame(
      cohort = sprintf("%s %s", cf$sex, cf$sport),
      test = cf$test, metric = m,
      have = m %in% names(cf$panel_metrics),
      stringsAsFactors = FALSE)
  }
}
miss <- do.call(rbind, miss_rows)
p_miss <- ggplot(miss, aes(metric, cohort, fill = have)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = BLUE_MED, `FALSE` = GREY_LIGHT),
                     labels = c("Missing", "Available"), name = NULL) +
  facet_wrap(~ test, scales = "free", nrow = 2) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)
JPG("paper/figures/fig_eda_missing.jpg", p_miss, w = 7.0, h = 5.0)

# 1e) JH vs age
keep_col <- "Jump.Height..Flight.Time._Both_CMJ"
ddd <- data_all[, c("gender","age_at_test", keep_col)]
ddd$y <- suppressWarnings(as.numeric(ddd[[keep_col]]))
ddd <- ddd[is.finite(ddd$y) & is.finite(ddd$age_at_test), ]
p_jh <- ggplot(ddd, aes(age_at_test, y, colour = gender)) +
  geom_point(alpha = 0.30, size = 0.7) +
  geom_smooth(formula = y ~ x, method = "loess", se = FALSE, linewidth = 1.0) +
  scale_colour_manual(values = SEX_PAL, name = NULL) +
  labs(x = "Age (years)", y = "CMJ jump height (cm)") +
  theme(legend.position = c(0.85, 0.15))
JPG("paper/figures/fig_eda_jh_age.jpg", p_jh, w = 4.5, h = 2.8)

# =============================================================================
# 2) Worked-example figures (Male Football CMJ) — large readable text
# =============================================================================
cat("Worked-example figures...\n")
ex <- cohort_fits[["Male_Football_CMJ"]]; stopifnot(!is.null(ex))

predict_quantile_grid <- function(obj, x_new, bm_new = NULL,
                                    ps = c(0.05, 0.25, 0.50, 0.75, 0.95)) {
  if (obj$type == "Linear") {
    res <- residuals(obj$fit)
    mu  <- if (!is.null(bm_new))
             predict(obj$fit, newdata = data.frame(x = x_new, bm = bm_new))
           else predict(obj$fit, newdata = data.frame(x = x_new))
    s   <- stats::sd(res)
    out <- sapply(ps, function(p) mu + qnorm(p) * s)
    if (length(x_new) == 1) out <- matrix(out, nrow = 1)
    colnames(out) <- paste0("q", sprintf("%02d", round(ps * 100)))
    return(out)
  }
  fam <- as.character(obj$fit$family)[1]
  qfun <- get(paste0("q", fam), envir = asNamespace("gamlss.dist"),
              mode = "function")
  newdat <- if (!is.null(bm_new)) data.frame(x = x_new, bm = bm_new)
            else data.frame(x = x_new)
  td <- attr(obj$fit, "training_data")
  pa <- suppressWarnings(predictAll(obj$fit, newdata = newdat,
                                      type = "response", data = td))
  pn <- intersect(c("mu","sigma","nu","tau"), names(pa))
  out <- sapply(ps, function(p) do.call(qfun, c(list(p = p), pa[pn])))
  if (length(x_new) == 1) out <- matrix(out, nrow = 1)
  colnames(out) <- paste0("q", sprintf("%02d", round(ps * 100)))
  out
}
predict_cdf_xs <- function(obj, y_new, x_new, bm_new = NULL) {
  if (obj$type == "Linear") {
    res <- residuals(obj$fit)
    mu  <- if (!is.null(bm_new))
             predict(obj$fit, newdata = data.frame(x = x_new, bm = bm_new))
           else predict(obj$fit, newdata = data.frame(x = x_new))
    return(as.numeric(stats::pnorm((y_new - mu) / stats::sd(res))))
  }
  fam <- as.character(obj$fit$family)[1]
  pfun <- get(paste0("p", fam), envir = asNamespace("gamlss.dist"),
              mode = "function")
  newdat <- if (!is.null(bm_new)) data.frame(x = x_new, bm = bm_new)
            else data.frame(x = x_new)
  td <- attr(obj$fit, "training_data")
  pa <- suppressWarnings(predictAll(obj$fit, newdata = newdat,
                                      type = "response", data = td))
  pn <- intersect(c("mu","sigma","nu","tau"), names(pa))
  as.numeric(do.call(pfun, c(list(q = y_new), pa[pn])))
}

ages_grid <- seq(min(ex$panel_metrics[[1]]$x),
                 max(ex$panel_metrics[[1]]$x), length.out = 50)
plots <- lapply(names(ex$panel_metrics), function(m) {
  pm <- ex$panel_metrics[[m]]
  bm_ref <- if (pm$mass_dep) rep(median(pm$bm, na.rm = TRUE), length(ages_grid)) else NULL
  q_mat <- predict_quantile_grid(pm$obj, ages_grid, bm_ref)
  bands <- data.frame(age = ages_grid,
                       q05 = q_mat[, "q05"], q25 = q_mat[, "q25"],
                       q50 = q_mat[, "q50"], q75 = q_mat[, "q75"],
                       q95 = q_mat[, "q95"])
  obs <- data.frame(age = pm$x, y = pm$y)
  ggplot() +
    geom_ribbon(data = bands, aes(age, ymin = q05, ymax = q95),
                fill = BLUE_LIGHT, alpha = 0.55) +
    geom_ribbon(data = bands, aes(age, ymin = q25, ymax = q75),
                fill = BLUE_MED,   alpha = 0.55) +
    geom_line  (data = bands, aes(age, q50),
                colour = BLUE_DARK, linewidth = 0.9) +
    geom_point (data = obs,  aes(age, y),
                alpha = 0.30, size = 0.7, colour = GREY_DARK) +
    labs(title = m, x = "Age (years)", y = pm$units) +
    theme(plot.title = element_text(size = 11, hjust = 0))
})
g_ex <- gridExtra::arrangeGrob(grobs = plots, nrow = 1)
JPG("paper/figures/fig_centiles_example.jpg", g_ex, w = 7.0, h = 2.6)

# 2c) Larger radar chart for one chosen athlete
ex_id <- "Male_Football_CMJ"
ex_cf <- cohort_fits[[ex_id]]
U <- ex_cf$joint_aligned$U
cm <- colMeans(U); cv <- cov(U)
d2 <- as.numeric(stats::mahalanobis(U, cm, cv + diag(1e-8, ncol(U))))
chosen_idx <- which.min(abs(d2 - quantile(d2, 0.9)))   # an interesting "stretched" athlete
chosen_pid <- ex_cf$joint_aligned$profile_id[chosen_idx]

make_radar <- function(athlete_id, cf, file) {
  rec <- data_all[data_all$profile_id == athlete_id & data_all$gender == cf$sex &
                    data_all$sport == cf$sport, ]
  panel_cols <- vapply(cf$panel_metrics, function(m) m$col, character(1))
  rec$cov <- rowSums(!is.na(rec[, panel_cols, drop = FALSE]))
  rec <- rec[order(-rec$cov, as.Date(rec$test_date), decreasing = c(FALSE, TRUE)), ][1, ]
  vals <- vapply(cf$panel_metrics, function(m) {
    yval <- suppressWarnings(as.numeric(rec[[m$col]]))
    if (!is.finite(yval)) return(NA_real_)
    100 * predict_cdf_xs(m$obj, yval, rec$age_at_test,
                          if (m$mass_dep) rec$bm else NULL)
  }, numeric(1))
  vals[is.na(vals)] <- 50
  pdat <- rbind(rep(100, length(vals)), rep(0, length(vals)), vals)
  pdat <- as.data.frame(pdat); colnames(pdat) <- names(cf$panel_metrics)
  jpeg(file, width = 6.5, height = 5.5, units = "in", res = 300, quality = 92)
  par(mar = c(2, 4, 3, 4))
  fmsb::radarchart(pdat, axistype = 1,
                    pcol = BLUE_DARK,
                    pfcol = scales::alpha(BLUE_MED, 0.45),
                    plwd = 2.5,
                    cglcol = "grey75", cglty = 1, cglwd = 0.6,
                    axislabcol = "grey40",
                    caxislabels = c("0", "25", "50", "75", "100"),
                    vlcex = 1.05, calcex = 0.95,
                    title = sprintf("Athlete radar (%s %s %s)",
                                     cf$sex, cf$sport, cf$test))
  dev.off()
  athlete_id
}
chosen_pid <- make_radar(chosen_pid, ex_cf, "paper/figures/fig_radar_example.jpg")

# 2d) PIT histograms for the worked example
pit_long <- list()
for (m in names(ex_cf$panel_metrics)) {
  pm <- ex_cf$panel_metrics[[m]]
  pit_vals <- predict_cdf_xs(pm$obj, pm$y, pm$x,
                               if (pm$mass_dep) pm$bm else NULL)
  pit_long[[m]] <- data.frame(metric = m, pit = pit_vals)
}
pit_df <- do.call(rbind, pit_long)
calp <- ggplot(pit_df, aes(pit)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1),
                  fill = BLUE_MED, colour = "white", linewidth = 0.4) +
  geom_hline(yintercept = nrow(pit_df) /
                            length(unique(pit_df$metric)) / 10,
             linetype = "22", colour = "grey40", linewidth = 0.5) +
  facet_wrap(~ metric, nrow = 1) +
  labs(x = "Probability Integral Transform (PIT)", y = "Count")
JPG("paper/figures/fig_calibration_example.jpg", calp, w = 7.0, h = 2.4)

# =============================================================================
# 3) (retired) Marginal CRPS heatmap and joint-winner barplot (R2).
# These two figures were dropped per author feedback: their information is
# now in tab_calibration_scores.tex (per-cohort C_M, C_J) and the prose
# text uses crps_df medians and joint_tab winner counts, both via
# numbers.tex macros. The crps_df summary is kept because numbers.tex
# consumes it for the CRPSmed macro.
# =============================================================================
crps_df <- margin_tab %>%
  rowwise() %>%
  mutate(best_score = min(c(Linear, NO, BCT), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(sex, sport, test) %>%
  summarise(mean_best = round(mean(best_score, na.rm = TRUE), 3),
             .groups = "drop")

# =============================================================================
# 5) Tables (LaTeX)
# Note: tab_cohorts, tab_marginal_competition, tab_joint_competition were
# retired (R4): the paper no longer references them. Their information is
# now in tab_calibration_scores, tab_fit_quality and tab_variable_selection.
# Only the helper used by the surviving tables is kept here.
# =============================================================================
cat("Tables...\n")

xtab_escape <- function(x) gsub("([%&#_])", "\\\\\\1", as.character(x))

# 5d) Reliability summary
rs <- rel_summary %>%
  rename(Test = test, Metric = metric,
          MedICC = median_ICC, MedMDC95 = median_MDC95, K = n_cohorts) %>%
  mutate(across(everything(), xtab_escape))
cat(c("\\begin{tabular}{llrrr}\n", "\\toprule\n",
      "Test & Metric & median ICC$_{(3,1)}$ & median MDC$_{95}$ & $k$ cohorts \\\\\n",
      "\\midrule\n",
      paste(apply(rs, 1, function(r) paste(r, collapse = " & ")), "\\\\\n"),
      "\\bottomrule\n", "\\end{tabular}\n"),
    file = "paper/tables/tab_reliability.tex", sep = "")

# 5e) score_athlete output table for Male Football CMJ — concrete utility
cat("score_athlete() output table...\n")
predict_cdf <- predict_cdf_xs

score_one <- function(athlete_id, cf) {
  rec <- data_all[data_all$profile_id == athlete_id & data_all$gender == cf$sex &
                    data_all$sport == cf$sport, ]
  panel_cols <- vapply(cf$panel_metrics, function(m) m$col, character(1))
  rec$cov <- rowSums(!is.na(rec[, panel_cols, drop = FALSE]))
  rec <- rec[order(-rec$cov, as.Date(rec$test_date), decreasing = c(FALSE, TRUE)), ][1, ]
  age <- rec$age_at_test; bm <- rec$bm
  rows <- list(); uvec <- numeric(0); names_u <- character(0)
  for (m in cf$panel_metrics) {
    yval <- suppressWarnings(as.numeric(rec[[m$col]]))
    if (!is.finite(yval)) next
    cdf <- predict_cdf(m$obj, yval, age, if (m$mass_dep) bm else NULL)
    rel <- rel_tab[rel_tab$sex == cf$sex & rel_tab$sport == cf$sport &
                    rel_tab$test == cf$test & rel_tab$metric == m$slug, ]
    mdc <- if (nrow(rel)) rel$MDC95[1] else NA_real_
    rows[[length(rows) + 1]] <- data.frame(
      metric = m$slug, value = signif(yval, 4),
      units  = m$units, model = m$best,
      pct    = round(100 * cdf, 1),
      MDC95  = signif(mdc, 3),
      stringsAsFactors = FALSE)
    if (m$slug %in% cf$joint_metrics) {
      uvec <- c(uvec, cdf); names_u <- c(names_u, m$slug)
    }
  }
  out <- do.call(rbind, rows)
  joint <- NA_real_; lo <- NA_real_; hi <- NA_real_
  jm <- cf$joint_model
  if (!is.null(jm) && length(uvec) >= 2) {
    common <- intersect(names_u, colnames(jm$U))
    if (length(common) >= 2) {
      U <- jm$U[, common, drop = FALSE]; ux <- uvec[match(common, names_u)]
      depth_pct <- function(U, x) {
        cm <- colMeans(U); cv <- stats::cov(U) + diag(1e-8, ncol(U))
        d2_a <- as.numeric(stats::mahalanobis(matrix(x, 1), cm, cv))
        d2_c <- as.numeric(stats::mahalanobis(U, cm, cv))
        100 * (1 - mean(d2_c <= d2_a))
      }
      joint <- depth_pct(U, ux)
      bs <- replicate(200, {
        idx <- sample.int(nrow(U), nrow(U), replace = TRUE)
        depth_pct(U[idx, , drop = FALSE], ux)
      })
      lo <- as.numeric(stats::quantile(bs, 0.025, na.rm = TRUE))
      hi <- as.numeric(stats::quantile(bs, 0.975, na.rm = TRUE))
    }
  }
  list(per_metric = out, joint_pct = joint, joint_lo = lo, joint_hi = hi,
       age = age, bm = bm, jm = jm$type, athlete_id = athlete_id)
}

scored <- score_one(chosen_pid, ex_cf)
band_text <- function(p) {
  if (is.na(p)) return("---")
  if (p < 25)  "Below cohort"
  else if (p < 50) "Below median"
  else if (p < 75) "Above median"
  else if (p < 90) "Above cohort"
  else "Elite tail"
}
rows <- scored$per_metric
rows$band <- vapply(rows$pct, band_text, character(1))
tex_safe <- function(s) gsub("_", "\\\\_", as.character(s))
tex_lines <- c(
  "\\begin{tabular}{lrlrrl}",
  "\\toprule",
  "Metric & Value & Units & Pct & MDC$_{95}$ & Band \\\\",
  "\\midrule",
  apply(rows, 1, function(r) {
    sprintf("\\texttt{%s} & %s & %s & %s & %s & %s \\\\",
            tex_safe(r["metric"]), r["value"], r["units"],
            r["pct"], r["MDC95"], r["band"])
  }),
  "\\midrule",
  sprintf("\\multicolumn{6}{l}{\\textbf{Joint percentile (%s):} %.1f\\%% \\quad 95\\%% bootstrap CI [%.1f, %.1f]} \\\\",
          tolower(scored$jm), scored$joint_pct, scored$joint_lo, scored$joint_hi),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(tex_lines, "paper/tables/tab_score_example.tex")

# Also save a brief textual context for the paper to reference
writeLines(c(
  sprintf("\\newcommand{\\AthleteAge}{%.1f}", scored$age),
  sprintf("\\newcommand{\\AthleteBM}{%.1f}",  scored$bm),
  sprintf("\\newcommand{\\AthleteJoint}{%.1f}", scored$joint_pct),
  sprintf("\\newcommand{\\AthleteJointLo}{%.1f}", scored$joint_lo),
  sprintf("\\newcommand{\\AthleteJointHi}{%.1f}", scored$joint_hi),
  sprintf("\\newcommand{\\AthleteJM}{%s}", scored$jm)
), "paper/numbers_athlete.tex", sep = "\n")

# =============================================================================
# 6) numbers.tex
# =============================================================================
n_athletes <- length(unique(data_all$profile_id))
n_tests    <- length(unique(data_all$test_id))
mean_age   <- round(mean(data_all$age_at_test, na.rm = TRUE), 1)
sd_age     <- round(stats::sd(data_all$age_at_test, na.rm = TRUE), 1)
n_female   <- length(unique(data_all$profile_id[data_all$gender == "Female"]))
n_male     <- length(unique(data_all$profile_id[data_all$gender == "Male"]))
n_sports   <- length(unique(data_all$sport))
n_cohorts  <- nrow(cohort_tab)
n_joint    <- sum(cohort_tab$joint_model %in% c("Vine","Gaussian","Indep"))
n_vine     <- sum(joint_tab$Best == "Vine")
n_gauss    <- sum(joint_tab$Best == "Gaussian")
n_indep    <- sum(joint_tab$Best == "Indep")
mean_gain  <- round(mean(joint_tab$GainOverIndep_pct, na.rm = TRUE), 1)
crps_med   <- round(median(crps_df$mean_best, na.rm = TRUE), 3)

best_table <- table(margin_tab$Best)
n_lin <- as.integer(best_table["Linear"]); n_lin <- if (is.na(n_lin)) 0 else n_lin
n_NO  <- as.integer(best_table["NO"]);     n_NO  <- if (is.na(n_NO))  0 else n_NO
n_BCT <- as.integer(best_table["BCT"]);    n_BCT <- if (is.na(n_BCT)) 0 else n_BCT

# ---- Reliability / variable-selection headline numbers (V1) ----------------
n_sports_cohort <- length(unique(cohort_tab$sport))   # sports with >=1 surviving cohort
icc_med  <- round(median(rel_tab$ICC,   na.rm = TRUE), 2)
icc_q1   <- round(quantile(rel_tab$ICC, 0.25, na.rm = TRUE), 2)
icc_q3   <- round(quantile(rel_tab$ICC, 0.75, na.rm = TRUE), 2)
icc_pass <- round(100 * mean(rel_tab$ICC >= 0.5, na.rm = TRUE), 0)

# variable-selection summary: median, range, fraction retired
nsel_per_cohort <- sapply(cohort_fits,
                            function(c) length(c$joint_metrics))
nsel_med    <- as.integer(median(nsel_per_cohort))
n_retired   <- sum(nsel_per_cohort == 0)
n_kept_full <- sum(nsel_per_cohort == 4)

# CMJ patterns (most-frequent metric retained)
cmj_keys <- grep("_CMJ$", names(cohort_fits), value = TRUE)
cmj_metrics <- unlist(lapply(cohort_fits[cmj_keys],
                              function(c) c$joint_metrics))
cmj_top <- if (length(cmj_metrics)) {
  tt <- sort(table(cmj_metrics), decreasing = TRUE); names(tt)[1]
} else "n/a"
cmj_top_freq <- if (length(cmj_metrics)) {
  as.integer(sort(table(cmj_metrics), decreasing = TRUE)[1])
} else 0L
cmj_n_cohorts <- length(cmj_keys)

# Best / worst joint cohort by winning C_J
joint_tab$Cwin <- pmin(joint_tab$Indep,
                        joint_tab$Gaussian,
                        joint_tab$Vine,
                        na.rm = TRUE)
.best  <- joint_tab[which.min(joint_tab$Cwin), ]
.worst <- joint_tab[which.max(joint_tab$Cwin), ]
.lab   <- function(r) sprintf("%s %s %s", r$sex, r$sport, r$test)
best_cohort_lab  <- .lab(.best)
worst_cohort_lab <- .lab(.worst)
best_cohort_cwin  <- round(.best$Cwin, 2)
worst_cohort_cwin <- round(.worst$Cwin, 2)

cat(sprintf(
"\\newcommand{\\Nathletes}{%d}
\\newcommand{\\Ntests}{%d}
\\newcommand{\\MeanAge}{%.1f}
\\newcommand{\\SdAge}{%.1f}
\\newcommand{\\Nfemale}{%d}
\\newcommand{\\Nmale}{%d}
\\newcommand{\\Nsports}{%d}
\\newcommand{\\NsportsCohort}{%d}
\\newcommand{\\Ncohorts}{%d}
\\newcommand{\\Njoint}{%d}
\\newcommand{\\Nvine}{%d}
\\newcommand{\\Ngauss}{%d}
\\newcommand{\\Nindep}{%d}
\\newcommand{\\MeanGain}{%.1f}
\\newcommand{\\CRPSmed}{%.3f}
\\newcommand{\\Nlin}{%d}
\\newcommand{\\NNO}{%d}
\\newcommand{\\NBCT}{%d}
\\newcommand{\\IccMed}{%.2f}
\\newcommand{\\IccQone}{%.2f}
\\newcommand{\\IccQthree}{%.2f}
\\newcommand{\\IccPass}{%d}
\\newcommand{\\NselMed}{%d}
\\newcommand{\\Nretired}{%d}
\\newcommand{\\NkeptFull}{%d}
\\newcommand{\\CmjTop}{%s}
\\newcommand{\\CmjTopFreq}{%d}
\\newcommand{\\CmjNcohorts}{%d}
\\newcommand{\\BestCohort}{%s}
\\newcommand{\\WorstCohort}{%s}
\\newcommand{\\BestCwin}{%.2f}
\\newcommand{\\WorstCwin}{%.2f}
",
n_athletes, n_tests, mean_age, sd_age,
n_female, n_male, n_sports, n_sports_cohort,
n_cohorts, n_joint,
n_vine, n_gauss, n_indep,
mean_gain, crps_med,
n_lin, n_NO, n_BCT,
icc_med, icc_q1, icc_q3, icc_pass,
nsel_med, n_retired, n_kept_full,
gsub("_", "\\\\_", cmj_top), cmj_top_freq, cmj_n_cohorts,
best_cohort_lab, worst_cohort_lab,
best_cohort_cwin, worst_cohort_cwin),
    file = "paper/numbers.tex")

cat("\nDONE — JPEG figures + tables in paper/.\n")

# =============================================================================
# 7) NEW (v3) — publication-quality coverage table, fit-quality table,
#    and calibration scatter. Idempotent; pure base R + ggplot2.
# =============================================================================
cat("v3 additions: coverage table, fit-quality table, calibration scatter...\n")

tex_safe_under <- function(x) gsub("_", "\\\\_", as.character(x))

# ---- 7a) tab_cohort_coverage.tex -------------------------------------------
# Wide table: rows = sport, columns = (Sex x Test); cell = n; blank if missing.
tests_present  <- unique(cohort_tab$test)
test_order     <- intersect(c("CMJ","DJ","SJ","SLJ","SLDJ","SLLAH","SLHAR"),
                              tests_present)
sex_order      <- c("Male","Female")
sport_order    <- sort(unique(cohort_tab$sport))

n_test <- length(test_order)
n_sex  <- length(sex_order)
align  <- paste0("l", paste(rep(strrep("c", n_test), n_sex), collapse = ""))

cov_lines <- c(
  "\\setlength{\\tabcolsep}{3pt}",
  paste0("\\begin{tabular}{", align, "}"),
  "\\toprule")

# multicolumn header row
hdr1 <- "Sport"
for (sx in sex_order)
  hdr1 <- paste0(hdr1, " & \\multicolumn{", n_test, "}{c}{", sx, "}")
hdr1 <- paste0(hdr1, " \\\\")

# cmidrule line
cmid <- ""
col0 <- 2L
for (i in seq_along(sex_order)) {
  a <- col0; b <- col0 + n_test - 1L
  cmid <- paste0(cmid, "\\cmidrule(lr){", a, "-", b, "}")
  col0 <- b + 1L
}

hdr2 <- ""
for (sx in sex_order) for (tt in test_order)
  hdr2 <- paste0(hdr2, " & ", tt)
hdr2 <- paste0(hdr2, " \\\\")

cov_lines <- c(cov_lines, hdr1, cmid, hdr2, "\\midrule")

for (sp in sport_order) {
  row <- tex_safe_under(sp)
  for (sx in sex_order) for (tt in test_order) {
    sub <- cohort_tab[cohort_tab$sex == sx & cohort_tab$sport == sp &
                       cohort_tab$test == tt, , drop = FALSE]
    cell <- if (nrow(sub) == 0) "" else as.character(sub$n[1])
    row <- paste0(row, " & ", cell)
  }
  cov_lines <- c(cov_lines, paste0(row, " \\\\"))
}
cov_lines <- c(cov_lines, "\\bottomrule", "\\end{tabular}",
                "\\setlength{\\tabcolsep}{6pt}")
writeLines(cov_lines, "paper/tables/tab_cohort_coverage.tex")

# ---- 7b) tab_fit_quality.tex -----------------------------------------------
# Cohort | n | p | C_M | C_J | mean GoF p | median 90% CI width
# C_M = mean best marginal CRPS/IQR per cohort (from margin_tab)
# C_J = winning joint Energy Score per cohort (from joint_tab)
# mean GoF p = mean across cohort metrics of KS-uniform(0,1) p-value on PIT
# median 90% CI width = median across cohort athletes of (q95 - q05) of
#   bootstrap depth percentile (B = N_BOOT_FQ); skipped where joint not fit.
N_BOOT_FQ      <- 200L
ATHLETES_PER   <- 30L  # cap per cohort for bootstrap, for runtime
RNG_FQ_SEED    <- 2026L

set.seed(RNG_FQ_SEED)

# C_M per cohort
cm_per <- margin_tab %>%
  rowwise() %>%
  mutate(best_score = min(c(Linear, NO, BCT), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(sex, sport, test) %>%
  summarise(C_M = mean(best_score, na.rm = TRUE), .groups = "drop")

# C_J per cohort (winning joint Energy Score; pulled from `Best` column)
cj_per <- joint_tab %>%
  rowwise() %>%
  mutate(C_J = c(Indep = Indep, Gaussian = Gaussian, Vine = Vine)[[Best]]) %>%
  ungroup() %>%
  select(sex, sport, test, C_J)

# GoF + bootstrap CI width — needs cohort_fits objects
gof_rows <- list()
for (key in names(cohort_fits)) {
  cf <- cohort_fits[[key]]
  ## GoF: KS uniform p across panel metrics (use the SAME predict_cdf_xs as above)
  ks_ps <- c()
  for (m in cf$panel_metrics) {
    pit_vals <- tryCatch(
      predict_cdf_xs(m$obj, m$y, m$x,
                       if (isTRUE(m$mass_dep)) m$bm else NULL),
      error = function(e) NA_real_)
    if (length(pit_vals) > 5 && all(is.finite(pit_vals))) {
      pit_vals <- pmin(pmax(pit_vals, 1e-6), 1 - 1e-6)
      ksp <- suppressWarnings(stats::ks.test(pit_vals, "punif")$p.value)
      ks_ps <- c(ks_ps, ksp)
    }
  }
  mean_gof <- if (length(ks_ps)) mean(ks_ps, na.rm = TRUE) else NA_real_

  ## Bootstrap CI width (only if joint model fit and U exists)
  med_w <- NA_real_
  ja <- cf$joint_aligned
  jm <- cf$joint_model
  if (!is.null(ja) && !is.null(ja$U) && ncol(ja$U) >= 2 && !is.null(jm)) {
    U <- as.matrix(ja$U); n_u <- nrow(U)
    if (n_u >= 10) {
      depth_pct <- function(Uref, x) {
        cm <- colMeans(Uref); cv <- stats::cov(Uref) + diag(1e-8, ncol(Uref))
        d2_a <- as.numeric(stats::mahalanobis(matrix(x, 1), cm, cv))
        d2_c <- as.numeric(stats::mahalanobis(Uref, cm, cv))
        100 * (1 - mean(d2_c <= d2_a))
      }
      take <- if (n_u <= ATHLETES_PER) seq_len(n_u)
              else sample.int(n_u, ATHLETES_PER)
      widths <- numeric(length(take))
      for (i in seq_along(take)) {
        ux <- U[take[i], ]
        bs <- replicate(N_BOOT_FQ, {
          idx <- sample.int(n_u, n_u, replace = TRUE)
          depth_pct(U[idx, , drop = FALSE], ux)
        })
        widths[i] <- as.numeric(quantile(bs, 0.95, na.rm = TRUE) -
                                  quantile(bs, 0.05, na.rm = TRUE))
      }
      med_w <- median(widths, na.rm = TRUE)
    }
  }
  gof_rows[[key]] <- data.frame(
    sex = cf$sex, sport = cf$sport, test = cf$test,
    mean_GoF_p = mean_gof, median_CI90 = med_w,
    stringsAsFactors = FALSE)
}
gof_df <- do.call(rbind, gof_rows)

# Merge into per-cohort fit-quality table; include binding cap and winning joint
# model so this single table replaces the older tab_cohorts.
fq <- cohort_tab %>%
  select(sex, sport, test, n, p_cohort, binding, joint_model) %>%
  left_join(cm_per, by = c("sex","sport","test")) %>%
  left_join(cj_per, by = c("sex","sport","test")) %>%
  left_join(gof_df, by = c("sex","sport","test")) %>%
  arrange(test, sex, sport)

fmt_num <- function(x, d = 2) {
  ifelse(is.na(x), "---", formatC(x, format = "f", digits = d))
}
fmt_pct <- function(x, d = 1) {
  ifelse(is.na(x), "---", formatC(x, format = "f", digits = d))
}

fq_lines <- c(
  "\\begin{tabular}{lllrrllrrrr}",
  "\\toprule",
  "Sex & Sport & Test & $n$ & $p$ & Binding & Joint & $C_M$ & $C_J$ & GoF $p$ & 90\\% CI \\\\",
  "\\midrule")
for (i in seq_len(nrow(fq))) {
  r <- fq[i, ]
  fq_lines <- c(fq_lines,
                 sprintf("%s & %s & %s & %d & %d & %s & %s & %s & %s & %s & %s \\\\",
                         tex_safe_under(r$sex), tex_safe_under(r$sport),
                         tex_safe_under(r$test), as.integer(r$n),
                         as.integer(r$p_cohort),
                         tex_safe_under(r$binding),
                         tex_safe_under(r$joint_model),
                         fmt_num(r$C_M, 2), fmt_num(r$C_J, 2),
                         fmt_num(r$mean_GoF_p, 2),
                         fmt_pct(r$median_CI90, 1)))
}
fq_lines <- c(fq_lines, "\\bottomrule", "\\end{tabular}")
writeLines(fq_lines, "paper/tables/tab_fit_quality.tex")

# ---- 7c) tab_calibration_scores.tex (R2 replacement of scatter figure) -----
# Compact per-cohort (C_M, C_J) numeric table + across-cohort summary, with
# the worked example flagged. Replaces the C_M/C_J scatter that the authors
# rejected as visually weak.
scores <- cm_per %>%
  inner_join(cj_per, by = c("sex", "sport", "test")) %>%
  arrange(test, sex, sport)
scores$is_example <- scores$sex == "Male" & scores$sport == "Football" &
                      scores$test == "CMJ"

cs_lines <- c(
  "\\begin{tabular}{lllrr}",
  "\\toprule",
  "Sex & Sport & Test & $C_M$ & $C_J$ \\\\",
  "\\midrule")
for (i in seq_len(nrow(scores))) {
  r <- scores[i, ]
  bold <- function(s) if (isTRUE(r$is_example)) sprintf("\\textbf{%s}", s) else s
  cs_lines <- c(cs_lines,
                 sprintf("%s & %s & %s & %s & %s \\\\",
                         bold(tex_safe_under(r$sex)),
                         bold(tex_safe_under(r$sport)),
                         bold(tex_safe_under(r$test)),
                         bold(fmt_num(r$C_M, 2)),
                         bold(fmt_num(r$C_J, 2))))
}
cs_lines <- c(cs_lines,
               "\\midrule",
               sprintf("\\multicolumn{3}{l}{Mean across cohorts} & %s & %s \\\\",
                       fmt_num(mean(scores$C_M, na.rm = TRUE), 2),
                       fmt_num(mean(scores$C_J, na.rm = TRUE), 2)),
               sprintf("\\multicolumn{3}{l}{Median across cohorts} & %s & %s \\\\",
                       fmt_num(median(scores$C_M, na.rm = TRUE), 2),
                       fmt_num(median(scores$C_J, na.rm = TRUE), 2)),
               "\\bottomrule",
               "\\end{tabular}")
writeLines(cs_lines, "paper/tables/tab_calibration_scores.tex")
wrote("paper/tables/tab_calibration_scores.tex")

# ---- 7d) tab_variable_selection.tex (R3 NEW) -------------------------------
# Per-cohort table of selected metrics: panel metrics that survived the
# reliability filter and the n >= 10 * choose(p,2) sample-size cap (i.e.
# the variables that actually entered the joint copula). joint_metrics
# is not in cohort_table.csv so we derive it from the saved cohort_fits
# objects (single source of truth).
.jm_rows <- list()
for (key in names(cohort_fits)) {
  cf <- cohort_fits[[key]]
  jm <- cf$joint_metrics
  jm_str <- if (length(jm)) paste(jm, collapse = ",") else ""
  .jm_rows[[key]] <- data.frame(
    sex = cf$sex, sport = cf$sport, test = cf$test,
    joint_metrics = jm_str,
    stringsAsFactors = FALSE)
}
joint_metrics_df <- do.call(rbind, .jm_rows)

varsel <- cohort_tab %>%
  select(sex, sport, test, n, joint_model) %>%
  left_join(joint_metrics_df, by = c("sex", "sport", "test")) %>%
  arrange(test, sex, sport)

fmt_metric_list <- function(metrics_str) {
  if (is.na(metrics_str) || nchar(metrics_str) == 0)
    return("\\emph{none retained}")
  ms <- strsplit(metrics_str, ",", fixed = TRUE)[[1]]
  ms <- gsub("^\\s+|\\s+$", "", ms)
  if (length(ms) > 6) ms <- c(ms[seq_len(6)], "\\ldots")
  paste(vapply(ms, tex_safe_under, character(1)), collapse = ", ")
}

vs_lines <- c(
  "\\begin{tabular}{lllrl}",
  "\\toprule",
  "Sex & Sport & Test & $n$ & Selected metrics (joint model input) \\\\",
  "\\midrule")
last_test <- ""
for (i in seq_len(nrow(varsel))) {
  r <- varsel[i, ]
  if (r$test != last_test && i > 1) {
    vs_lines <- c(vs_lines, "\\addlinespace[1pt]")
  }
  last_test <- r$test
  vs_lines <- c(vs_lines,
                 sprintf("%s & %s & %s & %d & %s \\\\",
                         tex_safe_under(r$sex), tex_safe_under(r$sport),
                         tex_safe_under(r$test), as.integer(r$n),
                         fmt_metric_list(r$joint_metrics)))
}
vs_lines <- c(vs_lines, "\\bottomrule", "\\end{tabular}")
writeLines(vs_lines, "paper/tables/tab_variable_selection.tex")
wrote("paper/tables/tab_variable_selection.tex")

cat("v3 additions written.\n")

# =============================================================================
# v4 (Z1+) — APPENDIX assets
# Reviewer-facing supplementary tables and figures. The STAI-X 2026 page rule
# is "8 pages for the body excluding references and appendices", so the
# appendix is unconstrained and is the right place for full per-cohort
# detail. Each artefact is wired into main.tex via \input{} or
# \includegraphics{} under the \appendix block.
# =============================================================================
cat("v4 additions: appendix tables and figures (Z1+)...\n")

suppressPackageStartupMessages({
  library(VineCopula)
})

# ---- A2.1) tab_reliability_full.tex --------------------------------------
# Per (cohort x metric) reliability row; longtable so it can run multi-page.
rt <- rel_tab[order(rel_tab$test, rel_tab$sex, rel_tab$sport, rel_tab$metric), ]
rf_lines <- c(
  "{\\small",
  "\\begin{longtable}{lllrlrrr}",
  "\\caption{\\textbf{Per-cohort reliability table (full):} every (sex, sport, test, metric) cell with its own ICC$_{(3,1)}$, SEM and MDC$_{95}$ from the test--retest pair construction in Methodology Stage~2.}\\label{tab:reliability_full}\\\\",
  "\\toprule",
  "Sex & Sport & Test & Metric & $n$ pairs & ICC$_{(3,1)}$ & SEM & MDC$_{95}$ \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{8}{l}{\\emph{(continued)}} \\\\",
  "\\toprule",
  "Sex & Sport & Test & Metric & $n$ pairs & ICC$_{(3,1)}$ & SEM & MDC$_{95}$ \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot")
last_test <- ""
for (i in seq_len(nrow(rt))) {
  r <- rt[i, ]
  if (r$test != last_test && i > 1) {
    rf_lines <- c(rf_lines, "\\addlinespace[1pt]")
  }
  last_test <- r$test
  rf_lines <- c(rf_lines,
    sprintf("%s & %s & %s & \\texttt{%s} & %d & %.2f & %.2f & %.2f \\\\",
      tex_safe_under(r$sex), tex_safe_under(r$sport), tex_safe_under(r$test),
      tex_safe_under(r$metric), as.integer(r$n_pairs),
      r$ICC, r$SEM, r$MDC95))
}
rf_lines <- c(rf_lines, "\\end{longtable}", "}")
writeLines(rf_lines, "paper/tables/tab_reliability_full.tex")
wrote("paper/tables/tab_reliability_full.tex")

# ---- A2.2) tab_marginal_competition_full.tex -----------------------------
# Per (cohort x metric) CRPS/IQR for Linear, NO, BCT, with winner bolded.
mt <- margin_tab[order(margin_tab$test, margin_tab$sex,
                        margin_tab$sport, margin_tab$metric), ]
fmt_or_dash <- function(x, d = 3) {
  ifelse(is.na(x) | !is.finite(x), "---",
         formatC(x, format = "f", digits = d))
}
bold_if <- function(s, is_winner) ifelse(is_winner, sprintf("\\textbf{%s}", s), s)
mf_lines <- c(
  "{\\small",
  "\\begin{longtable}{lllllrrrl}",
  "\\caption{\\textbf{Per-(cohort $\\times$ metric) marginal-competition (full):} 5-fold cross-validated CRPS$/$IQR (lower is sharper) for the three candidate marginals defined in Eqs.~(1)--(3); \\textbf{bold} = winner under Eq.~(4).}\\label{tab:marg_full}\\\\",
  "\\toprule",
  "Sex & Sport & Test & Metric & MassDep & Linear & NO & BCT & Winner \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{9}{l}{\\emph{(continued)}} \\\\",
  "\\toprule",
  "Sex & Sport & Test & Metric & MassDep & Linear & NO & BCT & Winner \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot")
last_test <- ""
for (i in seq_len(nrow(mt))) {
  r <- mt[i, ]
  if (r$test != last_test && i > 1) {
    mf_lines <- c(mf_lines, "\\addlinespace[1pt]")
  }
  last_test <- r$test
  is_lin <- r$Best == "Linear"; is_no <- r$Best == "NO"; is_bct <- r$Best == "BCT"
  mf_lines <- c(mf_lines,
    sprintf("%s & %s & %s & \\texttt{%s} & %s & %s & %s & %s & %s \\\\",
      tex_safe_under(r$sex), tex_safe_under(r$sport), tex_safe_under(r$test),
      tex_safe_under(r$metric),
      ifelse(isTRUE(as.logical(r$mass_dep)), "yes", "no"),
      bold_if(fmt_or_dash(r$Linear, 3), is_lin),
      bold_if(fmt_or_dash(r$NO, 3),     is_no),
      bold_if(fmt_or_dash(r$BCT, 3),    is_bct),
      r$Best))
}
mf_lines <- c(mf_lines, "\\end{longtable}", "}")
writeLines(mf_lines, "paper/tables/tab_marginal_competition_full.tex")
wrote("paper/tables/tab_marginal_competition_full.tex")

# ---- A2.3) tab_joint_competition_full.tex --------------------------------
jt <- joint_tab[order(joint_tab$test, joint_tab$sex, joint_tab$sport), ]
jf_lines <- c(
  "{\\small",
  "\\begin{longtable}{lllrrrrrlr}",
  "\\caption{\\textbf{Per-cohort joint-copula competition (full):} 5-fold cross-validated multivariate Energy Score for Independence, Gaussian copula and the regular vine of Eq.~(6); \\textbf{bold} = winner; ``Gain'' is the per-cent improvement of the winner over the Independence baseline.}\\label{tab:joint_full}\\\\",
  "\\toprule",
  "Sex & Sport & Test & $n$ & $p$ & ES Indep & ES Gaussian & ES Vine & Winner & Gain (\\%) \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{10}{l}{\\emph{(continued)}} \\\\",
  "\\toprule",
  "Sex & Sport & Test & $n$ & $p$ & ES Indep & ES Gaussian & ES Vine & Winner & Gain (\\%) \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot")
for (i in seq_len(nrow(jt))) {
  r <- jt[i, ]
  is_i <- r$Best == "Indep"; is_g <- r$Best == "Gaussian"; is_v <- r$Best == "Vine"
  jf_lines <- c(jf_lines,
    sprintf("%s & %s & %s & %d & %d & %s & %s & %s & %s & %s \\\\",
      tex_safe_under(r$sex), tex_safe_under(r$sport), tex_safe_under(r$test),
      as.integer(r$n), as.integer(r$p),
      bold_if(fmt_or_dash(r$Indep, 3),    is_i),
      bold_if(fmt_or_dash(r$Gaussian, 3), is_g),
      bold_if(fmt_or_dash(r$Vine, 3),     is_v),
      r$Best,
      ifelse(is.na(r$GainOverIndep_pct), "---",
             formatC(r$GainOverIndep_pct, format = "f", digits = 1))))
}
jf_lines <- c(jf_lines, "\\end{longtable}", "}")
writeLines(jf_lines, "paper/tables/tab_joint_competition_full.tex")
wrote("paper/tables/tab_joint_competition_full.tex")

# ---- A2.3b) tab_vine_structures.tex --------------------------------------
# One row per non-zero pair-copula edge, grouped by cohort. Extracted
# directly from the RVineMatrix slots (Matrix, family, par, tau) since
# VineCopula::summary is not exported.
extract_vine_edges <- function(rvm) {
  # VineCopula convention (RVineMatrix lower-triangular layout):
  #   pair-copula at (i, j) for i > j connects variable M[i, j] with
  #   variable M[j, j] (the column's "anchor"), conditional on
  #   M[(i+1):d, j]. Tree level is d - i + 1.
  d <- nrow(rvm$Matrix)
  if (is.null(d) || d < 2) return(NULL)
  rows <- list()
  for (i in 2:d) {
    for (j in 1:(i - 1L)) {
      fam_int <- rvm$family[i, j]
      tree_lev <- d - i + 1L
      v_a <- rvm$Matrix[i, j]
      v_b <- rvm$Matrix[j, j]
      cond_set <- if (i < d) rvm$Matrix[(i + 1L):d, j] else integer(0)
      lab <- if (length(cond_set) == 0L) {
        sprintf("%d,%d", v_a, v_b)
      } else {
        sprintf("%d,%d;%s", v_a, v_b, paste(cond_set, collapse = ","))
      }
      rows[[length(rows) + 1]] <- data.frame(
        tree = tree_lev,
        edge = lab,
        family_int = as.integer(fam_int),
        par   = as.numeric(rvm$par[i, j]),
        tau   = as.numeric(rvm$tau[i, j]),
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0L) NULL else do.call(rbind, rows)
}

vs_rows <- list()
for (key in names(cohort_fits)) {
  cf <- cohort_fits[[key]]
  jm <- cf$joint_model
  if (is.null(jm) || jm$type != "Vine" || is.null(jm$obj)) next
  rvm <- jm$obj
  edges <- extract_vine_edges(rvm)
  if (is.null(edges) || nrow(edges) == 0) next
  cohort_lab <- sprintf("%s %s %s", cf$sex, cf$sport, cf$test)
  vmap <- if (!is.null(rvm$names)) {
            paste(seq_along(rvm$names), rvm$names, sep = " = ", collapse = "; ")
          } else ""
  vine_type_lab <- if (!is.null(rvm$type)) as.character(rvm$type) else "vine"
  for (rr in seq_len(nrow(edges))) {
    er <- edges[rr, ]
    fam_int <- er$family_int
    fam_lab <- if (fam_int == 0L) "Indep"
               else tryCatch(VineCopula::BiCopName(fam_int, short = TRUE),
                              error = function(e) sprintf("fam=%d", fam_int))
    vs_rows[[length(vs_rows) + 1]] <- data.frame(
      cohort = cohort_lab,
      vine_type = vine_type_lab,
      tree = as.integer(er$tree),
      edge = er$edge,
      family = fam_lab,
      par = if (fam_int == 0L) NA_real_ else er$par,
      tau = er$tau,
      vmap = vmap,
      stringsAsFactors = FALSE)
  }
}
vs_df <- if (length(vs_rows)) do.call(rbind, vs_rows) else
         data.frame(cohort=character(0), vine_type=character(0),
                    tree=integer(0), edge=character(0), family=character(0),
                    par=numeric(0), tau=numeric(0), vmap=character(0))

vsx_lines <- c(
  "{\\small",
  "\\begin{longtable}{lcllrr}",
  "\\caption{\\textbf{Fitted R-vine structures} for cohorts whose vine wins the joint competition (one row per pair-copula edge). ``Edge'' uses VineCopula's notation: \\texttt{i,j} is an unconditional pair, \\texttt{i,j;k} is conditional on $u_k$. Family abbreviations: \\texttt{N}=Gaussian, \\texttt{t}=Student-$t$, \\texttt{C}=Clayton, \\texttt{G}=Gumbel, \\texttt{F}=Frank, \\texttt{J}=Joe, \\texttt{SC}=survival-Clayton, \\texttt{SG}=survival-Gumbel, \\texttt{Indep}=independence. $\\tau$ is the implied Kendall's $\\tau$.}\\label{tab:vine_structures}\\\\",
  "\\toprule",
  "Cohort & Tree & Edge & Family & Par & $\\tau$ \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{6}{l}{\\emph{(continued)}} \\\\",
  "\\toprule",
  "Cohort & Tree & Edge & Family & Par & $\\tau$ \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot")
last_cohort <- ""
for (i in seq_len(nrow(vs_df))) {
  r <- vs_df[i, ]
  if (r$cohort != last_cohort) {
    if (i > 1) vsx_lines <- c(vsx_lines, "\\midrule")
    vsx_lines <- c(vsx_lines,
      sprintf("\\multicolumn{6}{l}{\\textbf{%s} \\hfill {\\footnotesize \\textit{type:} %s; \\textit{var.\\ index:} %s}} \\\\",
              tex_safe_under(r$cohort), r$vine_type,
              tex_safe_under(r$vmap)))
  }
  last_cohort <- r$cohort
  vsx_lines <- c(vsx_lines,
    sprintf("        & %d & \\texttt{%s} & %s & %s & %s \\\\",
            as.integer(r$tree), r$edge, r$family,
            ifelse(is.na(r$par), "---", formatC(r$par, format = "f", digits = 3)),
            formatC(r$tau, format = "f", digits = 2)))
}
vsx_lines <- c(vsx_lines, "\\end{longtable}", "}")
writeLines(vsx_lines, "paper/tables/tab_vine_structures.tex")
wrote("paper/tables/tab_vine_structures.tex")

# ---- A2.4) fig_pit_grid.jpg ----------------------------------------------
# One panel per (cohort, metric) showing the PIT histogram of fitted
# residuals; flat under correct calibration.
pit_long <- list()
for (key in names(cohort_fits)) {
  cf <- cohort_fits[[key]]
  if (is.null(cf$joint_aligned) || is.null(cf$joint_aligned$U) ||
      ncol(cf$joint_aligned$U) < 1) next
  U <- cf$joint_aligned$U
  for (m in colnames(U)) {
    pit_long[[length(pit_long) + 1]] <- data.frame(
      cohort = sprintf("%s %s %s", cf$sex, cf$sport, cf$test),
      metric = m, pit = as.numeric(U[, m]),
      stringsAsFactors = FALSE)
  }
}
pit_df <- do.call(rbind, pit_long)
pit_df$panel <- paste(pit_df$cohort, pit_df$metric, sep = "\n")
n_panels <- length(unique(pit_df$panel))
n_col <- 5L
n_row <- ceiling(n_panels / n_col)
fig_h_pit <- max(7.0, 1.1 * n_row + 0.6)

# Per-panel reference line: expected count if PIT is uniform (n/10 bins).
ref_df <- aggregate(pit ~ panel, data = pit_df, FUN = function(x) length(x) / 10)
names(ref_df)[2] <- "ref_y"
p_pit <- ggplot(pit_df, aes(pit)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1),
                  fill = BLUE_MED, colour = "white", linewidth = 0.25) +
  geom_hline(data = ref_df, aes(yintercept = ref_y),
             linetype = "22", colour = "grey45", linewidth = 0.4) +
  facet_wrap(~ panel, ncol = n_col, scales = "free_y") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Probability Integral Transform (PIT)", y = "Count") +
  theme(strip.text = element_text(size = 6.5, lineheight = 0.85),
        axis.text  = element_text(size = 6.5),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey95"))
JPG("paper/figures/fig_pit_grid.jpg", p_pit, w = 7.0, h = fig_h_pit, dpi = 300)
wrote("paper/figures/fig_pit_grid.jpg")

# ---- A2.5) fig_centile_bands_grid.jpg ------------------------------------
# JumpHeight 5/25/50/75/95 bands per CMJ cohort, raw points overlaid.
band_panels <- list()
for (key in names(cohort_fits)) {
  cf <- cohort_fits[[key]]
  if (cf$test != "CMJ") next
  pm <- cf$panel_metrics$JumpHeight
  if (is.null(pm) || is.null(pm$obj)) next
  ages_g <- seq(min(pm$x, na.rm = TRUE), max(pm$x, na.rm = TRUE),
                length.out = 50)
  bm_ref <- if (isTRUE(pm$mass_dep))
              rep(median(pm$bm, na.rm = TRUE), length(ages_g)) else NULL
  q_mat <- tryCatch(predict_quantile_grid(pm$obj, ages_g, bm_ref),
                    error = function(e) NULL)
  if (is.null(q_mat)) next
  bands <- data.frame(age = ages_g,
                       q05 = q_mat[, "q05"], q25 = q_mat[, "q25"],
                       q50 = q_mat[, "q50"], q75 = q_mat[, "q75"],
                       q95 = q_mat[, "q95"],
                       cohort = sprintf("%s %s", cf$sex, cf$sport),
                       stringsAsFactors = FALSE)
  obs <- data.frame(age = pm$x, y = pm$y,
                     cohort = sprintf("%s %s", cf$sex, cf$sport),
                     stringsAsFactors = FALSE)
  band_panels[[length(band_panels) + 1]] <- list(bands = bands, obs = obs)
}
if (length(band_panels)) {
  bands_all <- do.call(rbind, lapply(band_panels, `[[`, "bands"))
  obs_all   <- do.call(rbind, lapply(band_panels, `[[`, "obs"))
  p_bands <- ggplot() +
    geom_ribbon(data = bands_all,
                 aes(age, ymin = q05, ymax = q95),
                 fill = BLUE_LIGHT, alpha = 0.55) +
    geom_ribbon(data = bands_all,
                 aes(age, ymin = q25, ymax = q75),
                 fill = BLUE_MED, alpha = 0.55) +
    geom_line  (data = bands_all,
                 aes(age, q50),
                 colour = BLUE_DARK, linewidth = 0.7) +
    geom_point (data = obs_all,
                 aes(age, y),
                 alpha = 0.30, size = 0.55, colour = GREY_DARK) +
    facet_wrap(~ cohort, ncol = 4, scales = "free_y") +
    labs(x = "Age (years)", y = "CMJ JumpHeight (cm)") +
    theme(strip.text = element_text(size = 9, face = "bold"),
          axis.text  = element_text(size = 7))
  JPG("paper/figures/fig_centile_bands_grid.jpg", p_bands,
      w = 7.0, h = 5.0, dpi = 300)
  wrote("paper/figures/fig_centile_bands_grid.jpg")
}

# ---- A2.6) tab_varsel_sensitivity.tex ------------------------------------
# Re-run Stage 3 selection at ICC thresholds {0.4, 0.5, 0.7}; show the
# resulting metric set per cohort. Sample-size cap is preserved (drop from
# the END of the panel-priority list until n >= 10*choose(p,2) holds).
do_selection <- function(test_type, n_size, rel_subset, icc_thr) {
  panel <- panel_slugs[[test_type]]
  if (is.null(panel) || length(panel) == 0) return(character(0))
  pass <- character(0)
  for (m in panel) {
    rmm <- rel_subset[rel_subset$metric == m, ]
    icc_val <- if (nrow(rmm) > 0) rmm$ICC[1] else NA_real_
    if (is.finite(icc_val) && icc_val >= icc_thr) pass <- c(pass, m)
  }
  while (length(pass) > 1 && n_size < 10 * choose(length(pass), 2)) {
    pass <- pass[seq_len(length(pass) - 1L)]
  }
  pass
}
sens_rows <- list()
for (key in names(cohort_fits)) {
  cf <- cohort_fits[[key]]
  rel_sub <- rel_tab[rel_tab$sex == cf$sex & rel_tab$sport == cf$sport &
                      rel_tab$test == cf$test, ]
  s04 <- do_selection(cf$test, cf$n, rel_sub, 0.4)
  s05 <- do_selection(cf$test, cf$n, rel_sub, 0.5)
  s07 <- do_selection(cf$test, cf$n, rel_sub, 0.7)
  fmt_set <- function(s) if (length(s) == 0) "\\emph{none}"
                          else paste(vapply(s, tex_safe_under,
                                              character(1)), collapse = ", ")
  sens_rows[[length(sens_rows) + 1]] <- data.frame(
    sex = cf$sex, sport = cf$sport, test = cf$test, n = cf$n,
    s04 = fmt_set(s04), s05 = fmt_set(s05), s07 = fmt_set(s07),
    stringsAsFactors = FALSE)
}
sens_df <- do.call(rbind, sens_rows)
sens_df <- sens_df[order(sens_df$test, sens_df$sex, sens_df$sport), ]
sens_lines <- c(
  "{\\footnotesize",
  "\\begin{longtable}{lllrp{3.2cm}p{3.2cm}p{3.2cm}}",
  "\\caption{\\textbf{Sensitivity of variable selection to the ICC$_{(3,1)}$ threshold.} For each cohort, the metrics that would have entered the joint copula at thresholds $0.4$, $0.5$ (default) and $0.7$, after the sample-size cap $n\\!\\geq\\!10\\binom{p}{2}$ is applied. The $0.5$ column reproduces Table~\\ref{tab:varsel}.}\\label{tab:varsel_sensitivity}\\\\",
  "\\toprule",
  "Sex & Sport & Test & $n$ & ICC $\\geq 0.4$ & ICC $\\geq 0.5$ (default) & ICC $\\geq 0.7$ \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{7}{l}{\\emph{(continued)}} \\\\",
  "\\toprule",
  "Sex & Sport & Test & $n$ & ICC $\\geq 0.4$ & ICC $\\geq 0.5$ (default) & ICC $\\geq 0.7$ \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot")
last_test <- ""
for (i in seq_len(nrow(sens_df))) {
  r <- sens_df[i, ]
  if (r$test != last_test && i > 1) {
    sens_lines <- c(sens_lines, "\\addlinespace[1pt]")
  }
  last_test <- r$test
  sens_lines <- c(sens_lines,
    sprintf("%s & %s & %s & %d & %s & %s & %s \\\\",
            tex_safe_under(r$sex), tex_safe_under(r$sport),
            tex_safe_under(r$test), as.integer(r$n),
            r$s04, r$s05, r$s07))
}
sens_lines <- c(sens_lines, "\\end{longtable}", "}")
writeLines(sens_lines, "paper/tables/tab_varsel_sensitivity.tex")
wrote("paper/tables/tab_varsel_sensitivity.tex")

# ---- A2.7) tab_score_example_extended.tex --------------------------------
# Same score_one() output for 5 randomly-selected athletes from the
# Male Football CMJ cohort (different from the worked-example athlete).
set.seed(42L)
mfcmj <- cohort_fits[["Male_Football_CMJ"]]
all_pids <- mfcmj$joint_aligned$profile_id
candidate_pids <- setdiff(unique(all_pids), chosen_pid)
ext_n <- min(5L, length(candidate_pids))
ext_pids <- sample(candidate_pids, ext_n)

ext_lines <- c(
  "{\\small",
  "\\begin{tabular}{rrrrrrrrr}",
  "\\toprule",
  "Athlete & Age & BM & JH pct & PP$_{\\text{BM}}$ pct & mRSI pct & Joint pct & 95\\% lo & 95\\% hi \\\\",
  "\\midrule")
for (i in seq_along(ext_pids)) {
  pid <- ext_pids[i]
  s <- tryCatch(score_one(pid, mfcmj), error = function(e) NULL)
  if (is.null(s)) next
  pm <- s$per_metric
  jh_p  <- pm$pct[pm$metric == "JumpHeight"][1]
  pp_p  <- pm$pct[pm$metric == "PeakPower_BM"][1]
  rsi_p <- pm$pct[pm$metric == "mRSI"][1]
  ext_lines <- c(ext_lines,
    sprintf("A%d & %.1f & %.1f & %s & %s & %s & %s & %s & %s \\\\",
            i, s$age, s$bm,
            ifelse(is.na(jh_p),  "---", formatC(jh_p,  format = "f", digits = 1)),
            ifelse(is.na(pp_p),  "---", formatC(pp_p,  format = "f", digits = 1)),
            ifelse(is.na(rsi_p), "---", formatC(rsi_p, format = "f", digits = 1)),
            ifelse(is.na(s$joint_pct), "---", formatC(s$joint_pct, format = "f", digits = 1)),
            ifelse(is.na(s$joint_lo),  "---", formatC(s$joint_lo,  format = "f", digits = 1)),
            ifelse(is.na(s$joint_hi),  "---", formatC(s$joint_hi,  format = "f", digits = 1))))
}
ext_lines <- c(ext_lines, "\\bottomrule", "\\end{tabular}", "}")
writeLines(ext_lines, "paper/tables/tab_score_example_extended.tex")
wrote("paper/tables/tab_score_example_extended.tex")

cat("v4 (appendix) additions written.\n")

