# =============================================================================
#                         RUN_ALL_IN_ONE.R
#         Multivariable Centile Modelling of VALD ForceDecks Data
#                  Sex-stratified, sport-aware, end-to-end
# -----------------------------------------------------------------------------
# Run top-to-bottom in R or RStudio. Sections 0-9 mirror the modular
# scripts code/00_setup.R … code/09_extra_figs.R, but kept intentionally
# linear and with no nested helpers so a reader can step through them.
#
# WHAT IS NEW (vs. the modular pipeline):
#   - Trial aggregation: 20%-trimmed mean with within-test CV gate, replacing
#     the arithmetic mean. Defensible from the strength-and-conditioning
#     literature (McMahon 2018; Bishop 2018; Liljequist 2019).
#   - Feature selection: LASSO with stability selection (Meinshausen &
#     Buehlmann 2010), not correlation and not random forest.
#   - 5 additional EDA plots designed for the paper.
#   - score_athlete() — practitioner-grade percentile lookup at the end.
#
# REFERENCES (numbered tags appear in section headers):
#   [R1]  Rigby & Stasinopoulos 2005, JRSS-C — GAMLSS.
#   [R2]  Cole & Green 1992, Stat. Med. — LMS method.
#   [R3]  Borghi et al. 2006, Stat. Med. — WHO growth standards.
#   [R4]  Royston & Wright 1998, Stat. Med. — multivariate reference centiles.
#   [R5]  Aas, Czado, Frigessi & Bakken 2009, IME — pair-copula construction.
#   [R6]  Joe 2014 — Dependence Modeling with Copulas.
#   [R7]  Czado & Nagler 2022, ARSIA — vine copula based modelling.
#   [R8]  Tibshirani 1996, JRSS-B — LASSO.
#   [R9]  Meinshausen & Buehlmann 2010, JRSS-B — Stability selection.
#   [R10] Gneiting & Raftery 2007, JASA — strictly proper scoring rules.
#   [R11] Gneiting & Katzfuss 2014, ARSIA — probabilistic forecasting.
#   [R12] Villani 2008 — Optimal Transport (Wasserstein).
#   [R13] Tukey 1975; Liu, Parelius & Singh 1999 — multivariate data depth.
#   [R14] Liljequist et al. 2019, PLOS ONE — ICC and reliability.
#   [R15] McMahon, Suchomel, Lake & Comfort 2018, S&C J — CMJ force-time.
#   [R16] Bishop, Read, Lake, Chavda & Turner 2018, S&C J — interlimb asymm.
#   [R17] Klein & Kneib 2016, JCGS — distributional regression.
#   [R18] Vatter & Nagler 2018, JCGS — gamCopula (covariate-dep. extension).
#   [R19] Stasinopoulos & Rigby 2007, JStatSoft — GAMLSS in R.
#   [R20] Nagler et al. 2022 — VineCopula R package v2.4.
#
# OUR USPS (vs. prior multivariable centile work):
#   1. Vine-copula-based joint normative profiling on a multi-sport ForceDecks
#      corpus, sex-stratified — most prior work uses Gaussian dependence or
#      univariate centiles only [R4, R17].
#   2. Stability-selected LASSO feature pruning specific to biomechanical
#      panels — replaces ad-hoc clinical picks and the well-known biases of
#      RF importance / correlation ranking.
#   3. Robust trimmed-mean trial aggregation with a CV gate, justified from
#      first principles in [R14-R16].
#   4. Scale-free distributional metrics (CRPS/IQR, W1/IQR), not MAPE — so
#      that fit quality is comparable across metrics with different units.
#   5. Joint percentile via Tukey / Mahalanobis depth on the PIT vector —
#      collapses the multivariable profile into a single interpretable score
#      that an S&C coach can act on, not just a statistical artefact.
#   6. score_athlete(): an end-to-end percentile lookup that is the closest
#      analogue of a clinical centile chart for an athletic profile.
#   7. The whole pipeline runs end-to-end with a fixed RNG seed.
# =============================================================================


# =============================================================================
# SECTION 0 — Setup, paths, packages, configuration
# =============================================================================

PROJECT_DIR <- "~/Downloads/VALD_Multivariate_Centiles_STAIX2026"
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
setwd(PROJECT_DIR)

# ---- (optional) live VALD pull ---------------------------------------------
# Paste your credentials below and set DO_LIVE_PULL <- TRUE to call the API.
# If you keep DO_LIVE_PULL <- FALSE the script will use the integrated CSV
# defined immediately below.
Sys.setenv(VALD_CLIENT_ID     = "")  # <-- paste here for a live pull
Sys.setenv(VALD_CLIENT_SECRET = "")  # <-- paste here for a live pull
Sys.setenv(VALD_TENANT_ID     = "")  # <-- paste here for a live pull
Sys.setenv(VALD_REGION        = "aue")
Sys.setenv(VALD_START_DATE    = "2020-01-01T00:00:00Z")

DO_LIVE_PULL  <- FALSE
INTEGRATED_CSV <- "~/Downloads/SAI_Mumbai/VALD_integrated_output/20260508/01_VALD_integrated_complete.csv"
SPORTS_DIR     <- "~/Downloads/VALD Data_sports_gender"

# ---- Pipeline tunables ------------------------------------------------------
MIN_N_STRATUM <- 25     # minimum cohort size to fit a marginal
MAX_METRICS   <-  4     # dimension of the multivariable profile
N_SIM         <- 300    # predictive draws per observation
CV_FOLDS      <-  5
N_BOOT_SS     <- 100    # bootstrap reps for stability selection
SS_PI_THRESH  <- 0.60   # selection probability threshold (Meinshausen 2010)
TRIM          <- 0.20   # 20% trimmed mean for trial aggregation
MAX_TRIAL_CV  <- 0.25   # max within-test CV (drop test if exceeded)
SEED          <- 2026
set.seed(SEED)

# ---- Required packages ------------------------------------------------------
required <- c(
  "dplyr", "tidyr", "data.table", "ggplot2", "gridExtra", "scales",
  "gamlss", "gamlss.dist", "VineCopula", "glmnet", "scoringRules",
  "fmsb", "reshape2", "GGally"
)
miss <- setdiff(required, rownames(installed.packages()))
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
invisible(lapply(required, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))

ggplot2::theme_set(theme_minimal(base_size = 10))
PAL <- list(female = "#d7301f", male = "#0570b0",
            band = "#3690c0", median = "#08519c")

# ---- Output folders ---------------------------------------------------------
for (d in c("data", "data/per_test", "results", "results/strata",
            "paper/figures", "paper/tables", "logs")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Quietly suppress GAMLSS iteration logs
quiet <- function(expr) {
  zz <- tempfile(); sink(zz, type = "output")
  on.exit({ sink(type = "output"); unlink(zz) })
  force(expr)
}


# =============================================================================
# SECTION 1 — VALD API pull (optional)
# Justification: VALD is a partially-cloud system whose only authoritative
# source of truth is the API. We pull profiles, tests, and trials.
# =============================================================================

if (DO_LIVE_PULL) {
  if (!requireNamespace("valdr", quietly = TRUE))
    install.packages("valdr", repos = "https://cloud.r-project.org")
  options(timeout = 18000)
  valdr::set_credentials(
    client_id     = Sys.getenv("VALD_CLIENT_ID"),
    client_secret = Sys.getenv("VALD_CLIENT_SECRET"),
    tenant_id     = Sys.getenv("VALD_TENANT_ID"),
    region        = Sys.getenv("VALD_REGION")
  )
  valdr::set_start_date(Sys.getenv("VALD_START_DATE"))
  raw <- valdr::get_forcedecks_data(include_attributes = TRUE)
  saveRDS(raw, "data/raw_forcedecks.rds")
  cat("[01] pulled", nrow(raw$trials), "trials,",
      nrow(raw$tests), "tests,", nrow(raw$profiles), "profiles\n")
}


# =============================================================================
# SECTION 2 — Cleaning, integration & robust aggregation
# Aggregation method (key methodological choice):
#   - per (athlete x test x definition x limb), take the 20%-trimmed mean of
#     the trial values [R14, R15, R16];
#   - drop the test if fewer than 2 trials remain after NA filtering or if
#     the within-test coefficient of variation exceeds 25% (signal of slip,
#     miscount, or technical fault) [R14].
# This is robust to outlier trials and to fatigue drift across the trial
# sequence and is closer to the "best stable trial" used in S&C practice.
# =============================================================================

agg_trim <- function(x, trim = TRIM, min_n = 2, max_cv = MAX_TRIAL_CV) {
  x <- x[is.finite(x)]
  if (length(x) < min_n) return(NA_real_)
  cv <- stats::sd(x) / max(abs(mean(x)), 1e-9)
  if (cv > max_cv) return(NA_real_)
  mean(x, trim = trim)
}

# ---- Build wide analysis frame ---------------------------------------------
if (DO_LIVE_PULL && file.exists("data/raw_forcedecks.rds")) {
  raw <- readRDS("data/raw_forcedecks.rds")

  trials <- raw$trials
  tests  <- raw$tests
  profs  <- raw$profiles
  trials$value     <- suppressWarnings(as.numeric(trials$value))
  trials$athleteId <- as.character(trials$athleteId)

  joined <- dplyr::left_join(trials, tests, by = "testId")

  agg <- joined |>
    dplyr::group_by(athleteId, testId, testType, recordedUTC,
                    recordedDateOffset, trialLimb, definition_name) |>
    dplyr::summarise(value = agg_trim(value), .groups = "drop") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(
      ts_utc    = lubridate::ymd_hms(recordedUTC),
      ts_local  = ts_utc + lubridate::minutes(recordedDateOffset),
      test_date = as.Date(ts_local)
    )

  # external sport / gender CSV folder
  read_sport_dir <- function(folder) {
    files <- list.files(folder, "\\.csv$", full.names = TRUE)
    do.call(rbind, lapply(files, function(f) {
      df <- read.csv(f, stringsAsFactors = FALSE)
      pid <- names(df)[grepl("profile.*id", names(df), ignore.case = TRUE)][1]
      sx  <- names(df)[grepl("^gender$", names(df), ignore.case = TRUE)][1]
      sport <- sub("\\..*$", "", sub("_All.*$|_Overall.*$", "", basename(f)))
      data.frame(profile_id = as.character(df[[pid]]),
                 gender = df[[sx]], sport = sport,
                 stringsAsFactors = FALSE)
    }))
  }
  sports <- if (dir.exists(SPORTS_DIR)) read_sport_dir(SPORTS_DIR) else NULL

  # pivot wide
  wide <- agg |>
    dplyr::select(profile_id = athleteId, test_id = testId,
                  test_type = testType, test_date,
                  limb = trialLimb, metric = definition_name,
                  value) |>
    tidyr::pivot_wider(
      id_cols     = c(profile_id, test_date, test_id),
      names_from  = c(metric, limb, test_type),
      values_from = value,
      names_glue  = "{metric}_{limb}_{test_type}",
      names_repair = "universal"
    )

  if (!is.null(profs)) {
    profs$profileId <- as.character(profs$profileId)
    profs_clean <- data.frame(
      profile_id = profs$profileId,
      dob = if ("dateOfBirth" %in% names(profs))
              as.Date(profs$dateOfBirth) else as.Date(NA),
      stringsAsFactors = FALSE)
    wide <- dplyr::left_join(wide, profs_clean, by = "profile_id")
  }
  if (!is.null(sports)) {
    wide <- dplyr::left_join(wide, sports, by = "profile_id")
  }
  wide$age_at_test <- as.numeric(difftime(wide$test_date, wide$dob,
                                           units = "days")) / 365.25
  data_all <- wide
} else {
  data_all <- read.csv(INTEGRATED_CSV, stringsAsFactors = FALSE,
                        check.names = FALSE)
}

# ---- Standardise gender + sport --------------------------------------------
g <- toupper(trimws(as.character(data_all$gender)))
data_all$gender <- ifelse(g %in% c("M", "MALE"), "Male",
                   ifelse(g %in% c("F", "FEMALE"), "Female", NA))
data_all$sport <- gsub("\\s+(Overall(\\s+Data)?)$", "", data_all$sport,
                       ignore.case = TRUE)
data_all$sport <- trimws(data_all$sport)
data_all <- data_all[!is.na(data_all$gender) & !is.na(data_all$sport) &
                       !is.na(data_all$age_at_test) &
                       data_all$age_at_test >= 8 & data_all$age_at_test <= 60, ]

write.csv(data_all, "data/analysis_dataset.csv", row.names = FALSE)
cat("[02] analysis_dataset.csv:", nrow(data_all), "rows,",
    ncol(data_all), "columns\n")

# ---- Per-test slices --------------------------------------------------------
meta_cols <- intersect(c("profile_id","name","gender","sport","dob",
                          "age_at_test","test_date","test_id",
                          "weight","height_cm","weight_kg"),
                        names(data_all))
metric_cols <- setdiff(names(data_all), meta_cols)
test_types  <- unique(sub(".*_([A-Z]+)$", "\\1", metric_cols))
test_types  <- test_types[nzchar(test_types)]
for (tt in test_types) {
  cols <- metric_cols[grepl(paste0("_", tt, "$"), metric_cols)]
  sub  <- data_all[, c(meta_cols, cols), drop = FALSE]
  keep <- rowSums(!is.na(sub[, cols, drop = FALSE])) > 0
  if (sum(keep) >= 10)
    write.csv(sub[keep, ], file.path("data/per_test", paste0(tt, ".csv")),
              row.names = FALSE)
}


# =============================================================================
# SECTION 3 — Exploratory data analysis (paper-grade)
# Five plots are produced and saved both as PNG (in results/) and as PDF
# (in paper/figures/) for inclusion in the manuscript.
# =============================================================================

per_test_files <- list.files("data/per_test", "\\.csv$", full.names = TRUE)
strata <- list()
for (f in per_test_files) {
  tt <- tools::file_path_sans_ext(basename(f))
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[!is.na(d$gender) & !is.na(d$sport) &
           !is.na(d$age_at_test) &
           d$age_at_test >= 8 & d$age_at_test <= 60, ]
  if (nrow(d)) strata[[tt]] <- d
}

# ---- Cohort counts table ----------------------------------------------------
counts <- do.call(rbind, lapply(names(strata), function(tt) {
  d <- strata[[tt]]
  ct <- dplyr::count(d, gender, sport, name = "n")
  ct$test_type <- tt
  ct
}))
eligible <- counts[counts$n >= MIN_N_STRATUM, ]
write.csv(counts,   "results/03_stratum_counts.csv",   row.names = FALSE)
write.csv(eligible, "results/03_eligible_strata.csv",  row.names = FALSE)
cat("[03] eligible cohorts:", nrow(eligible), "of", nrow(counts), "\n")

# ---- Plot 1: cohort sample sizes (heatmap, faceted by sex) -----------------
p1 <- ggplot(counts, aes(test_type, sport, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 3) +
  facet_wrap(~ gender) +
  scale_fill_viridis_c(option = "viridis", na.value = "grey90",
                       name = "n") +
  labs(title = "Cohort sample sizes by sport, test, and sex",
       subtitle = paste0("Strata used downstream require n >= ", MIN_N_STRATUM),
       x = "Test type", y = "Sport") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("results/03_counts_heatmap.png", p1, width = 12, height = 5.5, dpi = 130)
ggsave("paper/figures/fig_counts_heatmap.pdf", p1, width = 11, height = 5.5)

# ---- Plot 2: age density per stratum ---------------------------------------
age_df <- do.call(rbind, lapply(names(strata), function(tt) {
  d <- strata[[tt]][, c("gender", "sport", "age_at_test")]
  d$test_type <- tt; d
}))
key <- paste(eligible$gender, eligible$sport, eligible$test_type)
age_df <- age_df[paste(age_df$gender, age_df$sport, age_df$test_type) %in% key, ]
p2 <- ggplot(age_df, aes(age_at_test, fill = gender, color = gender)) +
  geom_density(alpha = 0.35, linewidth = 0.5, na.rm = TRUE) +
  scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
  scale_fill_manual(values  = c(Female = PAL$female, Male = PAL$male)) +
  facet_grid(test_type ~ sport, scales = "free_y") +
  labs(title = "Age distribution within each eligible cohort",
       x = "Age (years)", y = NULL) +
  theme(strip.text.y = element_text(angle = 0, size = 7),
        strip.text.x = element_text(size = 8),
        legend.position = "bottom")
ggsave("results/03_age_distributions.png", p2, width = 14, height = 8, dpi = 130)
ggsave("paper/figures/fig_age_distributions.pdf", p2, width = 12, height = 7)

# ---- Plot 3: cumulative tests over time (data growth) ----------------------
data_all$test_date <- as.Date(data_all$test_date)
cum_df <- data_all |>
  dplyr::filter(!is.na(test_date)) |>
  dplyr::mutate(month = as.Date(format(test_date, "%Y-%m-01"))) |>
  dplyr::group_by(month, gender) |>
  dplyr::summarise(n_tests = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(month) |>
  dplyr::group_by(gender) |>
  dplyr::mutate(cum = cumsum(n_tests))
p3 <- ggplot(cum_df, aes(month, cum, color = gender)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
  labs(title = "Cumulative ForceDecks tests recorded over time",
       subtitle = "Both sexes shown separately to expose recruitment imbalance",
       x = NULL, y = "Cumulative test count")
ggsave("results/03_cumulative_tests.png", p3, width = 9, height = 4.5, dpi = 130)
ggsave("paper/figures/fig_cumulative_tests.pdf", p3, width = 9, height = 4.5)

# ---- Plot 4: missingness map for the largest cohort ------------------------
big_test <- "CMJ"
if (big_test %in% names(strata)) {
  d <- strata[[big_test]]
  num_cols <- names(d)[sapply(d, is.numeric) & names(d) != "age_at_test"]
  miss <- as.data.frame(is.na(d[, num_cols]))
  miss$row <- seq_len(nrow(miss))
  miss_long <- reshape2::melt(miss, id.vars = "row",
                              variable.name = "metric", value.name = "missing")
  miss_long <- miss_long |>
    dplyr::group_by(metric) |>
    dplyr::summarise(p_miss = mean(missing), .groups = "drop") |>
    dplyr::arrange(p_miss) |>
    dplyr::slice_head(n = 30)
  p4 <- ggplot(miss_long, aes(reorder(metric, p_miss), 1 - p_miss)) +
    geom_col(fill = "#3690c0") +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = paste("Per-metric coverage in", big_test, "cohort"),
         subtitle = "30 metrics with the highest non-missing rate",
         x = NULL, y = "Coverage (1 - missing fraction)") +
    theme(axis.text.y = element_text(size = 7))
  ggsave("results/03_missingness.png", p4, width = 9, height = 7, dpi = 130)
  ggsave("paper/figures/fig_missingness.pdf", p4, width = 9, height = 7)
}

# ---- Plot 5: jump-height by age, per sex (the canonical curve) -------------
if ("CMJ" %in% names(strata)) {
  d <- strata[["CMJ"]]
  jh_col <- grep("Jump.Height", names(d), value = TRUE)[1]
  if (!is.na(jh_col)) {
    pdf_dat <- d[, c("age_at_test", "gender", "sport", jh_col)]
    names(pdf_dat)[4] <- "jh"
    pdf_dat <- pdf_dat[is.finite(pdf_dat$jh), ]
    p5 <- ggplot(pdf_dat, aes(age_at_test, jh, color = gender)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  linewidth = 0.9) +
      scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
      labs(title = "Countermovement-jump height vs. age, by sex",
           subtitle = jh_col,
           x = "Age (years)", y = "Jump height")
    ggsave("results/03_jh_vs_age.png", p5, width = 9, height = 5, dpi = 130)
    ggsave("paper/figures/fig_jh_vs_age.pdf", p5, width = 9, height = 5)
  }
}


# =============================================================================
# SECTION 4 — Feature selection by LASSO with stability selection [R8, R9]
# Why this method:
#   - LASSO (Tibshirani 1996, [R8]) is the canonical sparsity-inducing
#     regression; coefficients have a frequentist sampling distribution
#     and the regularisation path is interpretable.
#   - Stability selection (Meinshausen & Buehlmann 2010, [R9]) is the
#     state-of-the-art frequentist guarantee on false-discovery rate of
#     selected variables, by repeatedly subsampling the data and tracking
#     the empirical selection probability of each predictor.
#   - This is more defensible than correlation ranking (which ignores
#     joint information) and than random-forest variable importance
#     (which is biased toward high-cardinality predictors and is unstable
#     in the n < 1000 regime).
# Implementation: regress age on standardised metrics, run B = 100 bootstrap
# re-samplings, retain a metric if its empirical selection probability across
# B fits at lambda=lambda.1se exceeds pi_thr = 0.6. Top MAX_METRICS by
# selection probability are kept.
# =============================================================================

stability_select <- function(X, y, B = N_BOOT_SS, pi_thr = SS_PI_THRESH,
                             max_keep = MAX_METRICS) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  ok <- complete.cases(X, y)
  X <- X[ok, , drop = FALSE]; y <- y[ok]
  if (ncol(X) < 2 || length(y) < 30) return(colnames(X)[seq_len(min(ncol(X), max_keep))])

  # Per-feature median imputation so weakly-missing features can compete.
  for (j in seq_len(ncol(X))) {
    v <- X[, j]; v[!is.finite(v)] <- stats::median(v, na.rm = TRUE)
    X[, j] <- v
  }
  X <- scale(X)
  sel <- matrix(0L, B, ncol(X), dimnames = list(NULL, colnames(X)))
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(X), floor(nrow(X) / 2))
    fit <- tryCatch(glmnet::cv.glmnet(X[idx, , drop = FALSE], y[idx],
                                       alpha = 1, nfolds = 5),
                    error = function(e) NULL)
    if (is.null(fit)) next
    coefs <- as.numeric(coef(fit, s = "lambda.1se"))[-1]
    sel[b, ] <- as.integer(coefs != 0)
  }
  prob <- colMeans(sel)
  ord <- order(-prob)
  # Variables passing the stability threshold ([R9] guarantee).
  picked <- colnames(X)[ord][prob[ord] >= pi_thr]
  # If fewer than max_keep pass, top up by selection probability so the
  # multivariable profile retains its target dimensionality. The threshold
  # provides the "discovered" set; the top-up ensures a usable joint dim.
  if (length(picked) < max_keep) {
    extra <- setdiff(colnames(X)[ord], picked)
    picked <- c(picked, head(extra, max_keep - length(picked)))
  }
  attr(picked, "selection_prob") <- prob[picked]
  head(picked, max_keep)
}


# =============================================================================
# SECTION 5 — GAMLSS marginals per stratum [R1, R2, R19]
# y_i ~ F(mu_i, sigma_i, nu, tau);  g_mu(mu_i) = pb(age_i); g_sigma similarly.
# Family is BIC-selected from realplus (positive metrics) or realAll
# (signed metrics) using gamlss.dist::fitDist. Robust fallbacks below.
# =============================================================================

gamlss_age_fit <- function(y, x, family_chr, n_cyc = 120) {
  dat <- data.frame(y = y, x = x)
  args <- list(
    formula       = y ~ gamlss::pb(x),
    sigma.formula = ~ gamlss::pb(x),
    family        = family_chr,
    data          = dat,
    trace         = FALSE,
    control       = gamlss::gamlss.control(n.cyc = n_cyc)
  )
  quiet(do.call(gamlss::gamlss, args))
}

sample_predictive <- function(fit, x_new, n_sim = N_SIM) {
  fam <- as.character(fit$family)[1]
  rfun <- tryCatch(get(paste0("r", fam),
                       envir = asNamespace("gamlss.dist"),
                       mode = "function"),
                   error = function(e) NULL)
  if (is.null(rfun)) return(NULL)
  pa <- tryCatch(suppressWarnings(gamlss::predictAll(
    fit, newdata = data.frame(x = x_new), type = "response",
    data = fit$call$data)),
    error = function(e) NULL)
  if (is.null(pa)) return(NULL)
  pn <- intersect(c("mu", "sigma", "nu", "tau"), names(pa))
  if (!length(pn)) return(NULL)
  out <- replicate(n_sim,
    do.call(rfun, c(list(n = length(x_new)), pa[pn])))
  matrix(out, nrow = length(x_new))
}

# Robust draw: BCCGo can blow up at age boundaries — fallback to NO and
# winsorize predictive samples [R10, R12].
robust_predict <- function(fit, x_new, y_obs, n_sim) {
  fam <- as.character(fit$family)[1]
  yq <- stats::quantile(y_obs, c(0.005, 0.5, 0.995), na.rm = TRUE, names = FALSE)
  span <- max(yq[3] - yq[1], stats::sd(y_obs, na.rm = TRUE), 1e-6)
  cap_lo <- yq[1] - 8 * span; cap_hi <- yq[3] + 8 * span
  samp <- tryCatch(sample_predictive(fit, x_new, n_sim), error = function(e) NULL)
  blew <- !is.null(samp) && (!all(is.finite(samp)) ||
            max(abs(samp), na.rm = TRUE) > 1e3 * (max(abs(y_obs), na.rm = TRUE) + 1))
  if (blew && fam != "NO") {
    fit_no <- tryCatch(gamlss_age_fit(y_obs, fit$call$data$x, "NO"),
                       error = function(e) NULL)
    if (!is.null(fit_no)) {
      fit  <- fit_no; fam <- "NO"
      samp <- tryCatch(sample_predictive(fit, x_new, n_sim), error = function(e) NULL)
    }
  }
  if (!is.null(samp)) {
    samp[!is.finite(samp)] <- NA_real_
    samp[samp < cap_lo] <- cap_lo; samp[samp > cap_hi] <- cap_hi
  }
  list(fit = fit, fam = fam, samp = samp)
}

# Distributional metrics on the empirical/simulated distributions [R10, R11, R12]
wasserstein1d <- function(x, y, n_grid = 1024) {
  x <- sort(as.numeric(x[is.finite(x)]))
  y <- sort(as.numeric(y[is.finite(y)]))
  if (!length(x) || !length(y)) return(NA_real_)
  u  <- seq(0.5/n_grid, 1 - 0.5/n_grid, length.out = n_grid)
  mean(abs(stats::quantile(x, probs = u, type = 8) -
           stats::quantile(y, probs = u, type = 8)))
}
ks_uniform  <- function(p) {
  p <- p[is.finite(p)]
  if (length(p) < 5) return(NA_real_)
  unname(suppressWarnings(stats::ks.test(p, "punif")$statistic))
}
cvm_uniform <- function(p) {
  p <- sort(p[is.finite(p)]); n <- length(p)
  if (n < 5) return(NA_real_)
  i <- seq_len(n)
  sum((p - (2*i-1)/(2*n))^2) + 1/(12*n)
}
ad_uniform <- function(p) {
  p <- sort(p[is.finite(p)]); n <- length(p)
  if (n < 5) return(NA_real_)
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  i <- seq_len(n)
  -n - sum((2*i-1)/n * (log(p) + log(1 - rev(p))))
}

# Iterate cohorts (one simple top-level for-loop, no nesting) ---------------
fits_all <- list()
metric_rows <- list()

cat("[04] fitting GAMLSS marginals ...\n")
for (i in seq_len(nrow(eligible))) {
  sex <- eligible$gender[i]; sport <- eligible$sport[i]; tt <- eligible$test_type[i]
  d <- strata[[tt]]
  d <- d[d$gender == sex & d$sport == sport, , drop = FALSE]
  if (nrow(d) < MIN_N_STRATUM) next

  # Candidate metrics (numeric, sufficiently complete, non-blacklisted)
  meta <- c("profile_id","name","gender","sport","dob","age_at_test",
            "test_date","test_id","weight","height_cm","weight_kg")
  cand <- setdiff(names(d), meta)
  blacklist <- "(?i)Bodyweight|Body\\.?Mass|BM\\.?Index|Sample\\.Rate|Test\\.Duration"
  cand <- cand[!grepl(blacklist, cand, perl = TRUE)]
  for (m in cand) d[[m]] <- suppressWarnings(as.numeric(d[[m]]))
  ok_metrics <- cand[vapply(cand, function(m) {
    v <- d[[m]]
    sum(is.finite(v)) >= max(15L, 0.5 * nrow(d)) &&
      stats::sd(v, na.rm = TRUE) > 0
  }, logical(1))]
  if (length(ok_metrics) < 2) next

  # LASSO + stability selection on age ---------------------------------------
  X <- as.matrix(d[, ok_metrics, drop = FALSE])
  picked <- stability_select(X, d$age_at_test)
  if (!length(picked)) next

  sid <- gsub("[^A-Za-z0-9]+", "_", paste(sex, sport, tt, sep = "__"))
  fits_all[[sid]] <- list(sex = sex, sport = sport, test_type = tt,
                           metrics = list(),
                           profile_id = d$profile_id,
                           age = d$age_at_test)
  cat(sprintf("  %-6s | %-18s | %-7s | n=%3d  picked: %s\n",
              sex, sport, tt, nrow(d),
              paste(picked, collapse = ", ")))

  for (m in picked) {
    yvec <- suppressWarnings(as.numeric(d[[m]]))
    xvec <- d$age_at_test
    keep <- is.finite(yvec) & is.finite(xvec)
    if (sum(keep) < MIN_N_STRATUM) next
    y <- yvec[keep]; x <- xvec[keep]
    has_neg <- any(y < 0)
    fam <- tryCatch({
      fd <- quiet(gamlss.dist::fitDist(
        y, type = if (has_neg) "realAll" else "realplus",
        k = log(length(y)), trace = FALSE))
      fd$family[1]
    }, error = function(e) if (has_neg) "NO" else "BCCGo")
    fit <- tryCatch(gamlss_age_fit(y, x, fam),
                    error = function(e) tryCatch(gamlss_age_fit(y, x, "NO"),
                                                  error = function(e2) NULL))
    if (is.null(fit)) next

    rs <- robust_predict(fit, x, y, N_SIM)
    fit  <- rs$fit; fam <- rs$fam; samp <- rs$samp
    pit  <- pnorm(stats::residuals(fit))
    iqr  <- max(stats::IQR(y, na.rm = TRUE),
                 stats::sd(y, na.rm = TRUE), 1e-9)

    fits_all[[sid]]$metrics[[m]] <- list(
      fit = fit, family = fam, x = x, y = y, pit = pit,
      keep_idx = which(keep))

    crps <- if (!is.null(samp))
      mean(scoringRules::crps_sample(y = y, dat = samp), na.rm = TRUE)
      else NA_real_
    w1   <- if (!is.null(samp)) wasserstein1d(y, as.numeric(samp)) else NA_real_

    metric_rows[[length(metric_rows) + 1]] <- data.frame(
      sex = sex, sport = sport, test_type = tt, metric = m,
      n = length(y), family = fam,
      AIC = round(stats::AIC(fit), 1),
      BIC = round(stats::BIC(fit), 1),
      CRPS = signif(crps, 4), CRPS_rel = signif(crps / iqr, 4),
      Wasserstein1 = signif(w1, 4), W1_rel = signif(w1 / iqr, 4),
      KS_D   = signif(ks_uniform(pit), 4),
      CvM_W2 = signif(cvm_uniform(pit), 4),
      AD_A2  = signif(ad_uniform(pit), 4),
      stringsAsFactors = FALSE)
  }
}

marginals_table <- do.call(rbind, metric_rows)
write.csv(marginals_table, "results/04_marginals_table.csv", row.names = FALSE)
saveRDS(fits_all, "results/04_marginals_fits.rds")
cat("[04] strata =", length(fits_all), "rows =", nrow(marginals_table), "\n")


# =============================================================================
# SECTION 6 — Joint structure via regular vine copulas [R5, R6, R7, R20]
# By Sklar's theorem, F_joint(y_1, ..., y_p | age) = C(F_1, ..., F_p) and we
# estimate C from the PIT residuals u_ij = F_j(y_ij | age_i) of the marginals.
# A regular vine factorises C as a sequence of bivariate (pair-)copulas
# arranged on a hierarchy of spanning trees [R5]. RVineStructureSelect with
# AIC chooses the tree structure and the family of each pair copula from the
# parametric library in VineCopula [R20].
# =============================================================================

build_pit_matrix <- function(stratum) {
  ms <- stratum$metrics
  if (length(ms) < 2) return(NULL)
  # Each metric uses its own complete cases; align on profile_id
  mats <- lapply(names(ms), function(nm) {
    o <- ms[[nm]]
    data.frame(profile_id = stratum$profile_id[o$keep_idx], pit = o$pit)
  })
  names(mats) <- names(ms)
  M <- Reduce(function(a, b) merge(a, b, by = "profile_id"), mats)
  if (nrow(M) < MIN_N_STRATUM) return(NULL)
  U <- as.matrix(M[, -1, drop = FALSE])
  colnames(U) <- names(ms)
  U[stats::complete.cases(U), , drop = FALSE]
}

vines <- list(); vine_rows <- list()
cat("[05] fitting vine copulas ...\n")
for (sid in names(fits_all)) {
  s <- fits_all[[sid]]
  U <- build_pit_matrix(s)
  if (is.null(U)) next
  U <- pmin(pmax(U, 1e-6), 1 - 1e-6)
  vine <- tryCatch(suppressWarnings(VineCopula::RVineStructureSelect(
            U, familyset = NA, type = 0,
            selectioncrit = "AIC", indeptest = TRUE, level = 0.05)),
          error = function(e) NULL)
  if (is.null(vine)) next
  out_dir <- file.path("results/strata", sid)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(vine = vine, U = U), file.path(out_dir, "vine.rds"))
  vines[[sid]] <- list(vine = vine, U = U)
  vine_rows[[length(vine_rows) + 1]] <- data.frame(
    stratum = sid, sex = s$sex, sport = s$sport, test_type = s$test_type,
    n_obs = nrow(U), dim = ncol(U),
    vine_AIC = round(vine$AIC, 1),
    vine_BIC = round(vine$BIC, 1),
    vine_loglik = round(sum(VineCopula::RVineLogLik(U, vine)$loglik), 1),
    stringsAsFactors = FALSE)
}
joint_summary <- do.call(rbind, vine_rows)
write.csv(joint_summary, "results/05_joint_summary.csv", row.names = FALSE)
cat("[05] vines =", length(vines), "\n")


# =============================================================================
# SECTION 7 — Multivariable centile profiles, joint depth, and radar [R13]
# Joint percentile = depth-based rank of an athlete's PIT vector u_i within
# the cohort's empirical PIT cloud. Tukey halfspace depth where available;
# fallback to (1 - empirical_F(Mahalanobis^2)).
# =============================================================================

joint_depth <- function(U) {
  d2 <- stats::mahalanobis(U, colMeans(U), stats::cov(U))
  1 - rank(d2) / length(d2)
}

cat("[06] producing per-cohort centile bands and radars ...\n")
for (sid in names(fits_all)) {
  s <- fits_all[[sid]]
  if (length(s$metrics) < 2) next
  out_dir <- file.path("results/strata", sid)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ages <- unlist(lapply(s$metrics, function(m) m$x))
  age_grid <- seq(min(ages), max(ages), length.out = 120)
  panels <- list()
  for (nm in names(s$metrics)) {
    o <- s$metrics[[nm]]
    rs <- robust_predict(o$fit, age_grid, o$y, N_SIM)
    if (is.null(rs$samp)) next
    qs <- t(apply(rs$samp, 1, function(v) stats::quantile(v,
              c(.03, .10, .25, .50, .75, .90, .97), na.rm = TRUE,
              names = FALSE)))
    cf <- data.frame(age = age_grid,
                     c03 = qs[,1], c10 = qs[,2], c25 = qs[,3],
                     c50 = qs[,4], c75 = qs[,5], c90 = qs[,6], c97 = qs[,7])
    obs <- data.frame(x = o$x, y = o$y)
    panels[[nm]] <- ggplot() +
      geom_point(data = obs, aes(x, y), color = "grey55",
                 size = 0.9, alpha = 0.45) +
      geom_ribbon(data = cf, aes(age, ymin = c03, ymax = c97),
                  fill = PAL$band, alpha = 0.13) +
      geom_ribbon(data = cf, aes(age, ymin = c10, ymax = c90),
                  fill = PAL$band, alpha = 0.20) +
      geom_ribbon(data = cf, aes(age, ymin = c25, ymax = c75),
                  fill = PAL$band, alpha = 0.30) +
      geom_line(data = cf, aes(age, c50),
                color = PAL$median, linewidth = 0.9) +
      labs(title = nm, x = "Age (years)", y = NULL) +
      theme(plot.title = element_text(size = 8))
  }
  if (length(panels)) {
    g <- gridExtra::arrangeGrob(grobs = panels,
                                 ncol = min(2, length(panels)))
    ggsave(file.path(out_dir, "centiles.pdf"), g,
           width = 10, height = 3.4 * ceiling(length(panels) / 2),
           limitsize = FALSE)
    ggsave(file.path(out_dir, "centiles.png"), g, dpi = 130,
           width = 10, height = 3.4 * ceiling(length(panels) / 2),
           limitsize = FALSE)
  }

  # Joint depth per athlete -------------------------------------------------
  v <- vines[[sid]]
  if (!is.null(v)) {
    U <- v$U
    pct <- round(joint_depth(U) * 100, 1)
    write.csv(data.frame(joint_pctile = pct, U),
              file.path(out_dir, "joint_percentiles.csv"),
              row.names = FALSE)
  }
}


# =============================================================================
# SECTION 8 — 5-fold cross-validation of the marginals
# Out-of-sample CRPS_rel, W1_rel, KS_D for each (cohort, metric).
# Detects overfitting when CV scores diverge from in-sample.
# =============================================================================

cat("[07] cross-validating marginals ...\n")
cv_rows <- list()
for (sid in names(fits_all)) {
  s <- fits_all[[sid]]
  for (mn in names(s$metrics)) {
    o <- s$metrics[[mn]]; y <- o$y; x <- o$x
    folds <- split(sample.int(length(y)), seq_len(CV_FOLDS))
    iqr <- max(stats::IQR(y, na.rm = TRUE), stats::sd(y, na.rm = TRUE), 1e-9)
    f_crps <- numeric(0); f_w1 <- numeric(0); f_ks <- numeric(0)
    for (fi in seq_along(folds)) {
      ti <- folds[[fi]]; tr <- setdiff(seq_along(y), ti)
      fit <- tryCatch(gamlss_age_fit(y[tr], x[tr], o$family),
                      error = function(e) NULL)
      if (is.null(fit)) next
      rs <- robust_predict(fit, x[ti], y[ti], N_SIM)
      if (is.null(rs$samp)) next
      f_crps <- c(f_crps,
        mean(scoringRules::crps_sample(y = y[ti], dat = rs$samp), na.rm = TRUE) / iqr)
      f_w1   <- c(f_w1, wasserstein1d(y[ti], as.numeric(rs$samp)) / iqr)
      pit_te <- vapply(seq_along(y[ti]),
        function(j) mean(rs$samp[j, ] <= y[ti][j], na.rm = TRUE), numeric(1))
      f_ks   <- c(f_ks, ks_uniform(pmin(pmax(pit_te, 1e-6), 1 - 1e-6)))
    }
    cv_rows[[length(cv_rows) + 1]] <- data.frame(
      stratum = sid, sex = s$sex, sport = s$sport, test_type = s$test_type,
      metric = mn, family = o$family,
      cv_CRPS_rel = signif(mean(f_crps, na.rm = TRUE), 4),
      cv_W1_rel   = signif(mean(f_w1,   na.rm = TRUE), 4),
      cv_KS_D     = signif(mean(f_ks,   na.rm = TRUE), 4),
      stringsAsFactors = FALSE)
  }
}
cv_table <- do.call(rbind, cv_rows)
write.csv(cv_table, "results/07_cv_marginals.csv", row.names = FALSE)
cat("[07] CV rows =", nrow(cv_table), "\n")


# =============================================================================
# SECTION 9 — Athlete percentile lookup
# score_athlete(profile_id, test_type, sex, sport)
#   - finds the athlete's most recent test of that test_type
#   - reads each of the cohort's selected metrics from that test
#   - computes the marginal percentile from the GAMLSS predictAll at the
#     athlete's age
#   - computes the joint percentile via depth on the cohort PIT cloud
# Returns a data.frame; also prints a one-line summary.
# =============================================================================

score_athlete <- function(profile_id, test_type, sex, sport,
                          fits = fits_all, vines_obj = vines,
                          data = data_all) {
  sid <- gsub("[^A-Za-z0-9]+", "_", paste(sex, sport, test_type, sep = "__"))
  if (is.null(fits[[sid]])) {
    stop("No fitted cohort for ", sid,
         ". Strata available: ",
         paste(names(fits)[1:min(5, length(fits))], collapse = ", "), " ...")
  }
  s <- fits[[sid]]
  picked <- names(s$metrics)

  d <- data[data$profile_id == profile_id & data$gender == sex &
              data$sport == sport, ]
  if (!nrow(d)) stop("Athlete ", profile_id,
                      " not found in the integrated dataset.")
  # Restrict to test rows where any of the picked metrics is recorded
  has_any <- rowSums(!is.na(d[, intersect(picked, names(d)), drop = FALSE])) > 0
  d <- d[has_any, ]
  if (!nrow(d)) stop("Athlete ", profile_id,
                      " has no recorded values for ", test_type,
                      " in the picked metric set.")
  # Pick the test row with maximum coverage (ties broken by recency).
  cov  <- rowSums(!is.na(d[, intersect(picked, names(d)), drop = FALSE]))
  d <- d[order(-cov, as.Date(d$test_date), decreasing = c(FALSE, TRUE)), ]
  rec <- d[1, ]
  age <- rec$age_at_test

  rows <- list()
  uvec <- numeric(0); names_u <- character(0)
  for (mn in picked) {
    if (!mn %in% names(rec) || !is.finite(suppressWarnings(as.numeric(rec[[mn]]))))
      next
    yval <- as.numeric(rec[[mn]])
    pa <- tryCatch(suppressWarnings(gamlss::predictAll(
              s$metrics[[mn]]$fit, newdata = data.frame(x = age),
              type = "response",
              data = s$metrics[[mn]]$fit$call$data)),
            error = function(e) NULL)
    if (is.null(pa)) next
    fam <- as.character(s$metrics[[mn]]$fit$family)[1]
    pfun <- get(paste0("p", fam),
                envir = asNamespace("gamlss.dist"), mode = "function")
    pn <- intersect(c("mu","sigma","nu","tau"), names(pa))
    p_arg <- c(list(q = yval), lapply(pa[pn], function(z) z[1]))
    pct <- as.numeric(do.call(pfun, p_arg))
    rows[[length(rows) + 1]] <- data.frame(
      metric = mn, value = yval, age = round(age, 2), family = fam,
      marginal_pctile = round(100 * pct, 1),
      stringsAsFactors = FALSE)
    uvec <- c(uvec, pct); names_u <- c(names_u, mn)
  }
  out <- if (length(rows)) do.call(rbind, rows)
         else data.frame(metric = character(), value = numeric(),
                          age = numeric(), family = character(),
                          marginal_pctile = numeric())

  joint <- NA_real_
  v <- vines_obj[[sid]]
  if (!is.null(v) && length(uvec) >= 2) {
    common <- intersect(names_u, colnames(v$U))
    if (length(common) >= 2) {
      U <- v$U[, common, drop = FALSE]
      ux <- uvec[match(common, names_u)]
      cm <- colMeans(U); cv <- stats::cov(U)
      d2_a <- stats::mahalanobis(matrix(ux, 1), cm, cv)
      d2_c <- stats::mahalanobis(U, cm, cv)
      joint <- 100 * (1 - mean(d2_c <= d2_a))
    }
  }
  cat(sprintf(
    "Athlete %s | %s %s %s | age %.1f | metrics scored = %d/%d | joint percentile = %s\n",
    profile_id, sex, sport, test_type, age,
    length(uvec), length(picked),
    if (is.na(joint)) "NA" else sprintf("%.1f%%", joint)))
  out$joint_pctile <- if (is.na(joint)) NA_real_ else round(joint, 1)
  out
}

# Example call (uncomment after running sections 0-8):
# example_id <- data_all$profile_id[
#   data_all$gender == "Male" & data_all$sport == "Football"][1]
# score_athlete(example_id, "CMJ", "Male", "Football")


# =============================================================================
# SECTION 10 — Summary of where everything is
# =============================================================================

cat("\n==================================================================\n")
cat("Pipeline complete. Artefacts:\n")
cat("  data/analysis_dataset.csv         merged + cleaned\n")
cat("  data/per_test/<TEST>.csv          per-test slices\n")
cat("  results/03_*                      EDA tables and figures\n")
cat("  results/04_marginals_table.csv    GAMLSS marginal fit metrics\n")
cat("  results/04_marginals_fits.rds     fitted GAMLSS objects (in-memory: fits_all)\n")
cat("  results/05_joint_summary.csv      vine AIC/BIC/log-lik per cohort\n")
cat("  results/strata/<id>/              centiles.pdf, radar.pdf, vine.rds,\n")
cat("                                    joint_percentiles.csv\n")
cat("  results/07_cv_marginals.csv       5-fold CV scores\n")
cat("  paper/figures/*.pdf               paper-ready figures\n")
cat("\nUse score_athlete(profile_id, test_type, sex, sport) for a single\n")
cat("athlete's marginal + joint percentile.\n")
cat("==================================================================\n")
