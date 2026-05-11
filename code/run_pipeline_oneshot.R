# ===========================================================
# VALD Multivariable Centile Pipeline — one-shot script
# -----------------------------------------------------------
# Open in RStudio, place the cursor on the first executable
# line, and step through with Cmd+Enter (Mac) / Ctrl+Enter.
# Every intermediate object stays in the global environment
# so you can inspect it after any line.
#
# Sections (use the RStudio navigator, top-right of editor):
#   0  Setup, paths, packages, tunables
#   1  (optional) live VALD API pull
#   2  Integrate sport / gender → analysis_dataset.csv
#   3  EDA per (Sex × Sport × Test) cohort
#   4  GAMLSS marginals + distributional scoring
#   5  Vine copulas on PIT residuals
#   6  Multivariable centiles + radar + joint percentile
#   7  5-fold CV
#   8  Paper figures, tables, numbers.tex
#   9  Compile LaTeX paper (call from shell)
# ===========================================================


# ===========================================================
# 0  Setup ---------------------------------------------------
# ===========================================================

PROJECT_DIR <- "~/Downloads/VALD_Multivariate_Centiles_STAIX2026"
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
setwd(PROJECT_DIR)

CODE_DIR    <- file.path(PROJECT_DIR, "code")
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
PAPER_FIGS  <- file.path(PROJECT_DIR, "paper", "figures")
PAPER_TABS  <- file.path(PROJECT_DIR, "paper", "tables")
LOG_DIR     <- file.path(PROJECT_DIR, "logs")
for (d in c(DATA_DIR, RESULTS_DIR, PAPER_FIGS, PAPER_TABS, LOG_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

SEED          <- 2026L; set.seed(SEED)
MIN_N_STRATUM <- 25L     # min observations per (Sex × Sport × Test) cohort
MAX_METRICS   <- 4L      # multivariable dimension (top-k by RF importance)
N_SIM         <- 300L    # predictive draws per observation
CV_FOLDS      <- 5L

PAL <- list(female = "#d7301f", male = "#0570b0",
            band   = "#3690c0", median = "#08519c", ref = "firebrick")

# --- Required packages.  Set this once if you need to install:
# install.packages(c(
#   "dplyr","tidyr","data.table","ggplot2","gridExtra","scales",
#   "gamlss","gamlss.dist","VineCopula","randomForest",
#   "scoringRules","fmsb","caret","reshape2"
# ), repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(data.table)
  library(ggplot2); library(gridExtra); library(scales)
  library(gamlss); library(gamlss.dist); library(VineCopula)
  library(randomForest); library(scoringRules); library(fmsb)
  library(reshape2)
})
theme_set(theme_minimal(base_size = 10))
HAS_DEPTHPROC <- requireNamespace("DepthProc", quietly = TRUE)
HAS_VALDR     <- requireNamespace("valdr",     quietly = TRUE)

# --- Tiny helpers (inline so nothing depends on helpers.R) ---
log_step <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
                              ..., "\n", sep = "")
slugify  <- function(x) gsub("[^A-Za-z0-9]+", "_", trimws(x))
clean_sport  <- function(x) trimws(gsub("\\s+(Overall(\\s+Data)?)$", "", x,
                                         ignore.case = TRUE))
clean_gender <- function(x) {
  x <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("M", "MALE")   ~ "Male",
    x %in% c("F", "FEMALE") ~ "Female",
    TRUE                    ~ NA_character_
  )
}
quiet_call <- function(expr) {
  zz <- tempfile(); sink(zz, type = "output")
  on.exit({ sink(type = "output"); unlink(zz) })
  force(expr)
}

# GAMLSS fit using do.call so $call$family is a literal string
# (predict.gamlss(newdata=) breaks otherwise).
gamlss_age_fit <- function(y_vec, x_vec, family_chr,
                           n.cyc = 120L, quiet = TRUE) {
  dat <- data.frame(y = y_vec, x = x_vec)
  fit_call <- list(formula = y ~ gamlss::pb(x),
                   sigma.formula = ~ gamlss::pb(x),
                   family = family_chr, data = dat, trace = FALSE,
                   control = gamlss::gamlss.control(n.cyc = n.cyc))
  if (quiet) quiet_call(do.call(gamlss::gamlss, fit_call))
  else       do.call(gamlss::gamlss, fit_call)
}

sample_predictive <- function(fit, x_new, n_sim = 300) {
  fam  <- as.character(fit$family)[1]
  rfun <- tryCatch(get(paste0("r", fam),
                       envir = asNamespace("gamlss.dist"), mode = "function"),
                   error = function(e) NULL)
  if (is.null(rfun)) return(NULL)
  pa <- tryCatch(suppressWarnings(gamlss::predictAll(
          fit, newdata = data.frame(x = x_new),
          type = "response", data = fit$call$data)),
        error = function(e) NULL)
  if (is.null(pa)) return(NULL)
  pn <- intersect(c("mu", "sigma", "nu", "tau"), names(pa))
  if (!length(pn)) return(NULL)
  n  <- length(x_new)
  out <- matrix(NA_real_, n, n_sim)
  for (s in seq_len(n_sim)) {
    args <- c(list(n = n), lapply(pa[pn], identity))
    out[, s] <- tryCatch(do.call(rfun, args),
                         error = function(e) rep(NA_real_, n))
  }
  out
}

wasserstein1d <- function(x, y, n_grid = 1024L) {
  x <- sort(as.numeric(x[is.finite(x)]))
  y <- sort(as.numeric(y[is.finite(y)]))
  if (!length(x) || !length(y)) return(NA_real_)
  u  <- seq(0.5 / n_grid, 1 - 0.5 / n_grid, length.out = n_grid)
  qx <- stats::quantile(x, probs = u, names = FALSE, type = 8)
  qy <- stats::quantile(y, probs = u, names = FALSE, type = 8)
  mean(abs(qx - qy))
}
ks_uniform <- function(p) {
  p <- p[is.finite(p)]; if (length(p) < 5) return(NA_real_)
  unname(suppressWarnings(stats::ks.test(p, "punif")$statistic))
}
cvm_uniform <- function(p) {
  p <- sort(as.numeric(p[is.finite(p)])); n <- length(p)
  if (n < 5) return(NA_real_)
  i <- seq_len(n); sum((p - (2*i-1)/(2*n))^2) + 1/(12*n)
}
ad_uniform <- function(p) {
  p <- sort(as.numeric(p[is.finite(p)])); n <- length(p)
  if (n < 5) return(NA_real_)
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10); i <- seq_len(n)
  -n - sum((2*i-1)/n * (log(p) + log(1 - rev(p))))
}

# Robust predictive sampling with NO-fallback + winsorisation
robust_sample_predictive <- function(fit, x, y, n_sim, fam_used) {
  yq   <- stats::quantile(y, c(0.005, 0.5, 0.995), na.rm = TRUE, names = FALSE)
  span <- max(yq[3] - yq[1], stats::sd(y, na.rm = TRUE), 1e-6)
  cap_lo <- yq[1] - 8 * span; cap_hi <- yq[3] + 8 * span
  samp <- tryCatch(sample_predictive(fit, x, n_sim), error = function(e) NULL)
  blew <- !is.null(samp) && (
    !all(is.finite(samp)) ||
    max(abs(samp), na.rm = TRUE) > 1e3 * (max(abs(y), na.rm = TRUE) + 1))
  if (blew && fam_used != "NO") {
    fit_no <- tryCatch(suppressWarnings(gamlss_age_fit(y, x, "NO")),
                       error = function(e) NULL)
    if (!is.null(fit_no)) {
      fit <- fit_no; fam_used <- "NO"
      samp <- tryCatch(sample_predictive(fit, x, n_sim),
                       error = function(e) NULL)
    }
  }
  if (!is.null(samp)) {
    samp[!is.finite(samp)] <- NA_real_
    samp[samp < cap_lo] <- cap_lo; samp[samp > cap_hi] <- cap_hi
  }
  list(fit = fit, fam = fam_used, samp = samp)
}


# ===========================================================
# 1  (optional) Live VALD API pull  --------------------------
# ===========================================================
# Skip this block if you already have the integrated CSV cached.

DO_LIVE_PULL <- FALSE   # flip to TRUE and set env vars below to pull live

if (DO_LIVE_PULL && HAS_VALDR) {
  Sys.setenv(VALD_CLIENT_ID     = "PASTE_HERE")
  Sys.setenv(VALD_CLIENT_SECRET = "PASTE_HERE")
  Sys.setenv(VALD_TENANT_ID     = "PASTE_HERE")
  Sys.setenv(VALD_REGION        = "aue")
  Sys.setenv(VALD_START_DATE    = "2020-01-01T00:00:00Z")
  options(timeout = 18000)
  valdr::set_credentials(
    client_id     = Sys.getenv("VALD_CLIENT_ID"),
    client_secret = Sys.getenv("VALD_CLIENT_SECRET"),
    tenant_id     = Sys.getenv("VALD_TENANT_ID"),
    region        = Sys.getenv("VALD_REGION")
  )
  valdr::set_start_date(Sys.getenv("VALD_START_DATE"))
  raw <- valdr::get_forcedecks_data(include_attributes = TRUE)
  saveRDS(raw, file.path(DATA_DIR, "raw_forcedecks.rds"))
  log_step("Pulled VALD: profiles=", nrow(raw$profiles),
           " tests=", nrow(raw$tests), " trials=", nrow(raw$trials))
}


# ===========================================================
# 2  Integrate sport / gender  -------------------------------
# ===========================================================
# Falls back to a pre-integrated CSV if we don't have the raw API pull.

EXTERNAL_INTEGRATED_CSV <- "/Users/apple/Downloads/SAI_Mumbai/VALD_integrated_output/20260508/01_VALD_integrated_complete.csv"

if (file.exists(file.path(DATA_DIR, "raw_forcedecks.rds"))) {
  raw <- readRDS(file.path(DATA_DIR, "raw_forcedecks.rds"))
  trials <- raw$trials; tests <- raw$tests; profs <- raw$profiles
  trials$athleteId <- as.character(trials$athleteId)
  joined <- dplyr::left_join(trials, tests, by = "testId")
  agg <- joined |>
    dplyr::mutate(value  = suppressWarnings(as.numeric(value)),
                  weight = suppressWarnings(as.numeric(weight))) |>
    dplyr::group_by(athleteId, testId, testType, recordedUTC,
                    recordedDateOffset, trialLimb, definition_name) |>
    dplyr::summarise(mean_result = mean(value, na.rm = TRUE),
                     mean_weight = mean(weight, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(ts_utc   = lubridate::ymd_hms(recordedUTC),
                  ts_local = ts_utc + lubridate::minutes(recordedDateOffset),
                  test_date = as.Date(ts_local))
  integrated <- agg |>
    dplyr::select(profile_id = athleteId, test_id = testId,
                  test_type = testType, test_date,
                  limb = trialLimb, metric = definition_name,
                  value = mean_result, weight = mean_weight) |>
    tidyr::pivot_wider(
      id_cols = c(profile_id, test_date, test_id, weight),
      names_from = c(metric, limb, test_type), values_from = value,
      names_glue = "{metric}_{limb}_{test_type}",
      names_repair = "universal")
  # join external sport / gender
  sports_dir <- "/Users/apple/Downloads/VALD Data_sports_gender"
  sf <- list.files(sports_dir, pattern = "\\.csv$", full.names = TRUE)
  sports <- do.call(rbind, lapply(sf, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    pc <- names(df)[grepl("profile.*id", names(df), ignore.case = TRUE)][1]
    sc <- names(df)[grepl("^gender$",     names(df), ignore.case = TRUE)][1]
    if (is.na(pc) || is.na(sc)) return(NULL)
    sport <- sub("\\..*$", "",
                 sub("_All.*$", "", sub("_Overall.*$", "", basename(f))))
    data.frame(profile_id = as.character(df[[pc]]),
               gender = clean_gender(df[[sc]]),
               sport  = clean_sport(sport),
               stringsAsFactors = FALSE)
  })) |> dplyr::distinct()
  integrated <- dplyr::left_join(integrated, sports, by = "profile_id")
  integrated$age_at_test <- as.numeric(
    difftime(integrated$test_date, profs$dateOfBirth[
      match(integrated$profile_id, as.character(profs$profileId))],
      units = "days")) / 365.25
} else {
  log_step("No raw VALD pull found — using cached integrated CSV at\n  ",
           EXTERNAL_INTEGRATED_CSV)
  integrated <- read.csv(EXTERNAL_INTEGRATED_CSV,
                         stringsAsFactors = FALSE, check.names = FALSE)
  integrated$gender <- clean_gender(integrated$gender)
  integrated$sport  <- clean_sport(integrated$sport)
}
integrated <- integrated[!is.na(integrated$gender) &
                          !is.na(integrated$sport), , drop = FALSE]
write.csv(integrated, file.path(DATA_DIR, "analysis_dataset.csv"),
          row.names = FALSE)
log_step("analysis_dataset.csv  rows=", nrow(integrated),
         "  cols=", ncol(integrated))

# Per-test CSVs
meta <- intersect(c("profile_id","name","gender","sport","dob",
                    "age_at_test","test_date","test_id","weight",
                    "height_cm","weight_kg"), names(integrated))
metric_cols <- setdiff(names(integrated), meta)
test_types  <- unique(sub(".*_([A-Z]+)$", "\\1", metric_cols))
test_types  <- test_types[nzchar(test_types) & !is.na(test_types)]
per_test_dir <- file.path(DATA_DIR, "per_test")
dir.create(per_test_dir, recursive = TRUE, showWarnings = FALSE)
for (tt in test_types) {
  cols <- metric_cols[grepl(paste0("_", tt, "$"), metric_cols)]
  sub  <- integrated[, c(meta, cols), drop = FALSE]
  keep <- rowSums(!is.na(sub[, cols, drop = FALSE])) > 0
  if (sum(keep) >= 10)
    write.csv(sub[keep, , drop = FALSE],
              file.path(per_test_dir, sprintf("%s.csv", tt)), row.names = FALSE)
}
log_step("per-test CSVs: ", length(list.files(per_test_dir)))


# ===========================================================
# 3  EDA — eligible cohorts ----------------------------------
# ===========================================================

strata <- list()
for (f in list.files(per_test_dir, pattern = "\\.csv$", full.names = TRUE)) {
  tt <- tools::file_path_sans_ext(basename(f))
  d  <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d$gender <- clean_gender(d$gender); d$sport <- clean_sport(d$sport)
  d <- d[!is.na(d$gender) & !is.na(d$sport) &
          !is.na(d$age_at_test) & d$age_at_test >= 8 & d$age_at_test <= 60, ,
          drop = FALSE]
  strata[[tt]] <- d
}

ct_rows <- list()
for (tt in names(strata)) {
  d <- strata[[tt]]; if (!nrow(d)) next
  tab <- dplyr::count(d, gender, sport, name = "n"); tab$test_type <- tt
  ct_rows[[length(ct_rows) + 1]] <- tab
}
counts   <- do.call(rbind, ct_rows); counts <- counts[order(-counts$n), ]
eligible <- counts[counts$n >= MIN_N_STRATUM, , drop = FALSE]

write.csv(counts,   file.path(RESULTS_DIR, "03_stratum_counts.csv"),   row.names = FALSE)
write.csv(eligible, file.path(RESULTS_DIR, "03_eligible_strata.csv"), row.names = FALSE)
log_step("combos=", nrow(counts), "  eligible=", nrow(eligible))

p_counts <- ggplot(counts, aes(test_type, sport, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 3) +
  facet_wrap(~ gender) +
  scale_fill_viridis_c(option = "viridis", na.value = "grey90") +
  labs(title = "Eligible-cohort sample sizes",
       subtitle = paste0("Strata used downstream require n >= ", MIN_N_STRATUM),
       x = "Test type", y = "Sport") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(RESULTS_DIR, "03_counts_heatmap.png"),
       p_counts, width = 12, height = 5.5, dpi = 130)
ggsave(file.path(PAPER_FIGS, "fig_counts_heatmap.pdf"),
       p_counts, width = 11, height = 5.5)

# Inspect interactively:
# View(eligible); print(p_counts)


# ===========================================================
# 4  GAMLSS marginals  ---------------------------------------
# ===========================================================

select_metrics <- function(d, age_col = "age_at_test") {
  meta_drop <- c("profile_id","name","gender","sport","dob","age_at_test",
                 "test_date","test_id","weight","height_cm","weight_kg")
  blacklist <- "(?i)Bodyweight|Body\\.?Mass|BM\\.?Index|Sample\\.Rate|Test\\.Duration"
  cand <- setdiff(names(d), meta_drop)
  cand <- cand[!grepl(blacklist, cand, perl = TRUE)]
  for (m in cand) d[[m]] <- suppressWarnings(as.numeric(d[[m]]))
  cand <- cand[vapply(cand, function(m) {
    v <- d[[m]]; sum(is.finite(v)) >= max(15L, 0.5 * nrow(d)) &&
      stats::sd(v, na.rm = TRUE) > 0
  }, logical(1))]
  if (length(cand) < 2) return(character())
  X_imp <- as.data.frame(lapply(cand, function(m) {
    v <- d[[m]]; v[!is.finite(v)] <- stats::median(v, na.rm = TRUE); v
  })); names(X_imp) <- cand
  rf <- tryCatch(suppressWarnings(randomForest::randomForest(
    x = X_imp, y = d[[age_col]], ntree = 300, importance = TRUE)),
    error = function(e) NULL)
  if (is.null(rf)) return(head(cand, MAX_METRICS))
  imp <- randomForest::importance(rf, type = 1)
  rownames(imp)[order(-imp[, 1])][seq_len(min(MAX_METRICS, nrow(imp)))]
}

fit_one_metric <- function(y, x, n_sim = N_SIM) {
  has_neg <- any(y < 0, na.rm = TRUE)
  fam <- tryCatch(suppressWarnings({
    fd <- quiet_call(gamlss.dist::fitDist(
      y, type = if (has_neg) "realAll" else "realplus",
      k = log(length(y)), trace = FALSE))
    fd$family[1]
  }), error = function(e) if (has_neg) "NO" else "BCCGo")
  fit <- tryCatch(suppressWarnings(gamlss_age_fit(y, x, fam)),
                  error = function(e) tryCatch(
                    suppressWarnings(gamlss_age_fit(y, x, "NO")),
                    error = function(e2) NULL))
  if (is.null(fit)) return(NULL)
  fam_used <- as.character(fit$family)[1]
  rs <- robust_sample_predictive(fit, x, y, n_sim, fam_used)
  fit <- rs$fit; fam_used <- rs$fam; samp <- rs$samp

  pit <- pnorm(stats::residuals(fit))
  iqr_y <- max(stats::IQR(y, na.rm = TRUE), stats::sd(y, na.rm = TRUE), 1e-9)
  crps_v <- if (!is.null(samp)) tryCatch(
    mean(scoringRules::crps_sample(y = y, dat = samp), na.rm = TRUE),
    error = function(e) NA_real_) else NA_real_
  wass_v <- if (!is.null(samp)) wasserstein1d(y, as.numeric(samp)) else NA_real_

  list(fit = fit, family = fam_used, n = length(y), pit = pit,
       samp = samp, x = x, y = y,
       AIC = stats::AIC(fit), BIC = stats::BIC(fit),
       CRPS = crps_v, CRPS_rel = crps_v / iqr_y,
       Wasserstein1 = wass_v, W1_rel = wass_v / iqr_y,
       KS_D = ks_uniform(pit),
       CvM_W2 = cvm_uniform(pit), AD_A2 = ad_uniform(pit))
}

marg_rows <- list(); fits <- list()
for (i in seq_len(nrow(eligible))) {
  sex <- eligible$gender[i]; sport <- eligible$sport[i]; test <- eligible$test_type[i]
  d <- strata[[test]]; if (is.null(d)) next
  d <- d[d$gender == sex & d$sport == sport, , drop = FALSE]
  if (nrow(d) < MIN_N_STRATUM) next
  metrics <- select_metrics(d); if (!length(metrics)) next
  log_step(sprintf("%-6s | %-18s | %-7s | n=%d  metrics: %s",
                   sex, sport, test, nrow(d), paste(metrics, collapse = ", ")))
  sid <- sprintf("%s__%s__%s", slugify(sex), slugify(sport), slugify(test))
  fits[[sid]] <- list(sex = sex, sport = sport, test_type = test, metrics = list())
  for (m in metrics) {
    yv <- suppressWarnings(as.numeric(d[[m]])); xv <- d$age_at_test
    keep <- is.finite(yv) & is.finite(xv)
    if (sum(keep) < MIN_N_STRATUM) next
    fo <- tryCatch(fit_one_metric(yv[keep], xv[keep]), error = function(e) NULL)
    if (is.null(fo)) next
    fits[[sid]]$metrics[[m]] <- fo
    marg_rows[[length(marg_rows) + 1]] <- data.frame(
      sex = sex, sport = sport, test_type = test, metric = m,
      n = fo$n, family = fo$family,
      AIC = round(fo$AIC, 1), BIC = round(fo$BIC, 1),
      CRPS = signif(fo$CRPS, 4), CRPS_rel = signif(fo$CRPS_rel, 4),
      Wasserstein1 = signif(fo$Wasserstein1, 4), W1_rel = signif(fo$W1_rel, 4),
      KS_D = signif(fo$KS_D, 4),
      CvM_W2 = signif(fo$CvM_W2, 4), AD_A2 = signif(fo$AD_A2, 4),
      stringsAsFactors = FALSE)
  }
}

marg_table <- do.call(rbind, marg_rows); rownames(marg_table) <- NULL
write.csv(marg_table, file.path(RESULTS_DIR, "04_marginals_table.csv"), row.names = FALSE)
saveRDS(fits, file.path(RESULTS_DIR, "04_marginals_fits.rds"))
log_step("strata fit=", length(fits), "  marginal rows=", nrow(marg_table))


# ===========================================================
# 5  Vine copulas on PIT residuals  --------------------------
# ===========================================================

strata_dir <- file.path(RESULTS_DIR, "strata")
dir.create(strata_dir, recursive = TRUE, showWarnings = FALSE)

joint_rows <- list()
for (sid in names(fits)) {
  s <- fits[[sid]]; ms <- s$metrics
  if (length(ms) < 2) next
  n  <- min(vapply(ms, function(m) length(m$pit), integer(1)))
  if (n < MIN_N_STRATUM) next
  U  <- vapply(ms, function(m) m$pit[seq_len(n)], numeric(n))
  colnames(U) <- names(ms)
  U <- pmin(pmax(U, 1e-6), 1 - 1e-6)
  U <- U[stats::complete.cases(U), , drop = FALSE]
  vine <- tryCatch(suppressWarnings(VineCopula::RVineStructureSelect(
    U, familyset = NA, type = 0, selectioncrit = "AIC",
    indeptest = TRUE, level = 0.05)), error = function(e) NULL)
  if (is.null(vine)) { log_step("  vine failed: ", sid); next }
  sd <- file.path(strata_dir, sid); dir.create(sd, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(vine = vine, U = U), file.path(sd, "vine.rds"))
  joint_rows[[length(joint_rows) + 1]] <- data.frame(
    stratum = sid, sex = s$sex, sport = s$sport, test_type = s$test_type,
    n_obs = nrow(U), dim = ncol(U),
    vine_AIC = round(vine$AIC, 1), vine_BIC = round(vine$BIC, 1),
    vine_loglik = round(sum(VineCopula::RVineLogLik(U, vine)$loglik), 1),
    stringsAsFactors = FALSE)
  log_step(sprintf("  %s   AIC=%.1f  BIC=%.1f  n=%d  d=%d",
                   sid, vine$AIC, vine$BIC, nrow(U), ncol(U)))
}
joint_summary <- do.call(rbind, joint_rows)
write.csv(joint_summary, file.path(RESULTS_DIR, "05_joint_summary.csv"), row.names = FALSE)


# ===========================================================
# 6  Centile bands + radar + joint depth percentile ---------
# ===========================================================

centile_grid <- function(fit_obj, x_seq) {
  rs <- robust_sample_predictive(fit_obj$fit, x_seq, fit_obj$y,
                                  N_SIM, fit_obj$family)
  if (is.null(rs$samp)) return(NULL)
  qs <- t(apply(rs$samp, 1, function(v) stats::quantile(v,
        probs = c(.03,.10,.25,.50,.75,.90,.97), na.rm = TRUE, names = FALSE)))
  data.frame(age = x_seq,
             c03 = qs[,1], c10 = qs[,2], c25 = qs[,3], c50 = qs[,4],
             c75 = qs[,5], c90 = qs[,6], c97 = qs[,7])
}
joint_depth <- function(U) {
  if (HAS_DEPTHPROC) {
    dv <- tryCatch(DepthProc::depthTukey(U, U), error = function(e) NULL)
    if (!is.null(dv)) return(rank(dv) / length(dv))
  }
  md <- stats::mahalanobis(U, colMeans(U), stats::cov(U))
  1 - rank(md) / length(md)
}

for (sid in names(fits)) {
  s <- fits[[sid]]; if (length(s$metrics) < 2) next
  sd <- file.path(strata_dir, sid); dir.create(sd, recursive = TRUE, showWarnings = FALSE)
  age_lo <- min(unlist(lapply(s$metrics, function(m) min(m$x))), na.rm = TRUE)
  age_hi <- max(unlist(lapply(s$metrics, function(m) max(m$x))), na.rm = TRUE)
  x_seq <- seq(age_lo, age_hi, length.out = 120)
  panels <- lapply(names(s$metrics), function(nm) {
    o  <- s$metrics[[nm]]; cf <- centile_grid(o, x_seq)
    if (is.null(cf)) return(NULL)
    obs <- data.frame(x = o$x, y = o$y)
    ggplot() +
      geom_point(data = obs, aes(x, y), color = "grey50", size = 0.9, alpha = 0.45) +
      geom_ribbon(data = cf, aes(age, ymin = c03, ymax = c97), fill = PAL$band, alpha = 0.13) +
      geom_ribbon(data = cf, aes(age, ymin = c10, ymax = c90), fill = PAL$band, alpha = 0.20) +
      geom_ribbon(data = cf, aes(age, ymin = c25, ymax = c75), fill = PAL$band, alpha = 0.30) +
      geom_line(data = cf, aes(age, c50), color = PAL$median, linewidth = 0.9) +
      labs(title = nm, x = "Age (years)", y = NULL) +
      theme(plot.title = element_text(size = 8))
  })
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (length(panels)) {
    g <- gridExtra::arrangeGrob(grobs = panels, ncol = min(2, length(panels)))
    ggsave(file.path(sd, "centiles.pdf"), g,
           width = 10, height = 3.4 * ceiling(length(panels)/2), limitsize = FALSE)
    ggsave(file.path(sd, "centiles.png"), g,
           width = 10, height = 3.4 * ceiling(length(panels)/2),
           dpi = 130, limitsize = FALSE)
  }

  # joint depth percentile per athlete
  n  <- min(vapply(s$metrics, function(m) length(m$pit), integer(1)))
  U  <- vapply(s$metrics, function(m) m$pit[seq_len(n)], numeric(n))
  colnames(U) <- names(s$metrics)
  U <- U[stats::complete.cases(U), , drop = FALSE]
  if (nrow(U) >= MIN_N_STRATUM) {
    depth_pct <- round(joint_depth(U) * 100, 1)
    write.csv(data.frame(joint_pctile = depth_pct, U),
              file.path(sd, "joint_percentiles.csv"), row.names = FALSE)
    # 3-athlete radar (high / median / low joint depth)
    ord <- order(-depth_pct)
    picks <- c(ord[1], ord[ceiling(length(ord)/2)], ord[length(ord)])
    titles <- c("High joint depth", "Median joint depth", "Low joint depth")
    clr <- c("#238b45", PAL$median, "#cb181d")
    png(file.path(sd, "radar.png"), width = 1500, height = 500, res = 130)
    op <- par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))
    for (k in seq_along(picks)) {
      v <- as.numeric(U[picks[k], ]) * 100
      rdf <- data.frame(rbind(rep(100, ncol(U)), rep(0, ncol(U)), v))
      colnames(rdf) <- gsub("_Both_|_Right_|_Left_", "_", colnames(U))
      tryCatch(fmsb::radarchart(
        rdf, pcol = clr[k], pfcol = grDevices::adjustcolor(clr[k], 0.25),
        plwd = 2, plty = 1, cglcol = "grey75", cglty = 1,
        axislabcol = "grey40", vlcex = 0.5,
        title = sprintf("%s\n(joint pct %.0f%%)", titles[k], depth_pct[picks[k]])),
        error = function(e) NULL)
    }
    par(op); dev.off()
  }
}
log_step("Stage 6 done — per-cohort outputs in results/strata/")


# ===========================================================
# 7  5-fold CV  ----------------------------------------------
# ===========================================================
if (!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret", repos = "https://cloud.r-project.org", quiet = TRUE)
}

cv_one_metric <- function(y, x, family_chr, k = CV_FOLDS) {
  set.seed(SEED); folds <- caret::createFolds(x, k = k, returnTrain = FALSE)
  iqr_y <- max(stats::IQR(y, na.rm = TRUE), stats::sd(y, na.rm = TRUE), 1e-9)
  per <- list()
  for (fi in seq_along(folds)) {
    ti <- folds[[fi]]; tr <- setdiff(seq_along(y), ti)
    out <- tryCatch({
      fit <- suppressWarnings(gamlss_age_fit(y[tr], x[tr], family_chr))
      fam <- as.character(fit$family)[1]
      rs <- robust_sample_predictive(fit, x[ti], y[ti], N_SIM, fam)
      if (is.null(rs$samp)) return(NULL)
      fhat <- vapply(seq_along(y[ti]),
                     function(j) mean(rs$samp[j, ] <= y[ti][j], na.rm = TRUE),
                     numeric(1))
      pit  <- pmin(pmax(fhat, 1e-6), 1 - 1e-6)
      list(crps = mean(scoringRules::crps_sample(y = y[ti], dat = rs$samp),
                       na.rm = TRUE) / iqr_y,
           w1   = wasserstein1d(y[ti], as.numeric(rs$samp)) / iqr_y,
           ks   = ks_uniform(pit))
    }, error = function(e) NULL)
    if (!is.null(out)) per[[length(per) + 1]] <- out
  }
  if (!length(per)) return(c(crps = NA_real_, w1 = NA_real_, ks = NA_real_))
  c(crps = mean(vapply(per, `[[`, numeric(1), "crps"), na.rm = TRUE),
    w1   = mean(vapply(per, `[[`, numeric(1), "w1"),   na.rm = TRUE),
    ks   = mean(vapply(per, `[[`, numeric(1), "ks"),   na.rm = TRUE))
}

cv_rows <- list()
for (sid in names(fits)) {
  s <- fits[[sid]]
  for (mn in names(s$metrics)) {
    o <- s$metrics[[mn]]
    cv <- cv_one_metric(o$y, o$x, o$family)
    cv_rows[[length(cv_rows) + 1]] <- data.frame(
      stratum = sid, sex = s$sex, sport = s$sport, test_type = s$test_type,
      metric = mn, family = o$family,
      cv_CRPS_rel = signif(cv["crps"], 4),
      cv_W1_rel   = signif(cv["w1"],   4),
      cv_KS_D     = signif(cv["ks"],   4),
      stringsAsFactors = FALSE)
  }
}
cv_table <- do.call(rbind, cv_rows)
write.csv(cv_table, file.path(RESULTS_DIR, "07_cv_marginals.csv"), row.names = FALSE)
log_step("CV done. rows=", nrow(cv_table))


# ===========================================================
# 8  Paper figures, tables, numbers.tex ----------------------
# ===========================================================

mt <- marg_table; js <- joint_summary; cv <- cv_table

mean_tbl <- mt |>
  dplyr::group_by(sex, sport, test_type) |>
  dplyr::summarise(mean_CRPS_rel = mean(CRPS_rel, na.rm = TRUE),
                   mean_KS_D     = mean(KS_D,     na.rm = TRUE),
                   mean_AIC      = mean(AIC,      na.rm = TRUE),
                   n             = max(n), .groups = "drop")

p_heat <- ggplot(mean_tbl, aes(test_type, sport, fill = pmin(mean_CRPS_rel, 5))) +
  geom_tile(color = "white") +
  geom_text(aes(label = formatC(mean_CRPS_rel, format = "g", digits = 2)), size = 3) +
  facet_wrap(~ sex) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey90",
                        name = "Mean CRPS / IQR") +
  labs(title = "Marginal fit quality across cohorts",
       subtitle = "Each cell averages 4 metrics; lower is better.",
       x = "Test type", y = "Sport") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"))
ggsave(file.path(PAPER_FIGS, "fig_crps_heatmap.pdf"), p_heat, width = 10, height = 5.5)

n_test_levels <- length(unique(mt$test_type))
shape_pal <- (15:25)[seq_len(n_test_levels)]
mt_clip <- mt; mt_clip$CRPS_rel_c <- pmin(mt$CRPS_rel, 3)
p_cs <- ggplot(mt_clip, aes(KS_D, CRPS_rel_c, color = sex, shape = test_type)) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
  scale_shape_manual(values = shape_pal, name = "Test") +
  labs(title = "Calibration vs sharpness",
       subtitle = "Lower KS_D = better calibration; lower CRPS = sharper.",
       x = "KS_D on PIT", y = "CRPS / IQR (clipped at 3)")
ggsave(file.path(PAPER_FIGS, "fig_calibration_sharpness.pdf"), p_cs, width = 9, height = 5.5)

paired <- mt |>
  dplyr::group_by(sport, test_type, metric) |>
  dplyr::summarise(F = ifelse(any(sex == "Female"), CRPS_rel[sex == "Female"][1], NA_real_),
                   M = ifelse(any(sex == "Male"),   CRPS_rel[sex == "Male"][1],   NA_real_),
                   .groups = "drop") |>
  dplyr::filter(!is.na(F), !is.na(M))
if (nrow(paired) > 0) {
  p_pair <- ggplot(paired, aes(M, F, color = sport, shape = test_type)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
    geom_point(size = 2.4, alpha = 0.85) +
    scale_shape_manual(values = shape_pal, name = "Test") +
    coord_equal() +
    labs(title = "Female vs male CRPS / IQR (matched on metric)",
         x = "Male  CRPS / IQR", y = "Female  CRPS / IQR")
  ggsave(file.path(PAPER_FIGS, "fig_paired_sex.pdf"), p_pair, width = 7, height = 5.5)
}

# CV-vs-in-sample
mer <- merge(mt[, c("sex","sport","test_type","metric","CRPS_rel","KS_D")],
             cv[, c("sex","sport","test_type","metric","cv_CRPS_rel","cv_KS_D")],
             by = c("sex","sport","test_type","metric"))
p_cv <- ggplot(mer, aes(pmin(CRPS_rel, 3), pmin(cv_CRPS_rel, 3),
                        color = sex, shape = test_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
  scale_shape_manual(values = shape_pal, name = "Test") +
  labs(title = "Out-of-sample vs in-sample CRPS / IQR",
       x = "In-sample CRPS / IQR", y = "5-fold CV CRPS / IQR")
ggsave(file.path(PAPER_FIGS, "fig_cv_vs_insample.pdf"), p_cv, width = 9, height = 5.5)

# Vine dependence example (largest cohort)
vd <- list.dirs(strata_dir, recursive = FALSE)
big <- vd[which.max(vapply(vd, function(d) {
  fp <- file.path(d, "vine.rds")
  if (!file.exists(fp)) return(0L); nrow(readRDS(fp)$U)
}, integer(1)))]
if (length(big)) {
  obj <- readRDS(file.path(big, "vine.rds"))
  pdf(file.path(PAPER_FIGS, "fig_vine_pair_contours.pdf"), width = 8, height = 8)
  par(mar = c(2,2,2,2))
  tryCatch(VineCopula::contour(obj$vine, margins = "norm", cex = 0.6),
           error = function(e) plot.new())
  title(sprintf("Pair-copula contours: %s", basename(big)))
  dev.off()
  emp_long <- reshape2::melt(cor(obj$U, method = "spearman"))
  p_dep <- ggplot(emp_long, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 2.7) +
    scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#cb181d",
                          midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
    labs(title = paste("Empirical Spearman dependence:", basename(big)),
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))
  ggsave(file.path(PAPER_FIGS, "fig_dependence_heatmap.pdf"), p_dep, width = 8, height = 6)
}

# PIT density by sex (pooled)
pit_df <- do.call(rbind, lapply(names(fits), function(sid) {
  s <- fits[[sid]]
  do.call(rbind, lapply(s$metrics, function(m) data.frame(sex = s$sex, pit = m$pit)))
}))
p_pit <- ggplot(pit_df, aes(pit, fill = sex, color = sex)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
  scale_fill_manual(values  = c(Female = PAL$female, Male = PAL$male)) +
  labs(title = "PIT residuals across all marginals",
       subtitle = "Uniform reference shown as the dashed line at y = 1",
       x = "PIT", y = "Density")
ggsave(file.path(PAPER_FIGS, "fig_pit_density.pdf"), p_pit, width = 8, height = 4.5)

# Example multivariable centiles (largest stratum)
pick_score <- vapply(fits, function(x) {
  if (length(x$metrics) == 0) return(0L)
  yv <- x$metrics[[1]]$y
  if (is.null(yv)) 0L else as.integer(length(yv) * length(x$metrics))
}, integer(1))
pick <- names(fits)[which.max(pick_score)]
src  <- file.path(strata_dir, pick, "centiles.pdf")
if (file.exists(src)) file.copy(src, file.path(PAPER_FIGS, "fig_example_centiles.pdf"), overwrite = TRUE)

# LaTeX tables
top_strata   <- mean_tbl |> dplyr::arrange(mean_CRPS_rel)             |> dplyr::slice_head(n = 12)
worst_strata <- mean_tbl |> dplyr::arrange(dplyr::desc(mean_CRPS_rel)) |> dplyr::slice_head(n = 8)

write_latex_table <- function(df, path, caption, label) {
  cn <- gsub("_", "\\\\_", colnames(df))
  con <- file(path, "w"); on.exit(close(con))
  cat("\\begin{table}[t]\n\\centering\n\\small\n", file = con)
  cat(sprintf("\\caption{%s}\n\\label{%s}\n", caption, label), file = con)
  cat(sprintf("\\begin{tabular}{%s}\n\\toprule\n",
              paste(c("l", rep("r", ncol(df) - 1)), collapse = "")), file = con)
  cat(paste(cn, collapse = " & "), " \\\\\n\\midrule\n", file = con)
  for (i in seq_len(nrow(df))) {
    row <- vapply(df[i, ], function(v)
      if (is.numeric(v)) formatC(v, format = "g", digits = 3) else as.character(v),
      character(1))
    cat(paste(row, collapse = " & "), " \\\\\n", file = con)
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n", file = con)
}
write_latex_table(top_strata,   file.path(PAPER_TABS, "tab_top_strata.tex"),
                  "Cohorts with the best mean fit quality (lower CRPS/IQR is better).", "tab:top")
write_latex_table(worst_strata, file.path(PAPER_TABS, "tab_worst_strata.tex"),
                  "Cohorts where marginal fit is most challenged.", "tab:worst")

# numbers.tex (inline numbers used by the paper)
n_strata  <- length(fits); n_combos <- nrow(mean_tbl); n_metrics <- nrow(mt)
sports_n  <- length(unique(mt$sport)); tests_n <- length(unique(mt$test_type))
mean_crps <- round(mean(mt$CRPS_rel, na.rm = TRUE), 3)
median_ks <- round(median(mt$KS_D,    na.rm = TRUE), 3)
best_str  <- top_strata$sport[1];   best_test  <- top_strata$test_type[1]
best_sex  <- top_strata$sex[1];     best_val   <- round(top_strata$mean_CRPS_rel[1], 3)
worst_str <- worst_strata$sport[1]; worst_test <- worst_strata$test_type[1]
worst_sex <- worst_strata$sex[1];   worst_val  <- round(worst_strata$mean_CRPS_rel[1], 3)

inline_path <- file.path(PROJECT_DIR, "paper", "numbers.tex")
con <- file(inline_path, "w")
cat(sprintf("\\newcommand{\\NumStrata}{%d}\n",     n_strata),  file = con)
cat(sprintf("\\newcommand{\\NumCombos}{%d}\n",     n_combos),  file = con)
cat(sprintf("\\newcommand{\\NumMetricRows}{%d}\n", n_metrics), file = con)
cat(sprintf("\\newcommand{\\NumSports}{%d}\n",     sports_n),  file = con)
cat(sprintf("\\newcommand{\\NumTests}{%d}\n",      tests_n),   file = con)
cat(sprintf("\\newcommand{\\MeanCRPSrel}{%.3f}\n", mean_crps), file = con)
cat(sprintf("\\newcommand{\\MedianKS}{%.3f}\n",    median_ks), file = con)
cat(sprintf("\\newcommand{\\BestStratum}{%s %s (%s)}\n",  best_str,  best_test,  best_sex),  file = con)
cat(sprintf("\\newcommand{\\BestCRPSrel}{%.3f}\n",  best_val),  file = con)
cat(sprintf("\\newcommand{\\WorstStratum}{%s %s (%s)}\n", worst_str, worst_test, worst_sex), file = con)
cat(sprintf("\\newcommand{\\WorstCRPSrel}{%.3f}\n", worst_val), file = con)
if (!is.null(js) && nrow(js)) {
  cat(sprintf("\\newcommand{\\NumVines}{%d}\n",   nrow(js)), file = con)
  cat(sprintf("\\newcommand{\\MedVineAIC}{%.1f}\n",
              median(js$vine_AIC, na.rm = TRUE)), file = con)
}
close(con)
log_step("Paper assets written.")


# ===========================================================
# 9  Compile the paper (run in shell)  -----------------------
# ===========================================================
# system("cd paper && latexmk -pdf main.tex")
# OR from a terminal:
#   cd ~/Downloads/VALD_Multivariate_Centiles_STAIX2026/paper
#   latexmk -pdf main.tex
#   open main.pdf

log_step("ALL DONE.")
