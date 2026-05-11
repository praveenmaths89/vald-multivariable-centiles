# ------------------------------------------------------------
# 04_marginals_gamlss.R
# Fit GAMLSS marginals y ~ pb(age), sigma ~ pb(age) for each
# selected metric within each (Sex × Sport × Test) stratum.
# Distributional fit metrics are computed per metric:
#   CRPS, CRPS/IQR, Wasserstein-1, W1/IQR, KS_D, CvM_W2, AD_A2.
# Outputs:
#   results/04_marginals_table.csv
#   results/04_marginals_fits.rds
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

per_test_dir <- file.path(DATA_DIR, "per_test")
fits_path    <- file.path(RESULTS_DIR, "04_marginals_fits.rds")
table_path   <- file.path(RESULTS_DIR, "04_marginals_table.csv")

select_metrics <- function(d, age_col = "age_at_test") {
  meta_drop <- c("profile_id", "name", "gender", "sport", "dob",
                 "age_at_test", "test_date", "test_id", "weight",
                 "height_cm", "weight_kg")
  blacklist <- "(?i)Bodyweight|Body\\.?Mass|BM\\.?Index|Sample\\.Rate|Test\\.Duration"
  cand <- setdiff(names(d), meta_drop)
  cand <- cand[!grepl(blacklist, cand, perl = TRUE)]
  for (m in cand) d[[m]] <- suppressWarnings(as.numeric(d[[m]]))
  cand <- cand[vapply(cand, function(m) {
    v <- d[[m]]
    sum(is.finite(v)) >= max(15L, 0.5 * nrow(d)) &&
      stats::sd(v, na.rm = TRUE) > 0
  }, logical(1))]
  if (length(cand) < 2) return(character())

  X_imp <- as.data.frame(lapply(cand, function(m) {
    v <- d[[m]]; v[!is.finite(v)] <- stats::median(v, na.rm = TRUE); v
  }))
  names(X_imp) <- cand
  rf <- tryCatch(suppressWarnings(randomForest::randomForest(
    x = X_imp, y = d[[age_col]], ntree = 300, importance = TRUE
  )), error = function(e) NULL)
  if (is.null(rf)) return(head(cand, MAX_METRICS))
  imp <- randomForest::importance(rf, type = 1)
  rownames(imp)[order(-imp[, 1])][seq_len(min(MAX_METRICS, nrow(imp)))]
}

fit_one_metric <- function(y, x, n_sim = N_SIM) {
  has_neg <- any(y < 0, na.rm = TRUE)
  fam <- tryCatch({
    suppressWarnings({
      fd <- quiet_call(gamlss.dist::fitDist(
        y, type = if (has_neg) "realAll" else "realplus",
        k = log(length(y)), trace = FALSE
      ))
      fd$family[1]
    })
  }, error = function(e) if (has_neg) "NO" else "BCCGo")

  fit <- tryCatch(
    suppressWarnings(gamlss_age_fit(y, x, fam)),
    error = function(e) tryCatch(
      suppressWarnings(gamlss_age_fit(y, x, "NO")),
      error = function(e2) NULL
    )
  )
  if (is.null(fit)) return(NULL)
  fam_used <- as.character(fit$family)[1]

  rs <- robust_sample_predictive(fit, x, y, n_sim, fam_used)
  fit <- rs$fit; fam_used <- rs$fam; samp <- rs$samp

  pit <- pnorm(stats::residuals(fit))
  iqr_y <- stats::IQR(y, na.rm = TRUE)
  iqr_y <- if (is.finite(iqr_y) && iqr_y > 0) iqr_y else stats::sd(y, na.rm = TRUE) + 1e-9

  crps_v <- if (!is.null(samp)) tryCatch(
    mean(scoringRules::crps_sample(y = y, dat = samp), na.rm = TRUE),
    error = function(e) NA_real_) else NA_real_
  wass_v <- if (!is.null(samp)) wasserstein1d(y, as.numeric(samp)) else NA_real_

  list(
    fit = fit, family = fam_used, n = length(y), pit = pit,
    samp = samp, x = x, y = y,
    AIC = stats::AIC(fit), BIC = stats::BIC(fit),
    CRPS = crps_v, CRPS_rel = crps_v / iqr_y,
    Wasserstein1 = wass_v, W1_rel = wass_v / iqr_y,
    KS_D   = ks_uniform(pit),
    CvM_W2 = cvm_uniform(pit),
    AD_A2  = ad_uniform(pit)
  )
}

main <- function() {
  log_step("Stage 04 — GAMLSS marginals per stratum")
  env03 <- new.env()
  sys.source(file.path(CODE_DIR, "03_eda_strata.R"), envir = env03)
  prep <- env03$main()
  strata <- prep$strata
  eligible <- prep$eligible

  rows <- list()
  fits <- list()
  for (i in seq_len(nrow(eligible))) {
    sex   <- eligible$gender[i]
    sport <- eligible$sport[i]
    test  <- eligible$test_type[i]
    n_eli <- eligible$n[i]
    d <- strata[[test]]
    if (is.null(d)) next
    d <- d[d$gender == sex & d$sport == sport, , drop = FALSE]
    if (nrow(d) < MIN_N_STRATUM) next

    metrics <- select_metrics(d)
    if (!length(metrics)) next
    log_step(fmt_stratum(sex, sport, test, nrow(d)),
             "  metrics: ", paste(metrics, collapse = ", "))

    stratum_id <- sprintf("%s__%s__%s", slugify(sex), slugify(sport), slugify(test))
    fits[[stratum_id]] <- list(sex = sex, sport = sport, test_type = test,
                                metrics = list())

    for (m in metrics) {
      yvec <- suppressWarnings(as.numeric(d[[m]]))
      xvec <- d$age_at_test
      keep <- is.finite(yvec) & is.finite(xvec)
      if (sum(keep) < MIN_N_STRATUM) next
      fit_obj <- tryCatch(fit_one_metric(yvec[keep], xvec[keep]),
                          error = function(e) NULL)
      if (is.null(fit_obj)) next
      fits[[stratum_id]]$metrics[[m]] <- fit_obj
      rows[[length(rows) + 1]] <- data.frame(
        sex = sex, sport = sport, test_type = test, metric = m,
        n = fit_obj$n, family = fit_obj$family,
        AIC = round(fit_obj$AIC, 1), BIC = round(fit_obj$BIC, 1),
        CRPS         = signif(fit_obj$CRPS, 4),
        CRPS_rel     = signif(fit_obj$CRPS_rel, 4),
        Wasserstein1 = signif(fit_obj$Wasserstein1, 4),
        W1_rel       = signif(fit_obj$W1_rel, 4),
        KS_D         = signif(fit_obj$KS_D, 4),
        CvM_W2       = signif(fit_obj$CvM_W2, 4),
        AD_A2        = signif(fit_obj$AD_A2, 4),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) stop("No marginals fit. Check eligibility filter.")
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  utils::write.csv(tab, table_path, row.names = FALSE)
  saveRDS(fits, fits_path)
  log_step("Wrote ", table_path, " and ", fits_path)
  log_step(sprintf("Strata fit: %d  rows: %d", length(fits), nrow(tab)))
  invisible(list(table = tab, fits = fits, eligible = eligible))
}

if (sys.nframe() == 0) invisible(main())
