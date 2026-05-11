# ------------------------------------------------------------
# 07_validation.R
# Five-fold cross-validation of the marginal models per
# stratum, plus residual-PIT diagnostics for the joint vine.
# Distributional CV scores: out-of-fold CRPS_rel, W1_rel, KS_D.
# Output:
#   results/07_cv_marginals.csv
#   results/07_pit_rosenblatt.png  (per stratum, faceted)
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

fits_path <- file.path(RESULTS_DIR, "04_marginals_fits.rds")

cv_one_metric <- function(y, x, family_chr, k = CV_FOLDS) {
  set.seed(SEED)
  folds <- caret::createFolds(x, k = k, returnTrain = FALSE)
  iqr_y <- max(stats::IQR(y, na.rm = TRUE), stats::sd(y, na.rm = TRUE), 1e-9)
  per_fold <- list()
  for (fi in seq_along(folds)) {
    ti <- folds[[fi]]; tr <- setdiff(seq_along(y), ti)
    out <- tryCatch({
      fit <- suppressWarnings(gamlss_age_fit(y[tr], x[tr], family_chr))
      fam <- as.character(fit$family)[1]
      rs <- robust_sample_predictive(fit, x[ti], y[ti], N_SIM, fam)
      if (is.null(rs$samp)) return(NULL)
      pit <- {
        # Approximate PIT using gamlss::pdf -> not always available; use sample-based
        # F̂(y) = mean(samp <= y) per row
        fhat <- vapply(seq_along(y[ti]), function(j)
          mean(rs$samp[j, ] <= y[ti][j], na.rm = TRUE), numeric(1))
        pmin(pmax(fhat, 1e-6), 1 - 1e-6)
      }
      list(
        crps = mean(scoringRules::crps_sample(y = y[ti], dat = rs$samp), na.rm = TRUE) / iqr_y,
        w1   = wasserstein1d(y[ti], as.numeric(rs$samp)) / iqr_y,
        ks   = ks_uniform(pit)
      )
    }, error = function(e) NULL)
    if (!is.null(out)) per_fold[[length(per_fold) + 1]] <- out
  }
  if (!length(per_fold)) return(c(crps = NA_real_, w1 = NA_real_, ks = NA_real_))
  c(
    crps = mean(vapply(per_fold, `[[`, numeric(1), "crps"), na.rm = TRUE),
    w1   = mean(vapply(per_fold, `[[`, numeric(1), "w1"),   na.rm = TRUE),
    ks   = mean(vapply(per_fold, `[[`, numeric(1), "ks"),   na.rm = TRUE)
  )
}

main <- function() {
  log_step("Stage 07 — cross-validated distributional metrics")
  if (!file.exists(fits_path)) {
    env04 <- new.env(); sys.source(file.path(CODE_DIR, "04_marginals_gamlss.R"), envir = env04)
    invisible(env04$main())
  }
  fits <- readRDS(fits_path)
  if (!requireNamespace("caret", quietly = TRUE)) {
    install.packages("caret", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  rows <- list()
  for (sid in names(fits)) {
    s <- fits[[sid]]
    for (mn in names(s$metrics)) {
      o <- s$metrics[[mn]]
      cv <- cv_one_metric(o$y, o$x, o$family)
      rows[[length(rows) + 1]] <- data.frame(
        stratum = sid, sex = s$sex, sport = s$sport, test_type = s$test_type,
        metric = mn, family = o$family,
        cv_CRPS_rel = signif(cv["crps"], 4),
        cv_W1_rel   = signif(cv["w1"], 4),
        cv_KS_D     = signif(cv["ks"], 4),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows)) {
    tbl <- do.call(rbind, rows)
    utils::write.csv(tbl, file.path(RESULTS_DIR, "07_cv_marginals.csv"),
                     row.names = FALSE)
    log_step("Wrote 07_cv_marginals.csv (", nrow(tbl), " rows)")
  }
}

if (sys.nframe() == 0) main()
