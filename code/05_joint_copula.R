# ------------------------------------------------------------
# 05_joint_copula.R
# Fit a regular vine copula on the PIT-transformed residuals
# from stage-04 marginals, *within each stratum*. The pair-copula
# construction yields a fully specified joint predictive
# distribution from which joint centiles can be drawn.
#
# Outputs (per stratum):
#   results/strata/<id>/vine.rds            R-vine object
#   results/strata/<id>/vine_summary.csv    AIC, BIC, log-lik
# Master:
#   results/05_joint_summary.csv
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

fits_path  <- file.path(RESULTS_DIR, "04_marginals_fits.rds")
out_root   <- file.path(RESULTS_DIR, "strata")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

build_pit_matrix <- function(stratum) {
  ms <- stratum$metrics
  if (length(ms) < 2) return(NULL)
  n  <- min(vapply(ms, function(m) length(m$pit), integer(1)))
  if (n < MIN_N_STRATUM) return(NULL)
  mat <- vapply(ms, function(m) m$pit[seq_len(n)], numeric(n))
  colnames(mat) <- names(ms)
  mat <- pmin(pmax(mat, 1e-6), 1 - 1e-6)
  mat[stats::complete.cases(mat), , drop = FALSE]
}

fit_vine <- function(U) {
  if (is.null(U) || nrow(U) < MIN_N_STRATUM || ncol(U) < 2) return(NULL)
  tryCatch({
    suppressWarnings(VineCopula::RVineStructureSelect(
      U,
      familyset      = NA,         # all parametric pair-copulas
      type           = 0,          # R-vine
      selectioncrit  = "AIC",
      indeptest      = TRUE,
      level          = 0.05
    ))
  }, error = function(e) NULL)
}

main <- function() {
  log_step("Stage 05 — joint vine copulas")
  if (!file.exists(fits_path)) {
    env04 <- new.env(); sys.source(file.path(CODE_DIR, "04_marginals_gamlss.R"), envir = env04)
    invisible(env04$main())
  }
  fits <- readRDS(fits_path)
  rows <- list()
  for (sid in names(fits)) {
    s <- fits[[sid]]
    U <- build_pit_matrix(s)
    if (is.null(U)) {
      log_step("  skip ", sid, " (insufficient PIT matrix)")
      next
    }
    vine <- fit_vine(U)
    if (is.null(vine)) {
      log_step("  vine failed for ", sid)
      next
    }
    sd_dir <- file.path(out_root, sid)
    dir.create(sd_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(vine = vine, U = U), file.path(sd_dir, "vine.rds"))
    rows[[length(rows) + 1]] <- data.frame(
      stratum = sid, sex = s$sex, sport = s$sport, test_type = s$test_type,
      n_obs = nrow(U), dim = ncol(U),
      vine_AIC = round(vine$AIC, 1),
      vine_BIC = round(vine$BIC, 1),
      vine_loglik = round(VineCopula::RVineLogLik(U, vine)$loglik |> sum(), 1),
      stringsAsFactors = FALSE
    )
    log_step(sprintf("  %s   AIC=%.1f  BIC=%.1f  n=%d  d=%d",
                     sid, vine$AIC, vine$BIC, nrow(U), ncol(U)))
  }
  if (length(rows)) {
    tbl <- do.call(rbind, rows)
    utils::write.csv(tbl, file.path(RESULTS_DIR, "05_joint_summary.csv"),
                     row.names = FALSE)
    log_step("Wrote 05_joint_summary.csv")
  }
}

if (sys.nframe() == 0) main()
