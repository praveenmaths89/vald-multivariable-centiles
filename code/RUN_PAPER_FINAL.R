# =============================================================================
#                       RUN_PAPER_FINAL.R
#         Multivariable Centile Modelling of VALD ForceDecks Data
#                Final methodology (STAI-X 2026 submission)
#
# Pipeline (8 steps):
#   1.  Load the integrated dataset  (no API call here; see code/01_pull_vald.R
#       if you need a fresh API pull).
#   2.  Cohort definition + sample-size feasibility (Royston n_min).
#   3.  Reliability (ICC, SEM, MDC) per (cohort × panel-metric).
#   4.  Apply fixed clinical panel; per-cohort dimensionality:
#         p_cohort = panel size  ↓ drop ICC < 0.5  ↓ drop until n ≥ 10·C(p,2).
#   5.  Marginal model competition: Linear vs GAMLSS-NO vs GAMLSS-BCT(+bm).
#         Score = 5-fold CV CRPS / IQR. Per-metric winner kept.
#   6.  PIT residuals from winning marginals.
#   7.  Joint model competition: Independence vs Gaussian copula vs Vine.
#         Score = 5-fold CV multivariate Energy Score. Per-cohort winner kept.
#   8.  Athlete scoring: marginal percentiles (with MDC95) + joint percentile
#         (with bootstrap 95 % CI).
#
# Outputs (paper assets):
#   paper/figures/fig_methodology.pdf       <- TikZ block diagram (in main.tex)
#   paper/figures/fig_eda_*.pdf             <- 5 EDA figures
#   paper/figures/fig_centiles_example.pdf  <- worked-example multivariable centiles
#   paper/figures/fig_radar_example.pdf     <- worked-example radar
#   paper/tables/tab_cohorts.tex            <- cohort + sample-size table
#   paper/tables/tab_reliability.tex        <- ICC / SEM / MDC table
#   paper/tables/tab_marginal_competition.tex  <- Table A
#   paper/tables/tab_joint_competition.tex     <- Table B
#   paper/numbers.tex                       <- inline numbers macro file
#
# The eight modular code/0X_*.R scripts and code/RUN_ALL_IN_ONE.R from earlier
# iterations are kept untouched on disk for transparency.
# =============================================================================

# ---- Repo root anchor (portable; no absolute user paths) --------------------
.anchor_root <- function() {
  envroot <- Sys.getenv("VALD_PROJECT_DIR", unset = "")
  if (nzchar(envroot) && dir.exists(envroot)) return(normalizePath(envroot))
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(fa)) return(normalizePath(file.path(dirname(fa[1]), ".."),
                                         mustWork = FALSE))
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile))
    return(normalizePath(file.path(dirname(ofile), ".."), mustWork = FALSE))
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  if (basename(cwd) == "code") return(normalizePath(file.path(cwd, "..")))
  cwd
}
PROJECT_DIR <- .anchor_root()
setwd(PROJECT_DIR)

# ---- Tunables (kept minimal) ------------------------------------------------
MIN_N_STRATUM <- 25L
N_SIM         <- 300L
CV_FOLDS      <- 5L
N_BOOT        <- 200L     # bootstrap reps for CIs
SEED          <- 2026L
ICC_CUTOFF    <- 0.50     # below this, drop metric from JOINT model only
TAU_PRECISION <- 0.10     # Royston relative precision target
set.seed(SEED)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Packages ---------------------------------------------------------------
required <- c(
  "dplyr","tidyr","ggplot2","gridExtra","scales",
  "gamlss","gamlss.dist","VineCopula","scoringRules","reshape2","fmsb",
  "psych"   # ICC
)
miss <- setdiff(required, rownames(installed.packages()))
if (length(miss)) {
  stop(sprintf(
    "Missing required R packages for RUN_PAPER_FINAL.R: %s\n  ",
    paste(miss, collapse = ", ")),
    "Install with:\n",
    sprintf("    install.packages(c(%s))",
             paste0("\"", miss, "\"", collapse = ", ")),
    call. = FALSE)
}
invisible(lapply(required, function(p)
  suppressPackageStartupMessages(library(p, character.only=TRUE))))

# ---- Required input check ---------------------------------------------------
.input_csv <- "data/analysis_dataset.csv"
if (!file.exists(.input_csv)) {
  stop(sprintf(
    "Expected %s but it does not exist.\n  ", .input_csv),
    "See data/README.md for the schema. To produce it from a fresh VALD ",
    "pull, run:\n    Rscript code/01_pull_vald.R\n",
    "    Rscript code/02_integrate_sport_gender.R",
    call. = FALSE)
}

# ---- Per-artefact write logger ---------------------------------------------
wrote <- function(path) {
  sz <- tryCatch(file.info(path)$size, error = function(e) NA_real_)
  cat(sprintf("  -> wrote %s (%s)\n", path,
              if (is.na(sz)) "?" else
                if (sz > 1024L * 1024L) sprintf("%.1f MB", sz / 1024 / 1024)
                else sprintf("%.1f KB", sz / 1024)))
}

theme_set(theme_minimal(base_size = 10))
PAL <- list(female="#d7301f", male="#0570b0",
            band="#3690c0", median="#08519c",
            best="#238b45", flag="#fd8d3c")

for (d in c("results/final","paper/figures","paper/tables","logs"))
  if (!dir.exists(d)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

quiet <- function(expr) {
  zz <- tempfile()
  con <- file(zz, "w")
  sink(con, type="output")
  on.exit({
    if (sink.number() > 0) try(sink(NULL, type="output"), silent=TRUE)
    try(close(con), silent=TRUE)
    unlink(zz)
  }, add=TRUE)
  withCallingHandlers(force(expr),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
}
LOG_CON <- stderr()
log_step <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n",
      sep="", file=LOG_CON)
  flush(LOG_CON)
}

# =============================================================================
# STEP 1 — Load data
# =============================================================================
log_step("STEP 1 — load analysis_dataset.csv")
data_all <- read.csv(.input_csv,
                     stringsAsFactors=FALSE, check.names=FALSE)

# Standardise gender, sport, age range
g <- toupper(trimws(as.character(data_all$gender)))
data_all$gender <- ifelse(g %in% c("M","MALE"), "Male",
                   ifelse(g %in% c("F","FEMALE"), "Female", NA))
data_all$sport <- trimws(gsub("\\s+(Overall(\\s+Data)?)$", "",
                              data_all$sport, ignore.case=TRUE))
data_all <- data_all[!is.na(data_all$gender) & !is.na(data_all$sport) &
                       !is.na(data_all$age_at_test) &
                       data_all$age_at_test >= 8 &
                       data_all$age_at_test <= 60, ]
data_all$bm <- suppressWarnings(as.numeric(data_all$weight))
log_step(sprintf("  rows=%d  cols=%d  athletes=%d",
                 nrow(data_all), ncol(data_all),
                 length(unique(data_all$profile_id))))

# =============================================================================
# STEP 2 — Clinical panels per protocol (committed once)
# =============================================================================
# Each panel entry has:
#   slug:    short readable name used in tables/figures
#   regex:   pattern to find the actual VALD column name
#   units:   for the table caption
#   mass_dep: TRUE if metric units involve mass (force/impulse/power/stiffness)
#   ref:     short bibliographic justification for the metric

panels <- list(
  CMJ = list(
    list(slug="JumpHeight",   regex="^Jump\\.Height\\.\\.Flight\\.Time\\.\\_Both_CMJ$",
         units="cm",  mass_dep=FALSE, ref="Markovic 2007"),
    list(slug="PeakPower_BM", regex="^Concentric\\.Peak\\.Power\\.\\.\\.BM_Both_CMJ$",
         units="W/kg",mass_dep=FALSE, ref="Cormack 2008"),
    list(slug="mRSI",         regex="^RSI\\.modified_Both_CMJ$",
         units="m/s",mass_dep=FALSE, ref="Suchomel 2015"),
    list(slug="EccMeanForce", regex="^Eccentric\\.Mean\\.Force_Both_CMJ$",
         units="N",  mass_dep=TRUE,  ref="McMahon 2018")
  ),
  DJ = list(
    list(slug="RSI_FTCT",     regex="^RSI\\.\\.Flight\\.Time\\.Contact\\.Time\\.\\_Both_DJ$",
         units="m/s",mass_dep=FALSE, ref="Flanagan 2008"),
    list(slug="ContactTime",  regex="^Contact\\.Time_Both_DJ$",
         units="s",  mass_dep=FALSE, ref="Komi 2000"),
    list(slug="JumpHeight",   regex="^Jump\\.Height\\.\\.Flight\\.Time\\.\\_Both_DJ$",
         units="cm", mass_dep=FALSE, ref="Markovic 2007"),
    list(slug="PeakLandF_BW", regex="^Peak\\.Landing\\.Force\\.\\.\\.BW_Both_DJ$",
         units="N/kg",mass_dep=FALSE,ref="Bishop 2018")
  ),
  SJ = list(
    list(slug="JumpHeight",   regex="^Jump\\.Height\\.\\.Flight\\.Time\\.\\_Both_SJ$",
         units="cm",  mass_dep=FALSE, ref="Markovic 2007"),
    list(slug="PeakPower_BM", regex="^Concentric\\.Peak\\.Power\\.\\.\\.BM_Both_SJ$",
         units="W/kg",mass_dep=FALSE, ref="Cormack 2008"),
    list(slug="ConcMF_BM",    regex="^Concentric\\.Mean\\.Force\\.\\.\\.BM_Both_SJ$",
         units="N/kg",mass_dep=FALSE, ref="Samozino 2014"),
    list(slug="ForceAtPP",    regex="^Force\\.at\\.Peak\\.Power_Both_SJ$",
         units="N",  mass_dep=TRUE,  ref="Samozino 2014")
  ),
  SLJ = list(
    list(slug="ConcImp_R",    regex="^Concentric\\.Impulse_Right_SLJ$",
         units="N·s", mass_dep=TRUE,  ref="Bishop 2018"),
    list(slug="ConcImp_L",    regex="^Concentric\\.Impulse_Left_SLJ$",
         units="N·s", mass_dep=TRUE,  ref="Bishop 2018"),
    list(slug="ConcMF_R",     regex="^Concentric\\.Mean\\.Force_Right_SLJ$",
         units="N",  mass_dep=TRUE,  ref="—"),
    list(slug="ConcMF_L",     regex="^Concentric\\.Mean\\.Force_Left_SLJ$",
         units="N",  mass_dep=TRUE,  ref="—")
  ),
  SLDJ = list(
    list(slug="ConcImp_R",    regex="^Concentric\\.Impulse_Right_SLDJ$",
         units="N·s", mass_dep=TRUE,  ref="Bishop 2018"),
    list(slug="ConcImp_L",    regex="^Concentric\\.Impulse_Left_SLDJ$",
         units="N·s", mass_dep=TRUE,  ref="Bishop 2018"),
    list(slug="MeanLandPow_R",regex="^Mean\\.Landing\\.Power_Right_SLDJ$",
         units="W",  mass_dep=TRUE,  ref="—"),
    list(slug="MeanLandPow_L",regex="^Mean\\.Landing\\.Power_Left_SLDJ$",
         units="W",  mass_dep=TRUE,  ref="—")
  ),
  SLLAH = list(
    list(slug="DropLand_R",   regex="^Drop\\.Landing_Right_SLLAH$",
         units="N",  mass_dep=TRUE,  ref="Bishop 2018"),
    list(slug="DropLand_L",   regex="^Drop\\.Landing_Left_SLLAH$",
         units="N",  mass_dep=TRUE,  ref="Bishop 2018"),
    list(slug="PeakDLF_R",    regex="^Peak\\.Drop\\.Landing\\.Force_Right_SLLAH$",
         units="N",  mass_dep=TRUE,  ref="—")
  ),
  SLHAR = list(
    list(slug="PeakTakeoffF_R", regex="^Peak\\.Takeoff\\.Force_Right_SLHAR$",
         units="N", mass_dep=TRUE, ref="Bishop 2018"),
    list(slug="PeakTakeoffF_L", regex="^Peak\\.Takeoff\\.Force_Left_SLHAR$",
         units="N", mass_dep=TRUE, ref="Bishop 2018"),
    list(slug="PeakLandF_R",    regex="^Peak\\.First\\.Landing\\.Force_Right_SLHAR$",
         units="N", mass_dep=TRUE, ref="—"),
    list(slug="PeakLandF_L",    regex="^Peak\\.First\\.Landing\\.Force_Left_SLHAR$",
         units="N", mass_dep=TRUE, ref="—")
  )
)

resolve_panel <- function(panel, all_cols) {
  resolved <- list()
  for (p in panel) {
    hit <- grep(p$regex, all_cols, value=TRUE)
    if (length(hit)) {
      p$col <- hit[1]
      resolved[[length(resolved)+1]] <- p
    }
  }
  resolved
}

# =============================================================================
# STEP 3 — Cohort table + Royston sample-size diagnostic
# =============================================================================
log_step("STEP 2/3 — cohort table + Royston n_min")

royston_nmin <- function(y, p=0.97, tau=TAU_PRECISION) {
  y <- y[is.finite(y)]
  if (length(y) < 20) return(NA_integer_)
  iqr <- stats::IQR(y); if (!is.finite(iqr) || iqr <= 0) return(NA_integer_)
  qp  <- stats::quantile(y, probs=p, na.rm=TRUE, names=FALSE)
  d   <- stats::density(y, from=qp - iqr, to=qp + iqr, n=64)
  fhat <- approx(d$x, d$y, xout=qp, rule=2)$y
  if (!is.finite(fhat) || fhat <= 0) return(NA_integer_)
  ceiling(p*(1-p) / (tau^2 * iqr^2 * fhat^2))
}

cohort_rows <- list()
for (tt in names(panels)) {
  cols <- intersect(names(data_all), sapply(panels[[tt]], function(p) p$regex))
  resolved <- resolve_panel(panels[[tt]], names(data_all))
  if (!length(resolved)) next
  panel_cols <- sapply(resolved, function(p) p$col)
  d_tt <- data_all[, c("profile_id","gender","sport","age_at_test","bm",
                        "test_date", panel_cols), drop=FALSE]
  d_tt <- d_tt[rowSums(!is.na(d_tt[, panel_cols, drop=FALSE])) > 0, ]
  for (sx in c("Female","Male")) {
    dsx <- d_tt[d_tt$gender == sx, ]
    if (!nrow(dsx)) next
    by_sport <- split(dsx, dsx$sport)
    for (sp in names(by_sport)) {
      d <- by_sport[[sp]]
      n <- nrow(d)
      if (n < MIN_N_STRATUM) next
      n_min <- vapply(panel_cols, function(c) {
        v <- royston_nmin(d[[c]])
        if (is.na(v)) NA_integer_ else as.integer(v)
      }, integer(1))
      p_stat <- max(2L, floor((1 + sqrt(1 + 8*n/10)) / 2))   # n >= 10*C(p,2)
      cohort_rows[[length(cohort_rows)+1]] <- list(
        sex=sx, sport=sp, test=tt, n=n,
        panel_cols=panel_cols, n_min=n_min, p_stat=p_stat,
        d=d
      )
    }
  }
}
log_step(sprintf("  eligible cohorts: %d", length(cohort_rows)))

# =============================================================================
# STEP 4 — Reliability (ICC, SEM, MDC) per cohort × panel-metric
# =============================================================================
log_step("STEP 4 — reliability (ICC, SEM, MDC)")

icc_for_metric <- function(d, col_name, max_gap_days = 180) {
  ok <- d[is.finite(d[[col_name]]), ]
  if (!nrow(ok)) return(c(ICC=NA_real_, SEM=NA_real_, MDC95=NA_real_, n_pairs=0))
  ok$test_date <- as.Date(ok$test_date)
  ok <- ok[order(ok$profile_id, ok$test_date), ]
  pairs_df <- ok |>
    dplyr::group_by(profile_id) |>
    dplyr::filter(dplyr::n() >= 2) |>
    dplyr::slice_head(n=2) |>
    dplyr::ungroup()
  if (nrow(pairs_df) < 6) return(c(ICC=NA_real_, SEM=NA_real_, MDC95=NA_real_, n_pairs=0))
  spread <- tidyr::pivot_wider(pairs_df,
              id_cols=profile_id, names_from=test_date,
              values_from=tidyselect::all_of(col_name),
              values_fn = mean)
  m <- as.matrix(spread[, -1, drop=FALSE])
  good <- which(rowSums(!is.na(m)) >= 2)
  if (length(good) < 6) return(c(ICC=NA_real_, SEM=NA_real_, MDC95=NA_real_, n_pairs=0))
  M <- t(apply(m[good, , drop=FALSE], 1, function(r) r[which(!is.na(r))[1:2]]))
  res <- tryCatch(suppressMessages(suppressWarnings(quiet(psych::ICC(M, missing=FALSE)))),
                   error=function(e) NULL)
  if (is.null(res)) return(c(ICC=NA_real_, SEM=NA_real_, MDC95=NA_real_, n_pairs=0))
  icc31 <- res$results[res$results$type == "ICC3", "ICC"][1]
  sd_within <- sqrt(mean(apply(M, 1, var, na.rm=TRUE), na.rm=TRUE))
  sd_total  <- stats::sd(as.numeric(M), na.rm=TRUE)
  sem  <- sd_total * sqrt(max(0, 1 - icc31))
  mdc  <- 1.96 * sqrt(2) * sem
  c(ICC=icc31, SEM=sem, MDC95=mdc, n_pairs=length(good))
}

reliability_rows <- list()
for (cr in cohort_rows) {
  for (k in seq_along(cr$panel_cols)) {
    col <- cr$panel_cols[k]
    pe  <- panels[[cr$test]][[k]]
    r <- icc_for_metric(cr$d, col)
    reliability_rows[[length(reliability_rows)+1]] <- data.frame(
      sex=cr$sex, sport=cr$sport, test=cr$test, metric=pe$slug,
      n_pairs=as.integer(r["n_pairs"]),
      ICC=signif(r["ICC"], 3), SEM=signif(r["SEM"], 3),
      MDC95=signif(r["MDC95"], 3),
      stringsAsFactors=FALSE
    )
  }
}
reliability_table <- do.call(rbind, reliability_rows)
write.csv(reliability_table, "results/final/reliability_table.csv", row.names=FALSE)
wrote("results/final/reliability_table.csv")

# Per-protocol reliability summary (median ICC across cohorts that had pairs)
rel_summary <- reliability_table |>
  dplyr::filter(!is.na(ICC)) |>
  dplyr::group_by(test, metric) |>
  dplyr::summarise(median_ICC = round(median(ICC, na.rm=TRUE), 2),
                   median_MDC95 = signif(median(MDC95, na.rm=TRUE), 3),
                   n_cohorts = dplyr::n(), .groups="drop")
write.csv(rel_summary, "results/final/reliability_summary.csv", row.names=FALSE)
wrote("results/final/reliability_summary.csv")
log_step(sprintf("  reliability rows: %d", nrow(reliability_table)))

# =============================================================================
# STEP 5 — Marginal model competition per (cohort × metric)
#   M1: Linear              y = β0 + β1·age
#   M2: GAMLSS-NO           y ~ pb(age)
#   M3: GAMLSS-BCT(+bm)     y ~ pb(age) [+ pb(bm) if mass-dep]
# Score: 5-fold CV CRPS / IQR. Lower is better.
# =============================================================================

gamlss_fit_safe <- function(y, x, bm=NULL, family_chr, n_cyc=80) {
  dat <- if (is.null(bm)) data.frame(y=y, x=x) else data.frame(y=y, x=x, bm=bm)
  formula_mu <- if (is.null(bm)) y ~ gamlss::pb(x) else y ~ gamlss::pb(x) + gamlss::pb(bm)
  call <- bquote(gamlss::gamlss(
                   formula      = .(formula_mu),
                   sigma.formula= ~ gamlss::pb(x),
                   family       = .(family_chr),
                   data         = dat,
                   trace        = FALSE,
                   control      = gamlss::gamlss.control(n.cyc=.(n_cyc))))
  fit <- tryCatch(quiet(eval(call)), error=function(e) NULL)
  if (!is.null(fit)) attr(fit, "training_data") <- dat
  fit
}
predict_samples <- function(fit, x_new, bm_new=NULL, n_sim=N_SIM) {
  if (is.null(fit)) return(NULL)
  fam <- as.character(fit$family)[1]
  rfun <- tryCatch(get(paste0("r", fam),
                       envir=asNamespace("gamlss.dist"), mode="function"),
                   error=function(e) NULL)
  if (is.null(rfun)) return(NULL)
  newdat <- if (is.null(bm_new)) data.frame(x=x_new) else data.frame(x=x_new, bm=bm_new)
  td <- attr(fit, "training_data")
  pa <- tryCatch(suppressWarnings(gamlss::predictAll(
              fit, newdata=newdat, type="response", data=td)),
              error=function(e) NULL)
  if (is.null(pa)) return(NULL)
  pn <- intersect(c("mu","sigma","nu","tau"), names(pa))
  out <- tryCatch(
    replicate(n_sim,
       suppressWarnings(do.call(rfun, c(list(n=length(x_new)), pa[pn])))),
    error=function(e) NULL)
  if (is.null(out)) return(NULL)
  matrix(out, nrow=length(x_new))
}
winsorise <- function(samp, y) {
  yq <- stats::quantile(y, c(0.005, 0.995), na.rm=TRUE, names=FALSE)
  span <- max(yq[2] - yq[1], stats::sd(y, na.rm=TRUE), 1e-6)
  cap_lo <- yq[1] - 8*span; cap_hi <- yq[2] + 8*span
  samp[!is.finite(samp)] <- NA
  samp[samp < cap_lo] <- cap_lo; samp[samp > cap_hi] <- cap_hi
  samp
}

# Linear model "predictive samples" via residual bootstrap
linear_samples <- function(y_tr, x_tr, x_te, n_sim=N_SIM, bm_tr=NULL, bm_te=NULL) {
  if (!is.null(bm_tr) && !is.null(bm_te)) {
    dat <- data.frame(y=y_tr, x=x_tr, bm=bm_tr)
    fit <- lm(y ~ x + bm, data=dat)
    mu  <- predict(fit, newdata=data.frame(x=x_te, bm=bm_te))
  } else {
    dat <- data.frame(y=y_tr, x=x_tr)
    fit <- lm(y ~ x, data=dat)
    mu  <- predict(fit, newdata=data.frame(x=x_te))
  }
  res <- residuals(fit)
  matrix(rep(mu, n_sim) + sample(res, length(mu)*n_sim, replace=TRUE),
         nrow=length(x_te))
}

cv_score <- function(y, x, bm, model, k=CV_FOLDS) {
  set.seed(SEED)
  n <- length(y)
  if (n < 25) return(NA_real_)
  iqr <- max(stats::IQR(y, na.rm=TRUE), stats::sd(y, na.rm=TRUE), 1e-9)
  folds <- split(sample.int(n), rep(seq_len(k), length.out=n))
  scores <- numeric(0)
  for (fi in seq_along(folds)) {
    ti <- folds[[fi]]; tr <- setdiff(seq_len(n), ti)
    samp <- switch(model,
      "Linear" = linear_samples(y[tr], x[tr], x[ti], N_SIM,
                  if (!is.null(bm)) bm[tr] else NULL,
                  if (!is.null(bm)) bm[ti] else NULL),
      "NO"     = {
        fit <- gamlss_fit_safe(y[tr], x[tr], NULL, "NO")
        predict_samples(fit, x[ti])
      },
      "BCT"    = {
        bm_tr <- if (!is.null(bm)) bm[tr] else NULL
        bm_te <- if (!is.null(bm)) bm[ti] else NULL
        # BCT requires y > 0; if not, fall through to NO with bm covariate
        if (any(y[tr] <= 0, na.rm=TRUE)) {
          fit <- gamlss_fit_safe(y[tr], x[tr], bm_tr, "NO")
        } else {
          fit <- gamlss_fit_safe(y[tr], x[tr], bm_tr, "BCT")
          if (is.null(fit)) fit <- gamlss_fit_safe(y[tr], x[tr], bm_tr, "BCCGo")
          if (is.null(fit)) fit <- gamlss_fit_safe(y[tr], x[tr], bm_tr, "NO")
        }
        s <- predict_samples(fit, x[ti], bm_te)
        s
      })
    if (is.null(samp)) next
    samp <- winsorise(samp, y[tr])
    s <- tryCatch(mean(scoringRules::crps_sample(y=y[ti], dat=samp), na.rm=TRUE)/iqr,
                  error=function(e) NA_real_)
    if (is.finite(s)) scores <- c(scores, s)
  }
  if (!length(scores)) return(NA_real_)
  mean(scores)
}

log_step("STEP 5 — marginal competition (CV CRPS/IQR)")
marginal_rows <- list()
for (cr in cohort_rows) {
  d <- cr$d
  for (k in seq_along(cr$panel_cols)) {
    col <- cr$panel_cols[k]
    pe  <- panels[[cr$test]][[k]]
    y <- suppressWarnings(as.numeric(d[[col]]))
    x <- d$age_at_test
    bm <- if (pe$mass_dep) d$bm else NULL
    keep <- is.finite(y) & is.finite(x)
    if (!is.null(bm)) keep <- keep & is.finite(bm)
    if (sum(keep) < MIN_N_STRATUM) next
    ys <- y[keep]; xs <- x[keep]; bms <- if (!is.null(bm)) bm[keep] else NULL
    sLin <- cv_score(ys, xs, bms, "Linear")
    sNO  <- cv_score(ys, xs, bms, "NO")
    sBCT <- cv_score(ys, xs, bms, "BCT")
    scores <- c(Linear=sLin, NO=sNO, BCT=sBCT)
    best <- names(scores)[which.min(scores)]
    marginal_rows[[length(marginal_rows)+1]] <- data.frame(
      sex=cr$sex, sport=cr$sport, test=cr$test, metric=pe$slug,
      n=sum(keep), mass_dep=pe$mass_dep,
      Linear=signif(sLin, 3), NO=signif(sNO, 3),
      BCT=signif(sBCT, 3), Best=best,
      stringsAsFactors=FALSE)
  }
}
marginal_competition <- do.call(rbind, marginal_rows)
write.csv(marginal_competition, "results/final/marginal_competition.csv", row.names=FALSE)
wrote("results/final/marginal_competition.csv")
log_step(sprintf("  marginal rows: %d", nrow(marginal_competition)))
log_step(sprintf("  best counts: %s",
  paste(names(table(marginal_competition$Best)), table(marginal_competition$Best),
        sep="=", collapse=" ")))

# =============================================================================
# STEP 6 — Per-cohort dimensionality + fit winning marginals
# =============================================================================
fit_best_marginal <- function(y, x, bm, best) {
  switch(best,
    "Linear" = list(type="Linear",
                    fit = if (!is.null(bm)) lm(y ~ x + bm) else lm(y ~ x),
                    bm = bm, x=x, y=y),
    "NO"     = { f <- gamlss_fit_safe(y, x, NULL, "NO");
                  list(type="NO", fit=f, bm=NULL, x=x, y=y) },
    "BCT"    = { f <- gamlss_fit_safe(y, x, bm, "BCT");
                  if (is.null(f)) f <- gamlss_fit_safe(y, x, bm, "BCCGo")
                  list(type="BCT", fit=f, bm=bm, x=x, y=y) }
  )
}
predict_cdf <- function(obj, y_new, x_new, bm_new=NULL) {
  if (obj$type == "Linear") {
    res <- residuals(obj$fit)
    mu  <- if (!is.null(bm_new))
              predict(obj$fit, newdata=data.frame(x=x_new, bm=bm_new))
           else predict(obj$fit, newdata=data.frame(x=x_new))
    s   <- stats::sd(res)
    return(stats::pnorm((y_new - mu)/s))
  }
  if (is.null(obj$fit)) return(rep(NA_real_, length(y_new)))
  fam <- as.character(obj$fit$family)[1]
  pfun <- tryCatch(get(paste0("p", fam),
                       envir=asNamespace("gamlss.dist"), mode="function"),
                   error=function(e) NULL)
  if (is.null(pfun)) return(rep(NA_real_, length(y_new)))
  newdat <- if (!is.null(bm_new)) data.frame(x=x_new, bm=bm_new) else data.frame(x=x_new)
  td <- attr(obj$fit, "training_data")
  pa <- tryCatch(suppressWarnings(gamlss::predictAll(
            obj$fit, newdata=newdat, type="response", data=td)),
            error=function(e) NULL)
  if (is.null(pa)) return(rep(NA_real_, length(y_new)))
  pn <- intersect(c("mu","sigma","nu","tau"), names(pa))
  vapply(seq_along(y_new), function(i)
    do.call(pfun, c(list(q=y_new[i]), lapply(pa[pn], function(z) z[i]))),
    numeric(1))
}

# Build cohort fits
log_step("STEP 6 — fit winning marginals + per-cohort PIT matrices")
cohort_fits <- list()
for (cr in cohort_rows) {
  sid <- gsub("[^A-Za-z0-9]+","_", paste(cr$sex, cr$sport, cr$test, sep="_"))
  panel_metrics <- list()
  for (k in seq_along(cr$panel_cols)) {
    col <- cr$panel_cols[k]; pe <- panels[[cr$test]][[k]]
    y <- suppressWarnings(as.numeric(cr$d[[col]])); x <- cr$d$age_at_test
    bm <- if (pe$mass_dep) cr$d$bm else NULL
    keep <- is.finite(y) & is.finite(x)
    if (!is.null(bm)) keep <- keep & is.finite(bm)
    if (sum(keep) < MIN_N_STRATUM) next
    rowx <- which(marginal_competition$sex==cr$sex &
                  marginal_competition$sport==cr$sport &
                  marginal_competition$test==cr$test &
                  marginal_competition$metric==pe$slug)
    if (!length(rowx)) next
    best <- marginal_competition$Best[rowx[1]]
    obj  <- fit_best_marginal(y[keep], x[keep],
              if (!is.null(bm)) bm[keep] else NULL, best)
    if (is.null(obj$fit)) next
    panel_metrics[[pe$slug]] <- list(
      slug=pe$slug, col=col, units=pe$units, mass_dep=pe$mass_dep,
      best=best, obj=obj, x=x[keep], y=y[keep],
      bm=if (!is.null(bm)) bm[keep] else NULL,
      profile_id=cr$d$profile_id[keep])
  }
  if (length(panel_metrics) < 2) next

  # Apply ICC drop (drop metric only from JOINT model)
  rel_for_cohort <- reliability_table[reliability_table$sex==cr$sex &
                                       reliability_table$sport==cr$sport &
                                       reliability_table$test==cr$test, ]
  drop_for_joint <- rel_for_cohort$metric[!is.na(rel_for_cohort$ICC) &
                                          rel_for_cohort$ICC < ICC_CUTOFF]
  joint_metrics  <- setdiff(names(panel_metrics), drop_for_joint)

  # Build a JOINT-aligned sample: rows of cr$d that have ALL joint_metrics finite
  joint_cols <- vapply(joint_metrics, function(m) panel_metrics[[m]]$col, character(1))
  any_mass <- any(vapply(joint_metrics, function(m) panel_metrics[[m]]$mass_dep, logical(1)))
  align_keep <- rep(TRUE, nrow(cr$d))
  for (jc in joint_cols) align_keep <- align_keep & is.finite(suppressWarnings(as.numeric(cr$d[[jc]])))
  align_keep <- align_keep & is.finite(cr$d$age_at_test)
  if (any_mass) align_keep <- align_keep & is.finite(cr$d$bm)
  n_joint <- sum(align_keep)

  # Sample-size cap p such that n_joint >= 10 * C(p,2)
  p_stat_cap <- max(2L, floor((1 + sqrt(1 + 8*n_joint/10)) / 2))

  # Drop metrics until we satisfy the cap; drop the least informative
  # (smallest age-residual variance reduction) first.
  while (length(joint_metrics) > p_stat_cap && length(joint_metrics) > 2) {
    res_var <- vapply(joint_metrics, function(m) {
      pm <- panel_metrics[[m]]
      r <- residuals(lm(pm$y ~ pm$x))
      stats::var(r) / max(stats::var(pm$y), 1e-9)
    }, numeric(1))
    drop <- joint_metrics[which.max(res_var)]   # largest residual-to-variance = least informative
    joint_metrics <- setdiff(joint_metrics, drop)
    joint_cols <- joint_cols[joint_metrics]
    align_keep <- rep(TRUE, nrow(cr$d))
    for (jc in joint_cols) align_keep <- align_keep & is.finite(suppressWarnings(as.numeric(cr$d[[jc]])))
    align_keep <- align_keep & is.finite(cr$d$age_at_test)
    any_mass <- any(vapply(joint_metrics, function(m) panel_metrics[[m]]$mass_dep, logical(1)))
    if (any_mass) align_keep <- align_keep & is.finite(cr$d$bm)
    n_joint <- sum(align_keep)
    p_stat_cap <- max(2L, floor((1 + sqrt(1 + 8*n_joint/10)) / 2))
  }

  # Compute aligned PIT vectors for each joint metric
  if (n_joint < MIN_N_STRATUM || length(joint_metrics) < 2) {
    binding_cap <- "drop"
    joint_aligned <- list(profile_id=character(0), age=numeric(0),
                           bm=numeric(0), Y=matrix(numeric(0), 0, 0),
                           U=matrix(numeric(0), 0, 0))
  } else {
    sub <- cr$d[align_keep, , drop=FALSE]
    Y <- sapply(joint_metrics, function(m) {
      pm <- panel_metrics[[m]]
      suppressWarnings(as.numeric(sub[[pm$col]]))
    })
    U <- sapply(joint_metrics, function(m) {
      pm <- panel_metrics[[m]]
      cdf <- predict_cdf(pm$obj, Y[, m], sub$age_at_test,
                         if (pm$mass_dep) sub$bm else NULL)
      pmin(pmax(cdf, 1e-6), 1 - 1e-6)
    })
    if (is.null(dim(Y))) {Y <- matrix(Y, ncol=length(joint_metrics)); colnames(Y) <- joint_metrics}
    if (is.null(dim(U))) {U <- matrix(U, ncol=length(joint_metrics)); colnames(U) <- joint_metrics}
    binding_cap <- if (length(joint_metrics) == length(panel_metrics)) "panel"
                   else if (length(drop_for_joint)) "reliability" else "sample size"
    joint_aligned <- list(profile_id=sub$profile_id, age=sub$age_at_test,
                           bm=sub$bm, Y=Y, U=U)
  }

  cohort_fits[[sid]] <- list(
    sex=cr$sex, sport=cr$sport, test=cr$test, n=cr$n,
    panel_metrics=panel_metrics,
    joint_metrics=joint_metrics,
    n_joint=n_joint,
    p_stat=p_stat_cap,
    binding_cap=binding_cap,
    joint_aligned=joint_aligned
  )
}
log_step(sprintf("  fitted cohorts: %d", length(cohort_fits)))
saveRDS(cohort_fits, "results/final/cohort_fits.rds")
wrote("results/final/cohort_fits.rds")

# =============================================================================
# STEP 7 — Joint model competition per cohort
#   J1: Independence  -> ES of independent uniform marginals
#   J2: Gaussian copula
#   J3: Vine copula (RVineStructureSelect)
# Score: 5-fold CV multivariate Energy Score on the predictive sample.
# =============================================================================
log_step("STEP 7 — joint competition (CV Energy Score)")

# Build a per-cohort PIT matrix aligned by profile_id
build_PIT <- function(cf) {
  if (is.null(cf$joint_aligned) || nrow(cf$joint_aligned$U) < MIN_N_STRATUM)
    return(NULL)
  list(U = cf$joint_aligned$U, profile_id = cf$joint_aligned$profile_id)
}

# Multivariate energy score: ES(F, y) = E||X-y|| - 0.5 E||X-X'||
energy_score <- function(samples, y) {
  X <- samples; n_samp <- nrow(X)
  if (any(is.na(X))) return(NA_real_)
  diffs1 <- sqrt(rowSums((X - matrix(y, nrow=n_samp, ncol=length(y), byrow=TRUE))^2))
  idx2 <- sample.int(n_samp); X2 <- X[idx2, , drop=FALSE]
  diffs2 <- sqrt(rowSums((X - X2)^2))
  mean(diffs1) - 0.5*mean(diffs2)
}

# Per-cohort joint CV
joint_competition_rows <- list()
for (sid in names(cohort_fits)) {
  cf <- cohort_fits[[sid]]
  pit <- build_PIT(cf)
  if (is.null(pit) || ncol(pit$U) < 2) next
  U <- pit$U; n <- nrow(U); p <- ncol(U)

  set.seed(SEED)
  folds <- split(sample.int(n), rep(seq_len(CV_FOLDS), length.out=n))

  es_indep <- numeric(0); es_gauss <- numeric(0); es_vine <- numeric(0)
  for (fi in seq_along(folds)) {
    ti <- folds[[fi]]; tr <- setdiff(seq_len(n), ti)
    Utr <- U[tr, , drop=FALSE]; Ute <- U[ti, , drop=FALSE]

    # J1: Independence — sample uniform per dim
    samp_indep <- matrix(runif(N_SIM * p), nrow=N_SIM, ncol=p)
    es1 <- mean(vapply(seq_len(nrow(Ute)),
                  function(i) energy_score(samp_indep, Ute[i,]), numeric(1)),
                na.rm=TRUE)

    # J2: Gaussian copula via empirical correlation on Φ⁻¹(U)
    Z <- qnorm(Utr)
    R <- tryCatch(stats::cor(Z), error=function(e) diag(p))
    samp_g <- tryCatch({
      L <- chol(R + diag(1e-8, p))
      W <- matrix(rnorm(N_SIM * p), nrow=N_SIM) %*% L
      pnorm(W)
    }, error=function(e) NULL)
    es2 <- if (is.null(samp_g)) NA else
      mean(vapply(seq_len(nrow(Ute)),
                  function(i) energy_score(samp_g, Ute[i,]), numeric(1)),
           na.rm=TRUE)

    # J3: Vine copula
    vine <- tryCatch(suppressWarnings(VineCopula::RVineStructureSelect(
              Utr, familyset=c(0,1,2,3,4,5,6,13,14), # indep, Gauss, t, Clayton, Gumbel, Frank, Joe, surv-Clayton, surv-Gumbel
              selectioncrit="AIC", indeptest=TRUE, level=0.05)),
              error=function(e) NULL)
    samp_v <- if (is.null(vine)) NULL else
      tryCatch(VineCopula::RVineSim(N_SIM, vine), error=function(e) NULL)
    es3 <- if (is.null(samp_v)) NA else
      mean(vapply(seq_len(nrow(Ute)),
                  function(i) energy_score(samp_v, Ute[i,]), numeric(1)),
           na.rm=TRUE)

    es_indep <- c(es_indep, es1); es_gauss <- c(es_gauss, es2); es_vine <- c(es_vine, es3)
  }
  s_ind <- mean(es_indep, na.rm=TRUE)
  s_gau <- mean(es_gauss, na.rm=TRUE)
  s_vin <- mean(es_vine, na.rm=TRUE)
  scores <- c(Indep=s_ind, Gaussian=s_gau, Vine=s_vin)
  if (!any(is.finite(scores))) next
  best <- names(scores)[which.min(replace(scores, !is.finite(scores), Inf))]
  gain <- if (is.finite(s_ind)) (s_ind - scores[best]) / s_ind * 100 else NA_real_

  # Fit final winning joint model on all rows for downstream scoring
  if (best == "Vine") {
    final_vine <- tryCatch(suppressWarnings(VineCopula::RVineStructureSelect(
                  U, familyset=c(0,1,2,3,4,5,6,13,14),
                  selectioncrit="AIC", indeptest=TRUE, level=0.05)),
                  error=function(e) NULL)
    cohort_fits[[sid]]$joint_model <- list(type="Vine", obj=final_vine, U=U)
  } else if (best == "Gaussian") {
    Z <- qnorm(U); R <- cor(Z)
    cohort_fits[[sid]]$joint_model <- list(type="Gaussian", obj=R, U=U)
  } else {
    cohort_fits[[sid]]$joint_model <- list(type="Indep", obj=NULL, U=U)
  }

  joint_competition_rows[[length(joint_competition_rows)+1]] <- data.frame(
    sex=cf$sex, sport=cf$sport, test=cf$test, n=cf$n, p=ncol(U),
    binding=cf$binding_cap,
    Indep=signif(s_ind, 3), Gaussian=signif(s_gau, 3), Vine=signif(s_vin, 3),
    Best=best, GainOverIndep_pct=round(gain, 1),
    stringsAsFactors=FALSE)
}
joint_competition <- do.call(rbind, joint_competition_rows)
write.csv(joint_competition, "results/final/joint_competition.csv", row.names=FALSE)
wrote("results/final/joint_competition.csv")
saveRDS(cohort_fits, "results/final/cohort_fits.rds")
wrote("results/final/cohort_fits.rds")
log_step(sprintf("  joint rows: %d", nrow(joint_competition)))
log_step(sprintf("  best counts: %s",
  paste(names(table(joint_competition$Best)), table(joint_competition$Best),
        sep="=", collapse=" ")))

# =============================================================================
# STEP 8 — Athlete scoring + bootstrap CI on joint percentile
# =============================================================================
score_athlete <- function(profile_id, sex, sport, test_type,
                          fits = cohort_fits, data = data_all,
                          n_boot = N_BOOT) {
  sid <- gsub("[^A-Za-z0-9]+","_", paste(sex, sport, test_type, sep="_"))
  cf <- fits[[sid]]
  if (is.null(cf)) stop("No fitted cohort for ", sid)
  rec_all <- data[data$profile_id==profile_id & data$gender==sex & data$sport==sport, ]
  if (!nrow(rec_all)) stop("Athlete not found.")
  panel_cols <- vapply(cf$panel_metrics, function(m) m$col, character(1))
  rec_all$cov <- rowSums(!is.na(rec_all[, panel_cols, drop=FALSE]))
  rec_all <- rec_all[rec_all$cov > 0, ]
  if (!nrow(rec_all)) stop("Athlete has no panel values for ", test_type)
  rec <- rec_all[order(-rec_all$cov, as.Date(rec_all$test_date), decreasing=c(FALSE,TRUE)), ][1, ]
  age <- rec$age_at_test; bm  <- rec$bm
  rows <- list(); uvec <- numeric(0); names_u <- character(0)
  for (m in cf$panel_metrics) {
    yval <- suppressWarnings(as.numeric(rec[[m$col]]))
    if (!is.finite(yval)) next
    cdf <- predict_cdf(m$obj, yval, age, if (m$mass_dep) bm else NULL)
    rel <- reliability_table[reliability_table$sex==cf$sex &
                              reliability_table$sport==cf$sport &
                              reliability_table$test==cf$test &
                              reliability_table$metric==m$slug, ]
    mdc <- if (nrow(rel)) rel$MDC95[1] else NA_real_
    rows[[length(rows)+1]] <- data.frame(
      metric=m$slug, value=signif(yval, 4),
      units=m$units, model=m$best,
      marginal_pctile=round(100*cdf, 1),
      MDC95=signif(mdc, 3),
      stringsAsFactors=FALSE)
    if (m$slug %in% cf$joint_metrics) {
      uvec <- c(uvec, cdf); names_u <- c(names_u, m$slug)
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()

  joint <- NA_real_; joint_lo <- NA_real_; joint_hi <- NA_real_
  jm <- cf$joint_model
  if (!is.null(jm) && length(uvec) >= 2) {
    common <- intersect(names_u, colnames(jm$U))
    if (length(common) >= 2) {
      U <- jm$U[, common, drop=FALSE]; ux <- uvec[match(common, names_u)]
      depth_pct <- function(U, x) {
        cm <- colMeans(U); cv <- stats::cov(U) + diag(1e-8, ncol(U))
        d2_a <- as.numeric(stats::mahalanobis(matrix(x, 1), cm, cv))
        d2_c <- as.numeric(stats::mahalanobis(U, cm, cv))
        100 * (1 - mean(d2_c <= d2_a))
      }
      joint <- depth_pct(U, ux)
      bs <- replicate(n_boot, {
        idx <- sample.int(nrow(U), nrow(U), replace=TRUE)
        depth_pct(U[idx, , drop=FALSE], ux)
      })
      joint_lo <- as.numeric(stats::quantile(bs, 0.025, na.rm=TRUE))
      joint_hi <- as.numeric(stats::quantile(bs, 0.975, na.rm=TRUE))
    }
  }
  cat(sprintf("Athlete %s | %s %s %s | age %.1f bm %.1f kg | joint %s%s\n",
    profile_id, sex, sport, test_type, age, bm,
    if (is.na(joint)) "NA" else sprintf("%.1f%%", joint),
    if (is.na(joint)) "" else sprintf(" [%.1f, %.1f]", joint_lo, joint_hi)))
  out$joint_pctile <- if (is.na(joint)) NA else round(joint, 1)
  out$joint_lo <- if (is.na(joint)) NA else round(joint_lo, 1)
  out$joint_hi <- if (is.na(joint)) NA else round(joint_hi, 1)
  attr(out, "joint") <- list(point=joint, lo=joint_lo, hi=joint_hi,
                              joint_model=jm$type)
  out
}
saveRDS(cohort_fits, "results/final/cohort_fits.rds")
wrote("results/final/cohort_fits.rds")
log_step("STEP 8 — score_athlete() ready")

# =============================================================================
# Save a short cohort summary table for the paper
# =============================================================================
cohort_table <- do.call(rbind, lapply(names(cohort_fits), function(sid) {
  cf <- cohort_fits[[sid]]
  data.frame(
    sex=cf$sex, sport=cf$sport, test=cf$test, n=cf$n,
    panel_size=length(cf$panel_metrics),
    p_cohort=length(cf$joint_metrics),
    binding=cf$binding_cap,
    joint_model=if (is.null(cf$joint_model)) "—" else cf$joint_model$type,
    stringsAsFactors=FALSE)
}))
write.csv(cohort_table, "results/final/cohort_table.csv", row.names=FALSE)
wrote("results/final/cohort_table.csv")

cat("\n==== summary ====\n")
cat(sprintf("Eligible cohorts:        %d\n", length(cohort_fits)))
cat(sprintf("Marginals fit:           %d\n", nrow(marginal_competition)))
cat(sprintf("Marginal best counts:    %s\n",
   paste(names(table(marginal_competition$Best)),
         table(marginal_competition$Best), sep="=", collapse=" ")))
cat(sprintf("Joint best counts:       %s\n",
   paste(names(table(joint_competition$Best)),
         table(joint_competition$Best), sep="=", collapse=" ")))
cat("Outputs in results/final/.\n")
