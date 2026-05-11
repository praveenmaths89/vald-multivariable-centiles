# Run only Stages 8 (per-athlete scoring) and 9 (per-athlete radar plots).
# Reuses verification/results/cohort_fits.rds and reliability_table.csv from
# the previous VERIFY_ALL.R run so we don't repeat the 4-minute model fitting.
set.seed(20260508); options(stringsAsFactors = FALSE)
needed <- c("dplyr", "tidyr", "fmsb", "gamlss", "gamlss.dist", "VineCopula")
for (p in needed) suppressPackageStartupMessages(library(p, character.only = TRUE))

VERIFY_DIR <- "verification"; DATA_FILE <- "data/analysis_dataset.csv"
MIN_N      <- 25L
stopifnot(file.exists(file.path(VERIFY_DIR, "results", "cohort_fits.rds")))
stopifnot(file.exists(file.path(VERIFY_DIR, "results", "reliability_table.csv")))
dir.create(file.path(VERIFY_DIR, "figures", "athletes"),
           recursive = TRUE, showWarnings = FALSE)

# ---- helpers needed by predict_cdf (copied verbatim from VERIFY_ALL.R) -----
gamlss_dfun <- function(fam, prefix) {
  tryCatch(get(paste0(prefix, fam), envir = asNamespace("gamlss.dist"),
               mode = "function"), error = function(e) NULL)
}
gamlss_pa <- function(fit, x_new, bm_new) {
  nd <- if (!is.null(bm_new)) data.frame(x = x_new, bm = bm_new)
        else                 data.frame(x = x_new)
  tryCatch(suppressWarnings(gamlss::predictAll(
    fit, newdata = nd, type = "response",
    data = attr(fit, "training_data") %||% fit$call$data)),
    error = function(e) NULL)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
predict_cdf <- function(obj, y_new, x_new, bm_new = NULL) {
  if (obj$type == "Linear")
    return(stats::pnorm((y_new - predict(obj$fit, newdata =
      if (!is.null(bm_new)) data.frame(x = x_new, bm = bm_new)
      else                  data.frame(x = x_new))) / stats::sd(residuals(obj$fit))))
  if (is.null(obj$fit)) return(rep(NA_real_, length(y_new)))
  pfun <- gamlss_dfun(as.character(obj$fit$family)[1], "p")
  pa   <- if (!is.null(pfun)) gamlss_pa(obj$fit, x_new, bm_new) else NULL
  if (is.null(pa)) return(rep(NA_real_, length(y_new)))
  pn <- intersect(c("mu", "sigma", "nu", "tau"), names(pa))
  vapply(seq_along(y_new), function(i) do.call(pfun,
    c(list(q = y_new[i]), lapply(pa[pn], function(z) z[i]))), numeric(1))
}

# ---- load data + cleansing identical to Stage 1 ---------------------------
data_all <- read.csv(DATA_FILE, check.names = FALSE)
g <- toupper(trimws(as.character(data_all$gender)))
data_all$gender <- ifelse(g %in% c("M", "MALE"), "Male",
                   ifelse(g %in% c("F", "FEMALE"), "Female", NA))
data_all$sport <- trimws(gsub("\\s+(Overall(\\s+Data)?)$", "",
                              data_all$sport, ignore.case = TRUE))
data_all <- data_all[!is.na(data_all$gender) & !is.na(data_all$sport) &
                       !is.na(data_all$age_at_test) &
                       data_all$age_at_test >= 8 & data_all$age_at_test <= 60, ]
data_all$bm <- suppressWarnings(as.numeric(data_all$weight))

cohort_fits <- readRDS(file.path(VERIFY_DIR, "results", "cohort_fits.rds"))
reliability_table <- read.csv(file.path(VERIFY_DIR, "results", "reliability_table.csv"))

# =================== STAGE 8: score every athlete-test ==================
cat("\n--- STAGE 8: score every athlete-test ---\n"); t8 <- Sys.time()
band_marg  <- function(p) ifelse(is.na(p), "-",
              ifelse(p < 25, "Below cohort",
              ifelse(p < 50, "Below median",
              ifelse(p < 75, "Above median",
              ifelse(p < 90, "Above cohort", "Elite tail")))))
band_joint <- function(p) ifelse(is.na(p), "-",
              ifelse(p < 25, "Atypical multivariable",
              ifelse(p < 50, "Typical (below)",
              ifelse(p < 75, "Typical (above)",
              ifelse(p < 90, "Above-cohort multivariable", "Elite multivariable")))))
mahal_one <- function(Um, ux) {
  if (length(ux) < 2 || nrow(Um) < ncol(Um) + 2) return(NA_real_)
  cm <- colMeans(Um); cv <- cov(Um) + diag(1e-8, ncol(Um))
  dref <- as.numeric(stats::mahalanobis(Um, cm, cv))
  da   <- as.numeric(stats::mahalanobis(matrix(ux, 1), cm, cv))
  100 * (1 - mean(dref <= da))
}

long_rows <- list(); wide_rows <- list()
for (sid in names(cohort_fits)) {
  cf  <- cohort_fits[[sid]]; pms <- cf$panel_metrics; if (!length(pms)) next
  pcols  <- vapply(pms, function(m) m$col,   character(1))
  pslugs <- vapply(pms, function(m) m$slug,  character(1))
  punits <- vapply(pms, function(m) m$units, character(1))
  rel_c  <- reliability_table[reliability_table$sex == cf$sex &
            reliability_table$sport == cf$sport &
            reliability_table$test  == cf$test, , drop = FALSE]
  mdc95_lookup <- setNames(rel_c$MDC95, rel_c$metric)
  d <- data_all[data_all$gender == cf$sex & data_all$sport == cf$sport &
                  is.finite(data_all$age_at_test), , drop = FALSE]
  d <- d[rowSums(!is.na(d[, pcols, drop = FALSE])) > 0, , drop = FALSE]
  if (!nrow(d)) next
  pct_mat <- matrix(NA_real_, nrow = nrow(d), ncol = length(pms),
                     dimnames = list(NULL, pslugs))
  for (k in seq_along(pms)) {
    mt <- pms[[k]]
    yv <- suppressWarnings(as.numeric(d[[mt$col]])); xv <- d$age_at_test
    bv <- if (mt$mass_dep) d$bm else NULL
    keep <- is.finite(yv) & is.finite(xv) &
            (if (mt$mass_dep) is.finite(bv) else TRUE)
    if (!any(keep)) next
    cdf <- rep(NA_real_, length(yv))
    cdf[keep] <- tryCatch(predict_cdf(mt$obj, yv[keep], xv[keep],
                                       if (mt$mass_dep) bv[keep] else NULL),
                          error = function(e) rep(NA_real_, sum(keep)))
    pct_mat[, k] <- 100 * pmin(pmax(cdf, 0), 1)
  }
  jslugs <- cf$joint_metrics
  Uref   <- if (length(jslugs) >= 2) cf$joint_aligned$U else NULL
  joint_pct <- rep(NA_real_, nrow(d)); n_used <- rep(0L, nrow(d))
  jmodel    <- if (is.null(cf$joint_model)) "-" else cf$joint_model$type
  if (!is.null(Uref) && nrow(Uref) >= MIN_N) {
    j_idx <- match(jslugs, pslugs)
    Pj    <- pct_mat[, j_idx, drop = FALSE] / 100
    Pj[!is.finite(Pj)] <- NA
    Pj[Pj < 1e-6] <- 1e-6; Pj[Pj > 1 - 1e-6] <- 1 - 1e-6
    n_used <- rowSums(!is.na(Pj))
    for (i in which(n_used >= 2)) {
      avail <- which(!is.na(Pj[i, ]))
      Usub  <- Uref[, jslugs[avail], drop = FALSE]
      Usub  <- Usub[stats::complete.cases(Usub), , drop = FALSE]
      ux    <- as.numeric(Pj[i, avail])
      joint_pct[i] <- tryCatch(mahal_one(Usub, ux), error = function(e) NA_real_)
    }
  }
  for (k in seq_along(pms)) {
    yv <- suppressWarnings(as.numeric(d[[pcols[k]]]))
    long_rows[[length(long_rows) + 1]] <- data.frame(
      profile_id = d$profile_id, gender = cf$sex, sport = cf$sport,
      test_type = cf$test, test_date = as.character(d$test_date),
      age_at_test = round(d$age_at_test, 2), bm = round(d$bm, 1),
      metric = pslugs[k], units = punits[k], value = signif(yv, 4),
      marginal_pct = round(pct_mat[, k], 1), band = band_marg(pct_mat[, k]),
      MDC95 = signif(unname(mdc95_lookup[pslugs[k]]), 3),
      joint_pct = round(joint_pct, 1), joint_model = jmodel,
      n_metrics_used = n_used, stringsAsFactors = FALSE)
  }
  wide <- data.frame(profile_id = d$profile_id, gender = cf$sex,
    sport = cf$sport, test_type = cf$test,
    test_date = as.character(d$test_date),
    age_at_test = round(d$age_at_test, 2), bm = round(d$bm, 1),
    joint_pct = round(joint_pct, 1), joint_band = band_joint(joint_pct),
    joint_model = jmodel, n_metrics_used = n_used,
    stringsAsFactors = FALSE)
  for (k in seq_along(pms)) {
    s <- pslugs[k]
    wide[[paste0("pct_",   s)]] <- round(pct_mat[, k], 1)
    wide[[paste0("band_",  s)]] <- band_marg(pct_mat[, k])
    wide[[paste0("value_", s)]] <- signif(suppressWarnings(as.numeric(d[[pcols[k]]])), 4)
  }
  wide_rows[[length(wide_rows) + 1]] <- wide
  cat(sprintf("  scored cohort %-32s rows=%d\n", sid, nrow(d)))
}
athlete_long <- do.call(rbind, long_rows)
athlete_wide <- dplyr::bind_rows(wide_rows)
write.csv(athlete_long,
          file.path(VERIFY_DIR, "results/athlete_scores_long.csv"),
          row.names = FALSE)
write.csv(athlete_wide,
          file.path(VERIFY_DIR, "results/athlete_scores_wide.csv"),
          row.names = FALSE)
score_summary <- athlete_wide |>
  dplyr::group_by(gender, sport, test_type, joint_model) |>
  dplyr::summarise(n_test_rows      = dplyr::n(),
                   n_athletes       = dplyr::n_distinct(profile_id),
                   n_joint_scored   = sum(!is.na(joint_pct)),
                   median_joint_pct = round(median(joint_pct, na.rm = TRUE), 1),
                   pct_elite_joint  = round(100 * mean(joint_pct >= 90, na.rm = TRUE), 1),
                   .groups = "drop")
write.csv(score_summary,
          file.path(VERIFY_DIR, "results/athlete_scores_summary.csv"),
          row.names = FALSE)
cat(sprintf("  long=%d  wide=%d  athletes=%d  with joint pct=%d  (%.1fs)\n",
            nrow(athlete_long), nrow(athlete_wide),
            length(unique(athlete_wide$profile_id)),
            sum(!is.na(athlete_wide$joint_pct)),
            as.numeric(difftime(Sys.time(), t8, units = "secs"))))

# =================== STAGE 9: per-athlete radar plots ===================
cat("\n--- STAGE 9: per-athlete radar plots ---\n"); t9 <- Sys.time()
ath_root <- file.path(VERIFY_DIR, "figures", "athletes")
band_colour <- c("Atypical multivariable"     = "#cb181d",
                 "Typical (below)"            = "#fdae61",
                 "Typical (above)"            = "#fee08b",
                 "Above-cohort multivariable" = "#a6d96a",
                 "Elite multivariable"        = "#1a9850",
                 "-"                          = "#969696")
n_radar <- 0
for (sid in names(cohort_fits)) {
  cf  <- cohort_fits[[sid]]; pms <- cf$panel_metrics
  if (length(pms) < 3) { cat(sprintf("  skip %s (need >=3 axes)\n", sid)); next }
  pslugs <- vapply(pms, function(m) m$slug, character(1))
  this <- athlete_wide[athlete_wide$gender == cf$sex &
                         athlete_wide$sport  == cf$sport &
                         athlete_wide$test_type == cf$test, , drop = FALSE]
  if (!nrow(this)) next
  outdir <- file.path(ath_root, gsub("[^A-Za-z0-9]+", "_", cf$sex),
                       gsub("[^A-Za-z0-9]+", "_", cf$sport),
                       gsub("[^A-Za-z0-9]+", "_", cf$test))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pct_cols <- paste0("pct_", pslugs); vlabels <- gsub("_", " ", pslugs)
  for (i in seq_len(nrow(this))) {
    pcts <- as.numeric(this[i, pct_cols])
    pcts[!is.finite(pcts)] <- 50
    rdf <- as.data.frame(rbind(rep(100, length(pcts)),
                                rep(0,   length(pcts)), pcts))
    colnames(rdf) <- vlabels
    bnd <- this$joint_band[i]; clr <- unname(band_colour[bnd])
    if (is.na(clr)) clr <- "#737373"
    pid_safe  <- gsub("[^A-Za-z0-9]+", "_", as.character(this$profile_id[i]))
    date_safe <- gsub("[^A-Za-z0-9]+", "_", as.character(this$test_date[i]))
    fp <- file.path(outdir, sprintf("%s__%s.png", pid_safe, date_safe))
    png(fp, width = 1200, height = 1200, res = 200)
    op <- par(mar = c(2, 2, 4, 2))
    tryCatch(fmsb::radarchart(
      rdf, axistype = 1, seg = 4,
      caxislabels = c("0", "25", "50", "75", "100"),
      pcol = clr, pfcol = grDevices::adjustcolor(clr, 0.30),
      plwd = 2.5, plty = 1,
      cglcol = "grey75", cglty = 1, axislabcol = "grey40", vlcex = 0.70,
      title = sprintf("%s | %s | %s | %s\nage=%.1f  joint pct=%s  (%s; %s)",
                      this$gender[i], this$sport[i], this$test_type[i],
                      this$test_date[i], this$age_at_test[i],
                      ifelse(is.na(this$joint_pct[i]), "NA",
                             sprintf("%.0f", this$joint_pct[i])),
                      bnd, this$joint_model[i])),
      error = function(e) NULL)
    par(op); dev.off()
    n_radar <- n_radar + 1
  }
  cat(sprintf("  %-32s  PNGs=%d\n", sid, nrow(this)))
}
cat(sprintf("  total radar PNGs: %d  (%.1fs)\n  saved under %s\n",
            n_radar, as.numeric(difftime(Sys.time(), t9, units = "secs")),
            ath_root))
