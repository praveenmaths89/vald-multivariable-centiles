# ------------------------------------------------------------
# helpers.R
# Pure functions reused by multiple stages: I/O, GAMLSS wrappers,
# distributional metrics, and small plotting utilities.
# ------------------------------------------------------------

# Quiet a verbose call's stdout (e.g. GAMLSS RS iterations) -----
quiet_call <- function(expr) {
  zz <- tempfile()
  sink(zz, type = "output")
  on.exit({ sink(type = "output"); unlink(zz) })
  force(expr)
}

# Slugify free text into safe filenames ------------------------
slugify <- function(x) gsub("[^A-Za-z0-9]+", "_", trimws(x))

# Normalise sport name suffixes coming from the external CSV --
clean_sport <- function(x) {
  x <- gsub("\\s+(Overall(\\s+Data)?)$", "", x, ignore.case = TRUE)
  trimws(x)
}

# Standardise gender to {Male, Female, NA} ---------------------
clean_gender <- function(x) {
  x <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("M", "MALE")   ~ "Male",
    x %in% c("F", "FEMALE") ~ "Female",
    TRUE                    ~ NA_character_
  )
}

# GAMLSS fit using do.call so $call$family is a literal string
# (predict.gamlss(newdata=) breaks otherwise). -----------------
gamlss_age_fit <- function(y_vec, x_vec, family_chr,
                           n.cyc = 120L, quiet = TRUE) {
  dat <- data.frame(y = y_vec, x = x_vec)
  fit_call <- list(
    formula       = y ~ gamlss::pb(x),
    sigma.formula = ~ gamlss::pb(x),
    family        = family_chr,
    data          = dat,
    trace         = FALSE,
    control       = gamlss::gamlss.control(n.cyc = n.cyc)
  )
  if (quiet) {
    quiet_call(do.call(gamlss::gamlss, fit_call))
  } else {
    do.call(gamlss::gamlss, fit_call)
  }
}

# Sample n_sim predictive draws per row of x_new from the family
sample_predictive <- function(fit, x_new, n_sim = 300) {
  fam <- as.character(fit$family)[1]
  rfun <- tryCatch(
    get(paste0("r", fam), envir = asNamespace("gamlss.dist"), mode = "function"),
    error = function(e) NULL
  )
  if (is.null(rfun)) return(NULL)
  pa <- tryCatch(
    suppressWarnings(gamlss::predictAll(
      fit, newdata = data.frame(x = x_new), type = "response",
      data = fit$call$data
    )),
    error = function(e) NULL
  )
  if (is.null(pa)) return(NULL)
  pn <- intersect(c("mu", "sigma", "nu", "tau"), names(pa))
  if (!length(pn)) return(NULL)
  n  <- length(x_new)
  out <- matrix(NA_real_, n, n_sim)
  for (s in seq_len(n_sim)) {
    args <- c(list(n = n), lapply(pa[pn], identity))
    out[, s] <- tryCatch(do.call(rfun, args), error = function(e) rep(NA_real_, n))
  }
  out
}

# Wasserstein-1 distance via numeric quantile integration -----
wasserstein1d <- function(x, y, n_grid = 1024L) {
  x <- sort(as.numeric(x[is.finite(x)]))
  y <- sort(as.numeric(y[is.finite(y)]))
  if (!length(x) || !length(y)) return(NA_real_)
  u  <- seq(0.5 / n_grid, 1 - 0.5 / n_grid, length.out = n_grid)
  qx <- stats::quantile(x, probs = u, names = FALSE, type = 8)
  qy <- stats::quantile(y, probs = u, names = FALSE, type = 8)
  mean(abs(qx - qy))
}

# Cramér–von Mises and Anderson–Darling against Uniform(0,1) --
cvm_uniform <- function(p) {
  p <- sort(as.numeric(p[is.finite(p)]))
  n <- length(p)
  if (n < 5) return(NA_real_)
  i <- seq_len(n)
  sum((p - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n)
}
ad_uniform <- function(p) {
  p <- sort(as.numeric(p[is.finite(p)]))
  n <- length(p)
  if (n < 5) return(NA_real_)
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  i <- seq_len(n)
  -n - sum((2 * i - 1) / n * (log(p) + log(1 - rev(p))))
}

# Kolmogorov-Smirnov D against Uniform(0,1) -------------------
ks_uniform <- function(p) {
  p <- p[is.finite(p)]
  if (length(p) < 5) return(NA_real_)
  unname(suppressWarnings(stats::ks.test(p, "punif")$statistic))
}

# Robust predictive sampling: detect blow-up and refit Gaussian
robust_sample_predictive <- function(fit, x, y, n_sim, fam_used) {
  yq <- stats::quantile(y, c(0.005, 0.5, 0.995), na.rm = TRUE, names = FALSE)
  span <- max(yq[3] - yq[1], stats::sd(y, na.rm = TRUE), 1e-6)
  cap_lo <- yq[1] - 8 * span
  cap_hi <- yq[3] + 8 * span

  samp <- tryCatch(sample_predictive(fit, x, n_sim), error = function(e) NULL)
  blew <- !is.null(samp) && (
    !all(is.finite(samp)) ||
    max(abs(samp), na.rm = TRUE) > 1e3 * (max(abs(y), na.rm = TRUE) + 1)
  )
  if (blew && fam_used != "NO") {
    fit_no <- tryCatch(
      suppressWarnings(gamlss_age_fit(y, x, "NO")),
      error = function(e) NULL
    )
    if (!is.null(fit_no)) {
      fit <- fit_no
      fam_used <- "NO"
      samp <- tryCatch(sample_predictive(fit, x, n_sim), error = function(e) NULL)
    }
  }
  if (!is.null(samp)) {
    samp[!is.finite(samp)] <- NA_real_
    samp[samp < cap_lo] <- cap_lo
    samp[samp > cap_hi] <- cap_hi
  }
  list(fit = fit, fam = fam_used, samp = samp)
}

# Single-stratum compact summary line --------------------------
fmt_stratum <- function(sex, sport, test_type, n) {
  sprintf("%-6s | %-18s | %-7s | n=%d", sex, sport, test_type, n)
}

# Write a CSV that always exists, even when empty --------------
safe_write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}
