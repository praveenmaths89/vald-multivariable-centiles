# =============================================================================
# Quick demo: score one Male Football CMJ athlete with the saved cohort fits.
# Run AFTER RUN_PAPER_FINAL.R has produced results/final/cohort_fits.rds.
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

.required <- c("gamlss","gamlss.dist","VineCopula")
.miss <- setdiff(.required, rownames(installed.packages()))
if (length(.miss)) {
  stop(sprintf("Missing required R packages for RUN_DEMO_SCORE.R: %s\n  ",
                paste(.miss, collapse = ", ")),
        "Install with:\n",
        sprintf("    install.packages(c(%s))",
                 paste0("\"", .miss, "\"", collapse = ", ")),
        call. = FALSE)
}
suppressPackageStartupMessages({
  library(gamlss); library(gamlss.dist); library(VineCopula)
})

.need_inputs <- c("results/final/cohort_fits.rds",
                   "results/final/reliability_table.csv",
                   "data/analysis_dataset.csv")
.miss_in <- .need_inputs[!file.exists(.need_inputs)]
if (length(.miss_in)) {
  stop("RUN_DEMO_SCORE.R is missing required inputs:\n  ",
        paste(.miss_in, collapse = "\n  "),
        "\nRun the analysis pipeline first:\n    Rscript code/RUN_PAPER_FINAL.R",
        call. = FALSE)
}

cohort_fits <- readRDS("results/final/cohort_fits.rds")
rel_tab     <- read.csv("results/final/reliability_table.csv",
                          stringsAsFactors = FALSE)
data_all    <- read.csv("data/analysis_dataset.csv",
                          stringsAsFactors = FALSE, check.names = FALSE)
g <- toupper(trimws(as.character(data_all$gender)))
data_all$gender <- ifelse(g %in% c("M","MALE"), "Male",
                   ifelse(g %in% c("F","FEMALE"), "Female", NA))
data_all$sport <- trimws(gsub("\\s+(Overall(\\s+Data)?)$", "",
                              data_all$sport, ignore.case=TRUE))
data_all$bm <- suppressWarnings(as.numeric(data_all$weight))

predict_cdf <- function(obj, y_new, x_new, bm_new=NULL) {
  if (obj$type == "Linear") {
    res <- residuals(obj$fit)
    mu  <- if (!is.null(bm_new))
              predict(obj$fit, newdata=data.frame(x=x_new, bm=bm_new))
           else predict(obj$fit, newdata=data.frame(x=x_new))
    return(stats::pnorm((y_new - mu)/stats::sd(res)))
  }
  fam <- as.character(obj$fit$family)[1]
  pfun <- get(paste0("p", fam), envir=asNamespace("gamlss.dist"), mode="function")
  newdat <- if (!is.null(bm_new)) data.frame(x=x_new, bm=bm_new) else data.frame(x=x_new)
  td <- attr(obj$fit, "training_data")
  pa <- suppressWarnings(predictAll(obj$fit, newdata=newdat, type="response", data=td))
  pn <- intersect(c("mu","sigma","nu","tau"), names(pa))
  do.call(pfun, c(list(q=y_new), pa[pn]))
}

score_athlete <- function(profile_id, sex, sport, test_type, n_boot=200) {
  sid <- gsub("[^A-Za-z0-9]+","_", paste(sex, sport, test_type, sep="_"))
  cf <- cohort_fits[[sid]]
  if (is.null(cf)) stop("No cohort fit: ", sid)
  rec_all <- data_all[data_all$profile_id==profile_id & data_all$gender==sex &
                       data_all$sport==sport, ]
  if (!nrow(rec_all)) stop("Athlete not found")
  panel_cols <- vapply(cf$panel_metrics, function(m) m$col, character(1))
  rec_all$cov <- rowSums(!is.na(rec_all[, panel_cols, drop=FALSE]))
  rec_all <- rec_all[rec_all$cov > 0, ]
  if (!nrow(rec_all)) stop("No values for ", test_type)
  rec <- rec_all[order(-rec_all$cov,
                        as.Date(rec_all$test_date), decreasing=c(FALSE,TRUE)), ][1, ]
  age <- rec$age_at_test; bm <- rec$bm
  rows <- list(); uvec <- numeric(0); names_u <- character(0)
  for (m in cf$panel_metrics) {
    yval <- suppressWarnings(as.numeric(rec[[m$col]]))
    if (!is.finite(yval)) next
    cdf <- predict_cdf(m$obj, yval, age, if (m$mass_dep) bm else NULL)
    rel <- rel_tab[rel_tab$sex==cf$sex & rel_tab$sport==cf$sport &
                    rel_tab$test==cf$test & rel_tab$metric==m$slug, ]
    mdc <- if (nrow(rel)) rel$MDC95[1] else NA_real_
    rows[[length(rows)+1]] <- data.frame(
      metric=m$slug, value=signif(yval,4), units=m$units, model=m$best,
      pct=round(100*cdf,1), MDC95=signif(mdc,3))
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
      lo <- as.numeric(stats::quantile(bs, 0.025, na.rm=TRUE))
      hi <- as.numeric(stats::quantile(bs, 0.975, na.rm=TRUE))
    }
  }
  cat(sprintf("\n=== Athlete %s | %s %s %s | age %.1f, BM %.1f kg ===\n",
              profile_id, sex, sport, test_type, age, bm))
  print(out, row.names=FALSE)
  cat(sprintf("Joint percentile: %.1f%% [95%% boot CI: %.1f, %.1f]   joint model = %s\n\n",
              joint, lo, hi, jm$type %||% "—"))
  invisible(out)
}
`%||%` <- function(a,b) if (is.null(a)) b else a

# Pick the closest-to-50th-pct athlete in Male Football CMJ
ex_id <- "Male_Football_CMJ"; cf <- cohort_fits[[ex_id]]
U <- cf$joint_aligned$U
cm <- colMeans(U); cv <- cov(U)
d2 <- as.numeric(stats::mahalanobis(U, cm, cv + diag(1e-8, ncol(U))))
ord <- order(abs(d2 - median(d2)))
demo_pids <- cf$joint_aligned$profile_id[ord[1:3]]
for (pid in demo_pids) score_athlete(pid, "Male", "Football", "CMJ")
