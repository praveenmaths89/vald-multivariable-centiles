# VERIFY_ALL.R -- regenerate every body paper artefact (Chougale & Ananthakumar
# 2026, STAI-X 2026) from data/analysis_dataset.csv into verification/.  Run
# from repo root: Rscript code/VERIFY_ALL.R   (~1 min on a laptop).  Stages
# match paper Sec 3.1-3.6 (cleansing / reliability / marginal & joint
# competition / scoring); Stage 7 builds Sec 4-5 tables, figures and macros.
# Appendix tables + the cross-cohort PIT / centile-bands grids are omitted
# (the canonical pipeline produces them) to keep this script self-contained.
set.seed(20260508); options(stringsAsFactors = FALSE)
needed <- c("dplyr","tidyr","ggplot2","gridExtra","scales","fmsb",
            "gamlss","gamlss.dist","VineCopula","scoringRules","psych")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing R packages.  Install with:\n  ",
  sprintf("install.packages(c(%s))", paste0("\"", miss, "\"", collapse=", ")), call.=FALSE)
for (p in needed) suppressPackageStartupMessages(library(p, character.only=TRUE))
if (getRversion() < "4.2.0") stop("R >= 4.2.0 required.", call. = FALSE)
VERIFY_DIR <- "verification"; DATA_FILE <- "data/analysis_dataset.csv"
for (sub in c("results","tables","figures","numbers"))
  dir.create(file.path(VERIFY_DIR, sub), recursive=TRUE, showWarnings=FALSE)
if (!file.exists(DATA_FILE)) stop("Required input not found: ", DATA_FILE, call.=FALSE)
MIN_N <- 25L; ICC_CUTOFF <- 0.50    # min cohort size; ICC drop threshold (paper 3.1, 3.3)
N_SIM <- 300L; CV_FOLDS <- 5L       # predictive draws; CV folds            (paper 3.4-3.5)
N_BOOT <- 200L; SEED <- 2026L       # bootstrap reps; CV-fold seed          (paper 3.6)
TIMINGS <- list(); record <- function(stage, t0) TIMINGS[[stage]] <<- as.numeric(difftime(Sys.time(), t0, units="secs"))
# Clinical metric panels.  Each row = "slug|regex|units|mass_dep" (mass_dep="T"
# means body mass enters the GAMLSS location formula).  Parsed below into the
# nested list `panels` keyed by test code (CMJ, DJ, SJ, SLJ, SLDJ, SLLAH, SLHAR).
panel_lines <- list(
  CMJ = c("JumpHeight|^Jump\\.Height\\.\\.Flight\\.Time\\.\\_Both_CMJ$|cm|F", "PeakPower_BM|^Concentric\\.Peak\\.Power\\.\\.\\.BM_Both_CMJ$|W/kg|F",
          "mRSI|^RSI\\.modified_Both_CMJ$|m/s|F", "EccMeanForce|^Eccentric\\.Mean\\.Force_Both_CMJ$|N|T"),
  DJ = c("RSI_FTCT|^RSI\\.\\.Flight\\.Time\\.Contact\\.Time\\.\\_Both_DJ$|m/s|F", "ContactTime|^Contact\\.Time_Both_DJ$|s|F",
         "JumpHeight|^Jump\\.Height\\.\\.Flight\\.Time\\.\\_Both_DJ$|cm|F", "PeakLandF_BW|^Peak\\.Landing\\.Force\\.\\.\\.BW_Both_DJ$|N/kg|F"),
  SJ = c("JumpHeight|^Jump\\.Height\\.\\.Flight\\.Time\\.\\_Both_SJ$|cm|F", "PeakPower_BM|^Concentric\\.Peak\\.Power\\.\\.\\.BM_Both_SJ$|W/kg|F",
         "ConcMF_BM|^Concentric\\.Mean\\.Force\\.\\.\\.BM_Both_SJ$|N/kg|F", "ForceAtPP|^Force\\.at\\.Peak\\.Power_Both_SJ$|N|T"),
  SLJ = c("ConcImp_R|^Concentric\\.Impulse_Right_SLJ$|N.s|T", "ConcImp_L|^Concentric\\.Impulse_Left_SLJ$|N.s|T",
          "ConcMF_R|^Concentric\\.Mean\\.Force_Right_SLJ$|N|T", "ConcMF_L|^Concentric\\.Mean\\.Force_Left_SLJ$|N|T"),
  SLDJ = c("ConcImp_R|^Concentric\\.Impulse_Right_SLDJ$|N.s|T", "ConcImp_L|^Concentric\\.Impulse_Left_SLDJ$|N.s|T",
           "MeanLandPow_R|^Mean\\.Landing\\.Power_Right_SLDJ$|W|T", "MeanLandPow_L|^Mean\\.Landing\\.Power_Left_SLDJ$|W|T"),
  SLLAH = c("DropLand_R|^Drop\\.Landing_Right_SLLAH$|N|T", "DropLand_L|^Drop\\.Landing_Left_SLLAH$|N|T",
            "PeakDLF_R|^Peak\\.Drop\\.Landing\\.Force_Right_SLLAH$|N|T"),
  SLHAR = c("PeakTakeoffF_R|^Peak\\.Takeoff\\.Force_Right_SLHAR$|N|T", "PeakTakeoffF_L|^Peak\\.Takeoff\\.Force_Left_SLHAR$|N|T",
            "PeakLandF_R|^Peak\\.First\\.Landing\\.Force_Right_SLHAR$|N|T", "PeakLandF_L|^Peak\\.First\\.Landing\\.Force_Left_SLHAR$|N|T"))
panels <- lapply(panel_lines, function(rows) lapply(strsplit(rows, "\\|"), function(r)
  list(slug=r[1], regex=r[2], units=r[3], mass_dep=r[4]=="T")))
# Helpers (each used 3+ times).  GAMLSS per Rigby & Stasinopoulos (2005); BCT
# per Cole & Green (1992); ICC(3,1) per Shrout & Fleiss (1979); Koo & Li (2016).
gamlss_fit_safe <- function(y, x, bm = NULL, family_chr, n_cyc = 80) {
  dat <- if (is.null(bm)) data.frame(y=y, x=x) else data.frame(y=y, x=x, bm=bm)
  fmu <- if (is.null(bm)) y ~ gamlss::pb(x) else y ~ gamlss::pb(x) + gamlss::pb(bm)
  capture.output(fit <- tryCatch(suppressMessages(suppressWarnings(gamlss::gamlss(
      formula=fmu, sigma.formula=~gamlss::pb(x), family=family_chr, data=dat,
      trace=FALSE, control=gamlss::gamlss.control(n.cyc=n_cyc)))),
    error=function(e) NULL))
  if (!is.null(fit)) attr(fit, "training_data") <- dat
  fit
}
gamlss_pa <- function(fit, x_new, bm_new)                          # predictAll(...)
  tryCatch(suppressWarnings(gamlss::predictAll(fit, newdata=
      if (is.null(bm_new)) data.frame(x=x_new) else data.frame(x=x_new, bm=bm_new),
      type="response", data=attr(fit, "training_data"))), error=function(e) NULL)
gamlss_dfun <- function(fam, prefix) tryCatch(get(paste0(prefix, fam),
  envir=asNamespace("gamlss.dist"), mode="function"), error=function(e) NULL)
predict_samples <- function(fit, x_new, bm_new = NULL, n_sim = N_SIM) {
  if (is.null(fit)) return(NULL)
  rfun <- gamlss_dfun(as.character(fit$family)[1], "r")
  pa <- if (!is.null(rfun)) gamlss_pa(fit, x_new, bm_new) else NULL
  if (is.null(pa)) return(NULL)
  pn <- intersect(c("mu","sigma","nu","tau"), names(pa))
  out <- tryCatch(replicate(n_sim, suppressWarnings(do.call(rfun,
            c(list(n=length(x_new)), pa[pn])))), error=function(e) NULL)
  if (is.null(out)) NULL else matrix(out, nrow=length(x_new))
}
df_xb <- function(x, bm) if (!is.null(bm)) data.frame(x=x, bm=bm) else data.frame(x=x)
predict_cdf <- function(obj, y_new, x_new, bm_new = NULL) {
  if (obj$type == "Linear") return(stats::pnorm((y_new - predict(obj$fit,
    newdata=df_xb(x_new, bm_new))) / stats::sd(residuals(obj$fit))))
  if (is.null(obj$fit)) return(rep(NA_real_, length(y_new)))
  pfun <- gamlss_dfun(as.character(obj$fit$family)[1], "p")
  pa <- if (!is.null(pfun)) gamlss_pa(obj$fit, x_new, bm_new) else NULL
  if (is.null(pa)) return(rep(NA_real_, length(y_new)))
  pn <- intersect(c("mu","sigma","nu","tau"), names(pa))
  vapply(seq_along(y_new), function(i) do.call(pfun,
    c(list(q=y_new[i]), lapply(pa[pn], function(z) z[i]))), numeric(1))
}
icc_for_metric <- function(d, col_name) {                  # ICC(3,1) + SEM + MDC95
  ok <- d[is.finite(d[[col_name]]), ]
  na4 <- c(ICC=NA_real_, SEM=NA_real_, MDC95=NA_real_, n_pairs=0)
  if (!nrow(ok)) return(na4)
  ok$test_date <- as.Date(ok$test_date); ok <- ok[order(ok$profile_id, ok$test_date), ]
  pp <- ok |> dplyr::group_by(profile_id) |> dplyr::filter(dplyr::n() >= 2) |>
              dplyr::slice_head(n=2) |> dplyr::ungroup()
  if (nrow(pp) < 6) return(na4)
  M <- as.matrix(tidyr::pivot_wider(pp, id_cols=profile_id, names_from=test_date,
            values_from=tidyselect::all_of(col_name), values_fn=mean)[, -1, drop=FALSE])
  good <- which(rowSums(!is.na(M)) >= 2); if (length(good) < 6) return(na4)
  M <- t(apply(M[good, , drop=FALSE], 1, function(r) r[which(!is.na(r))[1:2]]))
  res <- tryCatch(suppressMessages(suppressWarnings(psych::ICC(M, missing=FALSE))),
                  error=function(e) NULL); if (is.null(res)) return(na4)
  icc <- res$results[res$results$type == "ICC3", "ICC"][1]
  sem <- stats::sd(as.numeric(M), na.rm=TRUE) * sqrt(max(0, 1 - icc))
  c(ICC=icc, SEM=sem, MDC95=1.96*sqrt(2)*sem, n_pairs=length(good))
}
save_jpg <- function(file, plt, w, h) {                # used 5+ times in Stage 7
  jpeg(file, width=w, height=h, units="in", res=300); print(plt); invisible(dev.off()) }
tex_safe <- function(s) gsub("([%&#_])", "\\\\\\1", as.character(s))
tex_table <- function(file, colspec, header, rows, footer = character(0))
  writeLines(c(sprintf("\\begin{tabular}{%s}", colspec), "\\toprule", header,
    "\\midrule", rows, if (length(footer)) c("\\midrule", footer) else NULL,
    "\\bottomrule", "\\end{tabular}"), file)
# STAGE 1 (paper Sec 3.1) -- Cleansing: gender/sport label normalisation,
# age restriction (8-60 yr), body-mass parse.
cat("\n--- STAGE 1: cleansing & loading ---\n"); t1 <- Sys.time()
data_all <- read.csv(DATA_FILE, check.names=FALSE)
g <- toupper(trimws(as.character(data_all$gender)))
data_all$gender <- ifelse(g %in% c("M","MALE"), "Male",
                   ifelse(g %in% c("F","FEMALE"), "Female", NA))
data_all$sport <- trimws(gsub("\\s+(Overall(\\s+Data)?)$", "", data_all$sport, ignore.case=TRUE))
data_all <- data_all[!is.na(data_all$gender) & !is.na(data_all$sport) &
                       !is.na(data_all$age_at_test) &
                       data_all$age_at_test >= 8 & data_all$age_at_test <= 60, ]
data_all$bm <- suppressWarnings(as.numeric(data_all$weight))
cat(sprintf("  rows=%d  athletes=%d  sports=%d\n", nrow(data_all),
            length(unique(data_all$profile_id)), length(unique(data_all$sport))))
record("Stage 1  cleansing", t1)
# STAGE 2 (paper Sec 3.2) -- One cohort per (sex x sport x test) with n>=MIN_N;
# per (cohort x metric) compute ICC(3,1), SEM, MDC95 from same-athlete repeats.
cat("\n--- STAGE 2: cohort table + reliability ---\n"); t2 <- Sys.time()
cohort_rows <- list()
for (tt in names(panels)) {
  resolved <- list()
  for (p in panels[[tt]]) { hit <- grep(p$regex, names(data_all), value=TRUE)
    if (length(hit)) { p$col <- hit[1]; resolved[[length(resolved)+1]] <- p } }
  if (!length(resolved)) next
  pcs <- vapply(resolved, function(p) p$col, character(1))
  d_tt <- data_all[rowSums(!is.na(data_all[, pcs, drop=FALSE])) > 0,
            c("profile_id","gender","sport","age_at_test","bm","test_date", pcs), drop=FALSE]
  for (sx in c("Female","Male")) for (sp in unique(d_tt$sport[d_tt$gender == sx])) {
    d <- d_tt[d_tt$gender == sx & d_tt$sport == sp, ]; if (nrow(d) < MIN_N) next
    cohort_rows[[length(cohort_rows)+1]] <- list(sex=sx, sport=sp, test=tt,
      n=nrow(d), resolved=resolved, panel_cols=pcs, d=d) }
}
rel_rows <- list()
for (cr in cohort_rows) for (k in seq_along(cr$panel_cols)) {
  r <- icc_for_metric(cr$d, cr$panel_cols[k])
  rel_rows[[length(rel_rows)+1]] <- data.frame(sex=cr$sex, sport=cr$sport, test=cr$test, metric=cr$resolved[[k]]$slug,
    n_pairs=as.integer(r["n_pairs"]), ICC=signif(r["ICC"],3), SEM=signif(r["SEM"],3), MDC95=signif(r["MDC95"],3))
}
reliability_table <- do.call(rbind, rel_rows)
write.csv(reliability_table, file.path(VERIFY_DIR,"results/reliability_table.csv"), row.names=FALSE)
reliability_summary <- reliability_table |> dplyr::filter(!is.na(ICC)) |>
  dplyr::group_by(test, metric) |>
  dplyr::summarise(median_ICC=round(median(ICC, na.rm=TRUE), 2),
                   median_MDC95=signif(median(MDC95, na.rm=TRUE), 3),
                   n_cohorts=dplyr::n(), .groups="drop")
write.csv(reliability_summary, file.path(VERIFY_DIR,"results/reliability_summary.csv"), row.names=FALSE)
cat(sprintf("  eligible cohorts: %d   reliability rows: %d\n",
            length(cohort_rows), nrow(reliability_table)))
record("Stage 2  reliability", t2)
# STAGE 3 (paper Sec 3.4) -- 5-fold CV CRPS/IQR (Gneiting & Raftery 2007) over
# Linear / GAMLSS-NO / GAMLSS-BCT marginals; lowest score wins for Stage 4.
cat("\n--- STAGE 3: marginal competition (CV CRPS / IQR) ---\n"); t3 <- Sys.time()
marginal_rows <- list()
for (cr in cohort_rows) for (k in seq_along(cr$panel_cols)) {
  pe <- cr$resolved[[k]]; col <- cr$panel_cols[k]
  y <- suppressWarnings(as.numeric(cr$d[[col]])); x <- cr$d$age_at_test
  bm <- if (pe$mass_dep) cr$d$bm else NULL
  keep <- is.finite(y) & is.finite(x); if (!is.null(bm)) keep <- keep & is.finite(bm)
  if (sum(keep) < MIN_N) next
  ys <- y[keep]; xs <- x[keep]; bms <- if (!is.null(bm)) bm[keep] else NULL
  iqr <- max(stats::IQR(ys, na.rm=TRUE), stats::sd(ys, na.rm=TRUE), 1e-9)
  scores <- c(Linear=NA_real_, NO=NA_real_, BCT=NA_real_)
  set.seed(SEED); folds <- split(sample.int(length(ys)),
                                  rep(seq_len(CV_FOLDS), length.out=length(ys)))
  for (mdl in names(scores)) {
    accum <- numeric(0)
    for (fi in seq_along(folds)) {
      ti <- folds[[fi]]; tr <- setdiff(seq_along(ys), ti)
      bm_tr <- if (!is.null(bms)) bms[tr] else NULL; bm_te <- if (!is.null(bms)) bms[ti] else NULL
      samp <- if (mdl == "Linear") {
        f <- lm(y ~ ., data=cbind(y=ys[tr], df_xb(xs[tr], bm_tr)))
        mu <- predict(f, newdata=df_xb(xs[ti], bm_te))
        matrix(rep(mu, N_SIM) + sample(residuals(f), length(mu)*N_SIM, replace=TRUE),
               nrow=length(xs[ti]))
      } else if (mdl == "NO") predict_samples(
                  gamlss_fit_safe(ys[tr], xs[tr], NULL, "NO"), xs[ti])
      else if (any(ys[tr] <= 0)) NULL else {
        f <- gamlss_fit_safe(ys[tr], xs[tr], bm_tr, "BCT")
        if (is.null(f)) f <- gamlss_fit_safe(ys[tr], xs[tr], bm_tr, "NO")
        predict_samples(f, xs[ti], bm_te)
      }
      if (is.null(samp)) next
      yq <- stats::quantile(ys[tr], c(0.005, 0.995), na.rm=TRUE, names=FALSE)
      span <- max(yq[2]-yq[1], stats::sd(ys[tr], na.rm=TRUE), 1e-6)
      samp[!is.finite(samp)] <- NA
      samp[samp < yq[1]-8*span] <- yq[1]-8*span; samp[samp > yq[2]+8*span] <- yq[2]+8*span
      s <- tryCatch(mean(scoringRules::crps_sample(y=ys[ti], dat=samp), na.rm=TRUE)/iqr,
                    error=function(e) NA_real_)
      if (is.finite(s)) accum <- c(accum, s)
    }
    if (length(accum)) scores[mdl] <- mean(accum)
  }
  marginal_rows[[length(marginal_rows)+1]] <- data.frame(sex=cr$sex, sport=cr$sport, test=cr$test,
    metric=pe$slug, n=sum(keep), mass_dep=pe$mass_dep, Linear=signif(scores["Linear"], 3),
    NO=signif(scores["NO"], 3), BCT=signif(scores["BCT"], 3), Best=names(scores)[which.min(scores)])
}
marginal_competition <- do.call(rbind, marginal_rows)
write.csv(marginal_competition, file.path(VERIFY_DIR,"results/marginal_competition.csv"), row.names=FALSE)
cat(sprintf("  best counts: %s\n", paste(names(table(marginal_competition$Best)),
            table(marginal_competition$Best), sep="=", collapse=" ")))
record("Stage 3  marginals", t3)
# STAGE 4 (paper Sec 3.3) -- ICC<0.5 drop -> joint metric set; refit winning
# marginal on full data; PIT matrix U feeds the joint copula in Stage 5.
cat("\n--- STAGE 4: final marginals + PIT ---\n"); t4 <- Sys.time()
cohort_fits <- list()
for (cr in cohort_rows) {
  sid <- gsub("[^A-Za-z0-9]+","_", paste(cr$sex, cr$sport, cr$test, sep="_"))
  pms <- list()
  for (k in seq_along(cr$panel_cols)) {
    pe <- cr$resolved[[k]]; col <- cr$panel_cols[k]
    y <- suppressWarnings(as.numeric(cr$d[[col]])); x <- cr$d$age_at_test
    bm <- if (pe$mass_dep) cr$d$bm else NULL
    keep <- is.finite(y) & is.finite(x); if (!is.null(bm)) keep <- keep & is.finite(bm)
    if (sum(keep) < MIN_N) next
    bm_k <- if (!is.null(bm)) bm[keep] else NULL
    rowx <- which(marginal_competition$sex==cr$sex & marginal_competition$sport==cr$sport &
                   marginal_competition$test==cr$test & marginal_competition$metric==pe$slug)
    if (!length(rowx)) next
    best <- marginal_competition$Best[rowx[1]]
    obj <- if (best == "Linear") list(type="Linear", fit=lm(y~., data=cbind(y=y[keep], df_xb(x[keep], bm_k))))
      else { f <- gamlss_fit_safe(y[keep], x[keep], if (best=="BCT") bm_k else NULL, if (best=="BCT") "BCT" else "NO")
             if (is.null(f) && best=="BCT") f <- gamlss_fit_safe(y[keep], x[keep], bm_k, "NO"); list(type=best, fit=f) }
    if (is.null(obj$fit)) next
    pms[[pe$slug]] <- list(slug=pe$slug, col=col, units=pe$units, mass_dep=pe$mass_dep,
      best=best, obj=obj, x=x[keep], y=y[keep], bm=bm_k, profile_id=cr$d$profile_id[keep])
  }
  if (length(pms) < 2) next
  rel_c <- reliability_table[reliability_table$sex==cr$sex &
              reliability_table$sport==cr$sport & reliability_table$test==cr$test, ]
  jm <- setdiff(names(pms), rel_c$metric[!is.na(rel_c$ICC) & rel_c$ICC < ICC_CUTOFF])
  if (length(jm) < 2) jm <- character(0)
  ja <- list(profile_id=character(0), age=numeric(0), bm=numeric(0), Y=NULL, U=NULL)
  if (length(jm) >= 2) {
    jcols <- vapply(jm, function(m) pms[[m]]$col, character(1))
    any_mass <- any(vapply(jm, function(m) pms[[m]]$mass_dep, logical(1)))
    al <- rowSums(sapply(jcols, function(jc) !is.finite(suppressWarnings(as.numeric(cr$d[[jc]]))))) == 0 &
          is.finite(cr$d$age_at_test) & (!any_mass | is.finite(cr$d$bm))
    if (sum(al) >= MIN_N) {
      sub <- cr$d[al, , drop=FALSE]
      Y <- sapply(jm, function(m) suppressWarnings(as.numeric(sub[[pms[[m]]$col]])))
      U <- sapply(jm, function(m) { pm <- pms[[m]]; pmin(pmax(predict_cdf(pm$obj,
        Y[, m], sub$age_at_test, if (pm$mass_dep) sub$bm else NULL), 1e-6), 1-1e-6) })
      if (is.null(dim(Y))) { Y <- matrix(Y, ncol=length(jm), dimnames=list(NULL, jm))
                             U <- matrix(U, ncol=length(jm), dimnames=list(NULL, jm)) }
      ja <- list(profile_id=sub$profile_id, age=sub$age_at_test, bm=sub$bm, Y=Y, U=U)
    }
  }
  cohort_fits[[sid]] <- list(sex=cr$sex, sport=cr$sport, test=cr$test, n=cr$n,
    panel_metrics=pms, joint_metrics=jm, joint_aligned=ja)
}
cat(sprintf("  fitted cohorts: %d\n", length(cohort_fits)))
record("Stage 4  final marginals", t4)
# STAGE 5 (paper Sec 3.5) -- Independence / Gaussian / R-vine scored by CV
# Energy Score; Vine via Dissmann et al. (2013) over the AIC pair-copula set.
cat("\n--- STAGE 5: joint copula competition (CV ES) ---\n"); t5 <- Sys.time()
fam_set <- c(0, 1, 2, 3, 4, 5, 6, 13, 14)
vine_fit <- function(U) tryCatch(suppressWarnings(VineCopula::RVineStructureSelect(
    U, familyset=fam_set, selectioncrit="AIC", indeptest=TRUE, level=0.05)),
    error=function(e) NULL)
es_one <- function(samp, Ute, p) {
  if (is.null(samp) || any(is.na(samp))) return(NA_real_)
  mean(vapply(seq_len(nrow(Ute)), function(i) {
    d1 <- sqrt(rowSums((samp - matrix(Ute[i,], N_SIM, p, byrow=TRUE))^2))
    d2 <- sqrt(rowSums((samp - samp[sample.int(N_SIM), , drop=FALSE])^2))
    mean(d1) - 0.5*mean(d2) }, numeric(1)), na.rm=TRUE)
}
joint_rows <- list()
for (sid in names(cohort_fits)) {
  cf <- cohort_fits[[sid]]; U <- cf$joint_aligned$U
  if (is.null(U) || nrow(U) < MIN_N || ncol(U) < 2) next
  n <- nrow(U); p <- ncol(U)
  set.seed(SEED); folds <- split(sample.int(n), rep(seq_len(CV_FOLDS), length.out=n))
  es_i <- es_g <- es_v <- numeric(0)
  for (fi in seq_along(folds)) {
    ti <- folds[[fi]]; Utr <- U[setdiff(seq_len(n), ti), , drop=FALSE]
    Ute <- U[ti, , drop=FALSE]
    samp_i <- matrix(runif(N_SIM*p), nrow=N_SIM)
    R <- tryCatch(stats::cor(qnorm(Utr)), error=function(e) diag(p))
    samp_g <- tryCatch(pnorm(matrix(rnorm(N_SIM*p), nrow=N_SIM) %*%
                              chol(R + diag(1e-8, p))), error=function(e) NULL)
    vfit <- vine_fit(Utr)
    samp_v <- if (is.null(vfit)) NULL else tryCatch(VineCopula::RVineSim(N_SIM, vfit),
                                                      error=function(e) NULL)
    es_i <- c(es_i, es_one(samp_i, Ute, p)); es_g <- c(es_g, es_one(samp_g, Ute, p))
    es_v <- c(es_v, es_one(samp_v, Ute, p))
  }
  scores <- c(Indep=mean(es_i, na.rm=TRUE), Gaussian=mean(es_g, na.rm=TRUE),
              Vine=mean(es_v, na.rm=TRUE))
  if (!any(is.finite(scores))) next
  best <- names(scores)[which.min(replace(scores, !is.finite(scores), Inf))]
  cohort_fits[[sid]]$joint_model <- if (best == "Vine") list(type="Vine", obj=vine_fit(U), U=U)
    else if (best == "Gaussian") list(type="Gaussian", obj=cor(qnorm(U)), U=U)
    else                         list(type="Indep",    obj=NULL,          U=U)
  joint_rows[[length(joint_rows)+1]] <- data.frame(sex=cf$sex, sport=cf$sport, test=cf$test, n=cf$n, p=p,
    Indep=signif(scores["Indep"],3), Gaussian=signif(scores["Gaussian"],3), Vine=signif(scores["Vine"],3),
    Best=best, GainOverIndep_pct=round((scores["Indep"]-scores[best])/scores["Indep"]*100, 1))
}
joint_competition <- do.call(rbind, joint_rows)
write.csv(joint_competition, file.path(VERIFY_DIR,"results/joint_competition.csv"), row.names=FALSE)
saveRDS(cohort_fits, file.path(VERIFY_DIR,"results/cohort_fits.rds"))
cohort_table <- do.call(rbind, lapply(cohort_fits, function(cf) data.frame(
  sex=cf$sex, sport=cf$sport, test=cf$test, n=cf$n,
  panel_size=length(cf$panel_metrics), p_cohort=length(cf$joint_metrics),
  joint_model=if (is.null(cf$joint_model)) "-" else cf$joint_model$type)))
write.csv(cohort_table, file.path(VERIFY_DIR,"results/cohort_table.csv"), row.names=FALSE)
cat(sprintf("  joint best counts: %s\n", paste(names(table(joint_competition$Best)),
            table(joint_competition$Best), sep="=", collapse=" ")))
record("Stage 5  joint copula", t5)
# STAGE 6 (paper Sec 3.6) -- Pick one Male Football CMJ athlete near the 90th
# Mahalanobis-depth percentile; per-metric centiles + 95% bootstrap CI on joint.
cat("\n--- STAGE 6: worked-example athlete scoring ---\n"); t6 <- Sys.time()
ex_cf <- cohort_fits[["Male_Football_CMJ"]]
stopifnot(!is.null(ex_cf), !is.null(ex_cf$joint_model))
mahal_pct <- function(Um, ux) {                        # used 3+ times below
  cm <- colMeans(Um); cv <- cov(Um) + diag(1e-8, ncol(Um))
  da <- as.numeric(stats::mahalanobis(matrix(ux, 1), cm, cv))
  list(d2=as.numeric(stats::mahalanobis(Um, cm, cv)),
       pct=100 * (1 - mean(stats::mahalanobis(Um, cm, cv) <= da)))
}
U <- ex_cf$joint_aligned$U
d2 <- mahal_pct(U, colMeans(U))$d2
ex_pid <- ex_cf$joint_aligned$profile_id[which.min(abs(d2 - quantile(d2, 0.9)))]
rec <- data_all[data_all$profile_id == ex_pid & data_all$gender == ex_cf$sex &
                  data_all$sport == ex_cf$sport, ]
pcols <- vapply(ex_cf$panel_metrics, function(m) m$col, character(1))
rec$cov <- rowSums(!is.na(rec[, pcols, drop=FALSE]))
rec <- rec[order(-rec$cov, as.Date(rec$test_date), decreasing=c(FALSE, TRUE)), ][1, ]
score_rows <- list(); uvec <- numeric(0); names_u <- character(0)
for (mt in ex_cf$panel_metrics) {
  yval <- suppressWarnings(as.numeric(rec[[mt$col]])); if (!is.finite(yval)) next
  cdf <- predict_cdf(mt$obj, yval, rec$age_at_test, if (mt$mass_dep) rec$bm else NULL)
  rel <- reliability_table[reliability_table$sex==ex_cf$sex & reliability_table$sport==ex_cf$sport &
            reliability_table$test==ex_cf$test & reliability_table$metric==mt$slug, ]
  band <- if (is.na(cdf)) "-" else c("Below cohort","Below median","Above median","Above cohort",
            "Elite tail")[findInterval(cdf, c(0.25,0.50,0.75,0.90))+1]
  score_rows[[length(score_rows)+1]] <- data.frame(metric=mt$slug, value=signif(yval, 4), units=mt$units,
    pct=round(100*cdf, 1), MDC95=if (nrow(rel)) signif(rel$MDC95[1], 3) else NA_real_, band=band)
  if (mt$slug %in% ex_cf$joint_metrics) { uvec <- c(uvec, cdf); names_u <- c(names_u, mt$slug) }
}
score_tbl <- do.call(rbind, score_rows)
common <- intersect(names_u, colnames(ex_cf$joint_model$U))
Ux <- ex_cf$joint_model$U[, common, drop=FALSE]; ux <- uvec[match(common, names_u)]
joint_pct <- mahal_pct(Ux, ux)$pct
boot <- replicate(N_BOOT, mahal_pct(
  Ux[sample.int(nrow(Ux), nrow(Ux), replace=TRUE), , drop=FALSE], ux)$pct)
joint_lo <- as.numeric(quantile(boot, 0.025)); joint_hi <- as.numeric(quantile(boot, 0.975))
cat(sprintf("  athlete=%s age=%.1f bm=%.1f kg joint=%.1f%% [%.1f, %.1f]\n",
            ex_pid, rec$age_at_test, rec$bm, joint_pct, joint_lo, joint_hi))
record("Stage 6  athlete scoring", t6)
# STAGE 7 (paper Sec 4-5) -- Asset generation: 4 EDA figs, worked-example
# centile/PIT/radar, body Tables 1-6, and inline number macros for the prose.
cat("\n--- STAGE 7: tables, figures, numbers ---\n"); t7 <- Sys.time()
fig_dir <- file.path(VERIFY_DIR,"figures"); tab_dir <- file.path(VERIFY_DIR,"tables")
num_dir <- file.path(VERIFY_DIR,"numbers")
# 7a) EDA figures (age histogram, cumulative coverage, JH-vs-age, missingness)
save_jpg(file.path(fig_dir,"fig_eda_age.jpg"), ggplot(data_all, aes(age_at_test)) +
  geom_histogram(bins=30, fill="#3F7AAD", colour="white", linewidth=0.3) +
  labs(x="Age (years)", y="Sessions") + theme_minimal(base_size=11), 4.5, 2.8)
cum_df <- data.frame(age=sort(data_all$age_at_test)); cum_df$cumf <- seq_along(cum_df$age)/nrow(cum_df)
save_jpg(file.path(fig_dir,"fig_eda_cumulative.jpg"),
  ggplot(cum_df, aes(age, cumf)) + geom_line(colour="#1F4F7F", linewidth=0.9) +
  scale_y_continuous(labels=scales::percent) + theme_minimal(base_size=11) +
  labs(x="Age (years)", y="Cumulative share of sessions"), 4.5, 2.8)
jh_df <- na.omit(data.frame(age=data_all$age_at_test, jh=data_all[["Jump.Height..Flight.Time._Both_CMJ"]]))
suppressMessages(save_jpg(file.path(fig_dir,"fig_eda_jh_age.jpg"),
  ggplot(jh_df, aes(age, jh)) + geom_point(alpha=0.25, colour="#444444", size=0.7) +
  geom_smooth(method="loess", se=FALSE, colour="#1F4F7F", linewidth=0.9) +
  labs(x="Age (years)", y="CMJ Jump Height (cm)") + theme_minimal(base_size=11), 4.5, 2.8))
miss_rows <- list()
for (cr in cohort_rows) for (k in seq_along(cr$panel_cols))
  miss_rows[[length(miss_rows)+1]] <- data.frame(test=cr$test, metric=cr$resolved[[k]]$slug,
    cohort=paste(cr$sex, cr$sport, sep="/"), pct_missing=mean(is.na(cr$d[[cr$panel_cols[k]]]))*100)
save_jpg(file.path(fig_dir,"fig_eda_missing.jpg"),
  ggplot(do.call(rbind, miss_rows), aes(metric, cohort, fill=pct_missing)) + geom_tile() +
  facet_grid(~ test, scales="free_x", space="free_x") +
  scale_fill_gradient(low="#EAF2F8", high="#1F4F7F", name="% missing") +
  theme_minimal(base_size=9) + labs(x=NULL, y=NULL) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)), 8.0, 5.0)
# 7b) Worked-example centile-bands triptych + per-metric PIT histograms + radar.
# For each metric: predict 5/25/50/75/95 centiles on a 50-point age grid (BM
# fixed at cohort median for mass-dependent metrics), then build the ribbon plot.
ages_grid <- seq(min(ex_cf$panel_metrics[[1]]$x), max(ex_cf$panel_metrics[[1]]$x), length.out=50)
qprobs <- c(0.05, 0.25, 0.50, 0.75, 0.95); band_plots <- list(); pit_long <- list()
for (mn in names(ex_cf$panel_metrics)) {
  pm <- ex_cf$panel_metrics[[mn]]
  bm_ref <- if (pm$mass_dep) rep(median(pm$bm, na.rm=TRUE), length(ages_grid)) else NULL
  if (pm$obj$type == "Linear") {
    mu_g <- predict(pm$obj$fit, newdata=df_xb(ages_grid, bm_ref))
    qmat <- sapply(qprobs, function(p) mu_g + qnorm(p)*stats::sd(residuals(pm$obj$fit)))
  } else {
    pa <- gamlss_pa(pm$obj$fit, ages_grid, bm_ref)
    pn <- intersect(c("mu","sigma","nu","tau"), names(pa))
    qfun <- gamlss_dfun(as.character(pm$obj$fit$family)[1], "q")
    qmat <- sapply(qprobs, function(p) do.call(qfun, c(list(p=p), pa[pn])))
  }
  bands <- data.frame(age=ages_grid, q05=qmat[,1], q25=qmat[,2], q50=qmat[,3], q75=qmat[,4], q95=qmat[,5])
  band_plots[[mn]] <- ggplot() +
    geom_ribbon(data=bands, aes(age, ymin=q05, ymax=q95), fill="#9DC3E6", alpha=0.55) +
    geom_ribbon(data=bands, aes(age, ymin=q25, ymax=q75), fill="#3F7AAD", alpha=0.55) +
    geom_line  (data=bands, aes(age, q50), colour="#1F4F7F", linewidth=0.9) +
    geom_point (data=data.frame(age=pm$x, y=pm$y), aes(age, y),
                alpha=0.3, size=0.7, colour="#444444") +
    labs(title=mn, x="Age (years)", y=pm$units) + theme_minimal(base_size=10)
  pit_long[[mn]] <- data.frame(metric=mn,
    pit=predict_cdf(pm$obj, pm$y, pm$x, if (pm$mass_dep) pm$bm else NULL))
}
jpeg(file.path(fig_dir,"fig_centiles_example.jpg"), width=7.0, height=2.6, units="in", res=300)
grid::grid.newpage(); grid::grid.draw(gridExtra::arrangeGrob(grobs=band_plots, nrow=1)); invisible(dev.off())
pit_df <- do.call(rbind, pit_long)
save_jpg(file.path(fig_dir,"fig_calibration_example.jpg"),
  ggplot(pit_df, aes(pit)) +
    geom_histogram(breaks=seq(0,1,0.1), fill="#3F7AAD", colour="white", linewidth=0.4) +
    geom_hline(yintercept=nrow(pit_df)/length(unique(pit_df$metric))/10,
               linetype="22", colour="grey40", linewidth=0.5) +
    facet_wrap(~ metric, nrow=1) + labs(x="PIT", y="Count") + theme_minimal(base_size=10),
  7.0, 2.4)
radar_pdat <- as.data.frame(rbind(rep(100, nrow(score_tbl)),
                                    rep(0, nrow(score_tbl)), score_tbl$pct))
colnames(radar_pdat) <- score_tbl$metric
jpeg(file.path(fig_dir,"fig_radar_example.jpg"), width=6.5, height=5.5, units="in", res=300)
par(mar=c(2,4,3,4))
fmsb::radarchart(radar_pdat, axistype=1, pcol="#1F4F7F", pfcol=scales::alpha("#3F7AAD", 0.45),
  plwd=2.5, cglcol="grey75", cglty=1, cglwd=0.6, axislabcol="grey40",
  caxislabels=c("0","25","50","75","100"), vlcex=1.05, calcex=0.95,
  title=sprintf("Athlete radar (%s %s %s)", ex_cf$sex, ex_cf$sport, ex_cf$test))
invisible(dev.off())
# 7c) Body Tables 1-6.  Vectorised sprintf builds the rows; tex_table assembles
# the booktabs tabular and writes the file.
tex_table(file.path(tab_dir,"tab_reliability.tex"), "llrrr",
  "Test & Metric & median ICC$_{(3,1)}$ & median MDC$_{95}$ & $k$ cohorts \\\\",
  with(reliability_summary, sprintf("%s & %s & %s & %s & %s \\\\",
       tex_safe(test), tex_safe(metric), median_ICC, median_MDC95, n_cohorts)))
vs_df <- cohort_table; vs_df$dropped <- vs_df$panel_size - vs_df$p_cohort
tex_table(file.path(tab_dir,"tab_variable_selection.tex"), "lllrrr",
  "Sex & Sport & Test & Panel & Joint $p$ & Dropped \\\\",
  with(vs_df, sprintf("%s & %s & %s & %s & %s & %s \\\\",
       sex, tex_safe(sport), test, panel_size, p_cohort, dropped)))
fq_df <- marginal_competition |> dplyr::group_by(test, metric) |>
  dplyr::summarise(med_lin=round(median(Linear, na.rm=TRUE), 3),
                   med_no =round(median(NO,     na.rm=TRUE), 3),
                   med_bct=round(median(BCT,    na.rm=TRUE), 3),
                   wins_BCT=sum(Best == "BCT"), .groups="drop")
tex_table(file.path(tab_dir,"tab_fit_quality.tex"), "llrrrr",
  "Test & Metric & median Linear & median NO & median BCT & BCT wins \\\\",
  with(fq_df, sprintf("%s & %s & %s & %s & %s & %s \\\\",
       tex_safe(test), tex_safe(metric), med_lin, med_no, med_bct, wins_BCT)))
cm_df <- marginal_competition |> dplyr::rowwise() |>
  dplyr::mutate(best_score = min(c(Linear, NO, BCT), na.rm=TRUE)) |> dplyr::ungroup() |>
  dplyr::group_by(sex, sport, test) |>
  dplyr::summarise(C_M=round(median(best_score, na.rm=TRUE), 3), .groups="drop")
cs_df <- merge(cm_df, joint_competition[, c("sex","sport","test","Best","Indep","Gaussian","Vine")],
               by=c("sex","sport","test"), all.x=TRUE)
cs_df$C_J <- with(cs_df, ifelse(Best=="Vine", Vine, ifelse(Best=="Gaussian", Gaussian, Indep)))
tex_table(file.path(tab_dir,"tab_calibration_scores.tex"), "lllrrl",
  "Sex & Sport & Test & $C_M$ & $C_J$ & Joint best \\\\",
  with(cs_df, sprintf("%s & %s & %s & %s & %s & %s \\\\", sex, tex_safe(sport), test,
       C_M, ifelse(is.na(C_J), "-", as.character(signif(C_J, 3))),
       ifelse(is.na(Best), "-", Best))))
tex_table(file.path(tab_dir,"tab_score_example.tex"), "lrlrrl",
  "Metric & Value & Units & Pct & MDC$_{95}$ & Band \\\\",
  with(score_tbl, sprintf("\\texttt{%s} & %s & %s & %s & %s & %s \\\\",
       tex_safe(metric), value, units, pct, MDC95, band)),
  sprintf("\\multicolumn{6}{l}{\\textbf{Joint percentile (%s):} %.1f\\%% \\quad 95\\%% bootstrap CI [%.1f, %.1f]} \\\\",
          tolower(ex_cf$joint_model$type), joint_pct, joint_lo, joint_hi))
cov_df <- cohort_table |> dplyr::group_by(sport, test) |>
  dplyr::summarise(F=sum(sex=="Female"), M=sum(sex=="Male"), .groups="drop")
tex_table(file.path(tab_dir,"tab_cohort_coverage.tex"), "llrr",
  "Sport & Test & Female & Male \\\\",
  with(cov_df, sprintf("%s & %s & %s & %s \\\\", tex_safe(sport), test, F, M)))
# 7d) Inline number macros (consumed by paper/main.tex prose)
crps_med <- median(apply(marginal_competition[,c("Linear","NO","BCT")], 1, min, na.rm=TRUE), na.rm=TRUE)
es_med   <- median(apply(joint_competition[,c("Indep","Gaussian","Vine")], 1, min, na.rm=TRUE), na.rm=TRUE)
writeLines(sprintf("\\newcommand{\\%s}{%s}",
  c("Nathletes","Ntests","Nsports","Ncohorts","CRPSmed","ESmed","NjointVineWins","NjointEligible"),
  c(length(unique(data_all$profile_id)), length(unique(data_all$test_id)), length(unique(data_all$sport)),
    nrow(cohort_table), sprintf("%.3f", crps_med), sprintf("%.3f", es_med),
    sum(joint_competition$Best == "Vine"), nrow(joint_competition))), file.path(num_dir,"numbers.tex"))
writeLines(sprintf("\\newcommand{\\%s}{%s}",
  c("AthleteAge","AthleteBM","AthleteJoint","AthleteJointLo","AthleteJointHi","AthleteJM"),
  c(sprintf("%.1f", c(rec$age_at_test, rec$bm, joint_pct, joint_lo, joint_hi)), ex_cf$joint_model$type)),
  file.path(num_dir,"numbers_athlete.tex"))
record("Stage 7  asset generation", t7)

# STAGE 8 -- Score EVERY athlete-test in EVERY (sex x sport x test) cohort.
# Marginal centile per panel metric uses the same predict_cdf() winner from
# Stage 4; joint percentile uses Mahalanobis depth on the cohort PIT matrix
# Ux restricted to whichever joint metrics the athlete actually has on that
# test occasion (so partially-missing rows are still scored).  Outputs:
#   verification/results/athlete_scores_long.csv     (1 row per athlete-test-metric)
#   verification/results/athlete_scores_wide.csv     (1 row per athlete-test occasion)
#   verification/results/athlete_scores_summary.csv  (per-cohort headline counts)
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
mahal_one <- function(Um, ux) {                    # depth percentile in [0,100]
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
    cdf[keep] <- tryCatch(
      predict_cdf(mt$obj, yv[keep], xv[keep], if (mt$mass_dep) bv[keep] else NULL),
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
}
athlete_long <- do.call(rbind, long_rows)
athlete_wide <- dplyr::bind_rows(wide_rows)
write.csv(athlete_long, file.path(VERIFY_DIR, "results/athlete_scores_long.csv"), row.names = FALSE)
write.csv(athlete_wide, file.path(VERIFY_DIR, "results/athlete_scores_wide.csv"), row.names = FALSE)
score_summary <- athlete_wide |>
  dplyr::group_by(gender, sport, test_type, joint_model) |>
  dplyr::summarise(n_test_rows      = dplyr::n(),
                   n_athletes       = dplyr::n_distinct(profile_id),
                   n_joint_scored   = sum(!is.na(joint_pct)),
                   median_joint_pct = round(median(joint_pct, na.rm = TRUE), 1),
                   pct_elite_joint  = round(100 * mean(joint_pct >= 90, na.rm = TRUE), 1),
                   .groups = "drop")
write.csv(score_summary, file.path(VERIFY_DIR, "results/athlete_scores_summary.csv"),
          row.names = FALSE)
cat(sprintf("  rows long=%d  rows wide=%d  athletes=%d  with joint pct=%d\n",
            nrow(athlete_long), nrow(athlete_wide),
            length(unique(athlete_wide$profile_id)),
            sum(!is.na(athlete_wide$joint_pct))))
record("Stage 8  athlete scoring", t8)

# STAGE 9 -- Per-athlete radar plot for EVERY athlete-test occasion.  One PNG
# per row of athlete_wide, organised as
#   verification/figures/athletes/<sex>/<sport>/<test>/<profile_id>__<date>.png
# Vertices = panel metrics; ring = marginal centile (0-100); colour = joint
# percentile band.
cat("\n--- STAGE 9: per-athlete radar plots ---\n"); t9 <- Sys.time()
ath_root <- file.path(VERIFY_DIR, "figures", "athletes")
dir.create(ath_root, recursive = TRUE, showWarnings = FALSE)
band_colour <- c("Atypical multivariable"      = "#cb181d",
                 "Typical (below)"             = "#fdae61",
                 "Typical (above)"             = "#fee08b",
                 "Above-cohort multivariable"  = "#a6d96a",
                 "Elite multivariable"         = "#1a9850",
                 "-"                           = "#969696")
n_radar <- 0
for (sid in names(cohort_fits)) {
  cf  <- cohort_fits[[sid]]; pms <- cf$panel_metrics
  if (length(pms) < 3) next                # radar needs >= 3 axes
  pslugs <- vapply(pms, function(m) m$slug, character(1))
  this <- athlete_wide[athlete_wide$gender == cf$sex &
                         athlete_wide$sport  == cf$sport &
                         athlete_wide$test_type == cf$test, , drop = FALSE]
  if (!nrow(this)) next
  outdir <- file.path(ath_root, gsub("[^A-Za-z0-9]+", "_", cf$sex),
                       gsub("[^A-Za-z0-9]+", "_", cf$sport),
                       gsub("[^A-Za-z0-9]+", "_", cf$test))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pct_cols <- paste0("pct_", pslugs)
  vlabels  <- gsub("_", " ", pslugs)
  for (i in seq_len(nrow(this))) {
    pcts <- as.numeric(this[i, pct_cols])
    pcts[!is.finite(pcts)] <- 50            # neutral fill so radar still draws
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
}
cat(sprintf("  radar PNGs written: %d  under %s\n", n_radar, ath_root))
record("Stage 9  athlete radars", t9)

# Final summary -- per-stage timings + every file written under verification/
cat("\n=========================================================\n  Verification complete\n=========================================================\n")
for (nm in names(TIMINGS)) cat(sprintf("  %-32s %7.1f s\n", nm, TIMINGS[[nm]]))
cat(sprintf("  %-32s %7.1f s\n", "TOTAL", sum(unlist(TIMINGS))))
fmt_size <- function(b) if (is.na(b)) "?" else if (b > 1048576) sprintf("%.1f MB", b/1048576) else
                        if (b > 1024) sprintf("%.1f KB", b/1024) else sprintf("%d B", as.integer(b))
ff <- sort(sub(paste0(VERIFY_DIR,"/?"), "", list.files(VERIFY_DIR, recursive=TRUE, full.names=TRUE)))
cat(sprintf("\nFiles written under %s/:\n", VERIFY_DIR)); maxw <- max(nchar(ff))
for (p in ff) cat(sprintf("  %-*s  %s\n", maxw, p, fmt_size(file.info(file.path(VERIFY_DIR, p))$size)))
invisible(NULL)
