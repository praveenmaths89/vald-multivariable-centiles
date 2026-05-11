# ------------------------------------------------------------
# 06_centile_profiles.R
# Translate marginals + vine into multivariable centile
# diagnostics:
#   - per-metric univariate centile bands (3/10/25/50/75/90/97)
#     conditional on age, faceted male vs female.
#   - joint percentile per athlete via Mahalanobis depth on
#     PIT residuals (or Tukey depth if DepthProc available),
#     summarising "how typical" the entire profile is.
#   - radar plots of representative athletes (high / median /
#     low joint depth) per stratum.
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
strata_dir <- file.path(RESULTS_DIR, "strata")

centile_grid <- function(fit, x_seq) {
  rs <- robust_sample_predictive(fit$fit, x_seq, fit$y, N_SIM, fit$family)
  if (is.null(rs$samp)) return(NULL)
  qs <- t(apply(rs$samp, 1, function(v) {
    stats::quantile(v, probs = c(.03, .10, .25, .50, .75, .90, .97),
                    na.rm = TRUE, names = FALSE)
  }))
  data.frame(
    age = x_seq,
    c03 = qs[, 1], c10 = qs[, 2], c25 = qs[, 3],
    c50 = qs[, 4], c75 = qs[, 5], c90 = qs[, 6], c97 = qs[, 7]
  )
}

plot_centiles <- function(stratum, out_path) {
  ms <- stratum$metrics
  if (!length(ms)) return(invisible())
  age_lo <- min(unlist(lapply(ms, function(m) min(m$x))), na.rm = TRUE)
  age_hi <- max(unlist(lapply(ms, function(m) max(m$x))), na.rm = TRUE)
  x_seq  <- seq(age_lo, age_hi, length.out = 120)

  panels <- lapply(names(ms), function(nm) {
    o <- ms[[nm]]
    cf <- centile_grid(o, x_seq)
    if (is.null(cf)) return(NULL)
    data.frame(metric = nm, x = o$x, y = o$y) -> obs
    ggplot() +
      geom_point(data = obs, aes(x, y), color = "grey50",
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
  })
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (!length(panels)) return(invisible())
  ggsave(out_path,
         gridExtra::arrangeGrob(grobs = panels,
                                 ncol = min(2, length(panels))),
         width = 10, height = 3.4 * ceiling(length(panels) / 2),
         dpi = 130, limitsize = FALSE)
}

joint_depth <- function(U) {
  if (HAS_DEPTHPROC) {
    dv <- tryCatch(DepthProc::depthTukey(U, U), error = function(e) NULL)
    if (!is.null(dv)) return(rank(dv) / length(dv))
  }
  md <- stats::mahalanobis(U, colMeans(U), stats::cov(U))
  1 - rank(md) / length(md)
}

radar_for_athletes <- function(U, ids, ttl) {
  if (length(ttl) != 3) return(invisible())
  pct <- joint_depth(U)
  ord <- order(-pct)
  picks <- c(ord[1], ord[ceiling(length(ord) / 2)], ord[length(ord)])
  rdf_list <- list()
  for (k in seq_along(picks)) {
    v <- as.numeric(U[picks[k], ]) * 100
    rdf_list[[k]] <- data.frame(rbind(rep(100, ncol(U)), rep(0, ncol(U)), v))
    colnames(rdf_list[[k]]) <- gsub("_Both_|_Right_|_Left_", "_",
                                    colnames(U), perl = TRUE)
  }
  list(rdf_list = rdf_list, picks = picks, pct = pct[picks], titles = ttl)
}

main <- function() {
  log_step("Stage 06 — multivariable centile profiles")
  if (!file.exists(fits_path)) {
    env04 <- new.env(); sys.source(file.path(CODE_DIR, "04_marginals_gamlss.R"), envir = env04)
    invisible(env04$main())
  }
  fits <- readRDS(fits_path)

  for (sid in names(fits)) {
    s <- fits[[sid]]
    if (length(s$metrics) < 2) next
    sd_dir <- file.path(strata_dir, sid)
    dir.create(sd_dir, recursive = TRUE, showWarnings = FALSE)

    plot_centiles(s, file.path(sd_dir, "centiles.png"))
    plot_centiles(s, file.path(sd_dir, "centiles.pdf"))

    # Joint depth percentile based on PIT residuals
    n  <- min(vapply(s$metrics, function(m) length(m$pit), integer(1)))
    U  <- vapply(s$metrics, function(m) m$pit[seq_len(n)], numeric(n))
    colnames(U) <- names(s$metrics)
    U <- U[stats::complete.cases(U), , drop = FALSE]
    if (nrow(U) >= MIN_N_STRATUM) {
      depth_pct <- round(joint_depth(U) * 100, 1)
      utils::write.csv(
        data.frame(joint_pctile = depth_pct, U),
        file.path(sd_dir, "joint_percentiles.csv"), row.names = FALSE
      )

      # radar
      r <- radar_for_athletes(U, NULL,
                              c("High joint depth", "Median joint depth", "Low joint depth"))
      if (!is.null(r)) {
        png(file.path(sd_dir, "radar.png"), width = 1500, height = 500, res = 130)
        op <- par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))
        clr <- c("#238b45", PAL$median, "#cb181d")
        for (k in seq_along(r$rdf_list)) {
          tryCatch(
            fmsb::radarchart(
              r$rdf_list[[k]], pcol = clr[k],
              pfcol = grDevices::adjustcolor(clr[k], 0.25),
              plwd = 2, plty = 1, cglcol = "grey75", cglty = 1,
              axislabcol = "grey40", vlcex = 0.5,
              title = sprintf("%s\n(joint pct %.0f%%)",
                              r$titles[k], r$pct[k] * 100)
            ),
            error = function(e) NULL
          )
        }
        par(op); dev.off()
        pdf(file.path(sd_dir, "radar.pdf"), width = 12, height = 4)
        op <- par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))
        for (k in seq_along(r$rdf_list)) {
          tryCatch(
            fmsb::radarchart(
              r$rdf_list[[k]], pcol = clr[k],
              pfcol = grDevices::adjustcolor(clr[k], 0.25),
              plwd = 2, plty = 1, cglcol = "grey75", cglty = 1,
              axislabcol = "grey40", vlcex = 0.5,
              title = sprintf("%s (joint pct %.0f%%)",
                              r$titles[k], r$pct[k] * 100)
            ),
            error = function(e) NULL
          )
        }
        par(op); dev.off()
      }
    }
  }
  log_step("Stage 06 done.")
}

if (sys.nframe() == 0) main()
