# ------------------------------------------------------------
# 03_eda_strata.R
# Per-stratum exploratory analysis: counts, age distribution,
# missingness map, and metric correlations within each
# (Sex × Sport × Test) cohort.
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

read_strata <- function() {
  files <- list.files(per_test_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No per-test CSVs in ", per_test_dir,
                            ". Run stage 02 first.")
  res <- list()
  for (f in files) {
    test_type <- tools::file_path_sans_ext(basename(f))
    d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    d$gender <- clean_gender(d$gender)
    d$sport  <- clean_sport(d$sport)
    d <- d[!is.na(d$gender) & !is.na(d$sport) &
            !is.na(d$age_at_test) & d$age_at_test >= 8 & d$age_at_test <= 60, ,
            drop = FALSE]
    res[[test_type]] <- d
  }
  res
}

stratum_counts <- function(strata) {
  rows <- list()
  for (tt in names(strata)) {
    d <- strata[[tt]]
    if (!nrow(d)) next
    tab <- dplyr::count(d, gender, sport, name = "n")
    tab$test_type <- tt
    rows[[length(rows) + 1]] <- tab
  }
  out <- do.call(rbind, rows)
  out[order(-out$n), ]
}

main <- function() {
  log_step("Stage 03 — EDA per (sex × sport × test)")
  strata <- read_strata()

  # ---- Counts ----
  ct <- stratum_counts(strata)
  utils::write.csv(ct, file.path(RESULTS_DIR, "03_stratum_counts.csv"),
                   row.names = FALSE)
  ct_eligible <- ct[ct$n >= MIN_N_STRATUM, , drop = FALSE]
  utils::write.csv(ct_eligible, file.path(RESULTS_DIR, "03_eligible_strata.csv"),
                   row.names = FALSE)
  log_step(sprintf("Total combos: %d  Eligible (n >= %d): %d",
                   nrow(ct), MIN_N_STRATUM, nrow(ct_eligible)))

  # ---- Heatmap of n by (sport × test), faceted by sex ----
  if (nrow(ct) > 0) {
    p <- ggplot(ct, aes(test_type, sport, fill = n)) +
      geom_tile(color = "white") +
      geom_text(aes(label = n), size = 3) +
      facet_wrap(~ gender) +
      scale_fill_viridis_c(option = "viridis", direction = 1, na.value = "grey90") +
      labs(title = "Eligible-cohort sample sizes",
           subtitle = paste0("Strata used downstream require n >= ",
                              MIN_N_STRATUM),
           x = "Test type", y = "Sport") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    ggsave(file.path(RESULTS_DIR, "03_counts_heatmap.png"),
           p, width = 12, height = 5.5, dpi = 130)
    ggsave(file.path(PAPER_FIGS, "fig_counts_heatmap.pdf"),
           p, width = 11, height = 5.5)
  }

  # ---- Age distribution per stratum ----
  ages <- list()
  for (tt in names(strata)) {
    d <- strata[[tt]]
    if (!nrow(d)) next
    d$test_type <- tt
    ages[[tt]]  <- d[, c("gender", "sport", "test_type", "age_at_test")]
  }
  age_df <- do.call(rbind, ages)
  if (nrow(age_df)) {
    age_df <- age_df[paste(age_df$gender, age_df$sport, age_df$test_type) %in%
                      paste(ct_eligible$gender, ct_eligible$sport, ct_eligible$test_type), ,
                      drop = FALSE]
    if (nrow(age_df)) {
      p <- ggplot(age_df, aes(age_at_test, fill = gender, color = gender)) +
        geom_density(alpha = 0.35, linewidth = 0.5, na.rm = TRUE) +
        scale_color_manual(values = c(Female = PAL$female, Male = PAL$male)) +
        scale_fill_manual(values  = c(Female = PAL$female, Male = PAL$male)) +
        facet_grid(test_type ~ sport, scales = "free_y", switch = "y") +
        labs(title = "Age distribution within each eligible cohort",
             x = "Age (years)", y = NULL) +
        theme(strip.text.y.left = element_text(angle = 0, size = 7),
              strip.text.x = element_text(size = 8),
              legend.position = "bottom")
      ggsave(file.path(RESULTS_DIR, "03_age_distributions.png"),
             p, width = 14, height = 8, dpi = 130, limitsize = FALSE)
      ggsave(file.path(PAPER_FIGS, "fig_age_distributions.pdf"),
             p, width = 12, height = 7, limitsize = FALSE)
    }
  }

  log_step("Stage 03 done.")
  invisible(list(strata = strata, counts = ct,
                 eligible = ct_eligible))
}

if (sys.nframe() == 0) invisible(main())
