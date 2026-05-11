# ------------------------------------------------------------
# 02_integrate_sport_gender.R
# Build the analysis-ready dataset by joining VALD trials with
# (a) the VALD profile metadata and (b) external sport / gender
# CSVs supplied by the federation.
#
# If the raw VALD pull from stage 01 is unavailable, this stage
# can use an already-integrated CSV via env var
# VALD_INTEGRATED_CSV.  Eligible per-test CSVs (02_<TEST>_dataset.csv)
# are written to data/per_test/.
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

raw_path     <- file.path(DATA_DIR, "raw_forcedecks.rds")
sports_dir   <- Sys.getenv("VALD_SPORTS_DIR", unset = "")
external_csv <- Sys.getenv("VALD_INTEGRATED_CSV", unset = "")
per_test_dir <- file.path(DATA_DIR, "per_test")
dir.create(per_test_dir, recursive = TRUE, showWarnings = FALSE)

load_external_sports <- function(folder) {
  if (!dir.exists(folder)) return(NULL)
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) return(NULL)
  rows <- list()
  for (f in files) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    pid_col <- names(df)[grepl("profile.*id", names(df), ignore.case = TRUE)][1]
    sex_col <- names(df)[grepl("^gender$", names(df), ignore.case = TRUE)][1]
    nm_col  <- names(df)[grepl("^name$", names(df), ignore.case = TRUE)][1]
    if (any(is.na(c(pid_col, sex_col, nm_col)))) next
    sport <- sub("\\..*$", "",
                 sub("_All.*$", "",
                     sub("_Overall.*$", "", basename(f))))
    rows[[length(rows) + 1]] <- data.frame(
      profile_id = as.character(df[[pid_col]]),
      gender     = clean_gender(df[[sex_col]]),
      name       = df[[nm_col]],
      sport      = clean_sport(sport),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(NULL)
  out <- dplyr::distinct(do.call(rbind, rows))
  out[!is.na(out$profile_id) & out$profile_id != "", , drop = FALSE]
}

# Merge raw API tables -> long test-level data -----------------
build_from_api <- function(raw, sports) {
  trials <- raw$trials
  tests  <- raw$tests
  profs  <- raw$profiles

  trials$athleteId <- as.character(trials$athleteId)
  joined <- dplyr::left_join(trials, tests, by = "testId")

  agg <- joined |>
    dplyr::mutate(
      value  = suppressWarnings(as.numeric(value)),
      weight = suppressWarnings(as.numeric(weight))
    ) |>
    dplyr::group_by(athleteId, testId, testType, recordedUTC,
                    recordedDateOffset, trialLimb, definition_name) |>
    dplyr::summarise(
      mean_result = mean(value, na.rm = TRUE),
      mean_weight = mean(weight, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      ts_utc   = lubridate::ymd_hms(recordedUTC),
      ts_local = ts_utc + lubridate::minutes(recordedDateOffset),
      test_date = as.Date(ts_local)
    )

  wide <- agg |>
    dplyr::select(profile_id = athleteId, test_id = testId,
                  test_type = testType, test_date,
                  limb = trialLimb, metric = definition_name,
                  value = mean_result, weight = mean_weight) |>
    tidyr::pivot_wider(
      id_cols     = c(profile_id, test_date, test_id, weight),
      names_from  = c(metric, limb, test_type),
      values_from = value,
      names_glue  = "{metric}_{limb}_{test_type}",
      names_repair = "universal"
    )

  if (!is.null(profs)) {
    profs$profileId <- as.character(profs$profileId)
    name_full <- if (all(c("givenName", "familyName") %in% names(profs))) {
      paste(profs$givenName, profs$familyName)
    } else if ("name" %in% names(profs)) profs$name else NA_character_
    profs_clean <- data.frame(
      profile_id = profs$profileId,
      name       = name_full,
      height_cm  = if ("heightInCm" %in% names(profs)) profs$heightInCm else NA_real_,
      weight_kg  = if ("weightInKg" %in% names(profs)) profs$weightInKg else NA_real_,
      dob        = if ("dateOfBirth" %in% names(profs)) as.Date(profs$dateOfBirth) else as.Date(NA),
      stringsAsFactors = FALSE
    )
    wide <- dplyr::left_join(wide, profs_clean, by = "profile_id")
  }

  if (!is.null(sports)) {
    sports_uniq <- sports |>
      dplyr::filter(!is.na(gender), !is.na(sport)) |>
      dplyr::group_by(profile_id) |>
      dplyr::slice(1) |>
      dplyr::ungroup()
    wide <- dplyr::left_join(wide, sports_uniq[, c("profile_id", "gender", "sport")],
                              by = "profile_id")
  }
  wide$age_at_test <- as.numeric(difftime(wide$test_date, wide$dob, units = "days")) / 365.25
  wide
}

write_per_test_files <- function(df, out_dir) {
  meta <- intersect(c("profile_id", "name", "gender", "sport", "dob",
                      "age_at_test", "test_date", "test_id", "weight",
                      "height_cm", "weight_kg"), names(df))
  metric_cols <- setdiff(names(df), meta)
  test_types  <- unique(sub(".*_([A-Z]+)$", "\\1", metric_cols))
  test_types  <- test_types[nzchar(test_types) & !is.na(test_types)]
  for (tt in test_types) {
    cols <- metric_cols[grepl(paste0("_", tt, "$"), metric_cols)]
    sub  <- df[, c(meta, cols), drop = FALSE]
    keep <- rowSums(!is.na(sub[, cols, drop = FALSE])) > 0
    if (sum(keep) >= 10) {
      out <- file.path(out_dir, sprintf("%s.csv", tt))
      utils::write.csv(sub[keep, , drop = FALSE], out, row.names = FALSE)
      log_step(sprintf("  per-test: %s  rows=%d", tt, sum(keep)))
    }
  }
}

main <- function() {
  log_step("Stage 02 — integrate sport / gender")
  sports <- load_external_sports(sports_dir)
  log_step("External sport/gender records: ",
           if (is.null(sports)) "<none>" else nrow(sports))

  if (file.exists(raw_path)) {
    raw <- readRDS(raw_path)
    integrated <- build_from_api(raw, sports)
    log_step(sprintf("Integrated from API: %d rows × %d cols",
                     nrow(integrated), ncol(integrated)))
  } else if (nzchar(external_csv) && file.exists(external_csv)) {
    log_step("Falling back to pre-integrated CSV: ", external_csv)
    integrated <- read.csv(external_csv, stringsAsFactors = FALSE,
                            check.names = FALSE)
    integrated$gender <- clean_gender(integrated$gender)
    integrated$sport  <- clean_sport(integrated$sport)
  } else {
    stop(
      "No integrated data source available. Either:\n",
      "  (a) point VALD_INTEGRATED_CSV at a pre-integrated wide CSV, or\n",
      "  (b) run code/01_pull_vald.R first (requires VALD_CLIENT_ID, ",
      "VALD_CLIENT_SECRET, VALD_TENANT_ID env vars).\n",
      "See data/README.md for the expected schema.",
      call. = FALSE)
  }

  integrated <- integrated[!is.na(integrated$gender) &
                              !is.na(integrated$sport), , drop = FALSE]
  utils::write.csv(integrated,
                   file.path(DATA_DIR, "analysis_dataset.csv"),
                   row.names = FALSE)
  log_step(sprintf("Wrote analysis_dataset.csv (%d rows, %d cols)",
                   nrow(integrated), ncol(integrated)))

  write_per_test_files(integrated, per_test_dir)
  log_step("Stage 02 done.")
}

if (sys.nframe() == 0) main()
