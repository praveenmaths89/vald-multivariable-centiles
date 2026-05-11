# ------------------------------------------------------------
# 01_pull_vald.R
# Pull profiles, tests, trials from the VALD ForceDecks API
# using `valdr`. Credentials must be supplied via environment
# variables — never paste them into source files.
#
#   VALD_CLIENT_ID
#   VALD_CLIENT_SECRET
#   VALD_TENANT_ID
#   VALD_REGION              (default: aue)
#   VALD_START_DATE          (default: 2020-01-01T00:00:00Z)
#
# If `valdr` is not installed or credentials are absent, the
# script falls back to an existing local CSV at
# ${VALD_FORCEDECKS_RAW} (a previously cached pull) or under
# data/raw_forcedecks_*.rds.
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

raw_path <- file.path(DATA_DIR, "raw_forcedecks.rds")
cache_csv <- Sys.getenv("VALD_FORCEDECKS_RAW", unset = "")

pull_with_valdr <- function() {
  client_id     <- Sys.getenv("VALD_CLIENT_ID")
  client_secret <- Sys.getenv("VALD_CLIENT_SECRET")
  tenant_id     <- Sys.getenv("VALD_TENANT_ID")
  region        <- Sys.getenv("VALD_REGION", unset = "aue")
  start_date    <- Sys.getenv("VALD_START_DATE",
                              unset = "2020-01-01T00:00:00Z")

  if (!nzchar(client_id) || !nzchar(client_secret) || !nzchar(tenant_id)) {
    return(NULL)
  }
  if (!HAS_VALDR) {
    log_step("valdr not installed; cannot pull from API.")
    return(NULL)
  }
  options(timeout = 18000)
  valdr::set_credentials(
    client_id = client_id, client_secret = client_secret,
    tenant_id = tenant_id, region = region
  )
  valdr::set_start_date(start_date)
  log_step("Calling valdr::get_forcedecks_data() …")
  data <- valdr::get_forcedecks_data(include_attributes = TRUE)
  list(
    profiles = as.data.frame(data$profiles),
    tests    = as.data.frame(data$tests),
    trials   = as.data.frame(data$trials)
  )
}

main <- function() {
  log_step("Stage 01 — pull VALD")

  data <- pull_with_valdr()
  if (is.null(data) && nzchar(cache_csv) && file.exists(cache_csv)) {
    log_step("Loading cached CSV from VALD_FORCEDECKS_RAW.")
    data <- readRDS(cache_csv)
  }
  if (is.null(data)) {
    log_step("No live API call possible and no cache configured. ",
             "Pipeline will fall back to integrated CSV in stage 02 ",
             "(VALD_INTEGRATED_CSV).")
    return(invisible(NULL))
  }

  saveRDS(data, raw_path)
  log_step("Saved raw pull to ", raw_path)
  log_step(sprintf("Profiles: %d  Tests: %d  Trials: %d",
                   nrow(data$profiles), nrow(data$tests), nrow(data$trials)))
}

if (sys.nframe() == 0) main()
