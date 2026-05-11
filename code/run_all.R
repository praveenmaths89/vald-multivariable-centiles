# ------------------------------------------------------------
# run_all.R — orchestrate the full pipeline.
# Usage:
#   Rscript code/run_all.R           # default: skip stage 01 if data already integrated
#   Rscript code/run_all.R --pull    # also try to pull from VALD API
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
do_pull <- "--pull" %in% args

.this_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(fa)) dirname(fa[1])
  else if (!is.null(sf <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL))) dirname(sf)
  else getwd()
}
source(file.path(.this_dir, "00_setup.R"))
log_step("Project: ", PROJECT_DIR)

run_stage <- function(file) {
  full <- file.path(CODE_DIR, file)
  log_step("--- ", file)
  env <- new.env()
  sys.source(full, envir = env)
  if (is.function(env$main)) {
    invisible(env$main())
  }
}

if (do_pull) run_stage("01_pull_vald.R")
run_stage("02_integrate_sport_gender.R")
run_stage("03_eda_strata.R")
run_stage("04_marginals_gamlss.R")
run_stage("05_joint_copula.R")
run_stage("06_centile_profiles.R")
run_stage("07_validation.R")
run_stage("08_paper_assets.R")

log_step("All stages complete.")
log_step("Compile paper:    cd paper && latexmk -pdf main.tex")
