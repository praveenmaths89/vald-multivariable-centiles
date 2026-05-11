# =============================================================================
#                       00_run_all.R  --  pipeline driver
# Runs the full STAI-X 2026 pipeline end-to-end and at the end calls the
# paper-asset audit (99_check_paper_assets.R).
#
# Usage
#   Rscript code/00_run_all.R                     # default: assumes
#                                                 #   data/analysis_dataset.csv
#                                                 #   already on disk
#   Rscript code/00_run_all.R --from-pull         # also try a fresh VALD API
#                                                 #   pull (requires VALD_*
#                                                 #   env vars)
#   Rscript code/00_run_all.R --skip-audit        # do not run the asset audit
#   Rscript code/00_run_all.R --no-mtime-check    # audit ignores mtimes
#
# Exit code: 0 on full success (incl. asset audit pass), 1 on any failure.
# =============================================================================

args        <- commandArgs(trailingOnly = TRUE)
DO_PULL     <- "--from-pull"     %in% args
SKIP_AUDIT  <- "--skip-audit"    %in% args
NO_MTIME    <- "--no-mtime-check" %in% args

# ---- Repo root anchor (portable; no absolute user paths) --------------------
.anchor_root <- function() {
  envroot <- Sys.getenv("VALD_PROJECT_DIR", unset = "")
  if (nzchar(envroot) && dir.exists(envroot)) return(normalizePath(envroot))
  .ca <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", .ca[grep("^--file=", .ca)])
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
CODE_DIR    <- file.path(PROJECT_DIR, "code")

banner <- function(title, ch = "=") {
  cat("\n", strrep(ch, 78), "\n", sep = "")
  cat("  ", title, "\n", sep = "")
  cat(strrep(ch, 78), "\n", sep = "")
}

fail <- function(msg) {
  cat("\nFATAL: ", msg, "\n", sep = "")
  quit(status = 1L)
}

human_dur <- function(secs) {
  if (secs < 60) sprintf("%.1fs", secs)
  else if (secs < 3600) sprintf("%dm %02ds", as.integer(secs / 60),
                                  as.integer(secs %% 60))
  else sprintf("%dh %02dm", as.integer(secs / 3600),
                as.integer((secs %% 3600) / 60))
}

# =============================================================================
# Pre-flight
# =============================================================================
banner("Pre-flight checks")

# 1. R version
rv <- getRversion()
cat(sprintf("  R version: %s\n", rv))
if (rv < "4.2.0")
  cat("  WARNING: R >= 4.2 is recommended (you have ", as.character(rv),
       ").\n", sep = "")

# 2. Package check (no install)
required <- c(
  # core
  "dplyr","tidyr","ggplot2","gridExtra","scales","reshape2",
  # statistical modelling
  "gamlss","gamlss.dist","VineCopula","scoringRules","psych",
  # plotting helpers
  "fmsb"
)
miss <- setdiff(required, rownames(installed.packages()))
if (length(miss)) {
  cat("\n  MISSING R packages:\n    ",
       paste(miss, collapse = ", "), "\n", sep = "")
  cat("  Install with:\n    install.packages(c(",
       paste0("\"", miss, "\"", collapse = ", "), "))\n", sep = "")
  fail(sprintf("install %d missing package(s) and re-run.", length(miss)))
} else {
  cat(sprintf("  All %d required CRAN packages installed.\n", length(required)))
}

# 3. Required input data
INPUT_CSV <- "data/analysis_dataset.csv"
if (!file.exists(INPUT_CSV)) {
  if (DO_PULL) {
    cat(sprintf("  %s not present; --from-pull set, will try to build it.\n",
                 INPUT_CSV))
  } else {
    cat("\n  MISSING ", INPUT_CSV, "\n", sep = "")
    cat("  Either:\n",
         "    (a) place the anonymised analysis CSV at that path",
         " (see data/README.md for the schema), or\n",
         "    (b) re-run with --from-pull and configure VALD_CLIENT_ID,",
         " VALD_CLIENT_SECRET, VALD_TENANT_ID env vars.\n", sep = "")
    fail("no input data, cannot proceed.")
  }
} else {
  sz <- file.info(INPUT_CSV)$size
  cat(sprintf("  Found %s (%.1f MB).\n", INPUT_CSV, sz / 1024 / 1024))
}

# =============================================================================
# Stage runner
# =============================================================================
timings <- list()

run_stage <- function(label, file, env_args = NULL) {
  banner(sprintf("%s  --  %s", label, file))
  full <- file.path(CODE_DIR, file)
  if (!file.exists(full)) fail(sprintf("script not found: %s", full))
  t0 <- Sys.time()
  env <- new.env()
  ok <- tryCatch({
    sys.source(full, envir = env)
    if (is.function(env$main)) {
      do.call(env$main, env_args %||% list())
    }
    TRUE
  }, error = function(e) {
    cat("\n  STAGE FAILED: ", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timings[[label]] <<- dt
  cat(sprintf("  -> %s in %s\n", if (ok) "completed" else "FAILED",
               human_dur(dt)))
  if (!ok) fail(sprintf("aborting at stage: %s", label))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# Pipeline
# =============================================================================
if (DO_PULL) {
  run_stage("Stage 0a", "01_pull_vald.R")
  run_stage("Stage 0b", "02_integrate_sport_gender.R")
  if (!file.exists(INPUT_CSV))
    fail(sprintf("after --from-pull, %s still missing.", INPUT_CSV))
}

run_stage("Stage 1 (analysis)", "RUN_PAPER_FINAL.R")
run_stage("Stage 2 (paper assets)", "MAKE_PAPER_ASSETS.R")

# =============================================================================
# Asset audit
# =============================================================================
if (!SKIP_AUDIT) {
  banner("Stage 3 (paper-asset audit)  --  99_check_paper_assets.R")
  audit_args <- if (NO_MTIME) c("--no-mtime-check") else character(0)
  t0 <- Sys.time()
  Sys.setenv(VALD_AUDIT_ARGS = paste(audit_args, collapse = " "))
  env <- new.env()
  ok <- tryCatch({
    sys.source(file.path(CODE_DIR, "99_check_paper_assets.R"), envir = env)
    isTRUE(env$.AUDIT_PASS)
  }, error = function(e) {
    cat("\n  AUDIT FAILED: ", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timings[["Stage 3 (audit)"]] <- dt
  cat(sprintf("  -> %s in %s\n", if (ok) "audit PASS" else "audit FAIL",
               human_dur(dt)))
  if (!ok) fail("paper-asset audit reported missing or stale assets.")
} else {
  cat("\n  --skip-audit set; not running 99_check_paper_assets.R\n")
}

# =============================================================================
# Final summary
# =============================================================================
banner("Pipeline finished", ch = "#")
total <- sum(unlist(timings))
for (lbl in names(timings))
  cat(sprintf("  %-30s  %s\n", lbl, human_dur(timings[[lbl]])))
cat(sprintf("  %-30s  %s\n", "TOTAL", human_dur(total)))
cat("\nNext step: rebuild the PDF with\n",
     "    cd paper && pdflatex -interaction=nonstopmode main.tex",
     " && bibtex main \\\n",
     "      && pdflatex -interaction=nonstopmode main.tex",
     " && pdflatex -interaction=nonstopmode main.tex\n", sep = "")

invisible(NULL)
