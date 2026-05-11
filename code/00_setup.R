# ------------------------------------------------------------
# 00_setup.R
# Common configuration, paths, and dependency management.
# Sourced by every numbered script and by run_all.R.
# ------------------------------------------------------------

# `or` fallback operator (defined first so paths can use it) --
`%||%` <- function(a, b) if (is.null(a)) b else a

# Paths -------------------------------------------------------
.detect_project_dir <- function() {
  env_root <- Sys.getenv("VALD_PROJECT_DIR", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) return(normalizePath(env_root))
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[grep("^--file=", args)])
  this_file <- if (length(fa)) fa[1] else
    tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(this_file) && nzchar(this_file)) {
    return(normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE))
  }
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  if (basename(cwd) == "code") return(normalizePath(file.path(cwd, "..")))
  cwd
}

PROJECT_DIR <- .detect_project_dir()
CODE_DIR    <- file.path(PROJECT_DIR, "code")
DATA_DIR    <- Sys.getenv("VALD_DATA_DIR",
  unset = file.path(PROJECT_DIR, "data"))
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
PAPER_FIGS  <- file.path(PROJECT_DIR, "paper", "figures")
PAPER_TABS  <- file.path(PROJECT_DIR, "paper", "tables")
LOG_DIR     <- file.path(PROJECT_DIR, "logs")

for (d in c(DATA_DIR, RESULTS_DIR, PAPER_FIGS, PAPER_TABS, LOG_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Reproducibility --------------------------------------------
SEED <- 2026L
set.seed(SEED)

# Tunables ---------------------------------------------------
MIN_N_STRATUM <- 25L     # min observations per (sex × sport × test) cohort
MAX_METRICS   <- 4L      # multivariable dimension (top-k by RF importance)
N_SIM         <- 300L    # predictive draws per observation
CV_FOLDS      <- 5L
QUIET_GAMLSS  <- TRUE

# Plot theme -------------------------------------------------
suppressPackageStartupMessages({
  if (!"package:ggplot2" %in% search()) library(ggplot2)
})
theme_set(theme_minimal(base_size = 10))
PAL <- list(
  female = "#d7301f",
  male   = "#0570b0",
  band   = "#3690c0",
  median = "#08519c",
  ref    = "firebrick"
)

# Dependency loader ------------------------------------------
.required_pkgs <- c(
  "dplyr", "tidyr", "data.table", "ggplot2", "gridExtra", "scales",
  "gamlss", "gamlss.dist", "VineCopula",
  "randomForest", "glmnet",
  "scoringRules", "moments",
  "fmsb", "knitr"
)
.optional_pkgs <- c("valdr", "DepthProc", "writexl")

ensure_pkgs <- function(pkgs, repo = "https://cloud.r-project.org") {
  miss <- setdiff(pkgs, rownames(installed.packages()))
  if (length(miss) && Sys.getenv("VALD_INSTALL_MISSING", "FALSE") == "TRUE") {
    message("VALD_INSTALL_MISSING=TRUE -> installing: ",
            paste(miss, collapse = ", "))
    try(install.packages(miss, repos = repo, quiet = TRUE), silent = TRUE)
    miss <- setdiff(pkgs, rownames(installed.packages()))
  }
  if (length(miss)) {
    stop(sprintf("Missing required R packages: %s\n  ",
                  paste(miss, collapse = ", ")),
          "Install with:\n",
          sprintf("    install.packages(c(%s))",
                   paste0("\"", miss, "\"", collapse = ", ")),
          call. = FALSE)
  }
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  invisible(ok)
}

ensure_pkgs(.required_pkgs)

suppressPackageStartupMessages({
  for (p in .required_pkgs) suppressWarnings(library(p, character.only = TRUE))
})
HAS_VALDR     <- requireNamespace("valdr", quietly = TRUE)
HAS_DEPTHPROC <- requireNamespace("DepthProc", quietly = TRUE)

# Logging helper --------------------------------------------
log_step <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
