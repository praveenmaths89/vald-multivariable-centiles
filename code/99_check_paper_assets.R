# =============================================================================
#                       99_check_paper_assets.R
# Parses paper/main.tex, extracts every \includegraphics{} and \input{tables/}
# reference, and verifies that the asset:
#   (a) exists on disk,
#   (b) was modified within the last MTIME_WINDOW seconds (so we know it was
#       just rebuilt by this run); pass --no-mtime-check to relax (b).
# Also walks results/final/*.csv and reports which CSVs the analysis writes.
#
# Exit code: 0 on success, 1 on any missing or stale asset.
# Sets `.AUDIT_PASS` in its calling environment so 00_run_all.R can read it.
#
# Usage
#   Rscript code/99_check_paper_assets.R
#   Rscript code/99_check_paper_assets.R --no-mtime-check
# =============================================================================

# Pick up flags either from the shell or via env var (set by 00_run_all.R)
.cli_args  <- commandArgs(trailingOnly = TRUE)
.env_args  <- strsplit(Sys.getenv("VALD_AUDIT_ARGS", unset = ""), "\\s+")[[1]]
.flags     <- c(.cli_args, .env_args)
NO_MTIME   <- "--no-mtime-check" %in% .flags

# ---- Repo root anchor -------------------------------------------------------
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

MTIME_WINDOW <- as.integer(Sys.getenv("VALD_AUDIT_MTIME_WINDOW",
                                       unset = "3600"))   # 1 hour

# ---- Parse paper/main.tex --------------------------------------------------
TEX <- "paper/main.tex"
if (!file.exists(TEX))
  stop(sprintf("Paper source not found: %s", TEX), call. = FALSE)
tex <- readLines(TEX, warn = FALSE)

# Strip simple % line comments (LaTeX comments) so we don't pick up commented refs.
tex <- sub("(?<!\\\\)%.*$", "", tex, perl = TRUE)

# All \includegraphics references, optional [opts]{path}.
inc_paths <- character()
inc_re <- "\\\\includegraphics\\s*(?:\\[[^]]*\\])?\\{([^}]+)\\}"
for (ln in tex) {
  m <- regmatches(ln, gregexpr(inc_re, ln, perl = TRUE))[[1]]
  for (full in m) {
    p <- sub(inc_re, "\\1", full, perl = TRUE)
    inc_paths <- c(inc_paths, p)
  }
}
inc_paths <- unique(inc_paths)

# All \input{...} references — keep only those that look like asset files
in_paths <- character()
in_re <- "\\\\input\\{([^}]+)\\}"
for (ln in tex) {
  m <- regmatches(ln, gregexpr(in_re, ln, perl = TRUE))[[1]]
  for (full in m) {
    p <- sub(in_re, "\\1", full, perl = TRUE)
    in_paths <- c(in_paths, p)
  }
}
in_paths <- unique(in_paths)

# Resolve includegraphics paths against \graphicspath{{figures/}} default
# (we don't parse multiple graphics paths; we just probe the obvious places).
resolve_fig <- function(p) {
  candidates <- c(
    file.path("paper", p),
    file.path("paper", "figures", p),
    file.path("paper", paste0(p, ".jpg")),
    file.path("paper", "figures", paste0(p, ".jpg")),
    file.path("paper", paste0(p, ".pdf")),
    file.path("paper", "figures", paste0(p, ".pdf")),
    file.path("paper", paste0(p, ".png")),
    file.path("paper", "figures", paste0(p, ".png"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[1] else file.path("paper", "figures", p)  # canonical guess
}
resolve_input <- function(p) {
  cand <- file.path("paper", p)
  if (file.exists(cand)) return(cand)
  cand2 <- file.path("paper", paste0(p, ".tex"))
  if (file.exists(cand2)) return(cand2)
  cand     # canonical, possibly missing
}

assets <- data.frame(
  ref  = c(inc_paths, in_paths),
  kind = c(rep("figure", length(inc_paths)),
            rep("table",  length(in_paths))),
  path = c(vapply(inc_paths, resolve_fig,   character(1)),
            vapply(in_paths,  resolve_input, character(1))),
  stringsAsFactors = FALSE
)
# Filter: ignore numbers.tex etc. that are not under tables/  -- include them
# only if they exist or are explicitly asset-like.

now <- Sys.time()
chk <- function(p) {
  if (!file.exists(p)) return(list(status = "MISSING",
                                     size_kb  = NA_real_,
                                     age_secs = NA_real_))
  inf <- file.info(p)
  age <- as.numeric(difftime(now, inf$mtime, units = "secs"))
  st  <- if (NO_MTIME || age <= MTIME_WINDOW) "OK" else "STALE"
  list(status = st, size_kb = inf$size / 1024, age_secs = age)
}
ck <- lapply(assets$path, chk)
assets$status   <- vapply(ck, `[[`, character(1), "status")
assets$size_kb  <- vapply(ck, `[[`, numeric(1),   "size_kb")
assets$age_secs <- vapply(ck, `[[`, numeric(1),   "age_secs")
assets <- assets[order(assets$status, assets$kind, assets$ref), ]

# ---- Pretty print -----------------------------------------------------------
cat("\n=== Paper-asset audit ===\n")
if (NO_MTIME) cat("  (mtime check disabled by --no-mtime-check)\n")
fmt_age <- function(a) {
  if (is.na(a)) "       --"
  else if (a < 60) sprintf("%5.1fs", a)
  else if (a < 3600) sprintf("%4.1fm", a / 60)
  else if (a < 86400) sprintf("%4.1fh", a / 3600)
  else sprintf("%4.1fd", a / 86400)
}
fmt_sz  <- function(s) if (is.na(s)) "    --" else sprintf("%6.1f KB", s)

w_path <- max(nchar(assets$path), 36L)
fmt_status <- function(s) switch(s,
  "OK"      = "[OK]   ",
  "STALE"   = "[STALE]",
  "MISSING" = "[FAIL] ",
  paste0("[", s, "]"))
for (i in seq_len(nrow(assets))) {
  r <- assets[i, ]
  cat(sprintf("  %s  %-7s  %-*s  %s  age %s\n",
               fmt_status(r$status), r$kind, w_path, r$path,
               fmt_sz(r$size_kb), fmt_age(r$age_secs)))
}

# ---- Results-CSV inventory --------------------------------------------------
cat("\n=== Pipeline outputs (results/final) ===\n")
expected_csvs <- c(
  "results/final/cohort_table.csv",
  "results/final/cohort_fits.rds",
  "results/final/marginal_competition.csv",
  "results/final/joint_competition.csv",
  "results/final/reliability_table.csv",
  "results/final/reliability_summary.csv"
)
res_status_ok <- TRUE
for (p in expected_csvs) {
  if (file.exists(p)) {
    inf <- file.info(p)
    cat(sprintf("  [OK]    %-44s %.1f KB  age %s\n",
                 p, inf$size / 1024,
                 fmt_age(as.numeric(difftime(now, inf$mtime, units = "secs")))))
  } else {
    res_status_ok <- FALSE
    cat(sprintf("  [FAIL]  %-44s MISSING\n", p))
  }
}

# ---- Summary + exit ---------------------------------------------------------
n_total   <- nrow(assets)
n_ok      <- sum(assets$status == "OK")
n_stale   <- sum(assets$status == "STALE")
n_missing <- sum(assets$status == "MISSING")

cat("\n=== Summary ===\n")
cat(sprintf("  Paper assets:    %d/%d OK", n_ok, n_total))
if (n_stale)   cat(sprintf(", %d STALE", n_stale))
if (n_missing) cat(sprintf(", %d MISSING", n_missing))
cat("\n")
cat(sprintf("  Pipeline CSVs:   %s\n",
             if (res_status_ok) "all present" else "INCOMPLETE"))

.AUDIT_PASS <- (n_missing == 0L) && (n_stale == 0L) && res_status_ok
if (sys.nframe() == 0) {
  if (!.AUDIT_PASS) quit(status = 1L) else quit(status = 0L)
}
invisible(.AUDIT_PASS)
