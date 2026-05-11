# Strips PII (name, dob) from the analysis dataset for public GitHub release.
# profile_id (UUID) and test_date are kept; all measurements are kept unchanged.
suppressPackageStartupMessages({})

in_dir  <- "data"
out_dir <- "data_anon"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "per_test"), showWarnings = FALSE)

strip_pii <- function(d) {
  d$name <- NULL
  d$dob  <- NULL
  d
}

cat("Anonymising analysis_dataset.csv ...\n")
d <- read.csv(file.path(in_dir, "analysis_dataset.csv"), check.names = FALSE)
write.csv(strip_pii(d), file.path(out_dir, "analysis_dataset.csv"), row.names = FALSE)

per_test_files <- list.files(file.path(in_dir, "per_test"), pattern = "\\.csv$")
for (f in per_test_files) {
  cat("Anonymising per_test/", f, " ...\n", sep = "")
  d <- read.csv(file.path(in_dir, "per_test", f), check.names = FALSE)
  write.csv(strip_pii(d), file.path(out_dir, "per_test", f), row.names = FALSE)
}
cat("Done. Output:", out_dir, "\n")
