# Data Availability

The underlying ForceDecks athlete performance corpus used in this study
is **not redistributed** in this repository.

The dataset is held by **PhysioQinesis / Performance Lab** (Thane,
Mumbai, India) in collaboration with the **Indian Institute of
Technology Bombay**. Researchers wishing to access the data for
replication or follow-on work should contact the corresponding author:

> Praveen D. Chougale — 22d1629@iitb.ac.in

Access is subject to the data sharing agreement governing the
collaboration between PhysioQinesis and IIT Bombay.

## What this repository does provide

- **Cohort-level results** under `results/final/` and `verification/`
  (sample sizes, model competition scores, reliability statistics).
- **A per-athlete summary table** at
  `verification/results/athlete_joint_scores.csv` containing only the
  joint percentile (and 95 % bootstrap CI where available) for each
  athlete-test occasion. Per-metric values and per-metric percentiles
  are not exposed in this repository.
- **All code** used to produce the paper (`code/`).
- **The compiled paper and its source** (`paper/`).

## Reproducing the analysis

Once data access has been granted, place the raw `analysis_dataset.csv`
under `data/` and run:

    Rscript code/anonymise_for_github.R   # strips PII for further sharing
    Rscript code/RUN_PAPER_FINAL.R        # ~5-15 minutes
    Rscript code/MAKE_PAPER_ASSETS.R      # ~1-3 minutes

To regenerate the per-athlete radar plots and the full per-metric
scoring tables, run:

    Rscript code/VERIFY_ALL.R             # ~10 minutes

These outputs are intentionally not committed to this repository.
