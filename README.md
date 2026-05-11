# Multivariable Centile Modelling of Force-Plate Athlete Profiles

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-pending--Zenodo-lightgrey.svg)](#)

A sex-stratified, reliability-aware pipeline for VALD ForceDecks data —
the companion code release for the STAI-X 2026 paper
*"Multivariable Centile Modelling of Force-Plate Athlete Profiles"*
by Praveen D. Chougale and Usha Ananthakumar (IIT Bombay).

## Abstract

Force-plate countermovement, drop-jump and single-leg protocols routinely
return dozens of biomechanical variables per athlete, but practitioners
benchmark each variable against a separate centile band and combine the
slices in their head. We present a fully reproducible pipeline that turns
a federation-scale corpus of VALD ForceDecks tests (529 athletes, 3,739
sessions, 14 sports) into *multivariable* centile profiles per
sex × sport × test cohort. The pipeline rests on three deliberately
conservative ingredients: (i) a fixed clinical metric panel per protocol
so cohorts are comparable; (ii) a 5-fold cross-validated *model
competition* between Linear, GAMLSS-NO and GAMLSS-BCT marginals under
CRPS/IQR; and (iii) a copula competition between Independence, Gaussian
and regular vine joint models under the multivariate Energy Score.
Per-cohort dimensionality is set by a reliability filter (drop metrics
with ICC₍₃,₁₎ < 0.5) and a sample-size cap (n ≥ 10·C(p,2)). Athlete
percentiles use Mahalanobis depth on probability-integral-transformed
residuals with a bootstrap interval.

## Quickstart

The whole pipeline runs end-to-end via a single driver
(`code/00_run_all.R`) that calls each stage in order and finishes with a
self-check that audits every figure and LaTeX table referenced by the
paper.

### Path A — verify on existing data

For a reviewer or collaborator who already has the anonymised
`data/analysis_dataset.csv` from the authors. No VALD credentials
needed. Approximate end-to-end runtime on a recent laptop:
**~1 minute for the R pipeline + ~5 seconds for the LaTeX build.**

```bash
git clone https://github.com/praveenmaths89/vald-multivariable-centiles.git
cd vald-multivariable-centiles

# 1. Place the anonymised analysis CSV at data/analysis_dataset.csv
#    (schema documented in data/README.md). The repo does NOT ship data.

# 2. Run the entire R pipeline + paper-asset audit
Rscript code/00_run_all.R
#   Stage 1 (analysis)              ~48s
#   Stage 2 (paper assets)          ~15s
#   Stage 3 (audit)                 <1s
#   TOTAL                           ~1m 03s
#   ...
#   === Summary ===
#     Paper assets:    18/18 OK
#     Pipeline CSVs:   all present

# 3. Compile the paper (pdflatex + bibtex)
cd paper
pdflatex -interaction=nonstopmode main.tex
bibtex   main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

### Path B — fresh pull from VALD

For a user with VALD ForceDecks API credentials who wants to rebuild
`data/analysis_dataset.csv` from scratch. Requires `valdr` to be
installed and the following environment variables:

| Variable | Required? | Description |
|----------|-----------|-------------|
| `VALD_CLIENT_ID`     | yes  | OAuth client id for the VALD Hub API |
| `VALD_CLIENT_SECRET` | yes  | OAuth client secret |
| `VALD_TENANT_ID`     | yes  | Tenant id of the federation |
| `VALD_REGION`        | no   | Region (default `aue`) |
| `VALD_START_DATE`    | no   | ISO-8601 start of the pull window (default `2020-01-01T00:00:00Z`) |
| `VALD_SPORTS_DIR`    | no   | Directory of per-sport rosters (CSV files with `profile_id`, `gender`, `name`) |
| `VALD_INTEGRATED_CSV`| no   | Bypass the API entirely by pointing at a pre-integrated wide CSV |

Approximate runtime: **API pull is network-bound (1–10 minutes
depending on cohort size); the rest matches Path A.**

```bash
git clone https://github.com/praveenmaths89/vald-multivariable-centiles.git
cd vald-multivariable-centiles

# 1. Configure credentials (zsh / bash):
export VALD_CLIENT_ID=...
export VALD_CLIENT_SECRET=...
export VALD_TENANT_ID=...
export VALD_REGION=aue                  # optional
export VALD_SPORTS_DIR=/path/to/rosters # optional

# 2. Run the full pipeline including the API pull and integration
Rscript code/00_run_all.R --from-pull
#   Stage 0a (API pull)             1-10m  (network-bound)
#   Stage 0b (integration)          ~10s
#   Stage 1 (analysis)              ~48s
#   Stage 2 (paper assets)          ~15s
#   Stage 3 (audit)                 <1s

# 3. Compile the paper (same as Path A)
cd paper
pdflatex -interaction=nonstopmode main.tex
bibtex   main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

### Useful flags on `00_run_all.R`

- `--from-pull`        — also run the VALD API pull and integration
  stages (`01_pull_vald.R` then `02_integrate_sport_gender.R`) before
  the analysis pipeline.
- `--skip-audit`       — skip the paper-asset audit at the end.
- `--no-mtime-check`   — relax the audit's "modified within the last
  hour" freshness check (useful when you only re-built one stage).

The audit script can also be run on its own:

```bash
Rscript code/99_check_paper_assets.R                 # strict check
Rscript code/99_check_paper_assets.R --no-mtime-check # existence only
```

It exits with status 1 if any asset referenced by `paper/main.tex` is
missing, so it can be chained into CI.

## Methodology (six stages)

1. **Load and clean** the integrated ForceDecks dataset; aggregate
   per-trial measurements to a test-level value via the trial median;
   drop tests whose within-test coefficient of variation exceeds 25 %.
2. **Define cohorts** (sex × sport × test type), apply the n ≥ 25
   eligibility filter, and compute Royston's closed-form n_min for the
   97th centile at relative precision τ = 0.10.
3. **Reliability**: ICC₍₃,₁₎, SEM and MDC₉₅ per cohort × metric, used to
   screen which metrics enter the joint copula.
4. **Marginal model competition** per (cohort, panel-metric): Linear,
   GAMLSS-NO and GAMLSS-BCT scored by 5-fold cross-validated CRPS/IQR;
   the winner is re-fit on the full cohort.
5. **Joint copula competition** per cohort on the PIT residual matrix:
   Independence, Gaussian and regular vine (Dißmann's algorithm) scored
   by 5-fold cross-validated multivariate Energy Score.
6. **Athlete scoring**: per-metric percentiles from the winning marginals,
   joint percentile from the Mahalanobis depth on the PIT vector, and
   a non-parametric 95 % bootstrap interval (B = 200).

## Outputs

After a clean run the following CSV/RDS artefacts are produced under
`results/final/`:

| File | What it contains |
|------|------------------|
| `cohort_table.csv` | One row per cohort: n, panel size, joint dimensionality p★, binding constraint, and winning joint model |
| `cohort_fits.rds` | Full fitted marginal objects, joint copula objects, and PIT matrices per cohort |
| `marginal_competition.csv` | 5-fold CV CRPS/IQR for Linear, GAMLSS-NO and GAMLSS-BCT per (cohort, metric) and the winner |
| `joint_competition.csv` | 5-fold CV multivariate Energy Score for Indep, Gaussian and Vine per cohort, plus the gain over Independence |
| `reliability_table.csv` | Per (cohort × metric) ICC₍₃,₁₎, SEM and MDC₉₅ |
| `reliability_summary.csv` | Per-protocol medians of ICC₍₃,₁₎ and MDC₉₅ across cohorts |

The corresponding figures and LaTeX tables used in the paper live in
`paper/figures/` and `paper/tables/`, and are regenerated end-to-end by
`code/MAKE_PAPER_ASSETS.R`.

## Data Availability

The underlying ForceDecks athlete performance corpus used in this study
is **not redistributed** in this repository. The dataset is held by
**PhysioQinesis / Performance Lab** (Thane, Mumbai, India) in
collaboration with the **Indian Institute of Technology Bombay**.
Researchers wishing to access the data for replication or follow-on
work should contact the corresponding author **Praveen D. Chougale —
22d1629@iitb.ac.in**. Access is subject to the data sharing agreement
governing the collaboration between PhysioQinesis and IIT Bombay. See
[`data/README.md`](data/README.md) for the full data availability
statement.

The only per-athlete file shipped with this repository is
`verification/results/athlete_joint_scores.csv` — one row per
athlete-test occasion with columns
`athlete_id` (UUID), `sex`, `sport`, `test_type`, `age_years`,
`joint_pct`, `joint_band`, `joint_model`, `n_metrics_used` (plus
`joint_lo` and `joint_hi` when bootstrap CIs are available, giving
either 9 or 11 columns). Per-metric values and per-metric percentiles
are intentionally **not** exposed here. See [`MANIFEST.md`](MANIFEST.md)
for the full file index.

## Citation

Please cite this work using the metadata in [`CITATION.cff`](CITATION.cff)
or the conference proceedings entry once the paper is published.

## Acknowledgements

The authors gratefully acknowledge **PhysioQinesis / Performance Lab —
Physiotherapy and Sports Performance Centre** for collecting and
providing the ForceDecks data used in this study, through their ongoing
collaboration with the Indian Institute of Technology Bombay. We also
thank the coaches and athletes whose performance data made this work
possible.

## License

Released under the [MIT License](LICENSE) — copyright © 2026
Praveen D. Chougale, Usha Ananthakumar.
