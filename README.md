# Disaggregating regional survey indicators to Polish counties

Reproducible R workflow accompanying the study *Disaggregating regional survey indicators to counties: a simulation study and empirical validation of small area methods in Poland*.

## Overview

The European Health Interview Survey (EHIS 2019) publishes health indicators for the 16 Polish voivodeships (NUTS-2). The Health Index for Poland (Indeks Zdrowia Polski, IZP) needs them for 380 counties (powiaty, LAU-1). With 16 regional values and 380 unknowns, every disaggregation method adds an assumption about how counties differ within a region. The study asks which methods recover county-level variation, how to validate them, and whether the resulting county estimates are fit for use in a composite index.

Fifteen estimators are compared: naive allocation (regional value), Tobler's pycnophylactic interpolation, macro-first ridge regression, mini-first bag-labelling learners (linear model, random forest), a Fay–Herriot synthetic estimator, Bayesian Gaussian disaggregation regression with a CAR effect (GDR), aggregated-likelihood GAMs with a Markov random field or Gaussian-process spatial effect (MRF-GAM, GP-GAM), **cross-level shrinkage** variants (GDR_XL, GP_XL), a multivariate extension with a shared latent field (MGDR), an ordinal (cumulative probit) variant for questions with ordered answers, and an infeasible oracle.

Evaluation uses three layers:

1. **Simulation** on the real geography and covariates of all 380 counties: seven scenarios representing distinct spatial structures of Poland (baseline socio-economic gradient, urban–rural gradient decaying with distance from cities, historical partitions, local hotspots, weak correlations, ecological fallacy, east–west gradient), calibrated to the within- to between-region variance ratio and regional predictability of real county variables, with exact and survey-error regional values.
2. **Pseudo-population validation** on eight real county variables (registered unemployment, turnout and vote shares 2023–2025) aggregated to regions and disaggregated back, with several covariate sets.
3. **Application** to eleven EHIS 2019 indicators with 90% prediction intervals and an indicator-level reliability assessment (signal-to-noise, agreement between methods, sensitivity to the cross-level assumption, agreement with external county data).

## Main findings

- Re-aggregation to regional values is not a test of a disaggregation: every benchmarked estimator reproduces the regional values exactly, while county-level error differs up to two-fold (*regional-validation paradox*).
- Standard covariate models rank counties well but overstate within-region contrasts. Cross-level shrinkage (GDR_XL) corrects the amplitude and was selected by a pre-specified minimax-regret rule: it beat naive allocation in 93% of evaluation cells, with a mean MSE reduction of 26% and a worst case of −3%.
- No method recovers structure the covariates cannot see (historical partitions, hotspots, cross-border gradients). Under an ecological fallacy every covariate method is worse than naive allocation.
- For EHIS 2019 the conditions are poor: covariates rarely explain regional differences, county estimates vary less than their own uncertainty, and the covariates are themselves IZP components. The recommended main version of the index therefore uses regional values for EHIS components; GDR_XL estimates are a sensitivity analysis, not a measurement.

These are ecological estimates conditional on an untestable cross-level assumption; they should not be used to rank individual counties.

A plain-language report in Polish for non-technical readers is available in [`docs/raport_dezagregacja_EHIS.docx`](docs/raport_dezagregacja_EHIS.docx).

## Repository structure

```text
R/
  00_config.R                   Paths, parameters, covariate sets, labels
  01_build_data.R               County database: boundaries, population weights, covariates,
                                historical partitions, contiguity, variance-ratio reference
  01b_build_izp.R               EHIS 2019 regional indicators and standard errors, full answer
                                distributions, IZP county covariates, cross-level coherence
  01c_build_extra_covars.R      Age structure, education, census population density, migrants,
                                fertility, medical staff; external county proxies for validation
  02_methods.R                  All estimators and county-level performance measures
  02b_ordinal.R                 Ordinal (cumulative probit) disaggregation
  03_scenarios.R                Data-generating processes S1–S7, calibration, survey-error model
  04_run_simulation.R           Main simulation
  05_signal_sweep.R             Covariate signal-strength sweep
  05b_sensitivity_gdr.R         Prior sensitivity (rho, variance prior, cross-level c)
  05c_multivariate_sim.R        Multivariate disaggregation (MGDR) simulation
  05d_ordinal_sim.R             Ordinal vs binary disaggregation simulation
  06_real_validation.R          Pseudo-population validation on real county variables
  06b_method_selection.R        Pre-specified method selection (minimax regret)
  06c_calibrate_xl_prior.R      Calibration of the cross-level shrinkage prior on real variables
  07_application_health.R       EHIS 2019 county estimates, MGDR, reliability assessment
  07b_ordinal_application.R     Ordinal variant for EHIS questions with ordered answers
  08_figures_tables.R           Figures (Lancet style) and tables
run_all.R                       Full analysis workflow
data/
  README_data.md                Data inputs, provenance and access requirements
docs/
  QUICK_START.md                How to prepare and run the analysis
  PIPELINE_DECISIONS.md         Pre-specified decisions and how they constrain later stages
  REPRODUCIBILITY_NOTES.md      Scope, run time, settings and known limitations
  raport_dezagregacja_EHIS.docx Plain-language report in Polish (methods, results, recommendations)
```

Input data and generated outputs are not tracked in this repository.

## Reproducing the analysis

1. Obtain the input files listed in [`data/README_data.md`](data/README_data.md) and place them in `data/raw/` (IZP and EHIS files in `data/raw/izp/`).
2. Open `regional-to-county-disaggregation-poland.Rproj` or set the R working directory to the repository root.
3. Run:

```r
Sys.setenv(SAE_CORES = 4)   # number of cores
source("run_all.R")
```

`run_all.R` installs missing packages and runs all steps in order. The full workflow is computationally intensive (several hours on four cores). A quick test run with few replicates is described in [`docs/QUICK_START.md`](docs/QUICK_START.md).

## Expected outputs

```text
output/results/   CSV, XLSX and RDS results, including
                  sim_metrics.csv, real_metrics.csv, tab_method_selection.csv,
                  selected_method.txt, xl_c_prior.csv,
                  app_health_estimates.xlsx (county estimates with 90% intervals,
                  reliability sheet "ocena_dezagregacji")
output/figures/   Manuscript figures 1–12 and supplementary figures S1–S10 (PNG, 400 dpi, and PDF)
```

## Software

Analyses were conducted in R 4.5. Main packages: sf, spdep, Matrix, mgcv, glmnet, ranger, dplyr, tidyr, purrr, readr, readxl, writexl, stringr, ggplot2, patchwork, ggrepel, scales, future, future.apply. Parallel computation uses `parallel::mclapply` on Linux/macOS and `future.apply` on Windows.

## Data availability

Public inputs (county and regional boundaries, election results, registered unemployment) come from Statistics Poland, the National Electoral Commission and the national boundary register. Project inputs (IZP calculator, EHIS 2019 regional tables with standard errors) are described in [`data/README_data.md`](data/README_data.md). No individual-level data are used.

## Citation

Bibliographic details will be added after publication.

## License

Code: MIT (see `LICENSE`). Data rights and licensing are separate; see [`data/README_data.md`](data/README_data.md).
