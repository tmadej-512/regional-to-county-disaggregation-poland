# Reproducibility notes

- **Seeds.** All random steps use fixed seeds derived from `SEED` in `R/00_config.R`.
- **Run time.** The full workflow takes several hours on four cores; the simulation (`04`), ordinal simulation (`05d`) and validation (`06`) dominate. Replicate counts are controlled by environment variables (see `docs/QUICK_START.md`).
- **Platform.** Tested on Windows (R 4.5, `future.apply`) and Linux (`mclapply`). Figures use Arial (`sans`) on Windows and Liberation Sans on Linux.
- **Package masking.** `R/00_config.R` binds `select`, `filter`, `mutate` and related verbs to dplyr/tidyr/purrr so that packages loaded earlier in the session (e.g. MASS) cannot mask them.
- **Scope of this release.** The results reported with the manuscript come from a run of 15 October 2026 of a working version that also included exploratory estimators (graph-Laplacian and fused regressions with covariates, stacking), additional block covariate sets and a reliability-rule calibration. These exploratory components did not change the selected method (GDR_XL) and are not part of this curated release.
- **Not independently rerun.** The curated workflow has not been rerun end-to-end from this release; minor numerical differences from the reported results are expected.
- **Outputs.** Everything is written to `output/` (ignored by git).
