# Quick start

1. Put the input files in `data/raw/` and `data/raw/izp/` (see `data/README_data.md`).
2. Open `regional-to-county-disaggregation-poland.Rproj` in RStudio, or set the working directory to the repository root.
3. Quick test run (a few minutes, few replicates):

```r
Sys.setenv(SAE_CORES = 4, SAE_NREP = 5, SAE_NREP_SWEEP = 3, SAE_NREP_REAL = 3,
           SAE_NREP_SENS = 3, SAE_NREP_MV = 3, SAE_NREP_ORD = 2, SAE_NREP_ORD_MV = 2,
           SAE_NREP_XLCAL = 2)
source("run_all.R")
```

4. Full run: start a fresh R session, set only `SAE_CORES`, and run `source("run_all.R")`. Run each step in a clean environment as `run_all.R` does (`source(s, local = new.env())`).
5. To redraw figures only: `source("R/08_figures_tables.R")`.

Useful settings (environment variables):

| Variable | Default | Meaning |
|---|---|---|
| `SAE_CORES` | cores − 1 | parallel workers |
| `SAE_NREP` | 200 | replicates per simulation cell |
| `SAE_XL_C0`, `SAE_XL_SD` | 0.5, 0.29 | cross-level prior used in simulation; the application uses the calibrated value from `06c` unless `SAE_XL_C0` is set |
| `APP_COVSET` | `izp` | covariate set in the application (`public`, `izp`, `core`, `wide`) |
| `APP_METHOD` | `auto` | main application method (`auto` = selected in `06b`) |
