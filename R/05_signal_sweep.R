## =============================================================================
## 05_signal_sweep.R — wpływ siły związku kowariaty–zmienna („rozmyte korelacje”)
## y = sqrt(r2) f_lin + sqrt((1-r2)/2) u_CAR + sqrt((1-r2)/2) eps,  r2 ∈ [0.05, 0.9]
## Wynik: output/results/sweep_metrics.csv
## =============================================================================
source("R/00_config.R"); source("R/02_methods.R"); source("R/03_scenarios.R")

pow <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |> arrange(pow)
W   <- readRDS(file.path(PATH_PROC, "W_queen.rds"))
P   <- make_problem(pow, W)
K   <- make_components(P, pow)

R2_GRID <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
NOISE <- c(exact = 0, survey = 0.4)
grid <- expand_grid(r2 = R2_GRID, noise = names(NOISE), rep = seq_len(N_REP_SWEEP))

one_run <- function(k) {
  r2 <- grid$r2[k]; nz <- grid$noise[k]; r <- grid$rep[k]
  set.seed(SEED + 7e5 + 1e4 * match(r2, R2_GRID) + 1e3 * match(nz, names(NOISE)) + r)
  y  <- simulate_truth("SWEEP", P, K, r2 = r2)
  zt <- agg(P, y)
  se <- if (NOISE[[nz]] > 0) survey_se(P, zt, NOISE[[nz]]) else rep(0, P$m)
  z  <- zt + rnorm(P$m, 0, se)
  fits <- run_methods(P, z, se, y_true = y, seed = r)
  evaluate_all(P, y, fits) |> mutate(r2 = r2, noise = nz, rep = r, .before = 1)
}
t0 <- Sys.time()
res <- par_lapply(seq_len(nrow(grid)), one_run)
message("Czas: ", format(Sys.time() - t0))
res <- res[!sapply(res, inherits, "try-error")]
write_csv(bind_rows(res), file.path(PATH_RES, "sweep_metrics.csv"))
