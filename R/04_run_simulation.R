## =============================================================================
## 04_run_simulation.R — badanie symulacyjne: 7 scenariuszy × 2 poziomy szumu
## Wynik: output/results/sim_metrics.csv (długi format), output/results/sim_dgp_diagnostics.csv,
##        output/results/sim_examples.rds (przykładowe replikacje do map)
## =============================================================================
source("R/00_config.R"); source("R/02_methods.R"); source("R/03_scenarios.R")

pow <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |> arrange(pow)
W   <- readRDS(file.path(PATH_PROC, "W_queen.rds"))
P   <- make_problem(pow, W)
K   <- make_components(P, pow)
lw  <- spdep::mat2listw(P$W, style = "W")

NOISE <- c(exact = 0, survey = 0.4)

## Trzy warianty projektu (cele kalibracji z prawdziwych zmiennych powiatowych, 01_build_data.R):
##  realistic      – stosunek wariancji wewnątrz/między ≈ mediana (~2) ORAZ Q2 regionalne ≈ mediana:
##                   różnice regionalne częściowo wyjaśnione kowariatami, jak w danych – ANALIZA GŁÓWNA
##  regional_noise – stosunek ≈ 2, ale różnice regionalne niewyjaśnione (psi = 0) – wariant pesymistyczny
##  heterogeneous  – bez komponentu regionalnego (duża zmienność wewnątrz woj.) – wariant optymistyczny
ref <- read_csv(file.path(PATH_RES, "tab_variance_ratio_reference.csv"), show_col_types = FALSE)
RATIO_TARGET <- median(ref$ratio)
Q2_TARGET    <- median(ref$q2_regional)
message(sprintf("Cele kalibracji: stosunek wariancji %.2f, Q2 regionalne %.2f", RATIO_TARGET, Q2_TARGET))

CAL <- bind_rows(lapply(names(SCEN_LABELS), function(sc) {
  a <- calibrate_design(sc, P, K, RATIO_TARGET, Q2_TARGET)
  b <- calibrate_design(sc, P, K, RATIO_TARGET, NULL)
  bind_rows(tibble(design = "realistic", scen = sc, lambda = a[["lambda"]], psi = a[["psi"]]),
            tibble(design = "regional_noise", scen = sc, lambda = b[["lambda"]], psi = 0),
            tibble(design = "heterogeneous", scen = sc, lambda = 0, psi = 0))
}))
write_csv(CAL, file.path(PATH_RES, "sim_lambda_calibration.csv"))
print(CAL, n = 30)

## warianty poboczne liczymy tylko tam, gdzie różnią się od głównego
grid <- bind_rows(
  expand_grid(design = "realistic", scen = names(SCEN_LABELS)),
  CAL |> filter(design != "realistic", lambda > 0) |> select(design, scen)) |>
  expand_grid(noise = names(NOISE), rep = seq_len(N_REP_SIM)) |>
  mutate(id = row_number())

one_run <- function(k) {
  s <- grid$scen[k]; nz <- grid$noise[k]; r <- grid$rep[k]; dsg <- grid$design[k]
  cal <- CAL[CAL$design == dsg & CAL$scen == s, ]
  set.seed(SEED + 1e5 * match(s, names(SCEN_LABELS)) + 1e4 * match(nz, names(NOISE)) + r +
             1e7 * match(dsg, c("realistic", "regional_noise", "heterogeneous")))
  y  <- simulate_truth(s, P, K, lambda = cal$lambda, psi = cal$psi)
  zt <- agg(P, y)
  se <- if (NOISE[[nz]] > 0) survey_se(P, zt, NOISE[[nz]]) else rep(0, P$m)
  z  <- zt + rnorm(P$m, 0, se)
  fits <- run_methods(P, z, se, y_true = y, seed = r)
  met <- evaluate_all(P, y, fits) |>
    mutate(design = dsg, scen = s, noise = nz, rep = r, q2_regional = regional_q2(z, P$Xbar), .before = 1)
  ## diagnostyka DGP (dla noise = exact)
  diag <- NULL
  if (nz == "exact") {
    yb <- zt[P$g]
    diag <- tibble(design = dsg, scen = s, rep = r,
                   ratio_within_between = sum(P$w * (y - yb)^2) / P$m / var(zt),
                   q2_regional = regional_q2(zt, P$Xbar),
                   R2_oracle = cor(y, fits$ORACLE$pred)^2,
                   moran_I = spdep::moran(y, lw, P$n, spdep::Szero(lw))$I)
  }
  ex <- if (r == 1) list(design = dsg, scen = s, noise = nz, y = y, z = z, se = se,
                         preds = map(fits, "pred"), sd_gdr = fits$GDR$sd) else NULL
  list(met = met, diag = diag, ex = ex)
}

message(sprintf("Symulacja: %d przebiegów na %d rdzeniach...", nrow(grid), N_CORES))
t0 <- Sys.time()
res <- par_lapply(seq_len(nrow(grid)), one_run)
message("Czas: ", format(Sys.time() - t0))
bad <- sapply(res, inherits, "try-error")
if (any(bad)) { warning(sum(bad), " przebiegów z błędem"); print(res[bad][1]) }
res <- res[!bad]

sim_metrics <- bind_rows(map(res, "met"))
sim_diag    <- bind_rows(map(res, "diag"))
sim_ex      <- compact(map(res, "ex"))

write_csv(sim_metrics, file.path(PATH_RES, "sim_metrics.csv"))
write_csv(sim_diag,    file.path(PATH_RES, "sim_dgp_diagnostics.csv"))
saveRDS(sim_ex,        file.path(PATH_RES, "sim_examples.rds"))

sim_diag |> group_by(design, scen) |> summarise(across(-rep, mean)) |> print()
