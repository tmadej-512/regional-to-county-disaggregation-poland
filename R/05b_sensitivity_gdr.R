## =============================================================================
## 05b_sensitivity_gdr.R — wrażliwość GDR, MRF-GAM i GP-GAM na ustawienia przyjęte z góry
##   rho       – autokorelacja przestrzenna efektu CAR/MRF            {0.80, 0.95, 0.99}
##   tau_ratio – położenie prioru wariancji wewnątrzwojewódzkiej      {1, 2, 4}
##   tau_sd    – szerokość tego prioru (SD na skali log)              {0.5, 1, 2}
##   c0        – GDR_XL: średnia prioru ułamka przenoszonej zależności {0.25, 0.5, 0.75}
##               (c0 = 1 przy SD = 0 to zwykły GDR – konfiguracja „default” metody GDR)
## Schemat „jeden czynnik naraz” wokół ustawień domyślnych (0.95, 2, 1) = 7 konfiguracji.
## Wszystkie konfiguracje liczone na TYCH SAMYCH danych w każdej replikacji (wspólne liczby losowe).
## Projekt: scenariusze w wariancie realistycznym, wartości woj. z błędem sondażowym.
## Wynik: output/results/sens_metrics.csv
## =============================================================================
source("R/00_config.R"); source("R/02_methods.R"); source("R/03_scenarios.R")

N_REP_SENS <- as.integer(Sys.getenv("SAE_NREP_SENS", "100"))

pow <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |> arrange(pow)
W   <- readRDS(file.path(PATH_PROC, "W_queen.rds"))

RHOS <- c(0.80, 0.95, 0.99)
PROB <- setNames(lapply(RHOS, function(r) make_problem(pow, W, rho = r)), as.character(RHOS))
P0   <- PROB[["0.95"]]
K    <- make_components(P0, pow)

## lambda i psi jak w analizie głównej (04_run_simulation.R)
lam_file <- file.path(PATH_RES, "sim_lambda_calibration.csv")
stopifnot(file.exists(lam_file))   # najpierw 04_run_simulation.R
CAL <- read_csv(lam_file, show_col_types = FALSE) |> filter(design == "realistic")

CONFIGS <- tribble(
  ~config,          ~rho, ~tau_ratio, ~tau_sd,
  "default",        0.95, 2,          1,
  "rho = 0.80",     0.80, 2,          1,
  "rho = 0.99",     0.99, 2,          1,
  "tau ratio = 1",  0.95, 1,          1,
  "tau ratio = 4",  0.95, 4,          1,
  "prior SD = 0.5", 0.95, 2,          0.5,
  "prior SD = 2",   0.95, 2,          2
)

grid <- expand_grid(scen = names(SCEN_LABELS), rep = seq_len(N_REP_SENS))

one_run <- function(k) {
  s <- grid$scen[k]; r <- grid$rep[k]
  set.seed(SEED + 9e6 + 1e5 * match(s, names(SCEN_LABELS)) + r)
  cal <- CAL[CAL$scen == s, ]
  y  <- simulate_truth(s, P0, K, lambda = cal$lambda, psi = cal$psi)
  zt <- agg(P0, y)
  se <- survey_se(P0, zt, 0.4)
  z  <- zt + rnorm(P0$m, 0, se)
  base <- metrics_pow(P0, y, z[P0$g]) |> mutate(method = "NAIVE", config = "default")
  out <- list(base)
  for (i in seq_len(nrow(CONFIGS))) {
    cf <- CONFIGS[i, ]
    P  <- PROB[[as.character(cf$rho)]]
    for (mm in c("GDR", "GDR_GAM", "GP_GAM")) {
      if (mm == "GP_GAM" && cf$rho != 0.95) next          # rho dotyczy tylko CAR/MRF
      f <- gdr_fit(P, z, se, tau_ratio = cf$tau_ratio, tau_sd = cf$tau_sd,
                   smooth = mm != "GDR", spatial = if (mm == "GP_GAM") "gp" else "car")
      out[[length(out) + 1]] <- metrics_pow(P, y, f$pred, f$sd) |>
        mutate(method = mm, config = cf$config)
    }
  }
  ## ściąganie międzypoziomowe: GDR_XL przy różnych średnich prioru c (SD jak domyślnie)
  for (c0 in c(0.25, 0.5, 0.75)) {
    f <- gdr_fit(P0, z, se, xl = list(c0 = c0, sd = XL_SD))
    out[[length(out) + 1]] <- metrics_pow(P0, y, f$pred, f$sd) |>
      mutate(method = "GDR_XL", config = paste0("c0 = ", c0))
  }
  bind_rows(out) |> mutate(scen = s, rep = r, .before = 1)
}

t0 <- Sys.time()
res <- par_lapply(seq_len(nrow(grid)), one_run)
message("Czas: ", format(Sys.time() - t0))
res <- res[!sapply(res, inherits, "try-error")]
write_csv(bind_rows(res), file.path(PATH_RES, "sens_metrics.csv"))
