## =============================================================================
## 05c_multivariate_sim.R — czy wspólna dezagregacja K wskaźników (MGDR) pomaga?
## DGP (K = 6 wskaźników, prawdziwa geografia i kowariaty):
##   y_k = sqrt(0.3) f_k(X) + sqrt(0.45) [ sqrt(c) l_k s + sqrt(1 - c) u_k ] + sqrt(0.25) r_k
##     f_k(X) – liniowy sygnał kowariat z losowymi współczynnikami (inny dla każdego wskaźnika)
##     s      – WSPÓLNE ukryte pole: gradient miejski (S2) + pole CAR (po połowie wariancji)
##     u_k    – pole specyficzne (CAR), r_k – efekt województwa (różnice regionalne)
##     l_k    – znak ładunku (±1): wskaźniki „korzystne” i „niekorzystne”
##     c      – siła wspólnego czynnika: 0 (brak), 0.5, 0.8
## Porównanie: NAIVE, GDR osobno dla każdego wskaźnika, MGDR łącznie; wartości woj. dokładne
## i z błędem sondażowym. Wynik: output/results/mv_metrics.csv
## =============================================================================
source("R/00_config.R"); source("R/02_methods.R"); source("R/03_scenarios.R")

N_REP_MV <- as.integer(Sys.getenv("SAE_NREP_MV", "100"))
K_MV <- 6
pow <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |> arrange(pow)
W   <- readRDS(file.path(PATH_PROC, "W_queen.rds"))
P   <- make_problem(pow, W)
K   <- make_components(P, pow)

grid <- expand_grid(c_common = c(0, 0.5, 0.8), noise = c("exact", "survey"), rep = seq_len(N_REP_MV))

one_run <- function(i) {
  cc <- grid$c_common[i]; nz <- grid$noise[i]; r <- grid$rep[i]
  set.seed(SEED + 8e6 + 1e5 * match(cc, c(0, 0.5, 0.8)) + 1e4 * (nz == "survey") + r)
  s <- zstd(sqrt(0.5) * K$f_dist + sqrt(0.5) * draw_car(P))
  sgn <- sample(c(-1, 1), K_MV, replace = TRUE)
  Y <- sapply(seq_len(K_MV), function(k) {
    fx <- zstd(as.numeric(P$X %*% rnorm(ncol(P$X))))
    lat <- sqrt(cc) * sgn[k] * s + sqrt(1 - cc) * draw_car(P)
    50 + 5 * zstd(sqrt(0.3) * fx + sqrt(0.45) * zstd(lat) + sqrt(0.25) * draw_regional(P))
  })
  colnames(Y) <- paste0("y", seq_len(K_MV))
  Zt <- apply(Y, 2, function(y) agg(P, y))
  SE <- if (nz == "survey") apply(Zt, 2, function(z) survey_se(P, z, 0.4)) else matrix(0, P$m, K_MV)
  Z  <- Zt + matrix(rnorm(length(Zt)), P$m) * SE
  mv <- mgdr_fit(P, Z, SE)
  bind_rows(lapply(seq_len(K_MV), function(k) bind_rows(
    metrics_pow(P, Y[, k], Z[P$g, k]) |> mutate(method = "NAIVE"),
    { g <- gdr_fit(P, Z[, k], SE[, k]); metrics_pow(P, Y[, k], g$pred, g$sd) |> mutate(method = "GDR") },
    metrics_pow(P, Y[, k], mv$pred[, k], mv$sd[, k]) |> mutate(method = "MGDR")) |>
      mutate(indicator = k))) |>
    mutate(c_common = cc, noise = nz, rep = r, .before = 1)
}

t0 <- Sys.time()
res <- par_lapply(seq_len(nrow(grid)), one_run)
message("Czas: ", format(Sys.time() - t0))
res <- res[!sapply(res, inherits, "try-error")]
write_csv(bind_rows(res), file.path(PATH_RES, "mv_metrics.csv"))
