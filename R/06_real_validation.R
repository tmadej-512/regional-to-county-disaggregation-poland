## =============================================================================
## 06_real_validation.R — walidacja na danych rzeczywistych („pseudo-populacja”)
## Prawdziwe zmienne powiatowe (PKW, BDL) agregujemy do województw (ważąc ludnością),
## dezagregujemy każdą metodą i porównujemy z prawdą powiatową.
##  (a) wariant dokładny: z = A y
##  (b) wariant sondażowy: z = A y + e, e ~ N(0, se^2) jak w EHIS, N_REP_REAL powtórzeń
## Wynik: output/results/real_metrics.csv, output/results/real_examples.rds
## =============================================================================
source("R/00_config.R"); source("R/02_methods.R"); source("R/03_scenarios.R")

pow <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |> arrange(pow)
W   <- readRDS(file.path(PATH_PROC, "W_queen.rds"))
izp <- read_csv(file.path(PATH_PROC, "izp_covars_2019.csv"), show_col_types = FALSE,
                col_types = cols(pow = "c", .default = "d"))
pow <- pow |> left_join(izp, by = "pow")
f_extra <- file.path(PATH_PROC, "extra_covars.csv")
if (file.exists(f_extra))
  pow <- pow |> left_join(read_csv(f_extra, show_col_types = FALSE, col_types = cols(pow = "c")), by = "pow") else
  message("Brak extra_covars.csv (01c_build_extra_covars.R) – bez zestawów core/wide z demografią.")

## zestawy kowariat (definicje: make_covsets() w 00_config.R): public, izp, core (krótka lista
## ustalona z góry), wide (wszystkie dostępne – test długiej listy ze ściąganiem współczynników)
COVSETS <- make_covsets(pow)
print(lengths(COVSETS))

## zmienne docelowe i kowariaty (bez opóźnionej wersji tej samej zmiennej)
REAL_TARGETS <- tribble(
  ~target,     ~label,                                  ~drop,
  "unemp23",   "Stopa bezrobocia 2023 (%)",              "unemp19",
  "turnout23", "Frekwencja – Sejm 2023 (%)",             "turnout19",
  "pis23",     "Poparcie PiS – Sejm 2023 (%)",           NA,
  "ko23",      "Poparcie KO – Sejm 2023 (%)",            NA,
  "konf23",    "Poparcie Konfederacji – Sejm 2023 (%)",  NA,
  "td23",      "Poparcie Trzeciej Drogi – Sejm 2023 (%)",NA,
  "trz25",     "Trzaskowski – II tura 2025 (%)",         NA,
  "right25",   "Mentzen + Braun – I tura 2025 (%)",      NA
)
write_csv(REAL_TARGETS, file.path(PATH_RES, "real_targets.csv"))

## Prior GDR bez błędnego koła: dla każdej zmiennej docelowej stosunek wariancji
## wewnątrz/między liczony jako mediana z POZOSTAŁYCH zmiennych referencyjnych, z pominięciem
## całej „rodziny” celu (np. dla pis23 pomijamy pis19 i pis23) – leave-one-variable-out.
ratio_tab <- read_csv(file.path(PATH_RES, "tab_variance_ratio_reference.csv"), show_col_types = FALSE)
family <- function(v) sub("\\d+$", "", v)
tau_lovo <- sapply(REAL_TARGETS$target, function(tg)
  median(ratio_tab$ratio[family(ratio_tab$variable) != family(tg)], na.rm = TRUE))
write_csv(tibble(target = names(tau_lovo), tau_ratio = tau_lovo), file.path(PATH_RES, "real_tau_lovo.csv"))
print(round(tau_lovo, 2))

all_met <- list(); examples <- list()
for (cs in names(COVSETS)) for (j in seq_len(nrow(REAL_TARGETS))) {
  tg <- REAL_TARGETS$target[j]
  ## bez opóźnionej wersji zmiennej docelowej (dla bezrobocia także wskaźnik IZP 2_2_3)
  drop <- REAL_TARGETS$drop[j]
  if (identical(drop, "unemp19")) drop <- c(drop, "izp_2_2_3")
  cv <- setdiff(COVSETS[[cs]], drop)
  P  <- make_problem(pow, W, covars = cv)
  y  <- pow[[tg]]
  zt <- agg(P, y)

  ## (a) dokładne
  set.seed(SEED + j + 100 * match(cs, names(COVSETS)))
  tr <- tau_lovo[[tg]]
  f <- run_methods(P, zt, rep(0, P$m), y_true = y, seed = 1, tau_ratio = tr)
  all_met[[length(all_met) + 1]] <- evaluate_all(P, y, f) |>
    mutate(covset = cs, target = tg, noise = "exact", rep = 0L, q2_regional = regional_q2(zt, P$Xbar),
           .before = 1)
  examples[[paste(cs, tg, sep = ":")]] <- list(y = y, z = zt, preds = map(f, "pred"), sd_gdr = f$GDR$sd)

  ## (b) z szumem sondażowym
  se <- survey_se(P, zt, 0.4)
  rr <- par_lapply(seq_len(N_REP_REAL), function(r) {
    set.seed(SEED + 1000 * j + 50000 * match(cs, names(COVSETS)) + r)
    z <- zt + rnorm(P$m, 0, se)
    ff <- run_methods(P, z, se, y_true = y, seed = r, tau_ratio = tr)
    evaluate_all(P, y, ff) |> mutate(covset = cs, target = tg, noise = "survey", rep = r,
                                     q2_regional = regional_q2(z, P$Xbar), .before = 1)
  })
  all_met[[length(all_met) + 1]] <- bind_rows(rr[!sapply(rr, inherits, "try-error")])
  message("Gotowe: ", cs, " / ", tg)
}
real_metrics <- bind_rows(all_met)
write_csv(real_metrics, file.path(PATH_RES, "real_metrics.csv"))
saveRDS(examples, file.path(PATH_RES, "real_examples.rds"))

## =============================================================================
## Wielowymiarowo: wszystkie 8 zmiennych naraz (MGDR) vs 8 osobnych GDR, ten sam zestaw kowariat
## (bez opóźnionych wersji celów: unemp19, turnout19, izp_2_2_3), wartości dokładne i z błędem.
## Wynik: output/results/real_mv_metrics.csv
## =============================================================================
mv_met <- list()
TG <- REAL_TARGETS$target
for (cs in names(COVSETS)) {
  cv <- setdiff(COVSETS[[cs]], c("unemp19", "turnout19", "izp_2_2_3"))
  P  <- make_problem(pow, W, covars = cv)
  Y  <- as.matrix(pow[, TG]); Zt <- apply(Y, 2, function(y) agg(P, y))
  SEt <- apply(Zt, 2, function(z) survey_se(P, z, 0.4))
  run_mv <- function(Z, SE, r, nz) {
    mv <- mgdr_fit(P, Z, SE, tau_ratio = median(tau_lovo))
    bind_rows(lapply(seq_along(TG), function(k) {
      g <- gdr_fit(P, Z[, k], SE[, k], tau_ratio = tau_lovo[[TG[k]]])
      bind_rows(metrics_pow(P, Y[, k], Z[P$g, k]) |> mutate(method = "NAIVE"),
                metrics_pow(P, Y[, k], g$pred, g$sd) |> mutate(method = "GDR"),
                metrics_pow(P, Y[, k], mv$pred[, k], mv$sd[, k]) |> mutate(method = "MGDR")) |>
        mutate(target = TG[k], loading = mv$loadings[[k]], common_share = mv$common_share[[k]])
    })) |> mutate(covset = cs, noise = nz, rep = r, .before = 1)
  }
  mv_met[[length(mv_met) + 1]] <- run_mv(Zt, matrix(0, P$m, length(TG)), 0L, "exact")
  rr <- par_lapply(seq_len(N_REP_REAL), function(r) {
    set.seed(SEED + 3e6 + 50000 * match(cs, names(COVSETS)) + r)
    run_mv(Zt + matrix(rnorm(length(Zt)), P$m) * SEt, SEt, r, "survey")
  })
  mv_met[[length(mv_met) + 1]] <- bind_rows(rr[!sapply(rr, inherits, "try-error")])
  message("Gotowe (wielowymiarowo): ", cs)
}
write_csv(bind_rows(mv_met), file.path(PATH_RES, "real_mv_metrics.csv"))

real_metrics |> filter(noise == "exact") |>
  group_by(covset, method) |> summarise(RMSE = mean(RMSE), CorWithin = mean(CorWithin)) |>
  print(n = 30)
