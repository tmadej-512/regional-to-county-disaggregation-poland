## =============================================================================
## 06c_calibrate_xl_prior.R — kalibracja prioru ściągania międzypoziomowego (c) na
## prawdziwych zmiennych powiatowych
##
## c = jaki ułamek zależności oszacowanej MIĘDZY województwami działa WEWNĄTRZ nich.
## Z danych wojewódzkich nie da się go oszacować, więc bierzemy go z danych zewnętrznych:
## 8 prawdziwych zmiennych powiatowych (06_real_validation.R) agregujemy do województw,
## dezagregujemy GDR_XL / GP_XL dla c z siatki i wybieramy c o najlepszym średnim zysku
## względem metody naiwnej. Wariant leave-one-target-out (LOTO) pokazuje, czy wybór jest
## stabilny, gdy zmienna celu nie bierze udziału w kalibracji.
## Symulacja (04) raportuje c ustalone z góry (XL_C0 = 0,5); aplikacja (07) używa c z tego kroku.
## Wynik: output/results/xl_c_calibration_grid.csv, xl_c_calibration_summary.csv,
##        xl_c_prior_loto.csv, xl_c_prior.csv
## Parametry: SAE_XL_CGRID (domyślnie "0,0.25,0.5,0.75,1"), SAE_NREP_XLCAL (domyślnie 20)
## =============================================================================
source("R/00_config.R"); source("R/02_methods.R"); source("R/03_scenarios.R")

C_GRID <- as.numeric(strsplit(Sys.getenv("SAE_XL_CGRID", "0,0.25,0.5,0.75,1"), ",")[[1]])
N_REP_XLCAL <- as.integer(Sys.getenv("SAE_NREP_XLCAL", "20"))

pow <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |> arrange(pow)
W   <- readRDS(file.path(PATH_PROC, "W_queen.rds"))
izp <- read_csv(file.path(PATH_PROC, "izp_covars_2019.csv"), show_col_types = FALSE,
                col_types = cols(pow = "c", .default = "d"))
pow <- pow |> left_join(izp, by = "pow")
f_extra <- file.path(PATH_PROC, "extra_covars.csv")
if (file.exists(f_extra))
  pow <- pow |> left_join(read_csv(f_extra, show_col_types = FALSE, col_types = cols(pow = "c")), by = "pow")
COVSETS <- make_covsets(pow)[c("public", "izp")]
TARGETS <- read_csv(file.path(PATH_RES, "real_targets.csv"), show_col_types = FALSE)
tau_lovo <- read_csv(file.path(PATH_RES, "real_tau_lovo.csv"), show_col_types = FALSE)

grid <- expand_grid(covset = names(COVSETS), target = TARGETS$target, rep = 0:N_REP_XLCAL)

one_run <- function(i) {
  cs <- grid$covset[i]; tg <- grid$target[i]; r <- grid$rep[i]
  drop <- TARGETS$drop[TARGETS$target == tg]
  if (identical(drop, "unemp19")) drop <- c(drop, "izp_2_2_3")
  P  <- make_problem(pow, W, covars = setdiff(COVSETS[[cs]], drop))
  y  <- pow[[tg]]; zt <- agg(P, y)
  if (r == 0) { z <- zt; se <- rep(0, P$m) } else {
    set.seed(SEED + 6e6 + 1e4 * match(tg, TARGETS$target) + 1e3 * match(cs, names(COVSETS)) + r)
    se <- survey_se(P, zt, 0.4); z <- zt + rnorm(P$m, 0, se)
  }
  tr <- tau_lovo$tau_ratio[tau_lovo$target == tg]
  base <- metrics_pow(P, y, z[P$g])
  bind_rows(lapply(C_GRID, function(cc) bind_rows(lapply(c("car", "gp"), function(sp) {
    f <- gdr_fit(P, z, se, tau_ratio = tr, spatial = sp, xl = list(c0 = cc, sd = XL_SD))
    metrics_pow(P, y, f$pred, f$sd) |>
      mutate(spatial = sp, c = cc, skill = 1 - RMSE^2 / base$RMSE^2)
  })))) |> mutate(covset = cs, target = tg, noise = ifelse(r == 0, "exact", "survey"), rep = r, .before = 1)
}

t0 <- Sys.time()
res <- par_lapply(seq_len(nrow(grid)), one_run)
message("Czas: ", format(Sys.time() - t0))
res <- bind_rows(res[!sapply(res, inherits, "try-error")])
write_csv(res, file.path(PATH_RES, "xl_c_calibration_grid.csv"))

## średnie po replikacjach, potem po celach i zestawach kowariat
cell <- res |> group_by(spatial, c, covset, target, noise) |>
  summarise(skill = mean(skill), cal_error = mean(abs(CalSlope - 1), na.rm = TRUE), .groups = "drop")
summ <- cell |> group_by(spatial, c) |>
  summarise(mean_skill = mean(skill), worst_skill = min(skill),
            mean_abs_cal_error = mean(cal_error), .groups = "drop")
write_csv(summ, file.path(PATH_RES, "xl_c_calibration_summary.csv"))

pick <- function(d) d |> arrange(desc(round(mean_skill, 3)), mean_abs_cal_error) |> slice(1)
prior <- summ |> group_by(spatial) |> group_modify(~ pick(.x)) |> ungroup()
write_csv(prior, file.path(PATH_RES, "xl_c_prior.csv"))

loto <- map_dfr(unique(cell$target), function(tg)
  cell |> filter(target != tg) |> group_by(spatial, c) |>
    summarise(mean_skill = mean(skill), mean_abs_cal_error = mean(cal_error), .groups = "drop") |>
    group_by(spatial) |> group_modify(~ pick(.x)) |> ungroup() |> mutate(heldout_target = tg))
write_csv(loto, file.path(PATH_RES, "xl_c_prior_loto.csv"))
print(summ |> mutate(across(where(is.numeric), ~ round(.x, 3))), n = 50)
print(prior)
