## =============================================================================
## 07b_ordinal_application.R — dezagregacja porządkowa wskaźników EHIS 2019 (x1–x5)
## Zamiast samego odsetka (np. „złe + bardzo złe”) używa całego rozkładu odpowiedzi na pytanie
## (02b_ordinal.R). x6–x8 mają tylko dwie istotne kategorie (+ „nie potrzebowali”), więc zostają
## binarne; x9–x11 nie mają rozkładów w tablicy SE.
## Wynik: output/results/app_health_ordinal.xlsx (estymacje, porównanie z wersją binarną, wyniki eta,
##        dopasowanie modelu progów), output/results/app_health_ordinal_long.rds
## =============================================================================
source("R/00_config.R"); source("R/02_methods.R"); source("R/02b_ordinal.R")

## definicje pytań (ORD_DEFS) w 02b_ordinal.R: x1 złe+bardzo złe, x2 nadwaga+otyłość,
## x3 warzywa rzadziej niż codziennie, x4 alkohol co najmniej raz w tygodniu, x5 palenie codziennie
ORD_METHODS <- c("PYCNO", "GDR", "GDR_XL", "GP_GAM")

pow <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |> arrange(pow)
W   <- readRDS(file.path(PATH_PROC, "W_queen.rds"))
izp <- read_csv(file.path(PATH_PROC, "izp_covars_2019.csv"), show_col_types = FALSE,
                col_types = cols(pow = "c", .default = "d"))
pow <- pow |> left_join(izp, by = "pow")
f_extra <- file.path(PATH_PROC, "extra_covars.csv")
if (file.exists(f_extra))
  pow <- pow |> left_join(read_csv(f_extra, show_col_types = FALSE, col_types = cols(pow = "c")), by = "pow")
APP_COVSET <- Sys.getenv("APP_COVSET", "izp")
SETS <- make_covsets(pow)
if (!APP_COVSET %in% names(SETS)) APP_COVSET <- "izp"
P <- make_problem(pow, W, covars = SETS[[APP_COVSET]])
woj_levels <- levels(factor(pow$woj))

cats <- read_csv(file.path(PATH_PROC, "ehis_categories_2019.csv"), show_col_types = FALSE,
                 col_types = cols(woj = "c"))

res <- list(); scores <- list(); thr <- list(); Z <- NULL; S <- NULL
for (j in seq_len(nrow(ORD_DEFS))) {
  tg <- ORD_DEFS$target[j]; num <- ORD_DEFS$num[[j]]
  cj <- cats |> filter(str_detect(variable, ORD_DEFS$var_pattern[j]))
  Pr <- cj |> select(woj, k, pct) |> pivot_wider(names_from = k, values_from = pct) |>
    arrange(match(woj, woj_levels))
  stopifnot(identical(Pr$woj, woj_levels))
  Pr <- as.matrix(Pr[, -1]) / 100
  neff <- cj |> distinct(woj, neff) |> arrange(match(woj, woj_levels)) |> pull(neff)
  od <- ord_disagg(P, Pr, neff, num, methods = ORD_METHODS, seed = 1)
  scores[[tg]] <- tibble(target = tg, woj = woj_levels, eta = od$scores$eta, se_eta = od$scores$se)
  thr[[tg]] <- tibble(target = tg, k = seq_along(od$scores$theta), theta = od$scores$theta,
                      r2_threshold = od$scores$r2_threshold,
                      categories = paste(unique(cj$category), collapse = " | "))
  Z <- cbind(Z, od$scores$eta); S <- cbind(S, od$scores$se)
  for (f in od$fits) res[[length(res) + 1]] <- tibble(
    !!!as.list(pow[c("pow", "woj", "nazwa", "w")]), target = tg, method = paste0("ORD_", f$method),
    est = f$est, lo90 = f$lo90, hi90 = f$hi90)
  message("Gotowe: ", tg, "  (R2 modelu progów = ", round(od$scores$r2_threshold, 3), ")")
}
colnames(Z) <- colnames(S) <- ORD_DEFS$target

## wszystkie wyniki porządkowe naraz: MGDR na eta (x1–x5 mają wspólne pole ukryte)
mv <- mgdr_fit(P, Z, S)
for (j in seq_len(nrow(ORD_DEFS))) {
  tg <- ORD_DEFS$target[j]; num <- ORD_DEFS$num[[j]]
  cj <- cats |> filter(str_detect(variable, ORD_DEFS$var_pattern[j]))
  Pr <- cj |> select(woj, k, pct) |> pivot_wider(names_from = k, values_from = pct) |>
    arrange(match(woj, woj_levels))
  Pr <- as.matrix(Pr[, -1]) / 100
  neff <- cj |> distinct(woj, neff) |> arrange(match(woj, woj_levels)) |> pull(neff)
  od <- ord_disagg(P, Pr, neff, num, fits = list(MGDR = list(pred = mv$pred[, j], sd = mv$sd[, j])))
  f <- od$fits$MGDR
  res[[length(res) + 1]] <- tibble(!!!as.list(pow[c("pow", "woj", "nazwa", "w")]), target = tg,
                                   method = "ORD_MGDR", est = f$est, lo90 = f$lo90, hi90 = f$hi90)
}
ord <- bind_rows(res)

## porównanie z wersją binarną (07_application_health.R), jeśli policzona
cmp <- tibble()
f_bin <- file.path(PATH_RES, "app_health_long.rds")
if (file.exists(f_bin)) {
  bin <- readRDS(f_bin) |> filter(target %in% ORD_DEFS$target) |>
    select(pow, target, method, est_bin = est, lo_bin = lo90, hi_bin = hi90)
  wd <- function(v, woj, w) v - ave(v * w, woj, FUN = sum) / ave(w, woj, FUN = sum)
  cmp <- ord |> mutate(base = sub("^ORD_", "", method)) |>
    inner_join(bin, by = c("pow", "target", "base" = "method")) |>
    group_by(target, method) |>
    summarise(cor_with_binary = cor(est, est_bin),
              cor_within_with_binary = cor(wd(est, woj, w), wd(est_bin, woj, w)),
              sd_within_ord = sd(wd(est, woj, w)), sd_within_bin = sd(wd(est_bin, woj, w)),
              width_ord = mean(hi90 - lo90), width_bin = mean(hi_bin - lo_bin), .groups = "drop")
  print(cmp |> mutate(across(where(is.numeric), ~ round(.x, 2))), n = 50)
}

write_xlsx(list(estymacje = ord, porownanie_z_binarna = cmp, wyniki_eta = bind_rows(scores),
                progi = bind_rows(thr),
                mgdr_ladunki = tibble(target = ORD_DEFS$target, loading = mv$loadings,
                                      common_share = mv$common_share)),
           file.path(PATH_RES, "app_health_ordinal.xlsx"))
saveRDS(ord, file.path(PATH_RES, "app_health_ordinal_long.rds"))
