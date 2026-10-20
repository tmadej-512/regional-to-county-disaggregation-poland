## =============================================================================
## 07_application_health.R — aplikacja: dezagregacja wskaźników zdrowotnych (EHIS 2019)
## z województw do powiatów, z przedziałami niepewności (GDR).
##
## Wejście (domyślnie): data/raw/najlepsze_modele_predykcja_woj.xlsx — kolumny woj, target, truth
##   (wartości wojewódzkie x1..x11 z arkusza „Do dezagregacji v2”).
## Opcjonalnie:
##   HEALTH_TARGETS_FILE – własny plik z wartościami woj. (kolumny: woj, target, value[, se])
##   HEALTH_SE_FILE      – błędy standardowe EHIS (kolumny: woj, target, se)
##   HEALTH_COVARS_FILE  – własne kowariaty powiatowe (np. z „IZP kalkulator”), kolumny: pow + zmienne
## Wynik: output/results/app_health_*.xlsx, output/figures/fig_app_*.png
## =============================================================================
source("R/00_config.R"); source("R/02_methods.R"); source("R/03_scenarios.R")

HEALTH_TARGETS_FILE <- Sys.getenv("HEALTH_TARGETS_FILE", "")
HEALTH_SE_FILE      <- Sys.getenv("HEALTH_SE_FILE", "")
HEALTH_COVARS_FILE  <- Sys.getenv("HEALTH_COVARS_FILE", "")
## metoda rekomendowana: "auto" = wybrana w 06b_method_selection.R (minimax regret, kryterium
## ustalone z góry); można wymusić, np. Sys.setenv(APP_METHOD = "GDR")
APP_METHOD          <- Sys.getenv("APP_METHOD", "auto")
APP_SCALE           <- Sys.getenv("APP_SCALE", "logit")    # "logit" dla wskaźników w %, "linear"
APP_METHODS         <- c("NAIVE", "PYCNO", "MACRO_B", "MINI_RF_B", "FH_B", "GDR", "GDR_GAM", "GP_GAM",
                         "GDR_XL", "GP_XL")
if (APP_METHOD == "auto") {
  f_sel <- file.path(PATH_RES, "selected_method.txt")
  APP_METHOD <- if (file.exists(f_sel)) readLines(f_sel)[1] else "GDR"
  message("Metoda główna aplikacji: ", APP_METHOD)
}
APP_METHODS <- union(APP_METHODS, APP_METHOD)

## prior ściągania międzypoziomowego skalibrowany na prawdziwych zmiennych (06c), o ile
## nie ustawiono go ręcznie zmienną SAE_XL_C0
f_xl <- file.path(PATH_RES, "xl_c_prior.csv")
if (file.exists(f_xl) && !nzchar(Sys.getenv("SAE_XL_C0"))) {
  xlp <- read_csv(f_xl, show_col_types = FALSE)
  if (any(xlp$spatial == "car")) XL_C0 <- xlp$c[xlp$spatial == "car"][1]
  if (any(xlp$spatial == "gp"))  XL_C0_GP <- xlp$c[xlp$spatial == "gp"][1]
  message("Prior c (06c): CAR = ", XL_C0, ", GP = ", if (exists("XL_C0_GP")) XL_C0_GP else XL_C0)
}

pow <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |> arrange(pow)
W   <- readRDS(file.path(PATH_PROC, "W_queen.rds"))
## domyślnie: zmienne strukturalne + kowariaty IZP 2019 (01b_build_izp.R)
izp <- read_csv(file.path(PATH_PROC, "izp_covars_2019.csv"), show_col_types = FALSE,
                col_types = cols(pow = "c", .default = "d"))
pow <- pow |> left_join(izp, by = "pow")
f_extra <- file.path(PATH_PROC, "extra_covars.csv")
if (file.exists(f_extra))
  pow <- pow |> left_join(read_csv(f_extra, show_col_types = FALSE, col_types = cols(pow = "c")), by = "pow")
## zestaw kowariat aplikacji: APP_COVSET = "izp" (domyślnie; w walidacji lepszy niż krótka lista "core"),
## "core" (krótka lista ustalona z góry), "wide" (wszystkie), "public"
APP_COVSET <- Sys.getenv("APP_COVSET", "izp")
SETS <- make_covsets(pow)
if (!APP_COVSET %in% names(SETS)) {
  message("Zestaw '", APP_COVSET, "' niedostępny – używam 'izp'."); APP_COVSET <- "izp"
}
covars <- SETS[[APP_COVSET]]
message("Kowariaty aplikacji (", APP_COVSET, "): ", paste(covars, collapse = ", "))

## --- opcjonalnie: własne kowariaty --------------------------------------------------
if (nzchar(HEALTH_COVARS_FILE)) {
  extra <- read_excel(HEALTH_COVARS_FILE) |> mutate(pow = sprintf("%04d", as.integer(pow)))
  new_cv <- setdiff(names(extra), "pow")
  pow <- pow |> left_join(extra, by = "pow")
  covars <- c(covars, new_cv)
  ## braki -> mediana (jak w dotychczasowym pipeline caret: medianImpute)
  for (v in new_cv) pow[[v]][is.na(pow[[v]])] <- median(pow[[v]], na.rm = TRUE)
}

## --- wartości wojewódzkie ------------------------------------------------------------
if (nzchar(HEALTH_TARGETS_FILE)) {
  tw <- read_excel(HEALTH_TARGETS_FILE)
} else {
  ## wartości woj. EHIS 2019 + SE wyprowadzone z tablic GUS (01b_build_izp.R)
  tw <- read_csv(file.path(PATH_PROC, "health_targets_2019.csv"), show_col_types = FALSE,
                 col_types = cols(woj = "c"))
}
labels_tab <- tw |> select(any_of(c("target", "id", "nazwa"))) |> distinct()
tw <- tw |> mutate(woj = sprintf("%02d", as.integer(woj)))
if (!"se" %in% names(tw)) tw$se <- 0
if (nzchar(HEALTH_SE_FILE)) {
  tw <- tw |> select(-se) |>
    left_join(read_excel(HEALTH_SE_FILE) |> mutate(woj = sprintf("%02d", as.integer(woj))),
              by = c("woj", "target"))
  tw$se[is.na(tw$se)] <- 0
}
targets <- unique(tw$target)
targets <- targets[order(readr::parse_number(targets))]

P <- make_problem(pow, W, covars = covars)
woj_levels <- levels(factor(pow$woj))

## Benchmarking na skali logitowej: dla każdego woj. przesunięcie delta_g takie, że
## sum_i w_i * 100 * expit(eta_i + delta_g) = z_g  (zachowuje zakres (0,100) i zgodność z woj.)
logit_shift_bench <- function(P, eta, z) {
  dl <- numeric(P$m)
  for (g in seq_len(P$m)) {
    ii <- which(P$g == g); wi <- P$w[ii]
    fn <- function(d) sum(wi * 100 * plogis(eta[ii] + d)) - z[g]
    dl[g] <- uniroot(fn, c(-30, 30), tol = 1e-10)$root
  }
  dl[P$g]
}

## --- dezagregacja ---------------------------------------------------------------------
long <- list(); gdr_par <- list(); diag_q2 <- list(); c_sens <- list()
for (tg in targets) {
  d <- tw |> filter(target == tg) |> arrange(match(woj, woj_levels))
  stopifnot(identical(d$woj, woj_levels))
  set.seed(SEED)
  use_logit <- APP_SCALE == "logit" && all(d$value > 0 & d$value < 100)
  ## diagnostyka: czy kowariaty przewidują wartości wojewódzkie (Q2 regionalne, LOO po woj.)
  zq <- if (use_logit) qlogis(d$value / 100) else d$value
  q2 <- regional_q2(zq, P$Xbar)
  diag_q2[[tg]] <- tibble(target = tg, q2_regional = q2, covariates_informative = q2 >= Q2_MIN,
    recommendation = if (q2 >= Q2_MIN)
      paste0("dezagregacja z kowariatami uzasadniona (", APP_METHOD, ")") else
      "kowariaty nie wyjaśniają różnic między woj. – raportować wariant bez kowariat (PYCNO/NAIVE) lub z wyraźnym zastrzeżeniem")
  if (use_logit) {
    ## skala logitowa: eta = logit(p); se metodą delta
    p  <- d$value / 100
    f  <- run_methods(P, qlogis(p), (d$se / 100) / (p * (1 - p)), methods = APP_METHODS, seed = 1)
    for (mth in names(f)) {
      dl <- logit_shift_bench(P, f[[mth]]$pred, d$value)
      sdv <- f[[mth]]$sd
      f[[mth]]$est  <- 100 * plogis(f[[mth]]$pred + dl)
      f[[mth]]$lo90 <- if (is.null(sdv)) NA_real_ else 100 * plogis(f[[mth]]$pred + dl - qnorm(.95) * sdv)
      f[[mth]]$hi90 <- if (is.null(sdv)) NA_real_ else 100 * plogis(f[[mth]]$pred + dl + qnorm(.95) * sdv)
      f[[mth]]$sd   <- if (is.null(sdv)) NA_real_ else (f[[mth]]$hi90 - f[[mth]]$lo90) / (2 * qnorm(.95))
    }
  } else {
    f <- run_methods(P, d$value, d$se, methods = APP_METHODS, seed = 1)
    for (mth in names(f)) {
      sdv <- f[[mth]]$sd; f[[mth]]$est <- f[[mth]]$pred
      f[[mth]]$lo90 <- if (is.null(sdv)) NA_real_ else f[[mth]]$pred - qnorm(.95) * sdv
      f[[mth]]$hi90 <- if (is.null(sdv)) NA_real_ else f[[mth]]$pred + qnorm(.95) * sdv
      if (is.null(sdv)) f[[mth]]$sd <- NA_real_
    }
  }
  ## wrażliwość na założenie o przenoszeniu zależności: GDR (c = 1) vs c = 0.25 (bez rozrzutu)
  zz  <- if (use_logit) qlogis(d$value / 100) else d$value
  sse <- if (use_logit) (d$se / 100) / ((d$value / 100) * (1 - d$value / 100)) else d$se
  f_lo <- gdr_fit(P, zz, sse, return_sd = FALSE, xl = list(c0 = 0.25, sd = 0))
  c_sens[[tg]] <- if (use_logit) 100 * plogis(f_lo$pred + logit_shift_bench(P, f_lo$pred, d$value)) else
    bench_add(P, f_lo$pred, d$value)
  gdr_par[[tg]] <- tibble(target = tg, scale = ifelse(use_logit, "logit", "liniowa"),
                          tau2 = f$GDR$tau2, phi = f$GDR$phi)
  for (mth in names(f)) {
    long[[length(long) + 1]] <- tibble(
      !!!as.list(pow[c("pow", "woj", "nazwa", "w")]), target = tg, method = mth,
      est = f[[mth]]$est, sd = f[[mth]]$sd, lo90 = f[[mth]]$lo90, hi90 = f[[mth]]$hi90)
  }
}
app <- bind_rows(long)

## --- wielowymiarowo: wszystkie wskaźniki naraz (MGDR), skala logitowa, benchmark logitowy ----------
if (APP_SCALE == "logit") {
  tw_w <- tw |> select(woj, target, value, se) |> arrange(match(woj, woj_levels))
  Zp <- sapply(targets, function(t) tw_w$value[tw_w$target == t])
  Sp <- sapply(targets, function(t) tw_w$se[tw_w$target == t])
  if (all(Zp > 0 & Zp < 100)) {
    pp <- Zp / 100
    mv <- mgdr_fit(P, qlogis(pp), (Sp / 100) / (pp * (1 - pp)))
    for (k in seq_along(targets)) {
      dl <- logit_shift_bench(P, mv$pred[, k], Zp[, k])
      long[[length(long) + 1]] <- tibble(
        !!!as.list(pow[c("pow", "woj", "nazwa", "w")]), target = targets[k], method = "MGDR",
        est  = 100 * plogis(mv$pred[, k] + dl),
        lo90 = 100 * plogis(mv$pred[, k] + dl - qnorm(.95) * mv$sd[, k]),
        hi90 = 100 * plogis(mv$pred[, k] + dl + qnorm(.95) * mv$sd[, k]),
        sd   = (hi90 - lo90) / (2 * qnorm(.95)))
    }
    mv_tab <- tibble(target = targets, loading = mv$loadings, common_share = mv$common_share)
    write_csv(mv_tab, file.path(PATH_RES, "app_health_mgdr_loadings.csv"))
    print(mv_tab |> mutate(across(where(is.numeric), ~ round(.x, 2))))
    app <- bind_rows(long)
  }
}

## =============================================================================
## OCENA DEZAGREGACJI – czy mapę danego wskaźnika można pokazać (prawdy powiatowej nie znamy)
## Kryteria ustalone z góry (EVAL_*), dla mapy modelowej EVAL_METHOD (domyślnie GDR):
##   q2            – Q2 regionalne (czy kowariaty tłumaczą różnice MIĘDZY województwami)
##   snr           – sygnał/szum: SD odchyleń wewnątrzwojewódzkich / średnia połowa przedziału 90%
##   agree_*       – korelacja odchyleń wewnątrzwojewódzkich z innymi metodami (PYCNO, GP-GAM, GDR_XL)
##   tercile_change – odsetek powiatów zmieniających tercyl, gdy zależność przenosi się w 25% (c = 0.25)
##                    zamiast w 100% (c = 1); top10_overlap – wspólna część górnego decyla
##   ext_*         – zgodność z danymi zewnętrznymi (external_proxies.csv, 01c): korelacja odchyleń
##                    wewnątrzwojewódzkich ze znakiem oczekiwanym; ext_calslope gdy ta sama skala
## Światło: czerwone = snr < EVAL_SNR_RED lub tercile_change > EVAL_TERC_RED lub zgodność zewn. < 0;
##          zielone  = q2 >= Q2_MIN i snr >= EVAL_SNR_GREEN i tercile_change <= EVAL_TERC_GREEN
##                     i (brak danych zewn. lub zgodność >= EVAL_EXT_GREEN); w pozostałych – żółte.
## =============================================================================
## metoda oceniana: wybrana w 06b, jeśli daje przedziały; inaczej GDR_XL
EVAL_METHOD <- Sys.getenv("EVAL_METHOD",
  if (APP_METHOD %in% c("GDR", "GDR_XL", "GP_XL", "GDR_GAM", "GP_GAM")) APP_METHOD else "GDR_XL")
## progi sygnału/szumu: na 16 prawdziwych zmiennych powiatowych (06) GDR_XL pobił metodę
## naiwną przy snr 0,14–0,38, więc próg czerwony 0,12 i zielony 0,30 (pierwotnie 0,25 i 0,5)
EVAL_SNR_GREEN <- 0.30; EVAL_SNR_RED <- 0.12
EVAL_TERC_GREEN <- 0.2; EVAL_TERC_RED <- 0.4; EVAL_EXT_GREEN <- 0.3
wdev <- function(v, woj, w) v - ave(v * w, woj, FUN = sum) / ave(w, woj, FUN = sum)
cor_s <- function(a, b) { ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 10 || sd(a[ok]) < 1e-10 || sd(b[ok]) < 1e-10) NA_real_ else cor(a[ok], b[ok]) }
f_ext <- file.path(PATH_PROC, "external_proxies.csv")
ext <- if (file.exists(f_ext)) read_csv(f_ext, show_col_types = FALSE, col_types = cols(pow = "c")) else NULL
terc <- function(v) cut(rank(v, ties.method = "first"), 3, labels = FALSE)
evals <- map_dfr(targets, function(tg) {
  a  <- app |> filter(target == tg)
  gm <- a |> filter(method == EVAL_METHOD)          # kolejność powiatów jak w P
  if (!nrow(gm)) return(NULL)
  dv <- wdev(gm$est, gm$woj, gm$w)
  other <- function(m) { o <- a |> filter(method == m); if (!nrow(o)) return(NA_real_)
    o <- o[match(gm$pow, o$pow), ]; cor_s(dv, wdev(o$est, o$woj, o$w)) }
  lo <- c_sens[[tg]]
  k <- ceiling(0.1 * nrow(gm))
  out <- tibble(target = tg, method = EVAL_METHOD,
    q2 = diag_q2[[tg]]$q2_regional,
    sd_within = sd(dv), mean_halfwidth = mean((gm$hi90 - gm$lo90) / 2, na.rm = TRUE),
    snr = sd(dv) / mean((gm$hi90 - gm$lo90) / 2, na.rm = TRUE),
    agree_pycno = other("PYCNO"), agree_gp_gam = other("GP_GAM"), agree_gdr_xl = other("GDR_XL"),
    tercile_change = mean(terc(gm$est) != terc(lo)),
    top10_overlap = length(intersect(order(-gm$est)[1:k], order(-lo)[1:k])) / k,
    ext_proxy = NA_character_, ext_cor_within = NA_real_, ext_cor_all = NA_real_, ext_calslope = NA_real_)
  if (!is.null(ext)) {
    e <- ext |> filter(target == tg)
    if (nrow(e)) {
      e1 <- e |> filter(proxy == proxy[1]); ev <- e1$value[match(gm$pow, e1$pow)]
      dvx <- wdev(ev, gm$woj, gm$w)
      out$ext_proxy <- e1$proxy[1]
      out$ext_cor_within <- e1$sign[1] * cor_s(dv, dvx)
      out$ext_cor_all <- e1$sign[1] * cor_s(gm$est, ev)
      if (isTRUE(e1$same_scale[1])) { ok <- is.finite(dvx)
        out$ext_calslope <- sum(dvx[ok] * dv[ok]) / sum(dv[ok]^2) }
    }
  }
  out
}) |>
  mutate(light = case_when(
    snr < EVAL_SNR_RED | tercile_change > EVAL_TERC_RED | (!is.na(ext_cor_within) & ext_cor_within < 0) ~ "czerwone",
    q2 >= Q2_MIN & snr >= EVAL_SNR_GREEN & tercile_change <= EVAL_TERC_GREEN &
      (is.na(ext_cor_within) | ext_cor_within >= EVAL_EXT_GREEN) ~ "zielone",
    TRUE ~ "żółte"),
    interpretation = case_when(
      light == "zielone"  ~ "mapę można pokazać (z przedziałami)",
      light == "żółte"    ~ "mapa jako scenariusz, z zastrzeżeniem; pokazać obok wariant ostrożny (PYCNO)",
      TRUE                ~ "nie pokazywać mapy powiatowej; raportować wartości wojewódzkie"))
write_csv(evals, file.path(PATH_RES, "app_health_evaluation.csv"))
print(evals |> select(target, q2, snr, tercile_change, ext_cor_within, light) |>
        mutate(across(where(is.numeric), ~ round(.x, 2))))

## kontrola spójności: agregaty = wartości woj. (dla metod z benchmarkingiem i se = 0)
chk <- app |> group_by(target, method, woj) |>
  summarise(agg = sum(est * w)) |>
  left_join(tw, by = c("woj", "target")) |>
  group_by(method) |> summarise(max_abs_diff = max(abs(agg - value)))
print(chk)

## --- porównanie z dotychczasowymi predykcjami (caret: cubist/glmnet/lars) -------------
norm_name <- function(x) x |> str_to_lower() |> str_remove("^powiat ") |>
  str_remove("^m\\. ?st\\. ") |> str_remove("^m\\. ") |> str_remove(" od \\d{4}$") |> str_trim()
F_PREV <- find_input("najlepsze_modele_predykcja_pow\\.xlsx$", required = FALSE)
if (is.na(F_PREV) && file.exists(file.path(PATH_RES, "najlepsze_modele_predykcja_pow.xlsx")))
  F_PREV <- file.path(PATH_RES, "najlepsze_modele_predykcja_pow.xlsx")
if (is.na(F_PREV)) {
  message("Brak pliku najlepsze_modele_predykcja_pow.xlsx – pomijam porównanie z poprzednimi predykcjami.")
  cmp <- tibble(target = character(), method = character())
} else {
prev <- read_excel(F_PREV) |>
  mutate(woj = sprintf("%02d", as.integer(woj)), key = norm_name(nazwa)) |>
  pivot_longer(starts_with("x"), names_to = "target", values_to = "prev")
keytab <- pow |> transmute(pow, woj, key = norm_name(nazwa))
prev <- prev |> left_join(keytab, by = c("woj", "key"))
message("Dopasowano poprzednie predykcje: ", sum(!is.na(prev$pow)), " / ", nrow(prev))

cmp <- app |> left_join(prev |> select(pow, target, prev), by = c("pow", "target")) |>
  group_by(target, method) |>
  mutate(dev_est = est - ave(est, woj), dev_prev = prev - ave(prev, woj)) |>
  summarise(cor_prev = suppressWarnings(cor(est, prev, use = "complete.obs")),
            cor_within_prev = suppressWarnings(cor(dev_est, dev_prev, use = "complete.obs")),
            sd_within_est = sd(dev_est), sd_within_prev = sd(dev_prev, na.rm = TRUE))
}
write_csv(cmp, file.path(PATH_RES, "app_health_compare_previous.csv"))

## --- zapis ---------------------------------------------------------------------------------
wide_primary <- app |> filter(method == APP_METHOD) |>
  select(pow, woj, nazwa, target, est, lo90, hi90) |>
  pivot_wider(id_cols = c(pow, woj, nazwa), names_from = target, values_from = c(est, lo90, hi90), names_glue = "{target}_{.value}")
write_xlsx(list(
  estymacje = wide_primary,
  wszystkie_metody = app,
  porownanie_z_poprzednimi = cmp,
  parametry_gdr = bind_rows(gdr_par),
  diagnostyka_q2 = bind_rows(diag_q2),
  ocena_dezagregacji = evals,
  mgdr_ladunki = if (exists("mv_tab")) mv_tab else tibble(),
  opis_wskaznikow = labels_tab,
  kowariaty = tibble(zestaw = APP_COVSET, kowariata = covars),
  spojnosc_z_woj = chk),
  file.path(PATH_RES, "app_health_estimates.xlsx"))
saveRDS(app, file.path(PATH_RES, "app_health_long.rds"))

## =============================================================================
## MODUŁY OPCJONALNE (wymagają danych spoza tego repozytorium)
## =============================================================================

## (1) Estymator syntetyczny demograficzny / standaryzacja pośrednia (Gonzalez 1973):
##     rates_woj: data.frame(woj, cell, rate)   – wskaźnik EHIS w komórkach wiek×płeć na poziomie woj.
##     pop_pow  : data.frame(pow, woj, cell, pop) – ludność powiatów w tych samych komórkach (GUS BDL)
##     Wynik: wskaźnik powiatowy = sum_c pop_share_{i,c} * rate_{g(i),c}, opcjonalnie z benchmarkiem.
synthetic_demographic <- function(rates_woj, pop_pow, z_woj = NULL, P = NULL) {
  est <- pop_pow |>
    group_by(pow) |> mutate(share = pop / sum(pop)) |> ungroup() |>
    left_join(rates_woj, by = c("woj", "cell")) |>
    group_by(pow, woj) |> summarise(est = sum(share * rate), .groups = "drop") |>
    arrange(pow)
  if (!is.null(z_woj) && !is.null(P)) est$est <- bench_add(P, est$est, z_woj)
  est
}

## (2) MRP (multilevel regression + poststratification) – szkic dla mikrodanych EHIS:
##     library(lme4)
##     fit <- glmer(y ~ wiek_kat * plec + wykszt + (1 | woj) + x_woj1 + x_woj2,
##                  family = binomial, data = ehis_micro, weights = waga)
##     cells$p <- predict(fit, newdata = cells, allow.new.levels = TRUE, type = "response")
##     est_pow <- cells |> group_by(pow) |> summarise(est = weighted.mean(p, pop))
##     Następnie benchmark do oficjalnych wartości woj. (bench_add / bench_mult).

## (3) Dezagregacja kompozycyjna (kategorie sumujące się do 100%):
##     transformacja ALR -> dezagregacja każdej współrzędnej (np. GDR) -> transformacja odwrotna
##     -> raking multiplikatywny do wartości woj. (iteracyjnie, IPF), np.:
disagg_compositional <- function(P, Zcat, method = m_gdr) {
  K <- ncol(Zcat); L <- log(Zcat[, -K, drop = FALSE] / Zcat[, K])
  Lhat <- sapply(seq_len(K - 1), function(k) method(P, L[, k], rep(0, P$m))$pred)
  E <- cbind(exp(Lhat), 1); Phat <- 100 * E / rowSums(E)
  for (it in 1:50) {                         # IPF: zgodność z woj. i suma = 100
    for (k in seq_len(K)) Phat[, k] <- bench_mult(P, Phat[, k], Zcat[, k])
    Phat <- 100 * Phat / rowSums(Phat)
  }
  Phat
}

write_csv(bind_rows(diag_q2), file.path(PATH_RES, "app_health_q2_diagnostics.csv"))
print(bind_rows(diag_q2) |> select(target, q2_regional, covariates_informative))
message("Aplikacja zakończona: ", length(targets), " wskaźników × ", length(APP_METHODS), " metod.")
