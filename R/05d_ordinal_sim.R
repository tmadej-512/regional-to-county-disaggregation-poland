## =============================================================================
## 05d_ordinal_sim.R — symulacja dla danych PORZĄDKOWYCH (pytania EHIS ze skalą odpowiedzi)
##
## Pytanie: czy dezagregacja całego rozkładu odpowiedzi (ORD_*, 02b_ordinal.R) daje lepsze
## szacunki powiatowe wskaźnika IZP niż dezagregacja samego odsetka (BIN_*, jak w 07)?
## I kiedy przestaje działać (złamane założenie wspólnych progów)?
##
## (A) JEDNO PYTANIE NARAZ
##   pytania: x1–x5 (liczba kategorii 5, 4, 5, 7, 3; progi z krajowego rozkładu EHIS 2019)
##   struktura przestrzenna eta_i: scenariusze z 03_scenarios.R (S1, S2, S3, S5, S7),
##     design realistic (kalibracja jak w 04) i no_regional (bez komponentu regionalnego)
##   skala eta: dobrana tak, by SD wartości wojewódzkich eta była jak w PRAWDZIWYCH danych EHIS
##     dla danego pytania (ord_scores na tablicy GUS); bez pliku kategorii: SAE_ORD_SETA (0,3)
##   kształt (shape): proportional (model poprawny), scale (różna rozpiętość odpowiedzi,
##     rośnie z urbanizacją), nonprop (nieproporcjonalne szanse)
##   wartości woj.: exact albo survey (losowanie wielomianowe z n_eff pytania × NEFF_MULT)
##   metody: PYCNO, GDR, GDR_XL w wersji BIN i ORD + NAIVE
## (B) WSZYSTKIE PYTANIA NARAZ (wielowymiarowo)
##   eta_k = sqrt(c) * l_k * s + sqrt(1 - c) * u_k (wspólne pole „profilu zdrowotnego” s),
##   c ∈ {0, 0.6}; porównanie: BIN_GDR, ORD_GDR, BIN_MGDR (MGDR na logitach odsetków),
##   ORD_MGDR (MGDR na wynikach porządkowych)
## Wynik: output/results/ordinal_sim_metrics.csv, output/results/ordinal_sim_mv_metrics.csv
## Liczba replikacji: SAE_NREP_ORD (domyślnie 30), SAE_NREP_ORD_MV (domyślnie 30); 9000 + 120 przebiegów
## =============================================================================
source("R/00_config.R"); source("R/02_methods.R"); source("R/02b_ordinal.R"); source("R/03_scenarios.R")

N_REP_ORD    <- as.integer(Sys.getenv("SAE_NREP_ORD", "30"))
N_REP_ORD_MV <- as.integer(Sys.getenv("SAE_NREP_ORD_MV", "30"))
S_ETA_DEFAULT <- as.numeric(Sys.getenv("SAE_ORD_SETA", "0.3"))
NEFF_MULT <- as.numeric(Sys.getenv("SAE_ORD_NEFF_MULT", "1"))   # >1: większa próba niż w EHIS
N_EFF_DEFAULT <- 800
ORD_G <- as.numeric(Sys.getenv("SAE_ORD_SCALE_G", "0.1"))     # siła złamania „scale”
ORD_H <- as.numeric(Sys.getenv("SAE_ORD_NONPROP_H", "0.8"))   # siła złamania „nonprop”
METH <- c("PYCNO", "GDR", "GDR_XL")

pow <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |> arrange(pow)
W   <- readRDS(file.path(PATH_PROC, "W_queen.rds"))
P   <- make_problem(pow, W)
K   <- make_components(P, pow)
CAL <- read_csv(file.path(PATH_RES, "sim_lambda_calibration.csv"), show_col_types = FALSE)
if (!"design" %in% names(CAL)) CAL$design <- "realistic"          # starszy format pliku
if (!"psi" %in% names(CAL)) CAL$psi <- 0
cal_row <- function(sc) { r <- CAL[CAL$design == "realistic" & CAL$scen == sc, ]
  if (nrow(r)) r[1, ] else tibble(lambda = 0, psi = 0) }
woj_levels <- levels(factor(pow$woj))

## --- parametry pytań: progi, licznik, n_eff i realna zmienność wojewódzka eta ------------------
f_cat <- file.path(PATH_PROC, "ehis_categories_2019.csv")
cats <- if (file.exists(f_cat)) read_csv(f_cat, show_col_types = FALSE, col_types = cols(woj = "c")) else NULL
QP <- lapply(seq_len(nrow(ORD_DEFS)), function(j) {
  pl <- ORD_DEFS$pl[[j]] / 100
  th <- qnorm(cumsum(pl / sum(pl)))[seq_len(length(pl) - 1)]
  neff <- rep(N_EFF_DEFAULT, P$m); sd_reg <- NA_real_
  if (!is.null(cats)) {
    cj <- cats |> filter(str_detect(variable, ORD_DEFS$var_pattern[j]))
    if (nrow(cj)) {
      Pr <- cj |> select(woj, k, pct) |> pivot_wider(names_from = k, values_from = pct) |>
        arrange(match(woj, woj_levels))
      ne <- cj |> distinct(woj, neff) |> arrange(match(woj, woj_levels)) |> pull(neff)
      if (length(ne) == P$m && all(is.finite(ne))) neff <- ne
      sd_reg <- sd(ord_scores(as.matrix(Pr[, -1]) / 100, neff)$eta)
    }
  }
  list(target = ORD_DEFS$target[j], theta = th, num = ORD_DEFS$num[[j]], neff = neff * NEFF_MULT,
       sd_reg = sd_reg)
})
names(QP) <- ORD_DEFS$target
print(tibble(target = names(QP), K = sapply(QP, function(q) length(q$theta) + 1),
             sd_reg_eta = sapply(QP, `[[`, "sd_reg"), neff_med = sapply(QP, function(q) median(q$neff))))

## skala eta: SD wartości wojewódzkich jak w danych (wariancja wewnątrz wynika ze scenariusza)
scale_eta <- function(y, q) {
  e <- zstd(y)
  if (is.finite(q$sd_reg)) e * q$sd_reg / sd(agg(P, e)) else e * S_ETA_DEFAULT
}
draw_regional_shares <- function(Pc, q, noise) {
  Pr_true <- as.matrix(P$Am %*% Pc)
  if (noise == "exact") return(list(Pr = Pr_true, neff = rep(1e6, P$m)))
  Pr <- t(vapply(seq_len(P$m), function(g) {
    n <- max(round(q$neff[g]), 20)
    as.numeric(rmultinom(1, n, Pr_true[g, ])) / n }, numeric(ncol(Pc))))
  Pr <- pmax(Pr, 1e-5)
  list(Pr = Pr / rowSums(Pr), neff = q$neff)
}
## wersja binarna: logit odsetka, SE dwumianowe, benchmark przesunięciem logitu (jak w 07)
bin_fit <- function(zb, neff, methods, seed, fits = NULL) {
  seb <- sqrt(zb * (1 - zb) / neff) / (zb * (1 - zb))
  if (is.null(fits)) fits <- run_methods(P, qlogis(zb), seb, methods = methods, seed = seed)
  lapply(fits, function(f) {
    dl <- vapply(seq_len(P$m), function(g) { ii <- which(P$g == g)
      uniroot(function(d) sum(P$w[ii] * plogis(f$pred[ii] + d)) - zb[g], c(-30, 30))$root }, 0)[P$g]
    est <- 100 * plogis(f$pred + dl)
    sd  <- if (is.null(f$sd)) NULL else
      100 * (plogis(f$pred + dl + 1.645 * f$sd) - plogis(f$pred + dl - 1.645 * f$sd)) / (2 * 1.645)
    list(est = est, sd = sd)
  })
}
ord_out <- function(od) lapply(od$fits, function(f)
  list(est = f$est, sd = if (all(is.na(f$lo90))) NULL else (f$hi90 - f$lo90) / (2 * 1.645)))

## =============================================================================
## (A) jedno pytanie naraz
## =============================================================================
grid <- expand_grid(target = names(QP), scen = c("S1_base", "S2_urban", "S3_zabory", "S5_fuzzy", "S7_gradient"),
                    design = c("realistic", "no_regional"), shape = c("proportional", "scale", "nonprop"),
                    noise = c("exact", "survey"), rep = seq_len(N_REP_ORD))

one_run <- function(i) {
  gr <- grid[i, ]; q <- QP[[gr$target]]
  set.seed(SEED + 7e6 + i)
  cal <- if (gr$design == "realistic") cal_row(gr$scen) else tibble(lambda = 0, psi = 0)
  y   <- simulate_truth(gr$scen, P, K, lambda = cal$lambda, psi = cal$psi)
  eta <- scale_eta(y, q)
  Pc  <- sim_ord_probs(eta, q$theta, gr$shape, aux = pow$urb, g = ORD_G, h = ORD_H)
  truth <- 100 * rowSums(Pc[, q$num, drop = FALSE])
  rs  <- draw_regional_shares(Pc, q, gr$noise)
  zb  <- rowSums(rs$Pr[, q$num, drop = FALSE])
  b   <- bin_fit(zb, rs$neff, METH, seed = gr$rep)
  o   <- ord_out(ord_disagg(P, rs$Pr, rs$neff, q$num, methods = METH, seed = gr$rep, ndraw = 200))
  bind_rows(
    metrics_pow(P, truth, (100 * zb)[P$g]) |> mutate(method = "NAIVE"),
    imap_dfr(b, ~ metrics_pow(P, truth, .x$est, .x$sd) |> mutate(method = paste0("BIN_", .y))),
    imap_dfr(o, ~ metrics_pow(P, truth, .x$est, .x$sd) |> mutate(method = paste0("ORD_", .y)))) |>
    mutate(!!!as.list(gr), .before = 1)
}

t0 <- Sys.time()
res <- par_lapply(seq_len(nrow(grid)), one_run)
message("(A) czas: ", format(Sys.time() - t0))
res <- res[!sapply(res, inherits, "try-error")]
write_csv(bind_rows(res), file.path(PATH_RES, "ordinal_sim_metrics.csv"))

## =============================================================================
## (B) wszystkie pytania naraz: wspólny „profil zdrowotny” powiatu
## =============================================================================
grid_mv <- expand_grid(c_common = c(0, 0.6), noise = c("exact", "survey"), rep = seq_len(N_REP_ORD_MV))
TG <- names(QP)

one_mv <- function(i) {
  gr <- grid_mv[i, ]
  set.seed(SEED + 7.5e6 + i)
  cal <- cal_row("S1_base")
  s <- zstd(simulate_truth("S2_urban", P, K))                    # wspólne pole: gradient miejski
  sgn <- sample(c(-1, 1), length(TG), replace = TRUE)
  sims <- lapply(seq_along(TG), function(k) {
    q <- QP[[TG[k]]]
    uk <- zstd(simulate_truth("S1_base", P, K, lambda = cal$lambda, psi = cal$psi))
    eta <- scale_eta(sqrt(gr$c_common) * sgn[k] * s + sqrt(1 - gr$c_common) * uk, q)
    Pc <- sim_ord_probs(eta, q$theta, "proportional")
    rs <- draw_regional_shares(Pc, q, gr$noise)
    list(q = q, truth = 100 * rowSums(Pc[, q$num, drop = FALSE]), rs = rs,
         zb = rowSums(rs$Pr[, q$num, drop = FALSE]), sc = ord_scores(rs$Pr, rs$neff))
  })
  ## MGDR na wynikach porządkowych i na logitach odsetków
  Zo <- sapply(sims, function(x) x$sc$eta); So <- sapply(sims, function(x) x$sc$se)
  Zb <- sapply(sims, function(x) qlogis(x$zb))
  Sb <- sapply(sims, function(x) sqrt(x$zb * (1 - x$zb) / x$rs$neff) / (x$zb * (1 - x$zb)))
  mo <- mgdr_fit(P, Zo, So); mb <- mgdr_fit(P, Zb, Sb)
  bind_rows(lapply(seq_along(TG), function(k) {
    x <- sims[[k]]
    gb <- bin_fit(x$zb, x$rs$neff, "GDR", seed = gr$rep)$GDR
    go <- ord_out(ord_disagg(P, x$rs$Pr, x$rs$neff, x$q$num, methods = "GDR", seed = gr$rep, ndraw = 200))$GDR
    bm <- bin_fit(x$zb, x$rs$neff, NULL, gr$rep, fits = list(M = list(pred = mb$pred[, k], sd = mb$sd[, k])))$M
    om <- ord_out(ord_disagg(P, x$rs$Pr, x$rs$neff, x$q$num,
                             fits = list(M = list(pred = mo$pred[, k], sd = mo$sd[, k]))))$M
    bind_rows(metrics_pow(P, x$truth, (100 * x$zb)[P$g]) |> mutate(method = "NAIVE"),
              metrics_pow(P, x$truth, gb$est, gb$sd) |> mutate(method = "BIN_GDR"),
              metrics_pow(P, x$truth, go$est, go$sd) |> mutate(method = "ORD_GDR"),
              metrics_pow(P, x$truth, bm$est, bm$sd) |> mutate(method = "BIN_MGDR"),
              metrics_pow(P, x$truth, om$est, om$sd) |> mutate(method = "ORD_MGDR")) |>
      mutate(target = TG[k])
  })) |> mutate(!!!as.list(gr), .before = 1)
}

t0 <- Sys.time()
res_mv <- par_lapply(seq_len(nrow(grid_mv)), one_mv)
message("(B) czas: ", format(Sys.time() - t0))
res_mv <- res_mv[!sapply(res_mv, inherits, "try-error")]
write_csv(bind_rows(res_mv), file.path(PATH_RES, "ordinal_sim_mv_metrics.csv"))
