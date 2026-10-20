## =============================================================================
## 03_scenarios.R — procesy generujące dane (DGP) na rzeczywistych 380 powiatach
## Każdy scenariusz zwraca wektor y (prawdziwe wartości powiatowe, skala ~ %),
## złożony ze standaryzowanych komponentów o zadanych udziałach wariancji.
## =============================================================================

zstd <- function(v) { s <- sd(v); if (s < 1e-12) v * 0 else (v - mean(v)) / s }

## Stałe komponenty (liczone raz)
make_components <- function(P, pow_df) {
  Xs <- P$X                                             # standaryzowane kowariaty obserwowane
  b  <- c(urb = 0.4, ldens = 0.3, city = 0.3, ldcap = -0.2, unemp19 = 0.5, turnout19 = -0.4)
  f_lin <- zstd(as.numeric(Xs[, names(b)] %*% b))

  ## miasto–wieś: nieliniowa funkcja urbanizacji + skok dla miast na prawach powiatu
  urb <- pow_df$urb
  f_urb <- zstd(1.2 * pow_df$city + plogis((urb - 60) / 8))
  ## suburbanizacja („obwarzanek”): powiaty ziemskie sąsiadujące z miastem > 150 tys.
  big_city <- pow_df$city == 1 & pow_df$pop > 120000
  suburb <- as.numeric(pow_df$city == 0 & (P$W %*% big_city) > 0)
  f_sub <- zstd(suburb)

  ## Gradient miasto–wieś (S2): potencjał miejski – suma po miastach na prawach powiatu,
  ## ważona ich wielkością i zanikająca wykładniczo z odległością (skala URBAN_BW_KM):
  ##   pot_i = sum_c (pop_c / 1e5) * exp(-d_ic / bw)
  ## Najwyżej w mieście, wysoko w powiatach przylegających, coraz niżej im dalej od miast.
  ## Kowariaty widzą to tylko częściowo (urbanizacja, gęstość, odległość od STOLICY woj.).
  cities <- which(pow_df$city == 1)
  xy <- cbind(pow_df$x_km, pow_df$y_km)
  d_city <- sqrt(outer(xy[, 1], xy[cities, 1], "-")^2 + outer(xy[, 2], xy[cities, 2], "-")^2)
  pot <- as.numeric(exp(-d_city / URBAN_BW_KM) %*% (pow_df$pop[cities] / 1e5))
  f_dist <- zstd(log(pot + 0.01))

  ## zabory
  zab_eff <- c(R = 1.0, A = 1.6, P = -0.6, ZO = -1.1)
  f_zab <- zstd(zab_eff[pow_df$zabor])

  ## gradient wschód–zachód (wsch. = wyżej)
  f_ew <- zstd(pow_df$x_km)

  ## komponenty do błędu ekologicznego: x̄_g (między) i x_i - x̄_g (wewnątrz)
  xb <- as.numeric(P$Xbar[P$g, "unemp19"]); xw <- Xs[, "unemp19"] - xb
  tb <- as.numeric(P$Xbar[P$g, "turnout19"]); tw <- Xs[, "turnout19"] - tb
  f_eco <- zstd(1.0 * zstd(xb) - 1.0 * zstd(xw) - 0.6 * zstd(tb) + 0.6 * zstd(tw))

  list(f_lin = f_lin, f_urb = f_urb, f_sub = f_sub, suburb = suburb, f_dist = f_dist,
       d_nearest_city = apply(d_city, 1, min), f_zab = f_zab, f_ew = f_ew,
       f_eco = f_eco, xy = cbind(pow_df$x_km, pow_df$y_km))
}

draw_car <- function(P) zstd(as.numeric(crossprod(P$Rchol, rnorm(P$n))))

draw_hotspots <- function(P, K, n_hot = 6, bw_km = 25) {
  ctr <- sample.int(P$n, n_hot)
  d2 <- sapply(ctr, function(c) colSums((t(K$xy) - K$xy[c, ])^2))
  amp <- runif(n_hot, 0.7, 1.3)
  h <- as.numeric(exp(-d2 / (2 * bw_km^2)) %*% amp)
  list(f = zstd(h), raw = h, centers = ctr)
}

## Komponent regionalny: niezależny efekt województwa (różnice regionalne niewyjaśnione
## kowariatami – polityka zdrowotna, kultura regionalna). Standaryzowany na poziomie powiatów.
draw_regional <- function(P) zstd(rnorm(P$m)[P$g])

## Stosunek wariancji wewnątrz- do międzywojewódzkiej (jak w 01_build_data.R)
variance_ratio <- function(P, y) {
  zt <- agg(P, y)
  (sum(P$w * (y - zt[P$g])^2) / P$m) / var(zt)
}

## Główna funkcja: zwraca y (średnia 50, odch. std. 5).
## lambda – udział wariancji komponentu regionalnego (stały w województwie);
## psi    – jaka część komponentu regionalnego jest wyjaśniona kowariatami: sqrt(psi) * średnia
##          wojewódzka sygnału kowariat (f_lin) + sqrt(1 - psi) * niezależny efekt województwa.
##          psi = 0: różnice regionalne całkowicie niewyjaśnione (wariant pesymistyczny).
simulate_truth <- function(scen, P, K, r2 = NULL, lambda = 0, psi = 0) {
  u <- draw_car(P); e <- rnorm(P$n)
  comp <- switch(scen,
    S1_base     = sqrt(.70) * K$f_lin + sqrt(.15) * u + sqrt(.15) * e,
    S2_urban    = sqrt(.65) * K$f_dist + sqrt(.20) * K$f_lin + sqrt(.15) * e,
    S3_zabory   = sqrt(.50) * K$f_zab + sqrt(.20) * K$f_lin + sqrt(.15) * u + sqrt(.15) * e,
    S4_hotspot  = { h <- draw_hotspots(P, K)
                    sqrt(.30) * h$f + sqrt(.50) * K$f_lin + sqrt(.10) * u + sqrt(.10) * e },
    S5_fuzzy    = sqrt(.15) * K$f_lin + sqrt(.15) * u + sqrt(.70) * e,
    S6_ecol     = sqrt(.60) * K$f_eco + sqrt(.15) * u + sqrt(.25) * e,
    S7_gradient = sqrt(.35) * K$f_ew + sqrt(.35) * u + sqrt(.15) * K$f_lin + sqrt(.15) * e,
    SWEEP       = sqrt(r2) * K$f_lin + sqrt((1 - r2) / 2) * u + sqrt((1 - r2) / 2) * e,
    stop("Nieznany scenariusz: ", scen))
  comp <- zstd(comp)
  if (lambda > 0) {
    reg_x <- zstd(agg(P, K$f_lin)[P$g])              # część regionalna związana z kowariatami
    reg <- zstd(sqrt(psi) * reg_x + sqrt(1 - psi) * draw_regional(P))
    comp <- sqrt(1 - lambda) * comp + sqrt(lambda) * reg
  }
  50 + 5 * zstd(comp)
}

## Q2 regionalne dla prawdziwych y (te same 3 składowe co w diagnostyce)
q2_of <- function(P, y) regional_q2(agg(P, y), P$Xbar)

## Kalibracja dwóch cech prawdziwych danych (01_build_data.R, tab_variance_ratio_reference.csv):
##   (1) stosunek wariancji wewnątrz/między ≈ ratio_target  -> lambda
##   (2) Q2 regionalne (jak dobrze kowariaty przewidują wartości woj.) ≈ q2_target -> psi
## Dla każdego psi z siatki szukamy najmniejszego lambda dającego (1), potem wybieramy psi,
## przy którym średnie Q2 jest najbliższe (2). psi = NULL -> tylko (1), psi = 0.
## Scenariusze, które już mają stosunek <= ratio_target, zostają bez zmian (lambda = psi = 0).
calibrate_design <- function(scen, P, K, ratio_target = 2, q2_target = NULL,
                             n_draws = 20, seed = 1) {
  stat_at <- function(lam, psi) {
    set.seed(seed)
    v <- replicate(n_draws, {
      y <- simulate_truth(scen, P, K, lambda = lam, psi = psi)
      c(variance_ratio(P, y), q2_of(P, y))
    })
    rowMeans(v)
  }
  if (stat_at(0, 0)[1] <= ratio_target) return(c(lambda = 0, psi = 0))
  lam_for <- function(psi) {
    for (lam in seq(0.02, 0.9, by = 0.02)) if (stat_at(lam, psi)[1] <= ratio_target) return(lam)
    0.9
  }
  if (is.null(q2_target)) return(c(lambda = lam_for(0), psi = 0))
  psis <- seq(0, 1, by = 0.1)
  res <- t(sapply(psis, function(ps) { l <- lam_for(ps); c(lambda = l, psi = ps, q2 = stat_at(l, ps)[2]) }))
  best <- which.min(abs(res[, "q2"] - q2_target))
  res[best, c("lambda", "psi")]
}
calibrate_lambda <- function(scen, P, K, target = 2, ...)          # zgodność wsteczna
  calibrate_design(scen, P, K, ratio_target = target)[["lambda"]]

## Szum sondażowy na poziomie woj. (jak w EHIS): SE ~ 1/sqrt(liczebność próby ~ ludność)
survey_se <- function(P, z_true, k = 0.4) {
  popg <- as.numeric(tapply(P$pop, P$g, sum))
  k * sd(z_true) * sqrt(mean(popg) / popg)
}
