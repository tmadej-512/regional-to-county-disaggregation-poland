## =============================================================================
## 02b_ordinal.R — dezagregacja PORZĄDKOWA (pytania EHIS ze skalą odpowiedzi)
##
## Pytania EHIS mają uporządkowane kategorie (samoocena zdrowia 1–5, BMI 4 klasy, częstość
## spożywania warzyw 5, picia alkoholu 7, palenie 3). Wskaźnik IZP to suma kilku kategorii
## (np. „złe + bardzo złe”), więc dotąd wyrzucaliśmy resztę rozkładu.
##
## Model (kumulatywny probit, wspólne progi):
##   osoba w powiecie i ma ukrytą skłonność  y* ~ N(eta_i, 1),  kategoria k gdy theta_{k-1} < y* <= theta_k
##   P(Y <= k | i) = Phi(theta_k - eta_i)
## Kroki:
##   (1) na poziomie województw: z udziałów wszystkich kategorii szacujemy JEDEN wynik eta_g
##       (uogólnione MNK na skali probitowej, kowariancja wielomianowa z n_eff) i progi theta_k
##       wspólne dla kraju; SE(eta_g) wykorzystuje informację ze wszystkich kategorii;
##   (2) eta_g dezagregujemy dowolną metodą z 02_methods.R (GDR, GDR_XL, PYCNO, ...) albo
##       wszystkie pytania naraz (MGDR);
##   (3) rozkład kategorii w powiecie: p_ik = Phi(theta_k - eta_i) - Phi(theta_{k-1} - eta_i);
##   (4) raking (IPF) do obserwowanych udziałów wojewódzkich WSZYSTKICH kategorii (dokładny
##       benchmark, wiersze sumują się do 1);
##   (5) wskaźnik = suma kategorii licznika; przedziały z losowań eta_i ~ N(pred, sd) z tą samą
##       poprawką rakingu (na skali logitowej) co estymacja punktowa.
## Przybliżenie: wartość wojewódzka eta_g ≈ suma w_i eta_i (dokładne, gdy zróżnicowanie eta między
## powiatami jest małe względem zmienności osób, sd = 1); błąd tego przybliżenia usuwa krok (4).
##
## Dlaczego nie regresja Dirichleta: model Dirichleta opisuje skład na poziomie jednostek, a my
## mamy tylko 16 składów wojewódzkich – dopasowany na nich byłby regresją „makro” (te same problemy
## z przenoszeniem zależności), nie zachowuje porządku kategorii, potrzebuje K−1 osobnych efektów
## i nie daje dokładnej zgodności z województwami (średnia Dirichleta nie agreguje się liniowo).
## Wariant kompozycyjny dla kategorii NIEuporządkowanych: disagg_compositional() w 07.
## =============================================================================

## (1) wyniki wojewódzkie: Pr (m × K, wiersze sumują się do 1), neff (m)
ord_scores <- function(Pr, neff, eps = 1e-4, maxit = 200, tol = 1e-10) {
  Pr <- Pr / rowSums(Pr)
  m <- nrow(Pr); K <- ncol(Pr); K1 <- K - 1
  C <- t(apply(Pr, 1, cumsum))[, seq_len(K1), drop = FALSE]
  C <- pmin(pmax(C, eps), 1 - eps)
  Y <- qnorm(C)                                           # y_gk = theta_k - eta_g
  Wi <- lapply(seq_len(m), function(g) {
    cg <- C[g, ]
    S  <- outer(seq_len(K1), seq_len(K1), function(a, b) cg[pmin(a, b)] * (1 - cg[pmax(a, b)])) / neff[g]
    J  <- diag(1 / dnorm(Y[g, ]), K1)
    solve(J %*% S %*% J + diag(1e-10, K1))
  })
  eta <- rep(0, m); theta <- colMeans(Y)
  SW <- Reduce(`+`, Wi)
  for (it in seq_len(maxit)) {
    eta_new <- vapply(seq_len(m), function(g) sum(Wi[[g]] %*% (theta - Y[g, ])) / sum(Wi[[g]]), 0)
    eta_new <- eta_new - mean(eta_new)                    # identyfikowalność: średnia eta = 0
    theta_new <- as.numeric(solve(SW, Reduce(`+`, lapply(seq_len(m), function(g)
      Wi[[g]] %*% (Y[g, ] + eta_new[g])))))
    done <- max(abs(eta_new - eta)) + max(abs(theta_new - theta)) < tol
    eta <- eta_new; theta <- theta_new
    if (done) break
  }
  if (any(diff(theta) <= 0)) warning("ord_scores: progi nie są rosnące – sprawdź kolejność kategorii")
  resid <- Y - (matrix(theta, m, K1, byrow = TRUE) - eta)   # reszty modelu wspólnych progów
  list(eta = eta, se = vapply(Wi, function(w) 1 / sqrt(sum(w)), 0), theta = theta,
       resid = resid,
       ## dopasowanie: odsetek zmienności probitów wyjaśniony przez jeden wynik na województwo
       r2_threshold = 1 - sum(resid^2) / sum((Y - matrix(colMeans(Y), m, K1, byrow = TRUE))^2))
}

## (3) rozkład kategorii przy danym eta (wektor n) i progach theta (K-1)
ord_probs <- function(eta, theta) {
  Cm <- pnorm(outer(-eta, theta, `+`))                   # P(Y <= k)
  cbind(Cm, 1) - cbind(0, Cm)
}

## (4) raking (IPF) do udziałów wojewódzkich Pr (m × K)
ipf_rake <- function(P, Pi, Pr, iter = 100, tol = 1e-12) {
  Pr <- Pr / rowSums(Pr)
  for (it in seq_len(iter)) {
    old <- Pi
    for (k in seq_len(ncol(Pi))) {
      a <- as.numeric(P$Am %*% Pi[, k])
      Pi[, k] <- Pi[, k] * (Pr[, k] / pmax(a, 1e-300))[P$g]
    }
    Pi <- Pi / rowSums(Pi)
    if (max(abs(Pi - old)) < tol) break
  }
  Pi
}

## (2)–(5) pełna ścieżka dla jednego pytania; num = indeksy kategorii licznika wskaźnika
## fits: wynik run_methods() na eta (opcjonalnie podany z zewnątrz, np. z MGDR)
ord_disagg <- function(P, Pr, neff, num, methods = c("PYCNO", "GDR", "GDR_XL"), fits = NULL,
                       ndraw = 400, seed = 1, tau_ratio = TAU_RATIO) {
  sc <- ord_scores(Pr, neff)
  if (is.null(fits)) fits <- run_methods(P, sc$eta, sc$se, methods = methods, seed = seed,
                                         tau_ratio = tau_ratio)
  Prn <- Pr / rowSums(Pr)
  out <- lapply(names(fits), function(mth) {
    f <- fits[[mth]]
    P0 <- ord_probs(f$pred, sc$theta)
    Pi <- ipf_rake(P, P0, Prn)
    s_raw <- rowSums(P0[, num, drop = FALSE]); s <- rowSums(Pi[, num, drop = FALSE])
    est <- 100 * s
    lo <- hi <- rep(NA_real_, P$n)
    if (!is.null(f$sd)) {
      set.seed(seed)
      adj <- qlogis(pmin(pmax(s, 1e-9), 1 - 1e-9)) - qlogis(pmin(pmax(s_raw, 1e-9), 1 - 1e-9))
      dr <- replicate(ndraw, {
        e <- f$pred + rnorm(P$n) * f$sd
        sr <- rowSums(ord_probs(e, sc$theta)[, num, drop = FALSE])
        plogis(qlogis(pmin(pmax(sr, 1e-9), 1 - 1e-9)) + adj)
      })
      lo <- 100 * apply(dr, 1, quantile, 0.05); hi <- 100 * apply(dr, 1, quantile, 0.95)
    }
    list(method = mth, est = est, lo90 = lo, hi90 = hi, probs = Pi)
  })
  names(out) <- names(fits)
  list(fits = out, scores = sc)
}

## =============================================================================
## Wspólne definicje pytań porządkowych (07b i symulacja 05d)
##   num    – kategorie licznika wskaźnika IZP (kolejność kategorii jak w tablicy GUS)
##   pl     – rozkład krajowy EHIS 2019 (%), z tablicy GUS; z niego progi w symulacji
## =============================================================================
ORD_DEFS <- tibble(
  target = c("x1", "x2", "x3", "x4", "x5"),
  var_pattern = c("^Ocena stanu zdrowia", "^Index masy ciała", "^Częstość spożywania warzyw",
                  "^Częstość picia alkoholu", "^Palenie tytoniu"),
  num = list(4:5, 3:4, 2:5, 1:4, 1),
  pl  = list(c(24.50, 42.43, 24.05, 7.46, 1.56),
             c(3.04, 40.40, 38.09, 18.47),
             c(51.74, 25.62, 18.89, 3.53, 0.22),
             c(1.56, 1.04, 3.80, 11.93, 15.18, 40.58, 25.91),
             c(19.41, 3.22, 77.37))
)

## Prawdziwe rozkłady kategorii w powiatach do SYMULACJI, z możliwością złamania założeń modelu:
##   shape = "proportional" – model poprawny: P(Y<=k|i) = Phi(theta_k - eta_i)
##   shape = "scale"        – różna rozpiętość odpowiedzi między powiatami (np. polaryzacja w
##                            miastach): P(Y<=k|i) = Phi((theta_k - eta_i) / s_i), s_i = exp(g * zstd(a_i));
##                            g = 0.1 daje s_i ok. 0,78–1,28
##   shape = "nonprop"      – nieproporcjonalne szanse: powiaty różnią się bardziej w jednym końcu
##                            skali: P(Y<=k|i) = Phi(theta_k - eta_i * (1 + h * (k - kbar) / (K-1)))
sim_ord_probs <- function(eta, theta, shape = c("proportional", "scale", "nonprop"),
                          aux = NULL, g = 0.1, h = 0.8) {
  shape <- match.arg(shape)
  K1 <- length(theta)
  Cm <- switch(shape,
    proportional = pnorm(outer(-eta, theta, `+`)),
    scale = { s <- exp(g * as.numeric(scale(aux))); pnorm(outer(-eta, theta, `+`) / s) },
    nonprop = { mult <- 1 + h * (seq_len(K1) - mean(seq_len(K1))) / max(K1 - 1, 1)
                pnorm(sweep(-outer(eta, mult), 2, theta, `+`)) })
  Cm <- t(apply(Cm, 1, cummax))                          # monotoniczność skumulowanych
  cbind(Cm, 1) - cbind(0, Cm)
}
