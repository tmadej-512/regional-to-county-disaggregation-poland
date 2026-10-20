## =============================================================================
## 02_methods.R — metody dezagregacji i miary jakości
##
## Notacja: i = 1..n powiaty, g = 1..m województwa, g(i) – województwo powiatu i.
##   w_i     – udział ludności powiatu w ludności województwa (sum_i∈g w_i = 1)
##   A       – macierz agregacji m×n, A[g,i] = w_i  (wartość woj. = A y)
##   X       – n×p macierz kowariat powiatowych (standaryzowanych)
##   Xbar    – A X : kowariaty wojewódzkie spójne z powiatowymi (średnie ważone)
##   z       – obserwowana wartość wojewódzka, z = A y + e,  e_g ~ N(0, se_g^2)
## Każda metoda zwraca listę: pred (n), opcjonalnie sd (n) – odch. std. predykcyjne.
## =============================================================================

## --- Przygotowanie struktury problemu (wywoływane raz) --------------------------
make_problem <- function(pow_df, W, covars = COVARS, rho = CAR_RHO, k_spline = 6) {
  n <- nrow(pow_df)
  g <- as.integer(factor(pow_df$woj))
  m <- max(g)
  A <- sparseMatrix(i = g, j = seq_len(n), x = pow_df$w, dims = c(m, n))
  Xraw <- as.matrix(pow_df[, covars])
  X <- scale(Xraw)
  Xbar <- as.matrix(A %*% X)

  ## Laplasjan grafu sąsiedztwa (do interpolacji piknofilaktycznej)
  Wd <- as.matrix(W)
  deg <- rowSums(Wd)
  L <- diag(deg) - Wd
  Qp <- L + diag(1e-3, n)                         # minimalna regularyzacja
  Qp_inv <- solve(Qp)
  At <- t(as.matrix(A))

  ## Właściwy CAR: Q = D - rho W ; R = Q^{-1} przeskalowane do średniej wariancji 1
  Qc <- diag(deg) - rho * Wd
  R <- solve(Qc)
  R <- R / mean(diag(R))
  Rchol <- chol(R)                                 # do losowania pól przestrzennych

  Am <- as.matrix(A)

  ## Gładkie efekty kowariat (MRF-GAM): dla każdej ciągłej kowariaty baza splajnu kubicznego
  ## (mgcv::smoothCon, bs = "cr", z warunkiem identyfikowalności) i jej macierz kary S.
  ## Postać modelu mieszanego GAM: część „falista” f_j = B_j b_j, b_j ~ N(0, S_j^+), więc
  ## Cov(f_j) = B_j S_j^+ B_j'. Część liniowa (przestrzeń zerowa kary) jest już w X.
  ## Sumujemy po kowariatach (wspólny parametr wygładzania) i skalujemy do średniej wariancji 1.
  Ssm <- matrix(0, n, n); n_smooth <- 0
  for (j in seq_len(ncol(X))) {
    xj <- X[, j]
    if (length(unique(xj)) < 2 * k_spline) next      # zmienne binarne/dyskretne – tylko liniowo
    sm <- mgcv::smoothCon(mgcv::s(xj, k = k_spline, bs = "cr"), data = data.frame(xj = xj),
                          absorb.cons = TRUE)[[1]]
    B <- sm$X; S <- sm$S[[1]]
    e <- eigen(S, symmetric = TRUE)
    pos <- e$values > max(e$values) * 1e-8
    Sinv <- e$vectors[, pos, drop = FALSE] %*% (t(e$vectors[, pos, drop = FALSE]) / e$values[pos])
    Ssm <- Ssm + B %*% Sinv %*% t(B)
    n_smooth <- n_smooth + 1
  }
  if (n_smooth > 0) Ssm <- Ssm / mean(diag(Ssm))

  ## Proces gaussowski na współrzędnych środków powiatów (GP-GAM): kowariancja Matérn 3/2
  ## (jak mgcv bs = "gp", m = 3), wariancja 1; zasięg wybierany z siatki GP_RANGES_KM
  ## przez wiarygodność brzegową w gdr_fit(). Jeśli brak współrzędnych – GP pomijany.
  gp <- list()
  if (all(c("x_km", "y_km") %in% names(pow_df))) {
    dd <- as.matrix(dist(cbind(pow_df$x_km, pow_df$y_km)))
    for (rg in GP_RANGES_KM) {
      a <- sqrt(3) * dd / rg
      Kg <- (1 + a) * exp(-a)
      gp[[as.character(rg)]] <- list(range = rg, RAt = Kg %*% At, ARAt = Am %*% Kg %*% At,
                                     diagR = rep(1, n))
    }
  }

  list(n = n, m = m, g = g, A = A, Am = Am, At = At, X = X, Xraw = Xraw, Xbar = Xbar,
       w = pow_df$w, pop = pow_df$pop, W = Wd, L = L,
       Qp_inv = Qp_inv, PQ_AT = Qp_inv %*% At, PQ_AQA = Am %*% Qp_inv %*% At,
       R = R, Rchol = Rchol, RAt = R %*% At, ARAt = Am %*% R %*% At, AAt = Am %*% At,
       diagR = diag(R), covars = covars, rho = rho,
       n_smooth = n_smooth, SsAt = Ssm %*% At, ASsAt = Am %*% Ssm %*% At, diagSs = diag(Ssm),
       gp = gp)
}

agg <- function(P, y) as.numeric(P$Am %*% y)
bench_add  <- function(P, yhat, target) yhat + (target - agg(P, yhat))[P$g]
bench_mult <- function(P, yhat, target) {
  a <- agg(P, yhat)
  if (any(yhat <= 0) || any(a <= 0) || any(target <= 0)) return(bench_add(P, yhat, target))
  yhat * (target / a)[P$g]
}

## --- M0: naiwna ---------------------------------------------------------------------
m_naive <- function(P, z, se) list(pred = z[P$g])

## --- M1: piknofilaktyczna (Tobler 1979), wersja zamknięta ------------------------------
## min y' L y  przy  A y = z   =>  y = Q^{-1} A' (A Q^{-1} A')^{-1} z
m_pycno <- function(P, z, se) {
  list(pred = as.numeric(P$PQ_AT %*% solve(P$PQ_AQA, z)))
}

## --- M2/M3: macro-first (regresja grzbietowa na 16 woj.) [+ benchmark] ----------------
fit_macro_ridge <- function(P, z) {
  ## LOOCV po lambdzie; glmnet standaryzuje wewnętrznie
  fit <- suppressWarnings(cv.glmnet(P$Xbar, z, alpha = 0, nfolds = P$m, grouped = FALSE,
                                    lambda = exp(seq(log(1e3), log(1e-3), length.out = 60)) * sd(z)))
  as.numeric(predict(fit, newx = P$X, s = "lambda.min"))
}
m_macro   <- function(P, z, se, yhat = NULL) list(pred = if (is.null(yhat)) fit_macro_ridge(P, z) else yhat)
m_macro_b <- function(P, z, se, yhat = NULL) {
  if (is.null(yhat)) yhat <- fit_macro_ridge(P, z)
  list(pred = bench_add(P, yhat, z))
}
m_macro_bm <- function(P, z, se, yhat = NULL) {
  if (is.null(yhat)) yhat <- fit_macro_ridge(P, z)
  list(pred = bench_mult(P, yhat, z))
}

## --- M4: mini-first LM („bag labelling”: y_i := z_g(i)) ---------------------------------
m_mini_lm <- function(P, z, se) {
  yl <- z[P$g]
  b <- qr.coef(qr(cbind(1, P$X)), yl)
  list(pred = as.numeric(cbind(1, P$X) %*% b))
}

## --- M5/M6: mini-first RF [+ benchmark] ---------------------------------------------------
## Zwraca predykcje in-sample (jak model końcowy caret na całym zbiorze) oraz out-of-bag.
fit_mini_rf <- function(P, z, seed = 1, tune = TUNE_RF) {
  df <- data.frame(y = z[P$g], P$X)
  p <- ncol(P$X)
  mtry <- max(1, floor(p / 3)); nodes <- 5
  if (tune) {
    ## strojenie leave-one-region-out: fold = całe województwo; kryterium = błąd zagregowanej
    ## predykcji dla pominiętego województwa względem jego wartości z (jedyna dostępna „prawda”)
    cfg <- expand_grid(mtry = unique(c(1, max(1, floor(p / 3)), p)), nodes = c(5, 25, 75))
    err <- sapply(seq_len(nrow(cfg)), function(i) {
      mean(sapply(seq_len(P$m), function(g) {
        tr <- P$g != g
        rf <- ranger(y ~ ., data = df[tr, ], num.trees = 150, mtry = cfg$mtry[i],
                     min.node.size = cfg$nodes[i], num.threads = 1, seed = seed)
        ph <- predict(rf, data = df[!tr, ], num.threads = 1)$predictions
        (sum(P$w[!tr] * ph) - z[g])^2
      }))
    })
    mtry <- cfg$mtry[which.min(err)]; nodes <- cfg$nodes[which.min(err)]
  }
  rf <- ranger(y ~ ., data = df, num.trees = RF_TREES, mtry = mtry, min.node.size = nodes,
               num.threads = 1, seed = seed)
  list(insample = predict(rf, data = df, num.threads = 1)$predictions,
       oob = rf$predictions, mtry = mtry, min.node.size = nodes)
}
m_mini_rf     <- function(P, z, se, fit = NULL) list(pred = (if (is.null(fit)) fit_mini_rf(P, z) else fit)$insample)
m_mini_rf_oob <- function(P, z, se, fit = NULL) list(pred = (if (is.null(fit)) fit_mini_rf(P, z) else fit)$oob)
m_mini_rf_b   <- function(P, z, se, fit = NULL) {
  if (is.null(fit)) fit <- fit_mini_rf(P, z)
  list(pred = bench_add(P, fit$insample, z))
}

## --- M7: Fay–Herriot + estymator syntetyczny + benchmark do EBLUP woj. -----------------------
## method = "REML" (estymacja) lub "ML" (porównanie modeli o różnych kowariatach przez AIC)
fh_fit <- function(Xd, z, se, method = c("REML", "ML")) {
  method <- match.arg(method)
  m <- length(z); p <- ncol(Xd); D <- se^2
  nll <- function(ls2) {
    V <- ls2 + D; Vi <- 1 / V
    XtViX <- crossprod(Xd, Xd * Vi)
    b <- solve(XtViX, crossprod(Xd, z * Vi))
    r <- z - Xd %*% b
    out <- 0.5 * (sum(log(V)) + sum(r^2 * Vi))
    if (method == "REML") out <- out + 0.5 * as.numeric(determinant(XtViX)$modulus)
    out
  }
  opt <- optimize(nll, c(1e-8, 10 * var(z) + 1e-6))
  s2 <- opt$minimum
  V <- s2 + D; Vi <- 1 / V
  b <- solve(crossprod(Xd, Xd * Vi), crossprod(Xd, z * Vi))
  gam <- s2 / V
  theta <- as.numeric(gam * z + (1 - gam) * Xd %*% b)
  list(beta = as.numeric(b), sigma2v = s2, theta = theta,
       loglik = -opt$objective - 0.5 * m * log(2 * pi), npar = p + 1)
}
m_fh_b <- function(P, z, se, k_max = 5, k_fixed = FH_K) {
  ## kowariaty jako składowe główne Xbar (16 obszarów nie uniesie wszystkich kowariat naraz).
  ## Domyślnie stałe k = 3 (FH_K); FH_K = "aic" – wybór k = 1..k_max kryterium AIC (ML).
  ## Uwaga: AIC na 16 obserwacjach okazał się niestabilny (gorsze wyniki w S4/S5 i na danych).
  pc <- prcomp(P$Xbar, center = TRUE, scale. = FALSE)
  se_ <- pmax(se, 1e-6)
  kmax <- min(k_max, ncol(pc$x), P$m - 3)
  if (identical(tolower(as.character(k_fixed)), "aic")) {
    ks <- seq_len(kmax)
    aic <- sapply(ks, function(k) {
      f <- fh_fit(cbind(1, pc$x[, 1:k, drop = FALSE]), z, se_, "ML")
      -2 * f$loglik + 2 * f$npar
    })
    k <- ks[which.min(aic)]
  } else k <- min(as.integer(k_fixed), kmax)
  Zw <- cbind(1, pc$x[, 1:k, drop = FALSE])
  Zp <- cbind(1, scale(P$X, center = pc$center, scale = FALSE) %*% pc$rotation[, 1:k, drop = FALSE])
  f <- fh_fit(Zw, z, se_, "REML")
  syn <- as.numeric(Zp %*% f$beta)
  list(pred = bench_add(P, syn, f$theta), k = k)
}

## --- M8: bayesowska dezagregacja gaussowska (GDR): kowariaty + CAR + nugget ----------------
## y = D beta + u,  beta ~ N(b0, B0),  u ~ N(0, tau2 [(1-phi) R + phi I])
## z = A y + e,     e ~ N(0, S)
## => posterior y | z jest gaussowski (kriging obszar-obszar z dryfem / przestrzenny Chow–Lin).
## Hiperparametry (tau2, phi, kappa) – MAP z marginalnej wiarygodności z (m wymiarów).
gdr_fit <- function(P, z, se, return_sd = TRUE, tau_ratio = TAU_RATIO, tau_sd = TAU_SD,
                    smooth = FALSE, spatial = c("car", "gp"), xl = NULL) {
  spatial <- match.arg(spatial)
  ## GP: dopasowanie dla każdego zasięgu z siatki, wybór najlepszej wartości funkcji celu (MAP)
  if (spatial == "gp") {
    if (!length(P$gp)) stop("Brak współrzędnych powiatów (x_km, y_km) – GP niedostępny.")
    fits <- lapply(P$gp, function(sp) gdr_core(P, z, se, return_sd, tau_ratio, tau_sd, smooth, sp, xl))
    best <- which.min(sapply(fits, `[[`, "objective"))
    out <- fits[[best]]; out$gp_range_km <- P$gp[[best]]$range
    return(out)
  }
  gdr_core(P, z, se, return_sd, tau_ratio, tau_sd, smooth,
           list(RAt = P$RAt, ARAt = P$ARAt, diagR = P$diagR), xl)
}

## Rdzeń: sp = składnik przestrzenny (RAt = R A', ARAt = A R A', diagR = diag(R))
gdr_core <- function(P, z, se, return_sd, tau_ratio, tau_sd, smooth, sp, xl = NULL) {
  ## smooth = TRUE -> MRF-GAM: dodatkowy składnik gładkich efektów kowariat
  ##   y = D beta + f(X) + u,  f ~ N(0, kappa_s var(z) Ss),  Ss = sum_j B_j S_j^+ B_j'
  ## (postać modelu mieszanego GAM z karą drugiej pochodnej; analogia mgcv:
  ##  y ~ x + s(x) + s(powiat, bs = "mrf"), ale dopasowana do 16 średnich wojewódzkich).
  ## xl = list(c0, sd) -> ściąganie MIĘDZYPOZIOMOWE (GDR_XL, GP_XL), zob. opis przy m_gdr_xl()
  smooth <- smooth && isTRUE(P$n_smooth > 0)
  m <- P$m; n <- P$n; p <- ncol(P$X)
  if (is.null(xl)) {
    D  <- cbind(1, P$X)
    AD <- cbind(1, P$Xbar)
  } else {
    Xb <- P$Xbar[P$g, , drop = FALSE]                  # część międzywojewódzka kowariat
    D  <- cbind(1, Xb, P$X - Xb)                       # + część wewnątrzwojewódzka (A (X - Xb) = 0)
    AD <- cbind(1, P$Xbar, matrix(0, m, p))
  }
  vz <- var(z)
  b0 <- c(mean(z), rep(0, ncol(D) - 1))
  psc <- if (isTRUE(KAPPA_PSCALE)) 1 / p else 1       # opcjonalnie: prior na łączną siłę sygnału
  S  <- diag(se^2, m)
  par0 <- c(log(tau_ratio * vz), 0, 0, if (smooth) 0)
  build <- function(th) {
    tau2 <- exp(th[1]); phi <- plogis(th[2]); kap <- exp(th[3])
    ks <- if (smooth) exp(th[4]) * vz else 0
    kv <- kap * vz * psc
    B0 <- if (is.null(xl)) diag(c(10 * vz + 1, rep(kv, p))) else {
      Sig <- matrix(c(1, xl$c0, xl$c0, xl$c0^2 + xl$sd^2), 2)   # Cov(beta_b, beta_w) na współczynnik
      as.matrix(Matrix::bdiag(10 * vz + 1, kronecker(Sig, diag(kv, p))))
    }
    V <- AD %*% B0 %*% t(AD) + tau2 * ((1 - phi) * sp$ARAt + phi * P$AAt) + S
    if (smooth) V <- V + ks * P$ASsAt
    list(tau2 = tau2, phi = phi, ks = ks, B0 = B0, V = (V + t(V)) / 2)
  }
  nlp <- function(th) {
    b <- build(th)
    cv <- tryCatch(chol(b$V), error = function(e) NULL)
    if (is.null(cv)) return(1e10)
    r <- z - AD %*% b0
    a <- backsolve(cv, r, transpose = TRUE)
    nll <- sum(log(diag(cv))) + 0.5 * sum(a^2)
    ## priory: log tau2 ~ N(log(tau_ratio * var(z)), tau_sd^2) – kalibracja na zmiennych
    ## referencyjnych (stosunek wariancji wewnątrz/między, zob. 01_build_data.R);
    ## logit phi ~ N(0, 1.5^2), log kappa ~ N(0, 1.5^2), log kappa_s ~ N(0, 1.5^2)
    nll + 0.5 * ((th[1] - log(tau_ratio * vz)) / tau_sd)^2 + 0.5 * (th[2] / 1.5)^2 +
      0.5 * (th[3] / 1.5)^2 + (if (smooth) 0.5 * (th[4] / 1.5)^2 else 0)
  }
  lo <- c(log(vz) - 8, -6, -8, if (smooth) -8); hi <- c(log(vz) + 8, 6, 6, if (smooth) 6)
  opt <- optim(par0, nlp, method = "L-BFGS-B", lower = lo, upper = hi)
  b <- build(opt$par)
  ## kowariancja y z z:  C A' = D B0 (AD)' + tau2[(1-phi) R A' + phi A'] (+ ks Ss A')
  CAt <- D %*% b$B0 %*% t(AD) + b$tau2 * ((1 - b$phi) * sp$RAt + b$phi * P$At)
  if (smooth) CAt <- CAt + b$ks * P$SsAt
  Vi_r <- solve(b$V, z - AD %*% b0)
  pred <- as.numeric(D %*% b0 + CAt %*% Vi_r)
  out <- list(pred = pred, tau2 = b$tau2, phi = b$phi, kappa_s = b$ks, objective = opt$value)
  if (return_sd) {
    diagC <- rowSums((D %*% b$B0) * D) + b$tau2 * ((1 - b$phi) * sp$diagR + b$phi)
    if (smooth) diagC <- diagC + b$ks * P$diagSs
    K <- t(solve(b$V, t(CAt)))
    v <- pmax(diagC - rowSums(K * CAt), 1e-12)
    if (!is.null(xl)) {
      ## przybliżenie gaussowskie liczy niepewność c względem PRIORU beta (kappa var(z));
      ## dodajemy składnik z oszacowanego nachylenia: Var(c) * (Xw beta_b_hat)^2 (metoda delta)
      bb <- as.numeric(b$B0 %*% t(AD) %*% Vi_r)[1 + seq_len(p)]
      v <- v + xl$sd^2 * as.numeric((P$X - P$Xbar[P$g, , drop = FALSE]) %*% bb)^2
    }
    out$sd <- sqrt(v)
  }
  out
}
m_gdr     <- function(P, z, se, tau_ratio = TAU_RATIO) gdr_fit(P, z, se, tau_ratio = tau_ratio)
m_gdr_gam <- function(P, z, se, tau_ratio = TAU_RATIO) gdr_fit(P, z, se, tau_ratio = tau_ratio, smooth = TRUE)
m_gp_gam  <- function(P, z, se, tau_ratio = TAU_RATIO)
  gdr_fit(P, z, se, tau_ratio = tau_ratio, smooth = TRUE, spatial = "gp")

## --- M8b: ściąganie międzypoziomowe (cross-level shrinkage) – GDR_XL, GP_XL -------------------
## Kowariaty dzielimy na część międzywojewódzką Xb = Xbar[g] i wewnątrzwojewódzką Xw = X - Xb:
##   y = b0 + Xb beta_b + Xw beta_w + u,   beta_w = c * beta_b + delta.
## 16 wartości woj. identyfikują tylko beta_b (A Xw = 0), więc o beta_w decyduje założenie o c –
## jaki ułamek zależności między województwami działa też wewnątrz nich. Zamiast przyjmować c = 1
## (jak GDR) dajemy mu rozkład o średniej XL_C0 i odch. std. XL_SD (domyślnie jak Uniform(0, 1):
## 0.5 i 0.29). W postaci gaussowskiej (dopasowanie momentów) to blokowy prior
##   Cov(beta_b, beta_w) = kappa var(z) [[1, c0], [c0, c0^2 + sd^2]] ⊗ I_p,
## więc wszystko zostaje w zamkniętej postaci GDR. Skutek: różnice wewnątrz województw są
## ściągnięte o czynnik c0, a niepewność co do przenoszenia zależności trafia do przedziałów.
## Wiarygodność brzegowa z jest taka sama jak w GDR (dane woj. nie mówią nic o c) – to uczciwe:
## c trzeba wziąć z wiedzy lub danych zewnętrznych (np. kalibracja na zmiennych referencyjnych).
XL_PRIOR <- function() list(c0 = XL_C0, sd = XL_SD)
m_gdr_xl <- function(P, z, se, tau_ratio = TAU_RATIO)
  gdr_fit(P, z, se, tau_ratio = tau_ratio, xl = XL_PRIOR())
m_gp_xl  <- function(P, z, se, tau_ratio = TAU_RATIO)
  gdr_fit(P, z, se, tau_ratio = tau_ratio, spatial = "gp",
          xl = list(c0 = if (exists("XL_C0_GP")) XL_C0_GP else XL_C0, sd = XL_SD))

## --- M9: wielowymiarowa dezagregacja gaussowska (MGDR) – K wskaźników jednocześnie ------------
## Każdy wskaźnik k (standaryzowany na poziomie woj.) na poziomie powiatów:
##   y_k = D beta_k + a_k s + u_k,   beta_k ~ N(0, B0),
##   s   ~ N(0, C)          – WSPÓLNE pole przestrzenne (ukryty czynnik, np. „profil zdrowotny”),
##   u_k ~ N(0, tau2 C)     – pole specyficzne dla wskaźnika,  C = (1-phi) R + phi I  (CAR + nugget)
##   z_k = A y_k + e_k,     e_k ~ N(0, diag(se_k^2))
## => Cov(y_k, y_l) = delta_kl D B0 D' + (a_k a_l + delta_kl tau2) C  (liniowy model koregionalizacji).
## Hiperparametry (tau2, phi, kappa, a_1..a_K) – MAP z łącznej wiarygodności K*m wartości woj.
## Zysk względem K osobnych GDR: (i) odfiltrowanie błędu sondażowego dzięki skorelowanym wskaźnikom,
## (ii) wspólny czynnik przestrzenny, (iii) hiperparametry z K*m zamiast m obserwacji.
## Przy braku błędu i identycznej strukturze pól zysk znika (autokrigeability) – raportujemy uczciwie.
mgdr_fit <- function(P, Z, SE = NULL, tau_ratio = TAU_RATIO, tau_sd = TAU_SD, return_sd = TRUE) {
  Z <- as.matrix(Z); K <- ncol(Z); m <- P$m; n <- P$n
  if (is.null(SE)) SE <- matrix(0, m, K)
  mu <- colMeans(Z); sdz <- apply(Z, 2, sd)
  Zs <- sweep(sweep(Z, 2, mu), 2, sdz, "/")
  Ss <- sweep(as.matrix(SE), 2, sdz, "/")
  D  <- cbind(1, P$X); AD <- cbind(1, P$Xbar); p1 <- ncol(D)
  r  <- as.numeric(Zs)                               # vec(Z): kolejne bloki = wskaźniki
  build <- function(th) {
    tau2 <- exp(th[1]); phi <- plogis(th[2]); kap <- exp(th[3]); a <- th[3 + seq_len(K)]
    B0 <- diag(c(10 + 1, rep(kap, p1 - 1)))
    ACA <- (1 - phi) * P$ARAt + phi * P$AAt
    V <- kronecker(diag(K), AD %*% B0 %*% t(AD)) + kronecker(tcrossprod(a) + tau2 * diag(K), ACA) +
      diag(as.numeric(Ss)^2)
    list(tau2 = tau2, phi = phi, a = a, B0 = B0, V = (V + t(V)) / 2)
  }
  nlp <- function(th) {
    b <- build(th)
    cv <- tryCatch(chol(b$V), error = function(e) NULL)
    if (is.null(cv)) return(1e10)
    x <- backsolve(cv, r, transpose = TRUE)
    sum(log(diag(cv))) + 0.5 * sum(x^2) +
      0.5 * ((th[1] - log(tau_ratio)) / tau_sd)^2 + 0.5 * (th[2] / 1.5)^2 + 0.5 * (th[3] / 1.5)^2 +
      0.5 * sum(th[3 + seq_len(K)]^2)               # a_k ~ N(0, 1) na skali standaryzowanej
  }
  ## start: ładunki z pierwszej składowej głównej reszt regionalnych
  pc1 <- prcomp(Zs)$rotation[, 1] * 0.5
  par0 <- c(log(tau_ratio), 0, 0, pc1)
  opt <- optim(par0, nlp, method = "L-BFGS-B",
               lower = c(-8, -6, -8, rep(-5, K)), upper = c(8, 6, 6, rep(5, K)))
  b <- build(opt$par)
  w <- matrix(solve(b$V, r), m, K)                   # V^{-1} r, kolumny = wskaźniki
  CAt <- (1 - b$phi) * P$RAt + b$phi * P$At          # C A'  (n × m)
  DBA <- D %*% b$B0 %*% t(AD)                        # D B0 (AD)'  (n × m)
  common <- as.numeric(w %*% b$a)                    # sum_l a_l w_l
  diagCC <- (1 - b$phi) * P$diagR + b$phi
  diagDB <- rowSums((D %*% b$B0) * D)
  pred <- sd <- matrix(NA_real_, n, K)
  Vi <- if (return_sd) chol2inv(chol(b$V)) else NULL
  for (k in seq_len(K)) {
    pred[, k] <- DBA %*% w[, k] + CAt %*% (b$a[k] * common + b$tau2 * w[, k])
    if (return_sd) {
      G <- do.call(cbind, lapply(seq_len(K), function(l)
        (k == l) * DBA + (b$a[k] * b$a[l] + (k == l) * b$tau2) * CAt))
      v <- diagDB + (b$a[k]^2 + b$tau2) * diagCC - rowSums((G %*% Vi) * G)
      sd[, k] <- sqrt(pmax(v, 1e-12)) * sdz[k]
    }
    pred[, k] <- mu[k] + sdz[k] * pred[, k]
  }
  colnames(pred) <- colnames(sd) <- colnames(Z)
  list(pred = pred, sd = if (return_sd) sd else NULL, loadings = setNames(b$a, colnames(Z)),
       tau2 = b$tau2, phi = b$phi,
       ## udział wspólnego czynnika w wariancji niewyjaśnionej kowariatami, dla każdego wskaźnika
       common_share = setNames(b$a^2 / (b$a^2 + b$tau2), colnames(Z)))
}

## --- Wyrocznia (niedostępna w praktyce): OLS na prawdziwych y --------------------------------
m_oracle <- function(P, y_true) {
  D <- cbind(1, P$X)
  list(pred = as.numeric(D %*% qr.coef(qr(D), y_true)))
}

## --- Uruchomienie wszystkich metod ------------------------------------------------------------
run_methods <- function(P, z, se, y_true = NULL, methods = METHOD_ORDER, seed = 1,
                        tau_ratio = TAU_RATIO) {
  out <- list()
  need_macro <- any(c("MACRO", "MACRO_B", "MACRO_BM") %in% methods)
  need_rf    <- any(c("MINI_RF", "MINI_RF_OOB", "MINI_RF_B") %in% methods)
  ym  <- if (need_macro) fit_macro_ridge(P, z) else NULL
  yrf <- if (need_rf) fit_mini_rf(P, z, seed) else NULL
  for (mth in methods) {
    t0 <- proc.time()[3]
    res <- switch(mth,
      NAIVE     = m_naive(P, z, se),
      PYCNO     = m_pycno(P, z, se),
      MACRO     = m_macro(P, z, se, ym),
      MACRO_B   = m_macro_b(P, z, se, ym),
      MACRO_BM  = m_macro_bm(P, z, se, ym),
      MINI_LM   = m_mini_lm(P, z, se),
      MINI_RF   = m_mini_rf(P, z, se, yrf),
      MINI_RF_OOB = m_mini_rf_oob(P, z, se, yrf),
      MINI_RF_B = m_mini_rf_b(P, z, se, yrf),
      FH_B      = m_fh_b(P, z, se),
      GDR       = m_gdr(P, z, se, tau_ratio),
      GDR_GAM   = m_gdr_gam(P, z, se, tau_ratio),
      GP_GAM    = m_gp_gam(P, z, se, tau_ratio),
      GDR_XL    = m_gdr_xl(P, z, se, tau_ratio),
      GP_XL     = m_gp_xl(P, z, se, tau_ratio),
      ORACLE    = if (!is.null(y_true)) m_oracle(P, y_true) else NULL)
    if (!is.null(res)) out[[mth]] <- res
  }
  out
}

## --- Miary jakości na poziomie powiatów (prawda znana) ----------------------------------------
metrics_pow <- function(P, y, pred, sd = NULL, z_obs = NULL, top_frac = 0.10) {
  e <- pred - y
  wm_y <- agg(P, y)[P$g]; wm_p <- agg(P, pred)[P$g]
  dy <- y - wm_y; dp <- pred - wm_p
  k <- ceiling(top_frac * P$n)
  top_t <- order(y, decreasing = TRUE)[1:k]
  top_p <- order(pred, decreasing = TRUE)[1:k]
  cw <- if (sd(dp) < 1e-10) 0 else cor(dy, dp)
  tibble(
    RMSE      = sqrt(mean(e^2)),
    MAE       = mean(abs(e)),
    Bias      = mean(e),
    MAPE      = 100 * mean(abs(e) / pmax(abs(y), 1e-6)),
    Cor       = if (sd(pred) < 1e-10) 0 else cor(y, pred),
    CorWithin = cw,                                     # korelacja odchyleń od średniej woj.
    Spearman  = if (sd(pred) < 1e-10) 0 else cor(y, pred, method = "spearman"),
    Top10     = length(intersect(top_t, top_p)) / k,    # trafność wskazania 10% powiatów o najwyższych wartościach
    ## nachylenie kalibracji odchyleń wewnątrzwojewódzkich: < 1 = estymator przesadza z różnicami
    ## w województwie (za duża amplituda), > 1 = je zaniża; dla metody naiwnej NA
    CalSlope  = if (sum(dp^2) < 1e-10) NA_real_ else sum(dy * dp) / sum(dp^2),
    WojRMSE   = sqrt(mean((agg(P, pred) - agg(P, y))^2)), # błąd po agregacji do woj. (dotychczasowa walidacja)
    Cover90   = if (is.null(sd)) NA_real_ else mean(abs(e) <= qnorm(0.95) * sd),
    Width90   = if (is.null(sd)) NA_real_ else mean(2 * qnorm(0.95) * sd)
  )
}

evaluate_all <- function(P, y, fits) {
  imap_dfr(fits, function(f, nm) metrics_pow(P, y, f$pred, f$sd) |> mutate(method = nm, .before = 1))
}
