## =============================================================================
## 00_config.R — wspólna konfiguracja projektu
## Dezagregacja wskaźników z poziomu województw (NUTS-2) do powiatów (LAU-1)
## =============================================================================

IS_WINDOWS <- .Platform$OS.type == "windows"
## Linux/macOS: wymuszamy UTF-8 (polskie znaki w nazwach kolumn PKW/GUS).
## Windows: R >= 4.2 działa natywnie w UTF-8 – nic nie ustawiamy.
if (!IS_WINDOWS) invisible(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")))
suppressPackageStartupMessages({
  library(sf)
  library(spdep)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(glmnet)
  library(ranger)
  library(mgcv)
  library(readxl)
  library(writexl)
  library(parallel)
})

options(dplyr.summarise.inform = FALSE, warn = 1)
sf_use_s2(FALSE)

## Ochrona przed maskowaniem: jeśli w sesji są załadowane pakiety z funkcjami o tych samych
## nazwach (np. MASS::select, raster::select, stats::filter, plyr::rename/summarise/mutate),
## skrypty i tak użyją wersji z dplyr/tidyr/purrr.
select    <- dplyr::select;    filter    <- dplyr::filter;    rename  <- dplyr::rename
mutate    <- dplyr::mutate;    summarise <- dplyr::summarise; arrange <- dplyr::arrange
count     <- dplyr::count;     first     <- dplyr::first;     lag     <- dplyr::lag
recode    <- dplyr::recode;    map       <- purrr::map;       extract <- tidyr::extract

## --- ścieżki ----------------------------------------------------------------
## Uruchamiaj skrypty z katalogu głównego projektu (tam, gdzie run_all.R).
PATH_ROOT  <- normalizePath(".", mustWork = TRUE)
PATH_RAW   <- file.path(PATH_ROOT, "data", "raw")
PATH_PROC  <- file.path(PATH_ROOT, "data", "processed")
PATH_OUT   <- file.path(PATH_ROOT, "output")
PATH_RES   <- file.path(PATH_OUT, "results")
PATH_FIG   <- file.path(PATH_OUT, "figures")
for (p in c(PATH_PROC, PATH_RES, PATH_FIG)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

## --- wyszukiwanie plików wejściowych ------------------------------------------
## Szuka pliku po wzorcu nazwy (bez względu na wielkość liter i ewentualne prefiksy
## w nazwie) kolejno w: data/raw/izp, data/raw, data/ (także podkatalogi).
## Dzięki temu działa zarówno układ z README, jak i dotychczasowy układ projektu IZP
## (np. data/IZP kalkulator_warianty.xlsx).
find_input <- function(pattern, required = TRUE) {
  dirs <- c(file.path(PATH_RAW, "izp"), PATH_RAW, file.path(PATH_ROOT, "data"))
  for (d in dirs) {
    if (!dir.exists(d)) next
    f <- list.files(d, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
    f <- f[!grepl("^~\\$", basename(f))]                 # pliki blokady Excela
    if (length(f)) return(f[1])
  }
  f <- list.files(file.path(PATH_ROOT, "data"), pattern = pattern, ignore.case = TRUE,
                  full.names = TRUE, recursive = TRUE)
  f <- f[!grepl("^~\\$", basename(f)) & !grepl("processed", f)]
  if (length(f)) return(f[1])
  if (required) stop("Nie znaleziono pliku pasującego do wzorca '", pattern, "'.\n",
                     "Szukano w: ", paste(dirs, collapse = "; "), " (oraz podkatalogach data/).\n",
                     "Katalog roboczy: ", PATH_ROOT, call. = FALSE)
  NA_character_
}

## --- parametry obliczeń -----------------------------------------------------
SEED        <- 20261007
N_CORES     <- as.integer(Sys.getenv("SAE_CORES", max(1L, parallel::detectCores() - 1L)))
N_REP_SIM   <- as.integer(Sys.getenv("SAE_NREP", "200"))   # replikacje na komórkę
N_REP_SWEEP <- as.integer(Sys.getenv("SAE_NREP_SWEEP", "100"))
N_REP_REAL  <- as.integer(Sys.getenv("SAE_NREP_REAL", "100"))
RF_TREES    <- 300L
## Prior dla wariancji wewnątrzwojewódzkiej w GDR: mediana stosunku wariancji
## wewnątrz/między województwami dla zmiennych referencyjnych (01_build_data.R, tab. ratio)
TAU_RATIO   <- 2
TAU_SD      <- 1
CAR_RHO     <- 0.95   # autokorelacja przestrzenna efektu CAR/MRF (analiza wrażliwości: 05b)
## Ściąganie międzypoziomowe (GDR_XL, GP_XL): jaki ułamek c zależności między województwami
## działa wewnątrz nich. Prior c ~ (średnia XL_C0, odch. std. XL_SD); domyślnie jak Uniform(0, 1).
## XL_C0 = 1 i XL_SD = 0 odtwarza zwykły GDR.
XL_C0       <- as.numeric(Sys.getenv("SAE_XL_C0", "0.5"))
XL_SD       <- as.numeric(Sys.getenv("SAE_XL_SD", format(sqrt(1 / 12))))
## Opcja: wariancja prioru współczynników dzielona przez liczbę kowariat (łączna siła sygnału
## nie rośnie z długością listy zmiennych). Domyślnie wyłączona, by nie zmieniać GDR.
KAPPA_PSCALE <- Sys.getenv("SAE_KAPPA_PSCALE", "0") == "1"
FH_K        <- Sys.getenv("SAE_FH_K", "3")   # liczba składowych głównych w FH: liczba lub "aic"
TUNE_RF     <- Sys.getenv("SAE_TUNE_RF", "0") == "1"   # strojenie RF (leave-one-region-out); wolne
## Diagnostyka „czy kowariaty wyjaśniają różnice MIĘDZY województwami”: Q2 regionalne
## (R^2 predykcyjne leave-one-region-out regresji z na 3 składowe główne Xbar).
## Q2_MIN – próg ustalony z góry: poniżej niego dezagregacja z kowariatami nie ma podstaw.
Q2_MIN      <- 0.2
GP_RANGES_KM <- c(50, 100, 200, 400)
URBAN_BW_KM  <- 30    # S2: skala zaniku „potencjału miejskiego” z odległością od miast (km)   # siatka zasięgów procesu gaussowskiego (GP-GAM)
CRS_PL      <- 2180   # PUWG 1992

## --- Q2 regionalne: R^2 predykcyjne (LOO po województwach) regresji z ~ PC1..PCk(Xbar) --------
regional_q2 <- function(z, Xbar, k = 3) {
  pc <- prcomp(Xbar, center = TRUE, scale. = TRUE)
  k <- min(k, ncol(pc$x), length(z) - 3)
  Z <- cbind(1, pc$x[, seq_len(k), drop = FALSE])
  H <- Z %*% solve(crossprod(Z), t(Z))
  e <- (z - H %*% z) / (1 - diag(H))                 # reszty leave-one-out (PRESS)
  1 - sum(e^2) / sum((z - mean(z))^2)
}

## --- obliczenia równoległe (Linux/macOS: fork; Windows: future.apply) -----------
## mclapply nie działa równolegle na Windows – tam używamy future.apply (sesje w tle),
## a jeśli pakietu brak, liczymy sekwencyjnie.
par_lapply <- function(X, FUN) {
  if (N_CORES <= 1) return(lapply(X, FUN))
  use_future <- IS_WINDOWS || Sys.getenv("SAE_FUTURE") == "1"
  if (!use_future) return(parallel::mclapply(X, FUN, mc.cores = N_CORES, mc.preschedule = TRUE))
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    message("Brak pakietu future.apply – liczę sekwencyjnie (install.packages('future.apply') przyspieszy).")
    return(lapply(X, FUN))
  }
  options(future.globals.maxSize = 4 * 1024^3)
  old <- future::plan(future::multisession, workers = N_CORES)
  on.exit(future::plan(old), add = TRUE)
  future.apply::future_lapply(X, FUN, future.seed = NULL,
    future.packages = c("Matrix", "dplyr", "tidyr", "purrr", "tibble", "glmnet", "ranger", "spdep"))
}

## --- kowariaty dostępne dla metod (poziom powiatu, agregowalne do woj.) ------
## urb      – odsetek wyborców w obwodach miejskich (proxy urbanizacji)
## ldens    – log gęstości zaludnienia (uprawnieni / km2)
## city     – miasto na prawach powiatu (0/1)
## ldcap    – log odległości od stolicy województwa (km)
## unemp19  – stopa bezrobocia rejestrowanego 2019 (BDL)
## turnout19– frekwencja w wyborach do Sejmu 2019 (proxy kapitału społecznego)
COVARS <- c("urb", "ldens", "city", "ldcap", "unemp19", "turnout19")

## --- zestawy kowariat (06_real_validation.R, 07_application_health.R) ----------------------
## public – 6 zmiennych publicznych (jak w symulacji)
## izp    – 4 strukturalne + 12 zmiennych IZP 2019 (01b_build_izp.R)
## core   – KRÓTKA LISTA USTALONA Z GÓRY (rekomendacja): urbanizacja, gęstość zaludnienia,
##          miasto na prawach powiatu, struktura wieku (65+, <20), wykształcenie wyższe i
##          co najwyżej podstawowe, wynagrodzenia, bezrobocie – mechanizmy działające tak samo
##          w powiecie i w województwie (skład ludności, rynek pracy, miejskość)
## wide   – wszystkie dostępne zmienne (test długiej listy ze ściąganiem współczynników)
## Dodatkowe zmienne (struktura wieku, wykształcenie, gęstość NSP, migranci, dzietność, kadry
## medyczne) buduje 01c_build_extra_covars.R -> data/processed/extra_covars.csv.
CORE_COVARS <- c("urb", "city", "ldens_pop", "demo_65plus", "demo_u20", "edu_high", "edu_low",
                 "izp_2_2_1", "izp_2_2_3")
make_covsets <- function(df) {
  nm  <- names(df)
  izp <- grep("^izp_", nm, value = TRUE)
  ext <- intersect(c("demo_65plus", "demo_u20", "ldens_pop", "edu_high", "migrants", "fertility",
                     "doctors", "nurses"), nm)              # edu_low pominięte: to samo co izp_2_1_1
  sets <- list(public = COVARS,
               izp    = c("urb", "ldens", "city", "ldcap", izp),
               core   = intersect(CORE_COVARS, nm),
               wide   = unique(c(COVARS, "city", izp, ext)))
  if (!all(c("demo_65plus", "edu_high") %in% nm)) sets$core <- NULL   # brak 01c -> bez core
  sets$wide <- setdiff(sets$wide, if (length(ext)) "ldens" else character())  # jedna gęstość
  sets
}
COVSET_LABELS <- c(public = "Public covariates (6)", izp = "IZP + structural (16)",
                   core = "Pre-specified short list (9)", wide = "All available covariates")

## --- etykiety (angielskie – rysunki i tabele do artykułu) -------------------
METHOD_LABELS <- c(
  NAIVE     = "Naive (regional value)",
  PYCNO     = "Pycnophylactic (Tobler)",
  MACRO     = "Macro-first ridge",
  MACRO_B   = "Macro-first + benchmarking",
  MACRO_BM  = "Macro-first + mult. raking",
  MINI_LM   = "Mini-first LM",
  MINI_RF   = "Mini-first RF",
  MINI_RF_OOB = "Mini-first RF (out-of-bag)",
  MINI_RF_B = "Mini-first RF + benchmarking",
  FH_B      = "Fay-Herriot synthetic + benchmarking",
  GDR       = "Bayesian disaggregation (CAR + X)",
  GDR_GAM   = "MRF-GAM disaggregation (smooth X + MRF)",
  GP_GAM    = "GP-GAM disaggregation (smooth X + Gaussian process)",
  GDR_XL    = "Bayesian disaggregation, cross-level shrinkage (CAR)",
  GP_XL     = "Bayesian disaggregation, cross-level shrinkage (GP)",
  MGDR      = "Multivariate disaggregation (shared spatial factor)",
  ORACLE    = "Oracle (OLS on true values)"
)
METHOD_ORDER <- names(METHOD_LABELS)

SCEN_LABELS <- c(
  S1_base     = "S1 Baseline",
  S2_urban    = "S2 Urban-rural",
  S3_zabory   = "S3 Historical partitions",
  S4_hotspot  = "S4 Local hotspots",
  S5_fuzzy    = "S5 Weak correlations",
  S6_ecol     = "S6 Ecological fallacy",
  S7_gradient = "S7 Spatial gradient"
)

HEALTH_LABELS_EN <- c(
  x1 = "Bad or very bad self-rated health", x2 = "Overweight or obesity",
  x3 = "Vegetables not eaten daily", x4 = "Alcohol at least weekly", x5 = "Daily smoking",
  x6 = "Delayed care (waiting times)", x7 = "Unaffordable medical care",
  x8 = "Unaffordable prescribed medicines", x9 = "Mammography within 3 years",
  x10 = "Cervical cytology within 3 years", x11 = "Influenza vaccination")

REAL_LABELS_EN <- c(
  unemp23 = "Registered unemployment 2023", turnout23 = "Turnout, Sejm 2023",
  pis23 = "PiS vote share 2023", ko23 = "KO vote share 2023", konf23 = "Konfederacja vote share 2023",
  td23 = "Third Way vote share 2023", trz25 = "Trzaskowski, 2025 run-off",
  right25 = "Mentzen + Braun, 2025 1st round")

theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid.minor = element_blank(),
                  strip.text = element_text(face = "bold"),
                  legend.position = "bottom"))
