## =============================================================================
## run_all.R — pełny pipeline: dane -> symulacje -> walidacja -> aplikacja -> rysunki
## Uruchomienie (z katalogu projektu):   Rscript run_all.R
## Szybki test (mało replikacji):         SAE_NREP=10 SAE_NREP_SWEEP=5 SAE_NREP_REAL=5 Rscript run_all.R
## =============================================================================
## --- pakiety: doinstaluj brakujące ---------------------------------------------
pkgs <- c("sf", "spdep", "Matrix", "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble",
          "ggplot2", "patchwork", "ggrepel", "scales", "glmnet", "ranger", "mgcv", "readxl", "writexl",
          "future", "future.apply")
miss <- setdiff(pkgs, rownames(installed.packages()))
if (length(miss)) {
  message("Instaluję brakujące pakiety: ", paste(miss, collapse = ", "))
  install.packages(miss)
}

steps <- c(
  "R/01_build_data.R",        # baza powiatowa: geometria, kowariaty, zabory, sąsiedztwo
  "R/01b_build_izp.R",        # dane IZP: wskaźniki EHIS 2019 + SE, kowariaty IZP, spójność
  "R/01c_build_extra_covars.R", # struktura wieku, wykształcenie, gęstość NSP 2021 + dane do walidacji zewn.
  "R/04_run_simulation.R",    # 7 scenariuszy × 2 poziomy szumu × N_REP_SIM
  "R/05_signal_sweep.R",      # siła związku kowariaty–zmienna
  "R/05b_sensitivity_gdr.R",  # wrażliwość GDR / MRF-GAM na rho i priory
  "R/05c_multivariate_sim.R", # dezagregacja wielowymiarowa (MGDR) vs osobne GDR
  "R/05d_ordinal_sim.R",      # dezagregacja porządkowa (cały rozkład odpowiedzi) vs binarna
  "R/06_real_validation.R",   # pseudo-populacje z prawdziwych zmiennych powiatowych
  "R/06b_method_selection.R", # wybór metody: minimax regret (kryterium ustalone z góry)
  "R/06c_calibrate_xl_prior.R", # kalibracja c (ściąganie międzypoziomowe) na prawdziwych zmiennych
  "R/07_application_health.R",# dezagregacja wskaźników zdrowotnych IZP/EHIS
  "R/07b_ordinal_application.R", # wariant porządkowy dla x1–x5
  "R/08_figures_tables.R"     # rysunki i tabele
)
for (s in steps) {
  message("\n==== ", s, " ====  ", format(Sys.time()))
  t0 <- Sys.time()
  source(s, local = new.env(), echo = FALSE)
  message("     zakończono w ", format(round(Sys.time() - t0, 1)))
}
