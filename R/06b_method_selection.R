## =============================================================================
## 06b_method_selection.R — wybór metody wg kryterium ustalonego Z GÓRY
##
## Reguła (zapisana przed obliczeniami, nie dopasowywana do wyników):
##   1. Dla każdej „komórki oceny” c liczymy skill_m(c) = 1 − MSE_m / MSE_naive
##      (średnia po replikacjach). Komórki:
##        – symulacja, wariant realistyczny: scenariusze z SELECTION_SCEN × {exact, survey}
##        – walidacja rzeczywista: {public, izp} × 8 zmiennych × {exact, survey}
##   2. Żal (regret) metody m w komórce c: best(c) − skill_m(c), gdzie best(c) to najlepsza
##      metoda WYKONALNA w tej komórce (bez wyroczni).
##   3. Wybieramy metodę o najmniejszym MAKSYMALNYM żalu (minimax regret); remis rozstrzyga
##      średni żal. Dodatkowo raportujemy średni skill i odsetek komórek lepszych od naiwnej.
## S6 (błąd ekologiczny) jest wyłączony z góry: żadna metoda nie identyfikuje zależności
## wewnątrzwojewódzkiej z danych wojewódzkich – raportujemy go osobno jako test granic.
## Wynik: output/results/tab_method_selection.csv, output/results/selected_method.txt
## =============================================================================
source("R/00_config.R")

SELECTION_SCEN <- c("S1_base", "S2_urban", "S3_zabory", "S4_hotspot", "S5_fuzzy", "S7_gradient")
CANDIDATES <- setdiff(METHOD_ORDER, "ORACLE")

skill_cells <- list()
f_sim <- file.path(PATH_RES, "sim_metrics.csv")
if (file.exists(f_sim)) {
  s <- read_csv(f_sim, show_col_types = FALSE, guess_max = 1e6)
  if (!"design" %in% names(s)) s$design <- "realistic"
  skill_cells$sim <- s |> filter(design == "realistic", scen %in% SELECTION_SCEN) |>
    group_by(scen, noise, rep) |> mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |>
    group_by(cell = paste("sim", scen, noise), method) |> summarise(skill = mean(skill), .groups = "drop")
}
f_real <- file.path(PATH_RES, "real_metrics.csv")
if (file.exists(f_real)) {
  r <- read_csv(f_real, show_col_types = FALSE, guess_max = 1e6)
  skill_cells$real <- r |> group_by(covset, target, noise, rep) |>
    mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |>
    group_by(cell = paste("real", covset, target, noise), method) |>
    summarise(skill = mean(skill), .groups = "drop")
}
stopifnot(length(skill_cells) > 0)

cells <- bind_rows(skill_cells, .id = "source") |> filter(method %in% CANDIDATES) |>
  group_by(cell) |> mutate(regret = max(skill) - skill) |> ungroup()

sel <- cells |> group_by(method) |>
  summarise(max_regret = max(regret), mean_regret = mean(regret),
            mean_skill = mean(skill), worst_skill = min(skill),
            share_better_than_naive = mean(skill > 0),
            max_regret_sim = max(regret[source == "sim"], -Inf),
            max_regret_real = max(regret[source == "real"], -Inf),
            n_cells = n(), .groups = "drop") |>
  arrange(max_regret, mean_regret) |>
  mutate(rank = row_number(), label = METHOD_LABELS[method], .before = 1)

write_csv(sel, file.path(PATH_RES, "tab_method_selection.csv"))
winner <- sel$method[1]
writeLines(winner, file.path(PATH_RES, "selected_method.txt"))
print(sel |> mutate(across(where(is.numeric), ~ round(.x, 3))), n = 20, width = 200)
message("Metoda wybrana wg minimax regret: ", winner, " (", METHOD_LABELS[[winner]], ")")
