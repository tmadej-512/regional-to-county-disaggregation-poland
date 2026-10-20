## =============================================================================
## 01b_build_izp.R — dane projektu IZP: wskaźniki EHIS 2019 (woj.) + kowariaty IZP (powiaty)
## Wejście (data/raw/izp):
##   ekstra_zmienne.xlsx         – arkusz „Do dezagregacji v2” (wartości woj.) i „Opis zmiennych”
##   tablice_woj_2019_se.xlsx    – GUS, EHIS 2019: % i SE wg województw (kolejność TERYT)
##   IZP_kalkulator_warianty.xlsx– arkusz „2019_dane”: zmienne IZP dla woj. i powiatów
## Wynik (data/processed):
##   health_targets_2019.csv     – woj, target, id, nazwa, value, se
##   izp_covars_2019.csv         – pow + kowariaty IZP
##   output/results/tab_izp_coherence.csv – spójność zmiennych między poziomami
## UWAGA: plik ehis_data_wide_2019.xlsx ma zamienione kody woj. 06/08/10 (sortowanie po
##   janitor::clean_names: „lodzkie” < „lubelskie”); dlatego SE czytamy z tablicy źródłowej.
## =============================================================================
source("R/00_config.R")
## pliki wyszukiwane po nazwie (data/raw/izp, data/raw, data/) – zob. find_input() w 00_config.R
F_EKSTRA <- find_input("ekstra_zmienne\\.xlsx$")
F_SE     <- find_input("tablice_woj_2019_se\\.xlsx$")
F_KALK   <- find_input("IZP[ _]?kalkulator_warianty\\.xlsx$")
message("Pliki IZP:\n  ", F_EKSTRA, "\n  ", F_SE, "\n  ", F_KALK)
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
WOJ_CODES <- sprintf("%02d", seq(2, 32, 2))

## --- 1. Wartości docelowe i opisy ---------------------------------------------------------
tv <- read_excel(F_EKSTRA, sheet = "Do dezagregacji v2") |>
  mutate(`Kod TERYT` = as.numeric(`Kod TERYT`)) |> filter(`Kod TERYT` > 0) |>
  transmute(woj = sprintf("%02d", as.integer(`Kod TERYT` %/% 100000)), id = `Nr wskaznika`,
            target = paste0("x", `ID wskaznika`), value = Wartość)
opis <- read_excel(F_EKSTRA, sheet = "Opis zmiennych",
                   col_names = FALSE, .name_repair = "minimal")
opis <- tibble(id = opis[[1]], nazwa = opis[[2]]) |> filter(!is.na(id))
tv <- tv |> left_join(opis, by = "id")

## KOREKTA x3 (4_2_1 „Zagrożenia wynikające z diety”): w arkuszu „Do dezagregacji v2” wartości x3
## są równe odsetkowi osób jedzących OWOCE co najmniej raz dziennie (zgodność z tablicą GUS do 0,05 pp),
## a nie – jak w opisie – odsetkowi osób NIE jedzących WARZYW codziennie. Domyślnie zastępujemy je
## właściwą definicją: 100 − odsetek jedzących warzywa ≥ 1 raz dziennie (arkusz „4_2_1”).
## Aby użyć wartości oryginalnych: Sys.setenv(SAE_X3_FIX = "0").
X3_FIX <- Sys.getenv("SAE_X3_FIX", "1") == "1"
if (X3_FIX) {
  veg <- read_excel(F_EKSTRA, sheet = "4_2_1", col_names = FALSE, .name_repair = "minimal")
  veg <- tibble(kod = as.character(veg[[1]]), daily = num(veg[[4]])) |>
    filter(grepl("^\\d{7}$", kod), kod != "0000000") |>
    transmute(woj = substr(kod, 1, 2), x3_fixed = 100 - daily)
  stopifnot(nrow(veg) == 16)
  tv <- tv |> left_join(veg, by = "woj") |>
    mutate(value = if_else(target == "x3", x3_fixed, value)) |> select(-x3_fixed)
  message("x3 zastąpiony: 100 − odsetek jedzących warzywa codziennie (arkusz 4_2_1).")
}

## --- 2. Tablica SE (GUS EHIS 2019) -----------------------------------------------------------
se_raw <- read_excel(F_SE, col_names = FALSE,
                     .name_repair = "minimal")
hdr <- as.character(unlist(se_raw[5, ]))
stopifnot(str_detect(hdr[5], "Dolno"), str_detect(hdr[20], "Zachodnio"))   # kolejność TERYT
se_tab <- se_raw[9:nrow(se_raw), ]
names(se_tab) <- c("variable", "category", "metric", "PL", WOJ_CODES)
se_tab <- se_tab |>
  mutate(variable = str_squish(variable), category = str_squish(category)) |>
  fill(variable, category) |>
  filter(metric %in% c("%", "SE")) |>
  pivot_longer(all_of(WOJ_CODES), names_to = "woj", values_to = "val") |>
  mutate(val = num(val)) |>
  select(-PL) |>
  pivot_wider(names_from = metric, values_from = val) |>
  rename(pct = `%`, se = SE)

## efektywna liczebność próby (z SE pojedynczych kategorii): n_eff = p(1-p)/SE^2
neff_var <- se_tab |> filter(pct > 0.5, se > 0) |>
  mutate(neff = (pct / 100) * (1 - pct / 100) / (se / 100)^2) |>
  group_by(variable, woj) |> summarise(neff = median(neff), .groups = "drop")
neff_all <- neff_var |> group_by(woj) |> summarise(neff = median(neff))

## pełne rozkłady kategorii pytań EHIS (kolejność kategorii jak w tablicy GUS) + n_eff –
## wejście dezagregacji porządkowej (02b_ordinal.R, 07b_ordinal_application.R)
cat_out <- se_tab |>
  group_by(variable, woj) |> mutate(k = row_number()) |> ungroup() |>
  left_join(neff_var, by = c("variable", "woj"))
write_csv(cat_out, file.path(PATH_PROC, "ehis_categories_2019.csv"))

## definicje wskaźników: zmienna EHIS, kategorie licznika, kategorie wyłączone z mianownika
defs <- tribble(
  ~target, ~var_pattern,                         ~num_pattern,                                 ~excl_pattern,       ~base_frac,
  "x1",   "^Ocena stanu zdrowia",                "^(złe|bardzo złe)$",                         NA,                  1,
  "x2",   "^Index masy ciała",                   "^(nadwaga|otyłość)$",                        NA,                  1,
  "x3",   "^Częstość spożywania warzyw",         "^(4-6|1-3|rzadziej|nigdy)",                  NA,                  1,
  "x4",   "^Częstość picia alkoholu",            "^(codziennie|5-6|3-4|1-2 dni)",              NA,                  1,
  "x5",   "^Palenie tytoniu",                    "^codziennie$",                               NA,                  1,
  "x6",   "^Opóźnienia w dostępie",              "^tak$",                                      "^nie potrzebowali$",1,
  "x7",   "^Możliwość sfinansowania opieki med", "^deklarowali",                               "^nie potrzebowali$",1,
  "x8",   "^Możliwość sfinansowania leków",      "^deklarowali",                               "^nie potrzebowali$",1,
  "x9",   NA, NA, NA, 0.5,      # mammografia – kobiety 15+ (brak w tablicy SE)
  "x10",  NA, NA, NA, 0.5,      # cytologia – kobiety 15+
  "x11",  NA, NA, NA, 1         # szczepienia przeciw grypie
)

se_out <- map_dfr(seq_len(nrow(defs)), function(k) {
  d <- defs[k, ]
  tvk <- tv |> filter(target == d$target)
  if (is.na(d$var_pattern)) {
    nf <- neff_all |> mutate(neff = neff * d$base_frac)
    return(tvk |> left_join(nf, by = "woj") |>
             mutate(recon = NA_real_, se = 100 * sqrt((value / 100) * (1 - value / 100) / neff),
                    se_method = "n_eff (mediana EHIS)"))
  }
  s <- se_tab |> filter(str_detect(variable, d$var_pattern))
  numr <- s |> filter(str_detect(category, d$num_pattern)) |> group_by(woj) |> summarise(p = sum(pct))
  excl <- if (is.na(d$excl_pattern)) tibble(woj = WOJ_CODES, q = 0) else
    s |> filter(str_detect(category, d$excl_pattern)) |> group_by(woj) |> summarise(q = sum(pct))
  nf <- neff_var |> filter(str_detect(variable, d$var_pattern))
  tvk |> left_join(numr, by = "woj") |> left_join(excl, by = "woj") |> left_join(nf, by = "woj") |>
    mutate(recon = 100 * p / (100 - q),                        # odtworzenie wartości z tablicy
           neff = neff * (1 - q / 100) * d$base_frac,
           se = 100 * sqrt((value / 100) * (1 - value / 100) / neff),
           se_method = "n_eff (kategorie zmiennej)") |>
    select(-p, -q, -variable)
})
chk <- se_out |> group_by(target) |> summarise(max_diff_recon = max(abs(recon - value)),
                                               mean_se = mean(se), mean_value = mean(value))
print(chk)
write_csv(chk, file.path(PATH_RES, "tab_health_targets_check.csv"))
health <- se_out |> select(woj, target, id, nazwa, value, se, se_method) |> arrange(target, woj)
write_csv(health, file.path(PATH_PROC, "health_targets_2019.csv"))

## --- 3. Kowariaty IZP 2019 (powiaty) ------------------------------------------------------------
raw <- read_excel(F_KALK, sheet = "2019_dane",
                  col_names = FALSE, .name_repair = "minimal")
ids <- as.character(unlist(raw[8, ])); nms <- as.character(unlist(raw[9, ]))
dat <- raw[10:nrow(raw), ]
names(dat) <- c("kod", "nazwa", "rodzaj", ids[-(1:3)])
dat <- dat |> mutate(across(-(1:3), num))
pw <- dat |> filter(rodzaj == "Powiat") |> mutate(pow = substr(kod, 1, 4))
wj <- dat |> filter(rodzaj == "Województwo") |> mutate(woj = substr(kod, 1, 2))

pow_base <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE)
stopifnot(setequal(pw$pow, pow_base$pow))

cand <- ids[-(1:3)]
avail <- map_dfr(cand, function(v) tibble(id = v,
  nazwa = nms[match(v, ids)],
  pow_ok = mean(!is.na(pw[[v]]) & pw[[v]] != 0),
  woj_ok = mean(!is.na(wj[[v]]) & wj[[v]] != 0)))
## kowariaty: dostępne dla >= 95% powiatów; bez zmiennych opartych na b. małych liczbach (zgony niemowląt)
IZP_COVARS <- avail |> filter(pow_ok >= 0.95, woj_ok == 1, id != "1_3_1") |> pull(id)

izp <- pw |> select(pow, all_of(IZP_COVARS))
for (v in IZP_COVARS) izp[[v]][is.na(izp[[v]]) | izp[[v]] == 0] <- median(izp[[v]], na.rm = TRUE)
names(izp) <- c("pow", paste0("izp_", IZP_COVARS))
write_csv(izp, file.path(PATH_PROC, "izp_covars_2019.csv"))

## spójność: oficjalna wartość woj. vs agregat powiatów ważony ludnością
coh <- map_dfr(IZP_COVARS, function(v) {
  a <- izp |> select(pow, x = paste0("izp_", v)) |> left_join(pow_base |> select(pow, woj, w), by = "pow") |>
    group_by(woj) |> summarise(agg = sum(w * x))
  o <- wj |> select(woj, off = all_of(v))
  j <- left_join(a, o, by = "woj")
  tibble(id = v, nazwa = nms[match(v, ids)], mean_woj_official = mean(j$off), mean_woj_aggregated = mean(j$agg),
         ratio = mean(j$off) / mean(j$agg), cor = cor(j$off, j$agg),
         max_rel_diff = max(abs(j$off - j$agg) / abs(j$agg)))
})
write_csv(coh, file.path(PATH_RES, "tab_izp_coherence.csv"))
print(coh |> mutate(across(where(is.numeric), ~ round(.x, 3))), width = 200)
message("Kowariaty IZP: ", paste(IZP_COVARS, collapse = ", "))
