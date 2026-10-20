## =============================================================================
## 01c_build_extra_covars.R — dodatkowe kowariaty powiatowe i zewnętrzne „kotwice”
## Źródło: arkusze z danymi źródłowymi w IZP_kalkulator_warianty.xlsx (GUS: NSP 2021, BDL, NFZ).
##
## (A) Kowariaty strukturalne (data/processed/extra_covars.csv):
##   demo_65plus  – odsetek ludności 65+ (NSP 2021; arkusz 2_1_1_old / 1_5)
##   demo_u20     – odsetek ludności poniżej 20 lat (NSP 2021)
##   ldens_pop    – log gęstości zaludnienia: ludność NSP 2021 / powierzchnia (zamiast wyborców/km2)
##   edu_high     – odsetek z wykształceniem wyższym wśród osób 20+ (NSP 2021)
##   edu_low      – odsetek z wykształceniem co najwyżej podstawowym (NSP 2021; = izp_2_1_1)
##   migrants     – odsetek migrantów (NSP 2021, arkusz 2_8_1)
##   fertility    – współczynnik dzietności (średnia 2021–2023, arkusz 1_7_1)
##   doctors      – lekarze na 10 tys. ludności, 2019 (arkusz 5_1_1)
##   nurses       – pielęgniarki i położne na 10 tys., 2019 (arkusz 5_1_2)
##   Uwaga: płeć – w kalkulatorze brak struktury wg płci na poziomie powiatów; można dołożyć
##   własny plik (HEALTH_COVARS_FILE w 07) albo eksport z BDL (kolumny: pow + zmienne).
##
## (B) Zmienne zdrowotne do WALIDACJI ZEWNĘTRZNEJ (NIE są kowariatami!)
##   (data/processed/external_proxies.csv: pow, target, proxy, value, sign, same_scale)
##   x1 zła samoocena zdrowia  <-> trwanie życia (średnia K i M, 2023), oczekiwany znak ujemny
##   x5 codzienne palenie       <-> hospitalizacje POChP J44 (2023), znak dodatni
##   x2 nadwaga/otyłość         <-> YLL cukrzyca i choroby nerek (2023), znak dodatni
##   Własne dane (np. NFZ: zgłaszalność na mammografię/cytologię, szczepienia przeciw grypie)
##   dołóż jako plik external_validation.csv w data/raw (kolumny: pow, target, value,
##   opcjonalnie sign, same_scale = TRUE gdy to ten sam wskaźnik w tej samej skali, np. %).
## =============================================================================
source("R/00_config.R")
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))

F_KALK <- find_input("IZP[ _]?kalkulator_warianty\\.xlsx$")
pow_base <- read_csv(file.path(PATH_PROC, "powiaty.csv"), show_col_types = FALSE) |>
  select(pow, woj, area_km2)

## Wczytanie arkusza: wiersze powiatów (kod TERYT 7-cyfrowy kończący się na 000, bez województw)
read_sheet <- function(sheet) {
  d <- as.data.frame(read_excel(F_KALK, sheet = sheet, col_names = FALSE, .name_repair = "minimal"))
  kod <- as.character(d[[1]])
  ok  <- grepl("^[0-9]{4}000$", kod) & !grepl("00000$", kod)
  out <- d[ok, , drop = FALSE]
  out$pow <- substr(kod[ok], 1, 4)
  out
}
col_num <- function(d, j) num(d[[j]])

## --- (A) kowariaty ---------------------------------------------------------------------------
ed  <- read_sheet("2_1_1_old")                     # NSP 2021: grupy wieku 20+ × wykształcenie
tot_cols  <- seq(3, 30, by = 3)                     # ogółem w grupach 20-24, …, 65+
high_cols <- tot_cols + 1                           # wyższe
adult <- rowSums(sapply(tot_cols, function(j) col_num(ed, j)))
high  <- rowSums(sapply(high_cols, function(j) col_num(ed, j)))
age65 <- col_num(ed, 30)
pop5  <- read_sheet("1_5")                          # liczba ludności NSP 2021 (kol. 4)
popt  <- col_num(pop5, 4)[match(ed$pow, pop5$pow)]
low   <- read_sheet("2_1_1")
mig   <- read_sheet("2_8_1")
fer   <- read_sheet("1_7_1")
doc   <- read_sheet("5_1_1")
nur   <- read_sheet("5_1_2")

extra <- tibble(pow = ed$pow,
                demo_65plus = 100 * age65 / popt,
                demo_u20    = 100 * (1 - adult / popt),
                edu_high    = 100 * high / adult,
                pop_nsp     = popt) |>
  left_join(tibble(pow = low$pow, edu_low = col_num(low, ncol(low) - 1)), by = "pow") |>
  left_join(tibble(pow = mig$pow, migrants = col_num(mig, 5)), by = "pow") |>
  left_join(tibble(pow = fer$pow, fertility = col_num(fer, 3)), by = "pow") |>
  left_join(tibble(pow = doc$pow, doctors = col_num(doc, 3)), by = "pow") |>
  left_join(tibble(pow = nur$pow, nurses = col_num(nur, 3)), by = "pow") |>
  right_join(pow_base, by = "pow") |>
  mutate(ldens_pop = log(pop_nsp / area_km2)) |>
  select(pow, demo_65plus, demo_u20, ldens_pop, edu_high, edu_low, migrants, fertility, doctors, nurses)

EXTRA_COVARS <- setdiff(names(extra), "pow")
miss <- sapply(extra[EXTRA_COVARS], function(v) sum(!is.finite(v)))
print(miss)
for (v in EXTRA_COVARS) extra[[v]][!is.finite(extra[[v]])] <- median(extra[[v]], na.rm = TRUE)
stopifnot(nrow(extra) == nrow(pow_base))
write_csv(extra, file.path(PATH_PROC, "extra_covars.csv"))
print(summary(extra[EXTRA_COVARS]))

## --- (B) zewnętrzne zmienne zdrowotne do walidacji -------------------------------------------
le_m <- read_sheet("1_3_2"); le_k <- read_sheet("1_3_3")
copd <- read_sheet("3_4_3"); yll_dm <- read_sheet("1_2_10")
lastnum <- function(d) col_num(d, max(which(sapply(d, function(x) any(is.finite(num(x)))) &
                                       seq_along(d) > 2 & names(d) != "pow")))
ext <- bind_rows(
  tibble(pow = le_m$pow, target = "x1", proxy = "Life expectancy at birth, 2023 (mean of men and women)",
         value = (lastnum(le_m) + lastnum(le_k)[match(le_m$pow, le_k$pow)]) / 2, sign = -1),
  tibble(pow = copd$pow, target = "x5", proxy = "COPD (J44) emergency hospitalisations per 10,000, 2023",
         value = lastnum(copd), sign = 1),
  tibble(pow = yll_dm$pow, target = "x2", proxy = "YLL, diabetes and kidney diseases, 2023",
         value = lastnum(yll_dm), sign = 1)) |>
  mutate(same_scale = FALSE, source = "IZP calculator")

F_EXT <- find_input("external_validation\\.csv$", required = FALSE)
if (!is.na(F_EXT)) {
  ue <- read_csv(F_EXT, show_col_types = FALSE, col_types = cols(pow = "c")) |>
    mutate(pow = sprintf("%04d", as.integer(pow)))
  if (!"sign" %in% names(ue)) ue$sign <- 1
  if (!"same_scale" %in% names(ue)) ue$same_scale <- FALSE
  if (!"proxy" %in% names(ue)) ue$proxy <- "user-supplied"
  ext <- bind_rows(ext, ue |> transmute(pow, target, proxy, value = as.numeric(value), sign,
                                        same_scale = as.logical(same_scale), source = basename(F_EXT)))
  message("Dołączono własne dane zewnętrzne: ", F_EXT)
}
ext <- ext |> filter(is.finite(value))
write_csv(ext, file.path(PATH_PROC, "external_proxies.csv"))
print(ext |> count(target, proxy))
