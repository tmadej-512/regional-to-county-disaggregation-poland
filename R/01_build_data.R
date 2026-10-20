## =============================================================================
## 01_build_data.R — baza powiatowa: geometria, wagi, kowariaty, zmienne testowe
## Wynik: data/processed/powiaty.gpkg, powiaty.csv, woj.gpkg, W_queen.rds
## Źródła (data/raw):
##   powiaty.json, wojewodztwa.json  – granice PRG (GeoJSON, kody TERYT JPT_KOD_JE)
##   2019-sejm.csv, 2023-sejm.csv    – PKW, wyniki wyborów do Sejmu wg gmin
##   prezydent2025_obwody.csv        – PKW, wybory prezydenckie 2025 wg obwodów
##   bezrobocie_2014_2025.csv        – GUS BDL, stopa bezrobocia rejestrowanego
##   2023-01-01.csv                  – TERC (GUS), nazwy jednostek
## =============================================================================
source("R/00_config.R")

num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", gsub("\\s", "", as.character(x)))))
pad6 <- function(x) sprintf("%06d", as.integer(num(x)))

## --- 1. Geometria ------------------------------------------------------------
pow_sf <- st_read(find_input("powiaty\\.json$"), quiet = TRUE) |>
  transmute(pow = JPT_KOD_JE) |>
  st_make_valid() |>
  st_transform(CRS_PL) |>
  group_by(pow) |> summarise(.groups = "drop") |>
  mutate(woj = substr(pow, 1, 2)) |>
  arrange(pow)
stopifnot(nrow(pow_sf) == 380)

woj_sf <- pow_sf |> group_by(woj) |> summarise(.groups = "drop")

terc <- read_delim(find_input("2023-01-01\\.csv$"), delim = ";",
                   col_types = cols(.default = "c"), locale = locale(encoding = "UTF-8"))
names(terc)[1] <- "WOJ"
pow_names <- terc |> filter(!is.na(POW), is.na(GMI)) |>
  transmute(pow = paste0(WOJ, POW), nazwa = NAZWA)
woj_names <- terc |> filter(is.na(POW)) |> transmute(woj = WOJ, woj_nazwa = str_to_lower(NAZWA))

## --- 2. Sejm 2019 i 2023 (gminy -> powiaty) ------------------------------------
read_sejm <- function(f, year) {
  d <- read_delim(f, delim = ";", col_types = cols(.default = "c"),
                  locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
  names(d)[1] <- "teryt"
  d <- d |> filter(!is.na(teryt), teryt != "")
  col <- function(pattern) {
    j <- grep(pattern, names(d), ignore.case = TRUE)[1]
    if (is.na(j)) return(rep(0, nrow(d)))
    v <- num(d[[j]]); v[is.na(v)] <- 0; v
  }
  out <- tibble(
    pow      = substr(pad6(d$teryt), 1, 4),
    elig     = col("^Liczba wyborców uprawnionych"),
    issued   = col("^Liczba wyborców, którym wydano karty do głosowania$"),
    valid    = col("^Liczba głosów ważnych oddanych łącznie"),
    pis      = col("PRAWO I SPRAWIEDLIWO"),
    ko       = col("KOALICJA OBYWATELSKA"),
    konf     = col("KONFEDERACJA"),
    td       = col("TRZECIA DROGA|POLSKIE STRONNICTWO LUDOWE"),
    lew      = col("LEWIC")
  )
  # w 2023 Warszawa (146501..146519) -> 1465; statki/zagranica już odrzucone
  out |> group_by(pow) |> summarise(across(everything(), sum)) |>
    transmute(pow,
              !!paste0("elig", year)    := elig,
              !!paste0("turnout", year) := 100 * issued / elig,
              !!paste0("pis", year)     := 100 * pis / valid,
              !!paste0("ko", year)      := 100 * ko / valid,
              !!paste0("konf", year)    := 100 * konf / valid,
              !!paste0("td", year)      := 100 * td / valid,
              !!paste0("lew", year)     := 100 * lew / valid)
}
s19 <- read_sejm(find_input("2019-sejm\\.csv$"), "19")
s23 <- read_sejm(find_input("2023-sejm\\.csv$"), "23")
## Wałbrzych (0265) w 2019 r. raportowany razem? — sprawdzamy kompletność niżej

## --- 3. Prezydent 2025 (obwody -> powiaty): urbanizacja i II tura -----------------
p25 <- read_delim(find_input("prezydent2025_obwody\\.csv$"), delim = ";",
                  col_types = cols(.default = "c"), show_col_types = FALSE)
cn <- names(p25)
gcol <- function(pattern) { j <- grep(pattern, cn)[1]; v <- num(p25[[j]]); v[is.na(v)] <- 0; v }
p25s <- tibble(
  teryt_g = p25$`Teryt Gminy`,
  typ_obw = p25$`Typ obwodu`,
  typ_obs = p25$`Typ obszaru`,
  elig_r1 = gcol("^Liczba wyborców uprawnionych.*_r1$"),
  valid_r2 = gcol("^Liczba głosów ważnych oddanych łącznie na obu"),
  trz_r2  = gcol("^TRZASKOWSKI.*_r2$"),
  valid_r1 = gcol("^Liczba głosów ważnych oddanych łącznie na wszystkich kandydatów"),
  mentzen = gcol("^MENTZEN"), braun = gcol("^BRAUN")
) |>
  filter(!is.na(teryt_g), typ_obs %in% c("miasto", "wieś", "miasto i wieś", "dzielnica w m.st. Warszawa")) |>
  mutate(pow = substr(pad6(teryt_g), 1, 4),
         urb_w = case_when(typ_obs %in% c("miasto", "dzielnica w m.st. Warszawa") ~ 1,
                           typ_obs == "miasto i wieś" ~ 0.5, TRUE ~ 0))
urb <- p25s |> filter(typ_obw == "stały") |>
  group_by(pow) |> summarise(urb = 100 * sum(urb_w * elig_r1) / sum(elig_r1))
pres <- p25s |> group_by(pow) |>
  summarise(trz25 = 100 * sum(trz_r2) / sum(valid_r2),
            right25 = 100 * sum(mentzen + braun) / sum(valid_r1))

## --- 4. Bezrobocie (BDL) ---------------------------------------------------------
bez <- read_delim(find_input("bezrobocie_2014_2025\\.csv$"), delim = ";",
                  col_types = cols(.default = "c"), show_col_types = FALSE)
names(bez)[1:2] <- c("kod", "nazwa")
ycol <- function(y) num(bez[[grep(paste0(";", y, ";"), names(bez))[1]]])
bez <- tibble(kod = bez$kod, unemp19 = ycol(2019), unemp23 = ycol(2023)) |>
  filter(nchar(kod) == 7, substr(kod, 5, 7) == "000", substr(kod, 3, 4) != "00") |>
  transmute(pow = substr(kod, 1, 4), unemp19, unemp23)

## --- 5. Zabory (przybliżenie stanu z 1914 r. + Ziemie Odzyskane po 1945) -------
## R = zabór rosyjski, A = austriacki, P = pruski (ziemie w granicach II RP),
## ZO = ziemie niemieckie do 1945 r. Przypisanie na poziomie całych powiatów
## (powiaty mieszane – wg przeważającej części); celowo upraszczające.
zabor_default <- c("02" = "ZO", "04" = "P", "06" = "R", "08" = "ZO", "10" = "R",
                   "12" = "A", "14" = "R", "16" = "ZO", "18" = "A", "20" = "R",
                   "22" = "P", "24" = "P", "26" = "R", "28" = "ZO", "30" = "P", "32" = "ZO")
zabor_override <- c(
  # kujawsko-pomorskie: Kujawy wsch. i ziemia dobrzyńska – zabór rosyjski
  "0401" = "R", "0405" = "R", "0408" = "R", "0411" = "R", "0412" = "R", "0418" = "R", "0464" = "R",
  # małopolskie: północ – zabór rosyjski
  "1208" = "R", "1212" = "R", "1214" = "R",
  # pomorskie: Pomorze Zach. i Powiśle, Wolne Miasto Gdańsk – ZO
  "2201" = "ZO", "2203" = "ZO", "2204" = "ZO", "2207" = "ZO", "2208" = "ZO", "2209" = "ZO",
  "2210" = "ZO", "2212" = "ZO", "2216" = "ZO", "2261" = "ZO", "2263" = "ZO", "2264" = "ZO",
  # śląskie: Zagłębie i Częstochowa – R; Śląsk Cieszyński, Żywiecczyzna, Jaworzno – A;
  # zachodnia część Górnego Śląska (w Niemczech do 1945) – ZO
  "2401" = "R", "2404" = "R", "2406" = "R", "2409" = "R", "2416" = "R", "2464" = "R",
  "2465" = "R", "2475" = "R",
  "2402" = "A", "2403" = "A", "2417" = "A", "2461" = "A", "2468" = "A",
  "2405" = "ZO", "2411" = "ZO", "2462" = "ZO", "2466" = "ZO", "2478" = "ZO",
  # warmińsko-mazurskie: ziemia działdowska i lubawska – zabór pruski w II RP
  "2803" = "P", "2812" = "P",
  # wielkopolskie: wschód (Kalisz, Konin) – R; Krajna/Piła – ZO
  "3007" = "R", "3009" = "R", "3010" = "R", "3023" = "R", "3027" = "R", "3061" = "R", "3062" = "R",
  "3019" = "ZO", "3031" = "ZO"
)

## --- 6. Złożenie bazy --------------------------------------------------------------
capitals <- c("0264", "0461", "0663", "0861", "1061", "1261", "1465", "1661",
              "1863", "2061", "2261", "2469", "2661", "2862", "3064", "3262")

cent <- st_centroid(st_geometry(pow_sf))
xy <- st_coordinates(cent) / 1000                                  # km
cap_xy <- xy[match(capitals, pow_sf$pow), ]
woj_cap <- setNames(capitals, substr(capitals, 1, 2))
cap_i <- match(woj_cap[pow_sf$woj], pow_sf$pow)
dcap <- sqrt(rowSums((xy - xy[cap_i, ])^2))

pow <- pow_sf |>
  mutate(area_km2 = as.numeric(st_area(pow_sf)) / 1e6,
         x_km = xy[, 1], y_km = xy[, 2],
         dcap_km = dcap,
         city = as.integer(as.integer(substr(pow, 3, 4)) >= 61),
         capital = as.integer(pow %in% capitals),
         zabor = unname(ifelse(pow %in% names(zabor_override), zabor_override[pow],
                               zabor_default[woj]))) |>
  left_join(pow_names, by = "pow") |>
  left_join(woj_names, by = "woj") |>
  left_join(s19, by = "pow") |>
  left_join(s23, by = "pow") |>
  left_join(urb, by = "pow") |>
  left_join(pres, by = "pow") |>
  left_join(bez, by = "pow")

## Wałbrzych (0265): jeżeli w którymś źródle brak – imputacja wartością z powiatu wałbrzyskiego
fill_from <- function(df, target = "0265", donor = "0221") {
  vars <- names(df)[sapply(st_drop_geometry(df), is.numeric)]
  for (v in vars) if (is.na(df[[v]][df$pow == target])) df[[v]][df$pow == target] <- df[[v]][df$pow == donor]
  df
}
pow <- fill_from(pow)

pow <- pow |>
  mutate(pop = elig23,                                  # waga: dorośli mieszkańcy (uprawnieni 2023)
         ldens = log(pop / area_km2),
         ldcap = log1p(dcap_km),
         w = pop / ave(pop, woj, FUN = sum))            # udział w ludności województwa

miss <- colSums(is.na(st_drop_geometry(pow)))
if (any(miss > 0)) print(miss[miss > 0])
stopifnot(all(miss[c(COVARS, "pop", "turnout23", "pis23", "unemp23", "trz25")] == 0))

## --- 7. Sąsiedztwo (queen) -----------------------------------------------------------
nb <- poly2nb(pow, queen = TRUE, snap = 50)   # tolerancja 50 m (PUWG92)
stopifnot(all(card(nb) > 0))
W <- as(nb2mat(nb, style = "B", zero.policy = TRUE), "CsparseMatrix")

## --- 8. Zapis ------------------------------------------------------------------------------
st_write(pow, file.path(PATH_PROC, "powiaty.gpkg"), delete_dsn = TRUE, quiet = TRUE)
st_write(woj_sf |> left_join(woj_names, by = "woj"), file.path(PATH_PROC, "woj.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
write_csv(st_drop_geometry(pow), file.path(PATH_PROC, "powiaty.csv"))
saveRDS(W, file.path(PATH_PROC, "W_queen.rds"))

message(sprintf("Baza: %d powiatów, %d województw; średnia liczba sąsiadów %.2f",
                nrow(pow), n_distinct(pow$woj), mean(card(nb))))
print(table(pow$zabor))
print(summary(st_drop_geometry(pow)[, c(COVARS, "unemp23", "turnout23", "pis23", "ko23", "konf23", "trz25")]))

## --- 9. Kalibracja: stosunek wariancji wewnątrz- do międzywojewódzkiej -----------------------
ref_vars <- c("turnout19", "turnout23", "pis19", "pis23", "ko19", "ko23", "konf19", "konf23",
              "td23", "lew19", "lew23", "unemp19", "unemp23", "trz25", "right25")
pd <- st_drop_geometry(pow)
Xs <- scale(as.matrix(pd[, COVARS]))
Xbar_ref <- apply(Xs, 2, function(x) tapply(pd$w * x, pd$woj, sum))
fam <- function(v) sub("\\d+$", "", v)
ratio_tab <- map_dfr(ref_vars, function(v) {
  y <- pd[[v]]; zg <- tapply(pd$w * y, pd$woj, sum)
  within <- sum(pd$w * (y - zg[pd$woj])^2) / 16
  ## Q2 regionalne: kowariaty bez opóźnionej wersji tej samej zmiennej (np. turnout19 dla turnout23)
  keep <- fam(COVARS) != fam(v)
  tibble(variable = v, var_between = var(zg), var_within = within, ratio = within / var(zg),
         q2_regional = regional_q2(as.numeric(zg), Xbar_ref[, keep, drop = FALSE]))
})
write_csv(ratio_tab, file.path(PATH_RES, "tab_variance_ratio_reference.csv"))
message(sprintf("Mediana stosunku wariancji wewnątrz/między: %.2f (TAU_RATIO w 00_config.R)",
                median(ratio_tab$ratio)))
message(sprintf("Mediana Q2 regionalnego (cel kalibracji symulacji): %.2f", median(ratio_tab$q2_regional)))
