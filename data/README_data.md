# Data access and provenance

No input data are included in this repository. Place the files below in `data/raw/` (IZP and EHIS files in `data/raw/izp/`). The scripts locate inputs by file name (`find_input()` in `R/00_config.R`), so the files may also sit directly in `data/` or its subfolders.

## Public inputs (`data/raw/`)

| File | Content | Source |
|---|---|---|
| `powiaty.json`, `wojewodztwa.json` | County and voivodeship boundaries (GeoJSON, TERYT codes in `JPT_KOD_JE`) | National register of boundaries (PRG), Head Office of Geodesy and Cartography; copy used: GitHub `waszkiewiczja/GeoJSON-Polska-Wojewodztwa-Powiaty-Gminy` |
| `2019-sejm.csv`, `2023-sejm.csv` | Sejm election results by municipality | National Electoral Commission (PKW); copy used: GitHub `DawidKK/electoral-drift` |
| `prezydent2025_obwody.csv` | 2025 presidential election results by polling district (urban/rural classification of districts) | PKW; copy used: GitHub `wybory2025/wybory2025` |
| `bezrobocie_2014_2025.csv` | Registered unemployment rate by county, 2014–2025 | Statistics Poland, Local Data Bank (BDL) export |
| `2023-01-01.csv` | TERC register of territorial units | Statistics Poland |

## Project inputs (`data/raw/izp/`)

| File | Content |
|---|---|
| `IZP_kalkulator_warianty.xlsx` (or `IZP kalkulator_warianty.xlsx`) | Health Index for Poland calculator: 2019 county indicators (sheet `2019_dane`) and source sheets (NSP 2021 education and age groups, population, migrants, fertility, medical staff, life expectancy, COPD hospitalisations, YLL) |
| `ekstra_zmienne.xlsx` | Regional EHIS 2019 values of the eleven IZP indicators (sheet `Do dezagregacji v2`), indicator descriptions, vegetable-consumption sheet `4_2_1` |
| `tablice_woj_2019_se.xlsx` | Statistics Poland EHIS 2019 regional tables with standard errors (full answer distributions) |

## Optional inputs

| File | Use |
|---|---|
| `external_validation.csv` | County-level external data for validation (columns `pow`, `target`, `value`; optional `sign`, `same_scale`), e.g. NFZ screening uptake or influenza vaccination |
| `najlepsze_modele_predykcja_pow.xlsx` | Earlier caret-based county predictions, used only for comparison |

## Notes

- `ehis_data_wide_2019.xlsx` from the earlier IZP workflow has swapped region codes 06/08/10; the pipeline reads standard errors directly from the Statistics Poland table instead.
- In `Do dezagregacji v2`, indicator x3 contains daily fruit consumption; by default it is replaced by 100 − daily vegetable consumption (`SAE_X3_FIX = "0"` keeps the original).
- Published regional values of some IZP covariates differ from the population-weighted aggregate of county values (wages up to 15%, unemployment up to 18%); the pipeline always aggregates county values.

**Before any public deposit:** verify redistribution rights for the IZP calculator and derived county estimates, add access dates and persistent identifiers, and cite Statistics Poland, PKW and GUGiK as original sources.
