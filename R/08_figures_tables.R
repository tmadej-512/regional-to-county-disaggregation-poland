## =============================================================================
## 08_figures_tables.R — figures (Lancet style, English) and tables for the paper
## Style: Arial-metric sans (Liberation Sans), 7-8 pt text, thin black axes, no
## background grid (light horizontal guides only where values are read off),
## bold panel letters A/B/C, no in-figure titles (captions live in the paper),
## Lancet palette (as in ggsci::pal_lancet).
## =============================================================================
source("R/00_config.R")

FONT <- if (.Platform$OS.type == "windows") "sans" else "Liberation Sans"   # Windows: "sans" = Arial
LANCET <- c(navy = "#00468B", red = "#ED0000", green = "#42B540", teal = "#0099B4",
            purple = "#925E9F", peach = "#FDAF91", crimson = "#AD002A", grey = "#ADB6B6",
            black = "#1B1919")
SEQ_BLUE <- c("#F2F6FB", "#C9DAEC", "#8FB2D8", "#4F84BC", "#1F5A9A", "#00468B", "#002A55")
SEQ_RED  <- c("#FFF3EE", "#FDD3C0", "#FDAF91", "#F07C5A", "#ED0000", "#AD002A")
DIV      <- c(neg = "#ED0000", mid = "#F7F7F7", pos = "#00468B")

FAMILY_EN <- c(NAIVE = "Reference", ORACLE = "Reference", PYCNO = "Spatial interpolation",
               MACRO = "Macro-first", MACRO_B = "Macro-first", MACRO_BM = "Macro-first",
               MINI_LM = "Mini-first", MINI_RF = "Mini-first", MINI_RF_OOB = "Mini-first", MINI_RF_B = "Mini-first",
               FH_B = "Model-based SAE", GDR = "Model-based SAE", GDR_GAM = "Model-based SAE",
               GP_GAM = "Model-based SAE", MGDR = "Model-based SAE",
               GDR_XL = "Model-based SAE", GP_XL = "Model-based SAE")
FAMILY_COL <- c("Reference" = "#6F7A7A", "Spatial interpolation" = LANCET[["green"]],
                "Macro-first" = LANCET[["navy"]], "Mini-first" = LANCET[["red"]],
                "Model-based SAE" = LANCET[["purple"]])
SHORT <- c(NAIVE = "Naive", PYCNO = "Pycnophylactic", MACRO = "Macro-first",
           MACRO_B = "Macro-first + B", MACRO_BM = "Macro + raking", MINI_LM = "Mini-first LM",
           MINI_RF = "Mini-first RF", MINI_RF_OOB = "Mini-first RF OOB", MINI_RF_B = "Mini-first RF + B", FH_B = "FH synthetic + B",
           GDR = "Bayesian (CAR + X)", GDR_GAM = "MRF-GAM", GP_GAM = "GP-GAM", MGDR = "Multivariate (MGDR)",
           GDR_XL = "Cross-level (CAR)", GP_XL = "Cross-level (GP)", ORACLE = "Oracle")
NOISE_EN <- c(exact = "Exact regional values", survey = "Regional values with survey error")
PART_COL <- c("Russian" = LANCET[["navy"]], "Austrian" = LANCET[["red"]],
              "Prussian" = LANCET[["green"]], "German until 1945" = LANCET[["purple"]])

## --- themes -----------------------------------------------------------------------
theme_lancet <- function(base = 8) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks.length = unit(1.5, "pt"),
          axis.text = element_text(colour = "black", size = base - 1),
          axis.title = element_text(size = base),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base, hjust = 0),
          legend.title = element_text(size = base - 0.5),
          legend.text = element_text(size = base - 1),
          legend.key.size = unit(8, "pt"),
          legend.background = element_blank(),
          legend.position = "bottom",
          plot.tag = element_text(face = "bold", size = base + 2, family = FONT),
          plot.title = element_text(size = base, face = "plain", hjust = 0),
          plot.margin = margin(4, 6, 4, 4))
}
theme_set(theme_lancet())
theme_map <- function(base = 8) {
  theme_void(base_size = base, base_family = FONT) +
    theme(legend.position = "bottom", legend.title = element_text(size = base - 0.5),
          legend.text = element_text(size = base - 1.5),
          legend.key.height = unit(4, "pt"), legend.key.width = unit(26, "pt"),
          plot.title = element_text(size = base - 0.5, hjust = 0.5, margin = margin(b = 2)),
          plot.tag = element_text(face = "bold", size = base + 2),
          plot.margin = margin(2, 2, 2, 2))
}
hgrid <- theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25))

save_fig <- function(p, name, w_mm, h_mm) {
  ggsave(file.path(PATH_FIG, paste0(name, ".png")), p, width = w_mm, height = h_mm, units = "mm",
         dpi = 400, bg = "white")
  ggsave(file.path(PATH_FIG, paste0(name, ".pdf")), p, width = w_mm, height = h_mm, units = "mm",
         bg = "white", device = cairo_pdf)
}
## etykiety z zabezpieczeniem: nieznany kod (np. metoda dodana w innej wersji kodu albo pusta nazwa)
## pokazujemy dosłownie, zamiast jako „NA”
lab_or_code <- function(x, dict) { x <- ifelse(is.na(x) | x == "", "<brak nazwy>", x)
  out <- unname(dict[x]); ifelse(is.na(out), x, out) }
mlab <- function(x) {
  l <- lab_or_code(x, METHOD_LABELS)
  factor(l, levels = unique(c(rev(METHOD_LABELS[METHOD_ORDER]), rev(sort(setdiff(l, METHOD_LABELS))))))
}
cslab <- function(x) { l <- lab_or_code(x, COVSET_LABELS)
  factor(l, levels = unique(c(COVSET_LABELS, sort(setdiff(l, COVSET_LABELS))))) }

## --- geometry ------------------------------------------------------------------------
pow_sf <- st_read(file.path(PATH_PROC, "powiaty.gpkg"), quiet = TRUE) |> arrange(pow)
woj_sf <- st_read(file.path(PATH_PROC, "woj.gpkg"), quiet = TRUE)
pow_sf <- st_simplify(pow_sf, dTolerance = 300, preserveTopology = TRUE)
woj_sf <- st_simplify(woj_sf, dTolerance = 300, preserveTopology = TRUE)

map_one <- function(vals, title = NULL, lims = NULL, pal = SEQ_BLUE, legend_title = NULL,
                    border_col = "white", border_lw = 0.3) {
  d <- pow_sf; d$v <- vals
  ggplot(d) + geom_sf(aes(fill = v), colour = NA) +
    geom_sf(data = woj_sf, fill = NA, colour = border_col, linewidth = border_lw) +
    scale_fill_gradientn(colours = pal, limits = lims, oob = scales::squish, name = legend_title) +
    labs(title = title) + theme_map()
}

## Grid of small maps as ONE ggplot (facet_grid): one shared legend, strip labels handled
## by ggplot (no overlapping tags), wrapped column/row titles.
## panels: list of list(row = "...", col = "...", v = numeric(n_powiat))
POW_GEOM <- st_sf(geometry = st_geometry(pow_sf))
map_grid <- function(panels, lims = NULL, pal = SEQ_BLUE, legend_title = NULL,
                     col_wrap = 11, row_wrap = 14, show_rows = TRUE, base = 7) {
  panels <- keep(panels, ~ !is.null(.x$v) && any(is.finite(.x$v)))
  rows <- unique(map_chr(panels, "row")); cols <- unique(map_chr(panels, "col"))
  d <- do.call(rbind, lapply(panels, function(p) {
    x <- POW_GEOM; x$v <- p$v; x$row <- p$row; x$col <- p$col; x }))
  d$row <- factor(str_wrap(d$row, row_wrap), levels = str_wrap(rows, row_wrap))
  d$col <- factor(str_wrap(d$col, col_wrap), levels = str_wrap(cols, col_wrap))
  ggplot(d) + geom_sf(aes(fill = v), colour = NA) +
    geom_sf(data = woj_sf, fill = NA, colour = "white", linewidth = 0.15) +
    facet_grid(if (show_rows) row ~ col else . ~ col, switch = "y") +
    scale_fill_gradientn(colours = pal, limits = lims, oob = scales::squish, name = legend_title) +
    theme_map(base) +
    theme(strip.text.x = element_text(size = base - 1, lineheight = 0.9, margin = margin(b = 2)),
          strip.text.y.left = element_text(size = base - 1, face = "bold", angle = 90,
                                           lineheight = 0.9, margin = margin(r = 2)),
          strip.placement = "outside",
          panel.spacing = unit(1.5, "pt"),
          legend.key.width = unit(30, "pt"), legend.key.height = unit(4, "pt"),
          legend.title = element_text(size = base - 0.5), legend.text = element_text(size = base - 1.5)) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
}
## compact legend settings used for dot/line plots with many keys
LEG_COMPACT <- theme(legend.text = element_text(size = 6.5), legend.key.size = unit(7, "pt"),
                     legend.spacing.x = unit(2, "pt"), legend.spacing.y = unit(1, "pt"),
                     legend.box.spacing = unit(3, "pt"), legend.margin = margin(0, 0, 0, 0))

## =============================================================================
## Figure 1 — study area
## =============================================================================
zab <- pow_sf |> mutate(zabor = factor(recode(zabor, R = "Russian", A = "Austrian", P = "Prussian",
                                              ZO = "German until 1945"), levels = names(PART_COL)))
f1a <- ggplot(zab) + geom_sf(aes(fill = zabor), colour = "white", linewidth = 0.05) +
  geom_sf(data = woj_sf, fill = NA, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = PART_COL, name = "Partition (1914)") +
  theme_map() + guides(fill = guide_legend(nrow = 2, title.position = "top"))
f1b <- map_one(pow_sf$urb, legend_title = "Urban share of electorate (%)") +
  guides(fill = guide_colourbar(title.position = "top"))
f1c <- map_one(pow_sf$unemp19, legend_title = "Registered unemployment 2019 (%)") +
  guides(fill = guide_colourbar(title.position = "top"))
f1d <- map_one(pow_sf$pis23, legend_title = "PiS vote share 2023 (%)") +
  guides(fill = guide_colourbar(title.position = "top"))
save_fig((f1a | f1b) / (f1c | f1d) + plot_annotation(tag_levels = "A"), "fig1_study_area", 170, 175)

## =============================================================================
## Figure 2 — simulation scenarios (one replicate)
## =============================================================================
ex <- readRDS(file.path(PATH_RES, "sim_examples.rds"))
ex_exact <- keep(ex, ~ .x$noise == "exact" && identical(.x$design, "realistic"))
names(ex_exact) <- map_chr(ex_exact, "scen")
f2 <- map(names(SCEN_LABELS), function(s)
  map_one(ex_exact[[s]]$y, SCEN_LABELS[[s]], lims = c(38, 62), legend_title = "Simulated county value"))
f2 <- wrap_plots(f2, ncol = 4) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom") & guides(fill = guide_colourbar(title.position = "top"))
save_fig(f2, "fig2_scenarios", 180, 112)

## =============================================================================
## Figure 3 — example disaggregations (S1, S3)
## =============================================================================
## all methods shown on example maps (those missing in a given run are dropped automatically)
SHOW_MAPS <- c("NAIVE", "PYCNO", "MACRO_B", "MINI_RF", "MINI_RF_OOB", "FH_B", "GDR", "GDR_GAM", "GP_GAM",
               "GDR_XL")
ex_panels <- function(e, row_lab) {
  c(list(list(row = row_lab, col = "Truth", v = e$y)),
    map(SHOW_MAPS, ~ list(row = row_lab, col = SHORT[[.x]], v = e$preds[[.x]])))
}
f3 <- map_grid(c(ex_panels(ex_exact[["S1_base"]], SCEN_LABELS[["S1_base"]]),
                 ex_panels(ex_exact[["S3_zabory"]], SCEN_LABELS[["S3_zabory"]])),
               lims = c(38, 62), legend_title = "County value")
save_fig(f3, "fig3_example_maps", 180, 68)

## =============================================================================
## Simulation results
## =============================================================================
DESIGN_EN <- c(realistic = "Main design (variance ratio and regional Q2 as in real data)",
               regional_noise = "Sensitivity: unexplained regional effects",
               heterogeneous = "Sensitivity: no regional component (high within-region heterogeneity)")
sim_all <- read_csv(file.path(PATH_RES, "sim_metrics.csv"), show_col_types = FALSE, guess_max = 1e6) |>
  group_by(design, scen, noise, rep) |>
  mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |> ungroup()
sim_sum_all <- sim_all |> group_by(design, scen, noise, method) |>
  summarise(RMSE_mean = mean(RMSE), RMSE_mcse = sd(RMSE) / sqrt(n()),
            skill = mean(skill), skill_q10 = quantile(skill, .1), skill_q90 = quantile(skill, .9),
            CorWithin = mean(CorWithin), Spearman = mean(Spearman), Top10 = mean(Top10),
            Bias = mean(Bias), WojRMSE = mean(WojRMSE), Cover90 = mean(Cover90),
            Width90 = mean(Width90), n_rep = n(), .groups = "drop")
write_csv(sim_sum_all, file.path(PATH_RES, "tab_sim_summary.csv"))
unk <- setdiff(unique(sim_all$method), names(METHOD_LABELS))
if (length(unk)) message("UWAGA: metody bez etykiety w sim_metrics.csv: ", paste(unk, collapse = ", "),
                         " (liczba wierszy: ", sum(sim_all$method %in% unk | is.na(sim_all$method)), ")")
## analiza główna = wariant realistyczny
sim     <- sim_all |> filter(design == "realistic")
sim_sum <- sim_sum_all |> filter(design == "realistic")

## Table 2 (paper): RMSE (MCSE) by scenario × method, survey-error variant + exact
tab2 <- sim_sum_all |> mutate(cell = sprintf("%.2f", RMSE_mean)) |>
  select(design, scen, noise, method, cell) |>
  pivot_wider(names_from = method, values_from = cell) |>
  mutate(scen = SCEN_LABELS[scen], noise = NOISE_EN[noise]) |>
  select(design, scen, noise, any_of(METHOD_ORDER))
write_csv(tab2, file.path(PATH_RES, "tab2_sim_rmse.csv"))
write_csv(read_csv(file.path(PATH_RES, "sim_dgp_diagnostics.csv"), show_col_types = FALSE) |>
            group_by(design, scen) |> summarise(across(-rep, mean)) |> mutate(scen = SCEN_LABELS[scen]),
          file.path(PATH_RES, "tab1_scenario_diagnostics.csv"))

## Figure 4 — MSE reduction relative to naive (main design); Figure S1 — sensitivity design
## ncol = NULL: one row when <= 5 scenarios, otherwise 4 columns; one-row version uses short
## method names on the y axis and wrapped panel titles so nothing is truncated
make_f4 <- function(ss, ncol = NULL) {
  n_s <- n_distinct(ss$scen)
  if (is.null(ncol)) ncol <- if (n_s <= 5) n_s else 4
  one_row <- ncol >= n_s
  sl <- str_wrap(SCEN_LABELS, if (one_row) 14 else 22)
  ylab_f <- if (one_row) function(x) factor(SHORT[x], levels = rev(SHORT[intersect(METHOD_ORDER, names(SHORT))]))
            else mlab
  d <- ss |> filter(method != "NAIVE") |>
    mutate(m = ylab_f(method), fam = factor(lab_or_code(method, FAMILY_EN), levels = c(names(FAMILY_COL), setdiff(unique(lab_or_code(method, FAMILY_EN)), names(FAMILY_COL)))),
           scen = factor(sl[scen], levels = sl), noise = NOISE_EN[noise],
           lo = pmax(skill_q10, -3), hi = pmin(skill_q90, 1), pt = pmax(skill, -3))
  ggplot(d, aes(x = pt, y = m, colour = fam, shape = noise)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.35, alpha = 0.55,
                   position = position_dodge(width = 0.65)) +
    geom_point(size = 1.4, position = position_dodge(width = 0.65), stroke = 0.5) +
    facet_wrap(~scen, ncol = ncol) +
    scale_colour_manual(values = FAMILY_COL, name = NULL) +
    scale_shape_manual(values = c(16, 21), name = NULL) +
    scale_x_continuous(limits = c(-3, 1), breaks = c(-3, -2, -1, 0, 1), labels = c("≤-3", "-2", "-1", "0", "1")) +
    labs(x = "Reduction in county-level MSE relative to naive allocation (1 - MSE/MSE naive)", y = NULL) +
    hgrid + LEG_COMPACT +
    theme(legend.box = "vertical", panel.spacing.x = unit(if (one_row) 4 else 8, "pt"),
          strip.text = element_text(size = 7, lineheight = 0.9, hjust = 0),
          axis.text.x = element_text(size = if (one_row) 5.5 else 7)) +
    guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(nrow = 1, order = 2))
}
## height: one row needs less
f4_h <- function(ss) if (n_distinct(ss$scen) <= 5) 100 else 140
save_fig(make_f4(sim_sum), "fig4_sim_skill", 180, f4_h(sim_sum))
## S1a: bez komponentu regionalnego (optymistyczny); S1b: różnice regionalne niewyjaśnione (pesymistyczny)
for (dn in c("heterogeneous", "regional_noise")) {
  ssd <- sim_sum_all |> filter(design == dn)
  if (nrow(ssd)) save_fig(make_f4(ssd),
                          c(heterogeneous = "figS1a_sim_skill_heterogeneous",
                            regional_noise = "figS1b_sim_skill_regional_noise")[[dn]], 180, f4_h(ssd))
}

## Figure 5 — within-region correlation (A) and top-decile hit rate (B)
f5d <- sim_sum |> filter(method != "NAIVE", noise == "survey") |>
  mutate(m = mlab(method), fam = factor(lab_or_code(method, FAMILY_EN), levels = c(names(FAMILY_COL), setdiff(unique(lab_or_code(method, FAMILY_EN)), names(FAMILY_COL)))),
         scen = factor(str_wrap(SCEN_LABELS[scen], 11), levels = str_wrap(SCEN_LABELS, 11)))
f5a <- ggplot(f5d, aes(CorWithin, m, colour = fam)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
  geom_point(size = 1.3) + facet_grid(~scen) +
  scale_colour_manual(values = FAMILY_COL, name = NULL) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5), labels = c("-.5", "0", ".5"), limits = c(-0.75, 1)) +
  labs(x = "Correlation of within-region deviations (estimate vs truth)", y = NULL) + hgrid
f5b <- ggplot(f5d, aes(Top10, m, colour = fam)) +
  geom_vline(xintercept = 0.10, colour = "grey40", linewidth = 0.3, linetype = 2) +
  geom_point(size = 1.3) + facet_grid(~scen) +
  scale_colour_manual(values = FAMILY_COL, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6), labels = c("0", ".3", ".6"), limits = c(0, 0.75)) +
  labs(x = "Share of true top-decile counties correctly identified (dashed: chance)", y = NULL) + hgrid
save_fig((f5a / f5b) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
           theme(legend.position = "bottom", strip.text = element_text(size = 6.5, hjust = 0.5, lineheight = 0.9),
                 axis.text.x = element_text(size = 5.5)) & LEG_COMPACT,
         "fig5_sim_within_top10", 180, 150)

## Figure 6 — the regional-validation paradox
f6d <- sim |> filter(scen %in% c("S1_base", "S3_zabory"), noise == "exact", method != "ORACLE") |>
  group_by(scen, method) |> summarise(WojRMSE = mean(WojRMSE), RMSE = mean(RMSE), .groups = "drop") |>
  mutate(fam = factor(lab_or_code(method, FAMILY_EN), levels = c(names(FAMILY_COL), setdiff(unique(lab_or_code(method, FAMILY_EN)), names(FAMILY_COL)))), lab = SHORT[method],
         scen = factor(SCEN_LABELS[scen], levels = SCEN_LABELS))
f6 <- ggplot(f6d, aes(WojRMSE, RMSE, colour = fam)) +
  geom_point(size = 1.8) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.3, family = FONT, colour = "black",
                           min.segment.length = 0, segment.colour = "grey60", segment.size = 0.25,
                           box.padding = 0.3, max.overlaps = 50, seed = 3) +
  facet_wrap(~scen, scales = "free_y") +
  scale_colour_manual(values = FAMILY_COL, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.15))) +
  labs(x = "RMSE after re-aggregation to regions (conventional validation)",
       y = "RMSE at county level (target of inference)") +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.25)) + LEG_COMPACT +
  guides(colour = guide_legend(nrow = 1))
save_fig(f6, "fig6_validation_paradox", 180, 85)

## Figure 7 — signal-strength sweep
sw <- read_csv(file.path(PATH_RES, "sweep_metrics.csv"), show_col_types = FALSE, guess_max = 1e6) |>
  group_by(r2, noise, rep) |> mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |> ungroup() |>
  group_by(r2, noise, method) |>
  summarise(skill = mean(skill), CorWithin = mean(CorWithin), Cover90 = mean(Cover90),
            .groups = "drop")
write_csv(sw, file.path(PATH_RES, "tab_sweep_summary.csv"))
sel <- c("PYCNO", "MACRO_B", "MINI_LM", "MINI_RF_B", "FH_B", "GDR", "GDR_GAM", "GP_GAM", "GDR_XL", "ORACLE")
SWEEP_COL <- c(PYCNO = LANCET[["green"]], MACRO_B = LANCET[["navy"]], MINI_LM = LANCET[["peach"]],
               MINI_RF_B = LANCET[["red"]], FH_B = LANCET[["teal"]], GDR = LANCET[["purple"]],
               GDR_GAM = LANCET[["crimson"]], GP_GAM = "#B07AA1", GDR_XL = "#E3A21A", ORACLE = LANCET[["black"]])
swd <- sw |> filter(method %in% sel) |>
  mutate(lab = factor(SHORT[method], levels = SHORT[sel]), noise = NOISE_EN[noise])
f7 <- ggplot(swd, aes(r2, pmax(skill, -1), colour = lab, linetype = lab)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
  geom_line(linewidth = 0.5) + geom_point(size = 0.9) +
  facet_wrap(~noise) +
  scale_colour_manual(values = setNames(SWEEP_COL, SHORT[sel]), name = NULL) +
  scale_linetype_manual(values = setNames(c(1, 1, 2, 2, 1, 1, 2, 4, 1, 3), SHORT[sel]), name = NULL) +
  scale_x_continuous(breaks = seq(0, 0.9, 0.3)) +
  coord_cartesian(ylim = c(-0.4, 1)) +
  labs(x = expression("Share of county-level variance explained by covariates (" * R^2 * ")"),
       y = "Reduction in MSE vs naive") +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
        legend.key.width = unit(16, "pt")) + LEG_COMPACT +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE))
save_fig(f7, "fig7_signal_sweep", 180, 85)

## Figure 8 — coverage of 90% intervals: GDR and MRF-GAM, main and sensitivity designs
f8d <- sim_sum_all |> filter(method %in% c("GDR", "GDR_GAM", "GP_GAM", "GDR_XL")) |>
  mutate(scen = factor(SCEN_LABELS[scen], levels = rev(SCEN_LABELS)), noise = NOISE_EN[noise],
         lab = factor(SHORT[method], levels = SHORT[c("GDR", "GDR_GAM", "GP_GAM", "GDR_XL")]),
         design = factor(c(realistic = "Main design", regional_noise = "Unexplained regional\neffects",
                           heterogeneous = "No regional\ncomponent")[design],
                         levels = c("Main design", "Unexplained regional\neffects", "No regional\ncomponent")))
f8 <- ggplot(f8d, aes(Cover90, scen, shape = noise, colour = lab)) +
  geom_vline(xintercept = 0.9, colour = "grey40", linewidth = 0.3, linetype = 2) +
  geom_point(size = 1.7, stroke = 0.6, position = position_dodge(width = 0.6)) +
  facet_wrap(~design) +
  scale_shape_manual(values = c(16, 21), name = NULL) +
  scale_colour_manual(values = setNames(c(LANCET[["purple"]], LANCET[["crimson"]], "#B07AA1", "#E3A21A"),
                                        SHORT[c("GDR", "GDR_GAM", "GP_GAM", "GDR_XL")]), name = NULL) +
  scale_x_continuous(limits = c(0.4, 1), labels = scales::percent) +
  labs(x = "Empirical coverage of nominal 90% prediction intervals", y = NULL) + hgrid +
  guides(shape = guide_legend(nrow = 1, order = 2), colour = guide_legend(nrow = 1, order = 1)) +
  theme(legend.box = "vertical", strip.text = element_text(lineheight = 0.9),
        axis.text.x = element_text(size = 6)) + LEG_COMPACT
save_fig(f8, "fig8_gdr_coverage", 180, 95)

## Figure S2 — example disaggregations for ALL scenarios (main design, exact values)
## one facet_grid: rows = scenarios (labels on the left), columns = truth + methods
fS2 <- map_grid(unlist(map(names(SCEN_LABELS), function(sc) {
  e <- ex_exact[[sc]]; if (is.null(e)) return(NULL)
  ex_panels(e, SCEN_LABELS[[sc]])
}), recursive = FALSE), lims = c(38, 62), legend_title = "County value", row_wrap = 12)
save_fig(fS2, "figS2_example_maps_all_scenarios", 180, 155)

## Figure S3 — the regional-validation paradox for ALL scenarios (main design, exact values)
fS3d <- sim |> filter(noise == "exact", method != "ORACLE") |>
  group_by(scen, method) |> summarise(WojRMSE = mean(WojRMSE), RMSE = mean(RMSE), .groups = "drop") |>
  mutate(fam = factor(lab_or_code(method, FAMILY_EN), levels = c(names(FAMILY_COL), setdiff(unique(lab_or_code(method, FAMILY_EN)), names(FAMILY_COL)))), lab = SHORT[method],
         scen = factor(SCEN_LABELS[scen], levels = SCEN_LABELS))
fS3 <- ggplot(fS3d, aes(WojRMSE, RMSE, colour = fam)) +
  geom_point(size = 1.5) +
  ggrepel::geom_text_repel(aes(label = lab), size = 1.9, family = FONT, colour = "black",
                           min.segment.length = 0, segment.colour = "grey60", segment.size = 0.2,
                           box.padding = 0.2, max.overlaps = 50, seed = 3) +
  facet_wrap(~scen, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = FAMILY_COL, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.15))) +
  labs(x = "RMSE after re-aggregation to regions (conventional validation)",
       y = "RMSE at county level (target of inference)") +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.25)) + LEG_COMPACT +
  guides(colour = guide_legend(nrow = 1))
save_fig(fS3, "figS3_validation_paradox_all_scenarios", 180, 200)

## Figure S5 — does regional Q2 (observable before disaggregation) predict the gain?
if ("q2_regional" %in% names(sim_all)) {
  q2m <- c("PYCNO", "MACRO_B", "MINI_RF_B", "FH_B", "GDR", "GDR_GAM", "GP_GAM", "GDR_XL")
  q2col <- c(PYCNO = LANCET[["green"]], MACRO_B = LANCET[["navy"]], MINI_RF_B = LANCET[["red"]],
             FH_B = LANCET[["teal"]], GDR = LANCET[["purple"]], GDR_GAM = LANCET[["crimson"]],
             GP_GAM = "#B07AA1", GDR_XL = "#E3A21A")
  q2d <- sim_all |> filter(method %in% q2m, scen != "S6_ecol") |>
    mutate(bin = cut(q2_regional, breaks = c(-Inf, seq(-0.2, 0.8, by = 0.1), Inf))) |>
    group_by(method, bin) |>
    summarise(q2 = median(q2_regional), skill = mean(skill), n = n(), .groups = "drop") |>
    filter(n >= 30) |>
    mutate(lab = factor(SHORT[method], levels = SHORT[q2m]))
  s5a <- ggplot(q2d, aes(q2, pmax(skill, -1.5), colour = lab)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
    geom_vline(xintercept = Q2_MIN, colour = "grey60", linewidth = 0.3, linetype = 3) +
    geom_line(linewidth = 0.5) + geom_point(size = 1) +
    scale_colour_manual(values = setNames(q2col[q2m], SHORT[q2m]), name = NULL) +
    labs(x = "Regional Q2 (leave-one-region-out)",
         y = "Reduction in MSE vs naive", subtitle = "Simulation (all designs, S6 excluded)") +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25)) + LEG_COMPACT +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))
  f_real <- file.path(PATH_RES, "real_metrics.csv")
  if (file.exists(f_real)) {
    rq <- read_csv(f_real, show_col_types = FALSE, guess_max = 1e6)
    if ("q2_regional" %in% names(rq)) {
      rq <- rq |> filter(noise == "exact") |> group_by(covset, target) |>
        mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |> ungroup() |>
        filter(method %in% c("GDR", "GDR_GAM", "GP_GAM", "GDR_XL", "FH_B", "MACRO_B")) |>
        mutate(lab = factor(SHORT[method], levels = SHORT[q2m]))
      s5b <- ggplot(rq, aes(q2_regional, pmax(skill, -1.5), colour = lab, shape = covset)) +
        geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
        geom_vline(xintercept = Q2_MIN, colour = "grey60", linewidth = 0.3, linetype = 3) +
        geom_point(size = 1.4) +
        scale_colour_manual(values = setNames(q2col[q2m], SHORT[q2m]), name = NULL) +
        scale_shape_manual(values = c(izp = 16, public = 21, core = 15, wide = 17), name = NULL,
                           labels = COVSET_LABELS) +
        guides(colour = "none") +   # kolory opisuje legenda panelu A (bez dublowania)
        labs(x = "Regional Q2", y = "Reduction in MSE vs naive",
             subtitle = "Real county variables (exact regional values)") +
        theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25))
      s5a <- (s5a | s5b) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
        theme(legend.position = "bottom", legend.box = "vertical") & LEG_COMPACT
    }
  }
  save_fig(s5a, "figS5_q2_diagnostic", 180, 100)
}

## Figure S6 — multivariate disaggregation (MGDR) vs separate GDR fits
f_mv <- file.path(PATH_RES, "mv_metrics.csv")
if (file.exists(f_mv)) {
  mvm <- read_csv(f_mv, show_col_types = FALSE, guess_max = 1e6) |>
    group_by(c_common, noise, rep, indicator) |>
    mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |> ungroup() |>
    filter(method != "NAIVE") |>
    group_by(c_common, noise, method) |>
    summarise(skill = mean(skill), lo = quantile(skill, .1), hi = quantile(skill, .9),
              Cover90 = mean(Cover90), .groups = "drop") |>
    mutate(noise = NOISE_EN[noise], lab = factor(SHORT[method], levels = SHORT[c("GDR", "MGDR")]),
           c_common = factor(c_common, labels = paste0("Shared factor ", c("absent", "moderate (0.5)", "strong (0.8)"))))
  write_csv(mvm, file.path(PATH_RES, "tab_mv_sim_summary.csv"))
  colmv <- setNames(c(LANCET[["purple"]], LANCET[["navy"]]), SHORT[c("GDR", "MGDR")])
  s6a <- ggplot(mvm, aes(skill, c_common, colour = lab)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
    geom_point(size = 1.6, position = position_dodge(width = 0.5)) + facet_wrap(~noise) +
    scale_colour_manual(values = colmv, name = NULL) +
    labs(x = "Reduction in MSE vs naive (simulation, 6 indicators)", y = NULL) + hgrid
  s6 <- s6a
  f_rmv <- file.path(PATH_RES, "real_mv_metrics.csv")
  if (file.exists(f_rmv)) {
    rmv <- read_csv(f_rmv, show_col_types = FALSE, guess_max = 1e6) |>
      group_by(covset, noise, rep, target) |>
      mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |> ungroup() |>
      filter(method != "NAIVE") |>
      group_by(covset, noise, target, method) |> summarise(skill = mean(skill), .groups = "drop") |>
      mutate(noise = NOISE_EN[noise], lab = factor(SHORT[method], levels = SHORT[c("GDR", "MGDR")]),
             label = factor(str_wrap(REAL_LABELS_EN[target], 22), levels = str_wrap(REAL_LABELS_EN, 22)),
             covset = cslab(covset))
    write_csv(rmv, file.path(PATH_RES, "tab_mv_real_summary.csv"))
    s6b <- ggplot(rmv, aes(pmax(skill, -1), label, colour = lab, shape = covset)) +
      geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
      geom_point(size = 1.4, position = position_dodge(width = 0.6)) + facet_wrap(~noise) +
      scale_colour_manual(values = colmv, name = NULL) +
      scale_shape_manual(values = c(21, 16, 15, 17), name = NULL, drop = TRUE) +
      labs(x = "Reduction in MSE vs naive (real county variables, fitted jointly)", y = NULL) + hgrid
    s6 <- (s6a / s6b) + plot_layout(guides = "collect", heights = c(1, 2)) +
      plot_annotation(tag_levels = "A") & theme(legend.position = "bottom") & LEG_COMPACT
  }
  save_fig(s6, "figS6_multivariate", 180, 150)
}

## Figure S7 — amplitude of within-region differences: calibration slope
## CalSlope = regression slope of true within-region deviations on estimated ones;
## 1 = right amplitude, < 1 = the method exaggerates differences inside regions.
if ("CalSlope" %in% names(sim_all)) {
  cs_m <- c("MACRO_B", "MINI_RF_B", "FH_B", "GDR", "GDR_GAM", "GP_GAM", "GDR_XL", "GP_XL")
  csd <- sim_all |> filter(design == "realistic", noise == "survey", method %in% cs_m) |>
    group_by(scen, method) |> summarise(cs = median(CalSlope, na.rm = TRUE), .groups = "drop") |>
    mutate(src = "Simulation, main design", row = SCEN_LABELS[scen])
  f_real <- file.path(PATH_RES, "real_metrics.csv")
  if (file.exists(f_real)) {
    rq <- read_csv(f_real, show_col_types = FALSE, guess_max = 1e6)
    if ("CalSlope" %in% names(rq))
      csd <- bind_rows(csd, rq |> filter(noise == "exact", covset == "izp", method %in% cs_m) |>
        group_by(target, method) |> summarise(cs = median(CalSlope, na.rm = TRUE), .groups = "drop") |>
        mutate(src = "Real county variables (IZP covariates)", row = REAL_LABELS_EN[target]))
  }
  csd <- csd |> filter(is.finite(cs)) |>
    mutate(lab = factor(SHORT[method], levels = SHORT[cs_m]),
           row = factor(str_wrap(row, 24), levels = rev(unique(str_wrap(row, 24)))),
           src = factor(src, levels = unique(src)))
  write_csv(csd, file.path(PATH_RES, "tab_calibration_slope.csv"))
  colCS <- setNames(c(LANCET[["navy"]], LANCET[["red"]], LANCET[["teal"]], LANCET[["purple"]],
                      LANCET[["crimson"]], "#B07AA1", "#E3A21A", "#8C6D1F"), SHORT[cs_m])
  s7 <- ggplot(csd, aes(pmin(pmax(cs, -0.5), 2), row, colour = lab)) +
    geom_vline(xintercept = 1, colour = "grey40", linewidth = 0.3, linetype = 2) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_point(size = 1.5, position = position_dodge(width = 0.6)) +
    facet_wrap(~src, scales = "free_y") +
    scale_colour_manual(values = colCS, name = NULL) +
    labs(x = "Calibration slope of within-region deviations (1 = right amplitude; <1 = exaggerated)", y = NULL) +
    hgrid + LEG_COMPACT + guides(colour = guide_legend(nrow = 2, byrow = TRUE))
  save_fig(s7, "figS7_calibration_slope", 180, 120)
}

## Figure S8 — ordinal (whole response distribution) vs binary disaggregation (05d, part A)
## Figure S9 — joint disaggregation of all ordinal questions (05d, part B)
f_ord <- file.path(PATH_RES, "ordinal_sim_metrics.csv")
if (file.exists(f_ord)) {
  om <- read_csv(f_ord, show_col_types = FALSE, guess_max = 1e6) |>
    group_by(target, scen, design, shape, noise, rep) |>
    mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |> ungroup() |>
    filter(method != "NAIVE") |>
    group_by(target, design, shape, noise, method) |>
    summarise(skill = mean(skill), Cover90 = mean(Cover90), CorWithin = mean(CorWithin),
              .groups = "drop") |>
    mutate(approach = ifelse(grepl("^ORD_", method), "Ordinal (all categories)", "Binary indicator"),
           base = factor(SHORT[sub("^(ORD|BIN)_", "", method)], levels = SHORT[c("PYCNO", "GDR", "GDR_XL")]),
           question = factor(HEALTH_LABELS_EN[target], levels = HEALTH_LABELS_EN[c("x1", "x2", "x3", "x4", "x5")]),
           shape = factor(recode(shape, proportional = "Model correct", scale = "Response spread varies",
                                 nonprop = "Non-proportional odds"),
                          levels = c("Model correct", "Response spread varies", "Non-proportional odds")),
           design = recode(design, realistic = "Main design", no_regional = "No regional component"),
           noise = NOISE_EN[noise])
  write_csv(om, file.path(PATH_RES, "tab_ordinal_sim_summary.csv"))
  s8 <- ggplot(om |> filter(noise == NOISE_EN[["survey"]]),
               aes(pmax(skill, -1), base, colour = approach, shape = design)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
    geom_point(size = 1.4, position = position_dodge(width = 0.7)) +
    facet_grid(shape ~ question, labeller = labeller(question = label_wrap_gen(14), shape = label_wrap_gen(12))) +
    scale_colour_manual(values = c("Binary indicator" = LANCET[["navy"]],
                                   "Ordinal (all categories)" = LANCET[["red"]]), name = NULL) +
    scale_shape_manual(values = c(16, 21), name = NULL) +
    labs(x = "Reduction in MSE vs naive (regional values with survey error; mean over scenarios)", y = NULL) +
    hgrid + LEG_COMPACT +
    theme(legend.box = "vertical", strip.text = element_text(size = 6.5, lineheight = 0.9, hjust = 0.5),
          axis.text.x = element_text(size = 5.5))
  save_fig(s8, "figS8_ordinal_vs_binary", 180, 130)
}
f_omv <- file.path(PATH_RES, "ordinal_sim_mv_metrics.csv")
if (file.exists(f_omv)) {
  omv <- read_csv(f_omv, show_col_types = FALSE, guess_max = 1e6) |>
    group_by(c_common, noise, rep, target) |>
    mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |> ungroup() |>
    filter(method != "NAIVE") |>
    group_by(c_common, noise, target, method) |> summarise(skill = mean(skill), .groups = "drop") |>
    mutate(question = factor(HEALTH_LABELS_EN[target], levels = HEALTH_LABELS_EN[c("x1", "x2", "x3", "x4", "x5")]),
           method = factor(method, levels = c("BIN_GDR", "ORD_GDR", "BIN_MGDR", "ORD_MGDR"),
                           labels = c("Binary, separate", "Ordinal, separate", "Binary, joint (MGDR)",
                                      "Ordinal, joint (MGDR)")),
           c_common = paste0("Shared factor ", ifelse(c_common > 0, "present", "absent")),
           noise = NOISE_EN[noise])
  write_csv(omv, file.path(PATH_RES, "tab_ordinal_sim_mv_summary.csv"))
  s9 <- ggplot(omv, aes(pmax(skill, -1), question, colour = method)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
    geom_point(size = 1.4, position = position_dodge(width = 0.7)) +
    facet_grid(c_common ~ noise) +
    scale_colour_manual(values = setNames(c(LANCET[["navy"]], LANCET[["red"]], LANCET[["teal"]], LANCET[["crimson"]]),
                                          c("Binary, separate", "Ordinal, separate", "Binary, joint (MGDR)",
                                            "Ordinal, joint (MGDR)")), name = NULL) +
    labs(x = "Reduction in MSE vs naive", y = NULL) + hgrid + LEG_COMPACT +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))
  save_fig(s9, "figS9_ordinal_multivariate", 180, 100)
}

## Figure S4 — sensitivity of GDR and MRF-GAM to rho and the variance prior
f_sens <- file.path(PATH_RES, "sens_metrics.csv")
if (file.exists(f_sens)) {
  sens <- read_csv(f_sens, show_col_types = FALSE, guess_max = 1e6) |>
    group_by(scen, rep) |>
    mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |> ungroup() |>
    filter(method != "NAIVE") |>
    group_by(scen, method, config) |>
    summarise(skill = mean(skill), Cover90 = mean(Cover90), Width90 = mean(Width90), .groups = "drop")
  write_csv(sens, file.path(PATH_RES, "tab_sensitivity_gdr.csv"))
  cfg_lv <- c("default", "rho = 0.80", "rho = 0.99", "tau ratio = 1", "tau ratio = 4",
              "prior SD = 0.5", "prior SD = 2", "c0 = 0.25", "c0 = 0.5", "c0 = 0.75")
  sd4 <- sens |> mutate(config = factor(config, levels = rev(cfg_lv)),
                        scen = factor(str_wrap(SCEN_LABELS[scen], 11), levels = str_wrap(SCEN_LABELS, 11)),
                        lab = factor(SHORT[method], levels = SHORT[c("GDR", "GDR_GAM", "GP_GAM", "GDR_XL")]))
  colS4 <- setNames(c(LANCET[["purple"]], LANCET[["crimson"]], "#B07AA1", "#E3A21A"),
                    SHORT[c("GDR", "GDR_GAM", "GP_GAM", "GDR_XL")])
  s4a <- ggplot(sd4, aes(skill, config, colour = lab)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
    geom_point(size = 1.3, position = position_dodge(width = 0.5)) + facet_grid(~scen) +
    scale_colour_manual(values = colS4, name = NULL) +
    labs(x = "Reduction in MSE vs naive", y = NULL) + hgrid
  s4b <- ggplot(sd4, aes(Cover90, config, colour = lab)) +
    geom_vline(xintercept = 0.9, colour = "grey40", linewidth = 0.3, linetype = 2) +
    geom_point(size = 1.3, position = position_dodge(width = 0.5)) + facet_grid(~scen) +
    scale_colour_manual(values = colS4, name = NULL) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Coverage of nominal 90% intervals", y = NULL) + hgrid
  save_fig((s4a / s4b) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
             theme(legend.position = "bottom", strip.text = element_text(size = 6.5, hjust = 0.5, lineheight = 0.9),
                   axis.text.x = element_text(size = 5.5)) & LEG_COMPACT,
           "figS4_sensitivity_gdr", 180, 140)
}

## =============================================================================
## Real-data validation
## =============================================================================
if (file.exists(file.path(PATH_RES, "real_metrics.csv"))) {
  rm <- read_csv(file.path(PATH_RES, "real_metrics.csv"), show_col_types = FALSE, guess_max = 1e6) |>
    group_by(covset, target, noise, rep) |>
    mutate(skill = 1 - RMSE^2 / RMSE[method == "NAIVE"]^2) |> ungroup()
  real_sum <- rm |> group_by(covset, target, noise, method) |>
    summarise(RMSE = mean(RMSE), skill = mean(skill), CorWithin = mean(CorWithin),
              Top10 = mean(Top10), Cover90 = mean(Cover90), WojRMSE = mean(WojRMSE), .groups = "drop") |>
    mutate(label = REAL_LABELS_EN[target])
  write_csv(real_sum, file.path(PATH_RES, "tab_real_summary.csv"))

  f9d <- real_sum |> filter(method != "NAIVE") |>
    mutate(m = mlab(method), label = factor(label, levels = REAL_LABELS_EN),
           covset = cslab(covset), noise = NOISE_EN[noise], v = pmax(pmin(skill, 1), -1))
  ## Figure 9: mean over the 8 target variables (range = min–max), by covariate set
  f9s <- f9d |> group_by(covset, noise, m, method) |>
    summarise(mean = mean(skill), lo = min(skill), hi = max(skill), .groups = "drop")
  covcol <- setNames(c(LANCET[["navy"]], LANCET[["red"]], LANCET[["green"]], LANCET[["purple"]],
                       LANCET[["grey"]], LANCET[["teal"]])[seq_along(levels(f9s$covset))], levels(f9s$covset))
  f9 <- ggplot(f9s, aes(pmax(mean, -1), m, colour = covset)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3, linetype = 2) +
    geom_linerange(aes(xmin = pmax(lo, -1), xmax = pmin(hi, 1)), linewidth = 0.35, alpha = 0.5,
                   position = position_dodge(width = 0.7)) +
    geom_point(size = 1.4, position = position_dodge(width = 0.7)) +
    facet_wrap(~noise) +
    scale_colour_manual(values = covcol, name = NULL) +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("≤-1", "-.5", "0", ".5", "1")) +
    labs(x = "Reduction in MSE vs naive, real county variables (point: mean of 8 targets; line: range)", y = NULL) +
    hgrid + LEG_COMPACT + guides(colour = guide_legend(nrow = 2, byrow = TRUE))
  save_fig(f9, "fig9_real_validation", 180, 150)
  ## Figure S10: full heatmap, one file per covariate set
  for (cs in levels(droplevels(f9d$covset))) {
    hm <- ggplot(f9d |> filter(covset == cs), aes(label, m, fill = v)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", skill), colour = abs(v) > 0.55), size = 2.1, family = FONT) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
      facet_wrap(~noise) +
      scale_fill_gradient2(low = DIV[["neg"]], mid = DIV[["mid"]], high = DIV[["pos"]], midpoint = 0,
                           limits = c(-1, 1), name = "Reduction in MSE vs naive") +
      scale_x_discrete(labels = function(x) str_wrap(x, 16)) +
      labs(x = NULL, y = NULL, subtitle = cs) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6, lineheight = 0.85), axis.line = element_blank(),
            axis.ticks = element_blank(), strip.text = element_text(hjust = 0.5),
            legend.key.width = unit(30, "pt"), legend.key.height = unit(5, "pt")) +
      guides(fill = guide_colourbar(title.position = "top"))
    save_fig(hm, paste0("figS10_real_validation_", gsub("[^A-Za-z0-9]+", "_", tolower(cs))), 180, 130)
  }

  rex <- readRDS(file.path(PATH_RES, "real_examples.rds"))
  ## each variable has its own colour scale -> one facet row per variable, legend on the right
  mk_real <- function(key, row_lab, unit_lab) {
    e <- rex[[key]]
    map_grid(ex_panels(e, row_lab), lims = quantile(e$y, c(.02, .98)), legend_title = unit_lab,
             row_wrap = 12) +
      theme(legend.position = "right", legend.key.height = unit(14, "pt"),
            legend.key.width = unit(4, "pt"), legend.title = element_text(size = 6),
            legend.text = element_text(size = 5.5)) +
      guides(fill = guide_colourbar(title.position = "top"))
  }
  f10 <- mk_real("izp:pis23", "PiS vote share 2023", "%") /
         mk_real("izp:unemp23", "Unemployment 2023", "%")
  save_fig(f10, "fig10_real_example_maps", 180, 72)
}

## =============================================================================
## Application — EHIS 2019 health indicators
## =============================================================================
if (file.exists(file.path(PATH_RES, "app_health_long.rds"))) {
  app <- readRDS(file.path(PATH_RES, "app_health_long.rds"))
  ht <- read_csv(file.path(PATH_PROC, "health_targets_2019.csv"), show_col_types = FALSE,
                 col_types = cols(woj = "c"))
  tg <- paste0("x", 1:11)
  ## metoda główna: wybrana w 06b (jeśli daje przedziały), inaczej GDR
  sel_m <- if (file.exists(file.path(PATH_RES, "selected_method.txt")))
    readLines(file.path(PATH_RES, "selected_method.txt"))[1] else "GDR"
  if (!sel_m %in% c("GDR", "GDR_GAM", "GP_GAM", "GDR_XL", "GP_XL")) sel_m <- "GDR"
  g <- app |> filter(method == sel_m)
  f11 <- map(tg, function(t) {
    d <- g |> filter(target == t); d <- d[match(pow_sf$pow, d$pow), ]
    map_one(d$est, str_wrap(HEALTH_LABELS_EN[[t]], 26), lims = quantile(d$est, c(.02, .98)), legend_title = "%") +
      theme(plot.title = element_text(size = 6.5, lineheight = 0.9), legend.key.width = unit(20, "pt"),
            legend.title = element_text(size = 6, vjust = 0.9))
  })
  save_fig(wrap_plots(f11, ncol = 4) + plot_annotation(tag_levels = "A"), "fig11_app_health_gdr", 180, 175)

  t0 <- "x6"
  d <- g |> filter(target == t0); d <- d[match(pow_sf$pow, d$pow), ]
  z <- ht |> filter(target == t0); zz <- z$value[match(pow_sf$woj, z$woj)]
  lims <- range(c(d$est, zz))
  f12 <- (map_one(zz, "Regional survey estimate\n(EHIS 2019)", lims, legend_title = "%") |
          map_one(d$est, paste0("County estimate\n(", SHORT[[sel_m]], ")"), lims, legend_title = "%") |
          map_one(d$hi90 - d$lo90, "Width of 90%\nprediction interval", NULL, pal = SEQ_RED,
                  legend_title = "Percentage points")) +
    plot_annotation(tag_levels = "A") & guides(fill = guide_colourbar(title.position = "top"))
  save_fig(f12, "fig12_app_uncertainty", 180, 72)

  ## Table: application summary per indicator
  tab_app <- g |> group_by(target) |>
    summarise(min = min(est), p10 = quantile(est, .1), median = median(est), p90 = quantile(est, .9),
              max = max(est), within_sd = sd(est - ave(est, woj)), mean_pi_width = mean(hi90 - lo90)) |>
    left_join(ht |> group_by(target) |> summarise(region_min = min(value), region_max = max(value),
                                                  mean_se = mean(se)), by = "target") |>
    mutate(indicator = HEALTH_LABELS_EN[target], o = readr::parse_number(target)) |> arrange(o) |> select(-o)
  write_csv(tab_app, file.path(PATH_RES, "tab_app_summary.csv"))
}
message("Figures written to ", PATH_FIG)
