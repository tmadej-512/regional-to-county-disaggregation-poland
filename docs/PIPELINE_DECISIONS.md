# Decision logic

Decisions are fixed before the final runs and passed downstream through files in `output/results/`.

| Stage | Evidence produced | Decision | Used by |
|---|---|---|---|
| 01 Data | coherence of published regional vs aggregated county covariates | regional covariates always computed as population-weighted means of county values (X̄ = AX) | all stages |
| 01 Data | within- to between-region variance ratio of reference county variables | prior location for within-region variance (`TAU_RATIO`); leave-one-variable-out in validation | 02, 04–07 |
| 03–04 Simulation | variance ratio and regional Q2 of real variables | calibration of the regional component (λ, ψ) for the main design | 04, 05b–05d |
| 04–06 Simulation and validation | county-level MSE reduction vs naive allocation in 44 cells (S6 excluded a priori as a limits test) | **method selection by minimax regret** → `selected_method.txt` | 07, 08 |
| 06c Calibration | GDR_XL / GP_XL skill over a grid of c on real county variables, with leave-one-target-out check | **cross-level prior c** → `xl_c_prior.csv` (simulation keeps the pre-specified c0 = 0.5) | 07 |
| 06 Validation | skill by covariate set | application covariate set `izp` (the pre-specified short list `core` did not outperform it) | 07 |
| 07 Application | Q2, signal-to-noise, agreement between methods, tercile change between c = 0.25 and c = 1, agreement with external county data | traffic-light reliability per indicator (`ocena_dezagregacji`); signal-to-noise thresholds 0.12 (red) and 0.30 (green), set from the real-variable validation | reporting |

## Rules

- **Validation target.** Methods are judged on county-level truth (simulation or real county variables), never on re-aggregation to regions.
- **Method selection.** Regret in a cell = gap in MSE reduction to the best feasible method in that cell. The selected method minimises the maximum regret; mean regret breaks ties.
- **Cross-level prior.** c is the share of the between-region covariate effect assumed to hold within regions; c = 1 reproduces standard GDR. It cannot be estimated from regional data and is taken from external real county variables.
- **Use in a composite index.** County EHIS estimates are reported as a scenario with intervals; the main version of the index uses regional values for EHIS components, because the disaggregation covariates are themselves index components.
