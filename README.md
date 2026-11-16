# Two-dose measles model with saturated treatment, calibrated to Ghana

Data and MATLAB code for the manuscript

> O.D. Makinde, I.Y. Seini, B. Seidu, S. Osman. *Modeling and Optimizing Control of Measles in West Africa: a two-dose model with saturated treatment, calibrated to Ghana.*

The repository holds every dataset used in the paper, the posterior sample from the Bayesian fit, and MATLAB/Octave code that reproduces the numerical results, tables and figures.

## Contents

| Folder | What it holds |
|---|---|
| `data/` | Input data for Ghana (monthly cases, demography, routine coverage, campaigns, case ages) and the Niamey 2003-04 outbreak. See `data/README.md` for sources and variable definitions. |
| `data/sources/` | Supporting data used only in the sensitivity analysis (DHS and WUENIC first-dose coverage, West Africa). |
| `published/` | Posterior draws and posterior summary from the fit reported in Section 8, and a reference file used to check the MATLAB estimation code against the original Python run. |
| `core/`, `sec*/` | MATLAB/Octave code, one folder per part of the paper (see below). |
| `results/`, `figures/` | Empty; the scripts write their output here. |

## Data in brief

| File | Content | Period | Source |
|---|---|---|---|
| `data/ghana_monthly_cases.csv` | Measles cases by final classification, suspected and discarded cases | Jan 2012 to May 2026 | WHO Immunization Data portal, provisional measles and rubella data |
| `data/ghana_demography.csv` | Population, crude birth rate, life expectancy | 2000 to 2027 (blank where not yet published) | World Bank, World Development Indicators |
| `data/ghana_coverage_wuenic.csv` | MCV1 and MCV2 coverage (%) | 1980 to 2025 | WHO/UNICEF estimates of national immunization coverage (WUENIC) |
| `data/ghana_sia.csv` | National measles campaigns: date, post-campaign survey coverage, target ages, model share parameter | 2013, 2018, 2024 | WHO supplementary immunization activities database |
| `data/ghana_outbreak_age_vaccination.csv` | Outbreak-reported cases by age group and vaccination status | 2018 to 2024 | WHO outbreak reporting (WIISE) |
| `data/ghana_annual_surveillance.csv` | Annual suspected, confirmed and discarded cases, incidence, discarded rate | 2012 to 2026 | WHO Immunization Data portal, surveillance performance file |
| `data/niamey_weekly.csv` | Weekly cases and deaths, Niamey (Niger) outbreak | Nov 2003 to Jun 2004 | Grais et al. (2008), via the R package `epimdr` |
| `data/sources/dhs_vs_wuenic_mcv1_west_africa.csv` | MCV1 from 63 DHS/MIS surveys matched to WUENIC | 1986 to 2024 | DHS Program API; WUENIC |

All WHO and World Bank data were retrieved on 6 November 2026. Provisional WHO counts are revised over time, so later downloads may differ.

## Running the code

The code needs no MATLAB toolboxes. It was tested in GNU Octave 8.4 and is written to run in current MATLAB releases.

```matlab
setup_paths            % add all folders to the path
run_all                % every section; set QUICK = true in run_all.m for a fast check
```

Each script can also be run on its own after `setup_paths`. With `QUICK = true` the whole run takes a few minutes. The full run takes several hours, most of it in the MCMC, profiles and sensitivity refits (Section 8); the optimal-control strategies (Section 10) take about an hour in Octave.

| Folder | Results in the manuscript |
|---|---|
| `sec03_thresholds/` | Control reproduction number and derived thresholds (Sections 3 and 7), elasticities (Table 3), Figure 5 |
| `sec05_equilibria/` | Numerical check of the equilibrium classification (Theorem 5.3), Figure 6 |
| `sec06_bifurcation/` | Hopf thresholds and first Lyapunov coefficients (Table 2), Findings 6.10 to 6.12, Figures 2 to 4 and 7 to 9 |
| `sec08_09_estimation/` | Data summaries, likelihood and priors, posterior mode, MCMC, profiles, validation, sensitivity refits and the Niamey fit (Section 8, Table 4, Figure 10) |
| `sec10_simulations/` | Numerical simulations (Section 9, Findings 9.1 to 9.5, Tables 5 and 6, Figures 11 to 16) |
| `sec11_optimal_control/` | Constant and optimal controls, cost-effectiveness (Section 10, Tables 7 and 8, Figures 17 to 19) |

The MCMC and the global sensitivity analysis use fixed random seeds, so repeated runs give the same numbers. Results from MATLAB and Octave agree with the original Python analysis to within Monte Carlo and optimiser tolerance; `published/python_reference_estimation.json` holds the posterior mode and its log-likelihood and log-posterior from the Python run, for comparison.

## Licences

- Code: MIT licence (`LICENSE`).
- Data: see `LICENSE-DATA.md`. Files derived from WHO sources are shared under CC BY-NC-SA 3.0 IGO terms, World Bank data under CC BY 4.0, and the Niamey data under GPL-3 as distributed in `epimdr`. Please cite the original sources listed in `data/README.md` as well as the paper.

## Citation

If you use these data or code, please cite the paper (details in `CITATION.cff`) and the original data sources.
