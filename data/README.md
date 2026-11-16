# Data dictionary and sources

All files are comma-separated, UTF-8, with one header row. Missing values are empty fields or `NaN`. Retrieval date for every online source: 6 November 2026.

## ghana_monthly_cases.csv

National monthly measles counts for Ghana, January 2012 to May 2026 (173 months).

| Column | Meaning |
|---|---|
| `year`, `month` | Calendar year and month (1 to 12) |
| `lab_confirmed` | Laboratory-confirmed measles cases |
| `epi_linked` | Cases confirmed by epidemiological link |
| `clinical` | Clinically compatible cases |
| `suspected` | Suspected measles cases notified |
| `discarded` | Suspected cases discarded as non-measles |

- Source: WHO Immunization Data portal, provisional measles and rubella data, file `404-table-web-epi-curve-data.xlsx` (https://immunizationdata.who.int/global?topic=Provisional-measles-and-rubella-data), rows for Ghana (ISO3 `GHA`).
- November and December 2022, November and December 2024, and December 2025 are absent from the WHO file and are left empty. They are treated as missing, not zero.
- The model observation in the paper is `lab_confirmed + epi_linked + clinical`. The fit uses January 2012 to October 2024 (152 observed months) and the hold-out period is January 2025 to May 2026 (16 observed months).
- May 2026 was the latest month in the file at retrieval and may be incomplete.
- These are provisional counts. They include clinically compatible cases and can differ from official annual totals (for 2024: 4,189 in the monthly file against 1,944 officially reported).

## ghana_demography.csv

| Column | Meaning |
|---|---|
| `year` | Calendar year |
| `population` | Total population, mid-year (World Bank `SP.POP.TOTL`) |
| `cbr_per1000` | Crude birth rate per 1,000 people (`SP.DYN.CBRT.IN`) |
| `life_expectancy` | Life expectancy at birth, years (`SP.DYN.LE00.IN`) |

- Source: World Bank, World Development Indicators, via https://api.worldbank.org/v2/.
- Years not yet published are `NaN`. The code extends the birth rate by the linear trend of 2019 to 2024 and the population by its last growth rate.
- In the model, the birth input is `cbr * population` and the exit rate is `cbr - g`, where `g = ln(P[y+1]/P[y])`, so net migration is absorbed in the exit rate.

## ghana_coverage_wuenic.csv

| Column | Meaning |
|---|---|
| `year` | Calendar year |
| `MCV1` | First-dose measles-containing vaccine coverage (%) |
| `MCV2` | Second-dose coverage (%); `NaN` before MCV2 was introduced |

- Source: WHO/UNICEF estimates of national immunization coverage (WUENIC), WHO WIISE table `MT_AD_COV_LONG`, coverage category `WUENIC` (https://xmart-api-public.who.int/WIISE/).
- The model uses `p1 = MCV1/100` and conditional second-dose uptake `p2c = min(MCV2/MCV1, 0.98)`.

## ghana_sia.csv

| Column | Meaning |
|---|---|
| `year`, `month`, `day` | Campaign date used in the model |
| `survey_coverage` | Post-campaign survey coverage (proportion) |
| `eta` | Share of model susceptibles in the target ages (modelling assumption, not data) |
| `age_group` | Target ages |

- Source of dates, target ages and survey coverage: WHO supplementary immunization activities database, WIISE table `V_SIA_MAIN_MR`.
- `eta` is set by the authors: 0.9 for the 2013 catch-up (9 months to 14 years) and 0.6 for the 2018 and 2024 follow-ups (9 to 59 months), based on the age distribution of cases (65% under five years).

## ghana_outbreak_age_vaccination.csv

Cases reported through WHO outbreak reporting, by year.

| Column | Meaning |
|---|---|
| `CODE`, `YEAR` | ISO3 code, year |
| `unvaccinated`, `vaccinated`, `unknown` | Cases by vaccination status |
| `lt_1_year`, `1_to_4_years`, `5_to_9_years`, `10_to_14_years`, `gteq_15_years`, `unknown_age` | Cases by age group |
| `all` | Total cases where reported |
| `deaths` | Deaths |

- Source: WHO WIISE table `AD_OUTBREAKS`. Negative sentinel codes in the original table were set to missing; 2021 has no record.
- The prior median of the reference force of infection, 1/5.2 per year, comes from these data: the 2,621 cases of known age in 2022 to 2024 have mean age 5.2 years when each group is placed at its midpoint (0.5, 2.5, 7, 12 and 22 years, the last for the open group of 15 and over).

## ghana_annual_surveillance.csv

| Column | Meaning |
|---|---|
| `year` | Calendar year |
| `suspected` | Suspected cases |
| `lab`, `epi`, `clin` | Laboratory-confirmed, epidemiologically linked, clinically compatible |
| `total` | Total measles cases |
| `incidence_per_million` | Measles incidence per million population |
| `discarded_rate_per_100k` | Non-measles discarded rate per 100,000 population |

- Source: WHO Immunization Data portal, surveillance performance file `403-table-web-reporting-data.xlsx`. Used for descriptive checks only.

## niamey_weekly.csv

| Column | Meaning |
|---|---|
| `week_index` | Week number from the start of the series (1 to 31) |
| `week` | Epidemiological week |
| `cases` | Reported cases, Niamey city (sum of the three districts) |
| `deaths` | Reported deaths |

- Source: Grais RF, Conlan AJK, Ferrari MJ, et al. Time is of the essence: exploring a measles outbreak response vaccination in Niamey, Niger. *J R Soc Interface* 2008;5:67-74. doi:10.1098/rsif.2007.1038. Data as distributed in the R package `epimdr` (O.N. Bjornstad), dataset `niamey`, licensed GPL-3; see `NIAMEY_DATA_LICENSE_GPL-3.txt`.
- Totals: 10,880 cases and 199 deaths over 31 weeks.

## sources/dhs_vs_wuenic_mcv1_west_africa.csv

| Column | Meaning |
|---|---|
| `CODE` | ISO3 country code |
| `SurveyYear`, `SurveyType` | Survey year and type (DHS or MIS) |
| `Value` | MCV1 coverage measured by the survey (%) |
| `mcv1_wuenic` | WUENIC MCV1 for the same country and year (%) |

- Sources: DHS Program API (https://api.dhsprogram.com/), indicator `CH_VACC_C_MSL` (measles vaccination, children aged 12 to 23 months, national total); WUENIC as above.
- The 63 matched surveys cover 13 West African countries from 1986 to 2024. WUENIC exceeds the survey value by 3.2 percentage points on average. The sensitivity analysis in Section 8.5 lowers MCV1 and MCV2 by this amount.
