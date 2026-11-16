# Data licences and attribution

The files in `data/` are extracts of public data, reformatted for this study. They remain subject to the terms of the original providers. When you reuse them, cite the provider as below and cite the paper.

| Files | Provider | Terms |
|---|---|---|
| `ghana_monthly_cases.csv`, `ghana_annual_surveillance.csv`, `ghana_coverage_wuenic.csv`, `ghana_sia.csv`, `ghana_outbreak_age_vaccination.csv`, WUENIC column of `sources/dhs_vs_wuenic_mcv1_west_africa.csv` | World Health Organization (Immunization Data portal and WIISE); WUENIC are joint WHO/UNICEF estimates | Reused with acknowledgement for non-commercial purposes under WHO's terms (CC BY-NC-SA 3.0 IGO, https://creativecommons.org/licenses/by-nc-sa/3.0/igo/). Adaptations must carry the same or equivalent licence. WHO does not endorse this work. |
| `ghana_demography.csv` | The World Bank, World Development Indicators | CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/) |
| `Value` column of `sources/dhs_vs_wuenic_mcv1_west_africa.csv` | The DHS Program (ICF), STATcompiler/API aggregate indicators | Cite The DHS Program; see https://dhsprogram.com/ for terms of use. |
| `niamey_weekly.csv` | Grais et al. (2008), as distributed in the R package `epimdr` by O.N. Bjornstad | GPL-3; full text in `data/NIAMEY_DATA_LICENSE_GPL-3.txt` |

The model-specific columns that are not source data (the `eta` column of `ghana_sia.csv`) and the files in `published/` are the authors' work and are released under CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/), which keeps them compatible with the WHO-derived inputs.

Suggested source citations:

- World Health Organization. Immunization Data portal: provisional measles and rubella data; WIISE data tables. Geneva: WHO. Accessed 6 November 2026. Licence: CC BY-NC-SA 3.0 IGO.
- WHO/UNICEF estimates of national immunization coverage (WUENIC), Ghana. Accessed 6 November 2026.
- The World Bank. World Development Indicators: SP.POP.TOTL, SP.DYN.CBRT.IN, SP.DYN.LE00.IN. Accessed 6 November 2026. Licence: CC BY 4.0.
- The DHS Program STATcompiler/API, indicator CH_VACC_C_MSL. ICF. Accessed 6 November 2026.
- Grais RF, Conlan AJK, Ferrari MJ, et al. Time is of the essence: exploring a measles outbreak response vaccination in Niamey, Niger. J R Soc Interface 2008;5:67-74. doi:10.1098/rsif.2007.1038.
- Bjornstad ON. epimdr: Functions and Data for "Epidemics: Models and Data in R". R package version 0.6-5, 2020. https://CRAN.R-project.org/package=epimdr
