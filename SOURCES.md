# Where the data came from

One row per file in `data-raw`. The dates are the dates I downloaded the
files.

## Invasive pneumococcal disease

Read by `01_rates.R` and `02_florida.R`.

| File | Source | Downloaded |
|---|---|---|
| `abcs_serotype_1998_2024.csv` | CDC Active Bacterial Core surveillance, 1998-2024 serotype data for invasive pneumococcal disease cases by age group and site. https://data.cdc.gov/Public-Health-Surveillance/1998-2024-Serotype-Data-for-Invasive-Pneumococcal-/qvzb-qs6p/about_data | 2026-11-04 |
| `abcs_aggregate_spn_2026-11-04.csv` | CDC Active Bacterial Core surveillance, *Streptococcus pneumoniae* rates. https://data.cdc.gov/Public-Health-Surveillance/Active-Bacterial-Core-surveillance-ABCs-Streptococ/en3s-hzsr/about_data | 2026-11-04 |

The serotype file gives case counts by serotype, age group, surveillance site
and year. The aggregate file gives case rates and death rates by age group and
year, and rates for serotype groups. The analysis uses the all-serotype rate
and the PCV13 rate from its Serotypes topic, and the case and death rates by
age.

## Florida population

Read by `02_florida.R`. No single Census file covers 2000-2024.

| File | Covers | Source | Downloaded |
|---|---|---|---|
| `st-est00int-agesex.csv` | 2000-2009 | Census Bureau intercensal estimates by single year of age and sex for states, 2000-2010. https://www.census.gov/data/datasets/time-series/demo/popest/intercensal-2000-2010-state.html | 2026-11-05 |
| `sc-est2020int-agesex-12.csv` | 2010-2019 | Census Bureau intercensal estimates for selected age groups by sex for Florida, 2010-2020. https://www.census.gov/data/tables/time-series/demo/popest/intercensal-2010-2020-state-detail.html | 2026-11-06 |
| `sc-est2025-syasex-12.csv` | 2020-2024 | Census Bureau Vintage 2025 estimates by single year of age and sex for Florida. https://www.census.gov/data/tables/time-series/demo/popest/2020s-state-detail.html | 2026-11-04 |

## PCV coverage in children

Read by `03_coverage.R`. These files are not in the repository because of
their size.

| File | Source | Downloaded |
|---|---|---|
| `nis/NISPUFyy.DAT` | CDC National Immunization Survey-Child public use data files, survey years 2005-2024 | 2026-11-05 |
| `nis/NISPUFyy.R` or `nis/NISPUFyy.sas` | CDC's program for reading each year's data file | 2026-11-05 |

To rerun `03_coverage.R`, download each survey year's data file and its R or
SAS program from CDC's National Immunization Surveys data pages
(https://www.cdc.gov/nis/php/datasets-child/index.html) and put them in
`data-raw/nis`. The public use files before 2005 have no pneumococcal
variables.

## Numbers entered in the scripts

`02_florida.R` holds three sets of counts that are not in a data file.

| Numbers | Source |
|---|---|
| 60 pneumonia and 28 invasive episodes | Bonten et al, N Engl J Med 2015, Table 2, placebo group, per-protocol, vaccine-type first episodes |
| Pneumonia and invasive disease cases by age group | LeBlanc et al, Vaccine 2022, Table 1 |
| 30-day deaths by age group, 7 of 244 and 31 of 439 | Provided by the LeBlanc et al investigators on 9 November 2026, reported in eMethods 4 of the Supplement |
