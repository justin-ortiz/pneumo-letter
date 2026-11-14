# Childhood pneumococcal vaccination and invasive disease in Florida adults

Code and results for a research letter on invasive pneumococcal disease (IPD)
averted in Florida adults aged 50 years or older by childhood pneumococcal
conjugate vaccine (PCV).

The analysis takes IPD rates for 1998-2024 from CDC's Active Bacterial Core
surveillance (ABCs), compares each year with the 1998-1999 rate before
childhood PCV, and applies the difference to the Florida population. Results
are averaged over 2011-2013 and over 2023-2024, with 95% confidence intervals
from a Monte Carlo simulation. The methods are in the letter and its
Supplement.

## What you need

R, version 4 or later. The scripts use base R only, so no packages need to be
installed. I ran them in RStudio on Windows.

## Folders

| Folder | What it holds |
|---|---|
| `R` | The four scripts, run in order |
| `data-raw` | The files downloaded from CDC and the Census Bureau, unchanged. `SOURCES.md` lists each one. |
| `output` | Everything the scripts write |

## How to rerun the analysis

1. Download this repository with the green Code button, then Download ZIP, and unzip it.
2. In RStudio, set the working directory to the unzipped folder, using
   Session, then Set Working Directory, then Choose Directory.
3. Run the scripts in order by typing each line in the console:

```
source("R/01_rates.R")
source("R/02_florida.R")
source("R/04_output.R")
```

`03_coverage.R` calculates childhood PCV coverage from the National
Immunization Survey files, which are too large for GitHub. The two files it
writes are already in `output`, so the other scripts run without it. To rerun
it, see `SOURCES.md`.

## What each script does

| Script | What it does | What it writes to `output` |
|---|---|---|
| `01_rates.R` | Calculates yearly IPD rates in the ABCs areas for the 13 PCV13 serotypes and 22 comparison serotypes, and the reduction from the 1998-1999 rate | `rates.csv` |
| `02_florida.R` | Applies the reductions to the Florida population, adds non-bacteremic pneumonia and deaths, and calculates 95% confidence intervals from 10,000 simulation draws. Also calculates cases averted by serotype group. | `florida_by_year.csv`, `florida_periods.csv`, `florida_deaths.csv`, `florida_deaths_by_year.csv`, `serotype_groups.csv` |
| `03_coverage.R` | Calculates PCV coverage in children in Florida and the US, pooled over adjacent survey years | `nis_florida_pooled.csv`, `nis_us_pooled.csv` |
| `04_output.R` | Makes the Table, the Figure and eFigure 1, and lists the numbers quoted in the text | `table1.txt`, `figure.png`, `figure.pdf`, `efigure1.png`, `efigure1.pdf`, `text_numbers.txt` |

`table1.txt` is tab-delimited, so it can be pasted into Word.

## Getting the same numbers

`02_florida.R` sets the random seed to 2026, so the confidence intervals come
out the same on every run. Changing the order of the random draws in that
script would change the intervals slightly.
