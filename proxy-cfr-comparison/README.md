# DoD2k speleothem overview

Single file: **`analysis.qmd`** — load DoD2k, overview plots. All R code is in the setup chunk at the top.

## Requirements

```bash
Rscript -e "install.packages(c('tidyverse', 'plotly'), repos='https://cloud.r-project.org')"
```

## Run

```bash
cd proxy-cfr-comparison
quarto render analysis.qmd
```

Or open `analysis.qmd` in RStudio and **Run All**.

Data: `../course/data/dod2k_v2.0/` (compact CSV). Cache: `output/cache/dod2k_timeseries.rds`.
