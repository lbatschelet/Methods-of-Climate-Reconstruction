# DoD2k speleothem overview

Single file: **`analysis.qmd`** — load DoD2k, overview plots, proxy–climate correlations. Code is grouped by section in the notebook.

## Requirements

```bash
Rscript -e "install.packages(c('tidyverse', 'plotly', 'terra'), repos='https://cloud.r-project.org')"
```

## Run

```bash
cd proxy-cfr-comparison
quarto render analysis.qmd
```

Or open `analysis.qmd` in RStudio and **Run All**.

**Data**

- DoD2k: `../course/data/dod2k_v2.0/` (compact CSV)
- CRU TS v4.07: downloaded on first run to `output/cache/cru_ts/` (~400 MB compressed; needs network + `gunzip`)

**Cache**

- `output/cache/dod2k_timeseries.rds` — long-format proxy catalog
- `output/cache/cru_site_climate.rds` — extracted grid-box climate at speleothem sites
- `output/cache/speleo_instrument_correlations.rds` — correlation results
