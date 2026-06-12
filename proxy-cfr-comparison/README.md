# Proxy CFR Comparison

Single Quarto notebook (`analysis.qmd`): native **`{python}`** and **`{r}`** chunks.

| Language | Role |
|----------|------|
| **Python** | Load/filter DoD2k ([official toolkit](https://lluecke.github.io/dod2k/)) |
| **R** (tidyverse) | CFR, validation, ggplot |

Python writes `output/cache/dod2k_timeseries.parquet`; R reads it with **arrow**. No `reticulate` bridge inside the R modules.

## Setup

```bash
cd proxy-cfr-comparison
bash scripts/setup_python.sh          # .venv + Jupyter kernel dod2k-cfr
bash scripts/setup_dod2k_repo.sh      # optional: dod2k_utilities
Rscript -e "install.packages('arrow', repos='https://cloud.r-project.org')"
```

R packages: **tidyverse**, **arrow**, **ggplot2**, **ncdf4**, **glue**, plus CFR deps (`foreach`, `doParallel`, `mvnfast`).

Quarto uses kernel **`dod2k-cfr`** and `QUARTO_PYTHON` from `_environment` — no manual Python picker needed.

```r
source("R/load_project.R")
load_project()   # loads tidyverse via analysis.qmd setup chunk
```

## Render

```bash
bash scripts/render.sh          # sets QUARTO_PYTHON to .venv (recommended)
quarto preview analysis.qmd     # or preview in IDE after setup_python.sh
```

Run chunks top-to-bottom: **setup-r** → **dod2k-load** (Python) → R sections.
