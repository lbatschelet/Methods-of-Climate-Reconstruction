# Proxy CFR Comparison (Quarto)

Quarto project comparing climate field reconstructions across proxy networks, with DoD2k v2.0 speleothems.

## Structure

```
proxy-cfr-comparison/
├── _quarto.yml          # Website
├── index.qmd            # Landing page
├── R/                   # Analysis pipeline
├── config/              # Experiment definitions
└── analysis/
    ├── report.qmd
    └── includes/        # Shared code chunks
```

Course materials live in **`../course/`** (data + exercises).

## Requirements

- R: `ggplot2`, `ncdf4`, `foreach`, `doParallel`, `mvnfast`
- [Quarto](https://quarto.org/)

## Quick start

```r
setwd("proxy-cfr-comparison")
source("R/load_project.R")

summarize_dod2k_archives("v2")
run_cfr_experiment("dod2k_speleothems", "course_split")
```

## Render website

```bash
cd proxy-cfr-comparison
quarto render
quarto preview
```

Set `quick_run: false` in the report YAML params for the full experiment batch.

## DoD2k metadata

Filter speleothems with **`archiveType = Speleothem`** and **`paleoData_proxy = d18O`** — see [DoD2k docs](https://lluecke.github.io/dod2k/).

Data: `../course/data/dod2k_v2.0/` (official compact CSV bundle).
