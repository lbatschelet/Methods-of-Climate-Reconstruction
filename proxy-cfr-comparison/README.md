# Proxy CFR Comparison

Single **Quarto notebook** (`analysis.qmd`) for exploring and comparing CFR experiments.  
Course data and exercises live in **`../course/`**.

## Workflow (like Jupyter)

1. Open `analysis.qmd` in RStudio, VS Code, or Positron.
2. Run chunks interactively as you work.
3. Render when you want the HTML report:

```bash
cd proxy-cfr-comparison
quarto preview analysis.qmd   # live preview while editing
quarto render analysis.qmd    # write output/analysis.html
```

In RStudio: use **Run** on chunks and the **Render** button.

## Parameters

Edit the YAML header in `analysis.qmd`:

| Param | Default | Meaning |
|-------|---------|---------|
| `quick_run` | `true` | 3 experiments; `false` = full batch |
| `overwrite_cfr` | `false` | Re-run CFR even if cached |

## R modules

Reusable code stays in `R/` and `config/` — the notebook sources them via `load_project()`.

```r
source("R/load_project.R")
load_project()
run_cfr_experiment("dod2k_speleothems", "course_split")
```

## Requirements

R: `ggplot2`, `ncdf4`, `foreach`, `doParallel`, `mvnfast`  
[Quarto](https://quarto.org/)
