# Methods of Climate Reconstruction

Repository for a climate reconstruction course and an independent **proxy CFR comparison** Quarto project.

## Layout

```
Methods-of-Climate-Reconstruction/
├── course/                      # Course data & exercises
│   ├── data/
│   ├── exercises/
│   └── materials/
└── proxy-cfr-comparison/        # Quarto analysis project
    ├── _quarto.yml
    ├── index.qmd
    ├── R/
    └── analysis/
```

| Part | Description |
|------|-------------|
| `course/` | Course data, exercises, `Master_CFR.R` |
| `proxy-cfr-comparison/` | DoD2k speleothem vs. tree-ring CFR comparison (Quarto) |

## Quick start

**Analysis project (R):**

```r
setwd("proxy-cfr-comparison")
source("R/load_project.R")
run_cfr_experiment("dod2k_speleothems", "course_split")
```

**Quarto website:**

```bash
cd proxy-cfr-comparison
quarto render
```

## Data

- Course datasets: `course/data/`
- DoD2k v2.0 compact CSV: `course/data/dod2k_v2.0/` ([documentation](https://lluecke.github.io/dod2k/))

## License

Course materials and data may be subject to their original licenses (CRUTEM, PAGES, DoD2k, etc.). Check source attributions in the course files.
