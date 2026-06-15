# Speleothem proxy screening for climate reconstruction

Course project: explore **DoD2k v2.0** speleothems, screen against **CRU TS** instrumental grid data, then (Step 2) attempt climate field reconstruction.

## Roadmap

| Step | Goal | Status |
|------|------|--------|
| **1** | Load DoD2k; map four proxy types (δ18O, δ13C, growth rate, Mg/Ca); correlate with local grid-box T and precipitation | `analysis.qmd` |
| **2** | Detrend proxies; regional composites; cal/val split (1901–1980 / 1981–2000) | `analysis.qmd` §2 |
| **3** | ENSO screening with lags, tuned Niño index, local summer rainfall | `analysis.qmd` §3 |

## Data

- **Proxies:** [DoD2k v2.0](https://essd.copernicus.org/articles/16/1933/2024) — `../course/data/dod2k_v2.0/`
- **Instrumental:** [CRU TS v4.07](https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.07/) — auto-download to `output/cache/cru_ts/`

## Run

```bash
cd proxy-cfr-comparison
Rscript -e "install.packages(c('tidyverse','plotly','terra'), repos='https://cloud.r-project.org')"
quarto render analysis.qmd
```

Output: `output/analysis.html`

## References

Bibliography: `references.bib` (cited in `analysis.qmd` via Quarto).

- @evans2026 — DoD2k database
- @kaushal2024 — SISALv3 speleothem compilation
- @fairchild2006; @fairchild2009; @fairchild2012 — proxy interpretation
- @tadros2016 — ENSO and drip-water geochemistry
- @kost2023 — La Vallina cave monitoring
- @bernal2024 — Brazilian drought and P−PET
- @harris2020 — CRU TS instrumental grid
