# Python (DoD2k)

Used by `{python}` chunks in `analysis.qmd`.

## One-command setup

From `proxy-cfr-comparison/`:

```bash
bash scripts/setup_python.sh
```

This creates `.venv`, installs `requirements.txt`, and registers the Jupyter kernel **`dod2k-cfr`**. Quarto picks up `.venv/bin/python` via `_environment` (`QUARTO_PYTHON`).

## Modules

| File | Purpose |
|------|---------|
| `project_env.py` | Project root & data paths for Quarto chunks |
| `dod2k_io.py` | Load compact CSV, filter, Parquet cache |
| `speleothem_viz.py` | Interactive Plotly map & timeline |

## Optional: dod2k_utilities

```bash
bash scripts/setup_dod2k_repo.sh
```

Clones `external/dod2k` for the official loader; compact CSV in `course/data/dod2k_v2.0/` works without it.
