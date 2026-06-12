"""
DoD2k v2.0 helpers — used from Quarto {python} chunks.

Docs: https://lluecke.github.io/dod2k/
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover
    pa = None
    pq = None

YEAR_MIN, YEAR_MAX = 0, 2000


def _candidate_dod2k_repos() -> list[Path]:
    env = os.environ.get("DOD2K_REPO")
    roots: list[Path] = []
    if env:
        roots.append(Path(env))
    here = Path(__file__).resolve()
    roots.extend([
        here.parent.parent / "external" / "dod2k",
        here.parent.parent.parent / "external" / "dod2k",
    ])
    return [p for p in roots if (p / "dod2k_utilities").is_dir()]


def ensure_dod2k_utilities() -> bool:
    for repo in _candidate_dod2k_repos():
        repo_str = str(repo.resolve())
        if repo_str not in sys.path:
            sys.path.insert(0, repo_str)
        try:
            import dod2k_utilities.ut_functions  # noqa: F401
            return True
        except ImportError:
            continue
    return False


def _parse_array_string(text: str) -> np.ndarray:
    if not text or (isinstance(text, float) and np.isnan(text)):
        return np.array([], dtype=np.float32)
    return np.array([float(x) for x in str(text).split(",") if x.strip()], dtype=np.float32)


def _as_float_array(x) -> np.ndarray:
    """Normalize a compact-CSV cell to a 1-D float array."""
    if isinstance(x, np.ndarray):
        return np.asarray(x, dtype=np.float32).ravel()
    if isinstance(x, (list, tuple)):
        if len(x) == 1 and isinstance(x[0], np.ndarray):
            return np.asarray(x[0], dtype=np.float32).ravel()
        return np.asarray(x, dtype=np.float32).ravel()
    return _parse_array_string(x)


def _read_compact_array_column(path: Path, key: str) -> pd.DataFrame:
    ids: list[str] = []
    data: list[np.ndarray] = []
    with path.open(newline="") as f:
        reader = csv.reader(f)
        for ii, row in enumerate(reader):
            if ii == 0:
                continue
            ids.append(row[0])
            data.append(_parse_array_string(",".join(row[1:])))
    return pd.DataFrame({key: data}, index=ids)


def _load_compact_pandas(data_dir: Path) -> pd.DataFrame:
    prefix = data_dir.name
    meta = pd.read_csv(
        data_dir / f"{prefix}_compact_metadata.csv",
        index_col=0,
        keep_default_na=False,
    )
    values = _read_compact_array_column(
        data_dir / f"{prefix}_compact_paleoData_values.csv",
        "paleoData_values",
    )
    years = _read_compact_array_column(
        data_dir / f"{prefix}_compact_year.csv",
        "year",
    )
    df = meta.join(values).join(years)
    df["datasetId"] = df.index.astype(str)
    df.index = range(len(df))
    df["year"] = df["year"].map(_as_float_array)
    df["paleoData_values"] = df["paleoData_values"].map(_as_float_array)
    df.name = prefix
    return df


def load_dod2k(data_dir: str | Path) -> pd.DataFrame:
    """Load DoD2k v2.0 compact CSV bundle."""
    data_dir = Path(data_dir).resolve()
    if not data_dir.is_dir():
        raise FileNotFoundError(f"DoD2k data directory not found: {data_dir}")

    if ensure_dod2k_utilities():
        from dod2k_utilities.ut_functions import load_compact_dataframe_from_csv

        parent = data_dir.parent
        rel = "/" + data_dir.name
        template = f"{data_dir.name}_compact_%s"
        old = os.getcwd()
        os.chdir(parent)
        try:
            return load_compact_dataframe_from_csv(
                data_dir.name,
                readfrom=(rel, template),
            )
        finally:
            os.chdir(old)

    return _load_compact_pandas(data_dir)


def summarize_archives(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["archiveType", "paleoData_proxy"], observed=True)
        .size()
        .reset_index(name="n")
        .sort_values("n", ascending=False)
    )


def filter_dod2k(
    df: pd.DataFrame,
    archive_types: Optional[Iterable[str]] = None,
    proxy_types: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    out = df.copy()
    if archive_types is not None:
        out = out[out["archiveType"].isin(list(archive_types))]
    if proxy_types is not None:
        out = out[out["paleoData_proxy"].isin(list(proxy_types))]
    return out.reset_index(drop=True)


def to_timeseries_tibble(df: pd.DataFrame) -> pd.DataFrame:
    """
    Long-format catalog for R (tidyverse): one row per (dataset_id, year).

    Columns: dataset_id, year, value, lon, lat, archive_type, proxy_type,
    interpretation
    """
    rows: list[dict] = []
    for _, rec in df.iterrows():
        years = np.asarray(rec["year"], dtype=float)
        values = np.asarray(rec["paleoData_values"], dtype=float)
        mask = (years >= YEAR_MIN) & (years <= YEAR_MAX)
        if not np.any(mask):
            continue
        base = {
            "dataset_id": str(rec["datasetId"]),
            "lon": float(rec["geo_meanLon"]),
            "lat": float(rec["geo_meanLat"]),
            "archive_type": str(rec["archiveType"]),
            "proxy_type": str(rec["paleoData_proxy"]),
            "interpretation": str(rec["interpretation_variable"]),
        }
        for y, v in zip(years[mask].astype(int), values[mask]):
            rows.append({**base, "year": int(y), "value": float(v)})

    return pd.DataFrame(rows)


def _write_parquet_bytes(tbl: pd.DataFrame, out: Path) -> None:
    """
    Write Parquet via in-memory buffer.

    Avoids pyarrow LocalFileSystem registration issues in RStudio/reticulate
    (ArrowKeyError: scheme 'file' already registered).
    """
    import io

    buf = io.BytesIO()
    tbl.to_parquet(buf, index=False, engine="pyarrow", compression="snappy")
    out.write_bytes(buf.getvalue())


def write_dod2k_cache(df: pd.DataFrame, cache_dir: str | Path) -> Path:
    """
    Export DoD2k to Parquet for R chunks (Quarto language bridge via files).

    Writes ``dod2k_timeseries.parquet`` under *cache_dir* (CSV.gz fallback).
    """
    if pq is None:
        raise ImportError("pyarrow is required: pip install pyarrow")

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    parquet_out = cache_dir / "dod2k_timeseries.parquet"
    csv_out = cache_dir / "dod2k_timeseries.csv.gz"
    tbl = to_timeseries_tibble(df)

    try:
        _write_parquet_bytes(tbl, parquet_out)
        return parquet_out
    except Exception as parquet_err:
        tbl.to_csv(csv_out, index=False, compression="gzip")
        raise RuntimeError(
            f"Parquet export failed ({parquet_err!r}); wrote CSV fallback to {csv_out}"
        ) from parquet_err
