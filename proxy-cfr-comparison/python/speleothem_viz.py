"""
Interactive DoD2k speleothem maps (Plotly) — matches DoD2k archive × proxy legend.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

YEAR_MIN, YEAR_MAX = 0, 2000

# DoD2k-style symbols (see dod2k plot notebooks / documentation)
SPELEO_PROXY_STYLE = {
    "d18O": {"symbol": "triangle-up", "color": "#b2182b", "label": "δ18O"},
    "d13C": {"symbol": "triangle-down", "color": "#d6604d", "label": "δ13C"},
    "growth rate": {"symbol": "square", "color": "#f4a582", "label": "growth rate"},
    "Mg/Ca": {"symbol": "circle", "color": "#92c5de", "label": "Mg/Ca"},
}

DEFAULT_STYLE = {"symbol": "diamond", "color": "#969696", "label": "other"}


def _series_label(site: str, proxy_type: str, dataset_id: str) -> str:
    style = _style_for(proxy_type)
    name = (site or "").strip()
    if not name or name == dataset_id or name == "nan":
        name = dataset_id.replace("sisal_", "").replace("iso2k_", "").replace("pages2k_", "")
    return f"{name} · {style['label']}"


def speleothem_records(df: pd.DataFrame) -> pd.DataFrame:
    """One row per speleothem series with coverage stats."""
    sp = df.loc[df["archiveType"] == "Speleothem"].copy()
    rows: list[dict] = []

    for _, rec in sp.iterrows():
        years = np.asarray(rec["year"], dtype=float)
        values = np.asarray(rec["paleoData_values"], dtype=float)
        mask = (years >= YEAR_MIN) & (years <= YEAR_MAX) & np.isfinite(values)
        if not np.any(mask):
            continue
        yr = years[mask].astype(int)
        rows.append(
            {
                "dataset_id": str(rec["datasetId"]),
                "site": str(rec.get("geo_siteName", "") or rec["datasetId"]),
                "proxy_type": str(rec["paleoData_proxy"]),
                "interpretation": str(rec.get("interpretation_variable", "")),
                "lon": float(rec["geo_meanLon"]),
                "lat": float(rec["geo_meanLat"]),
                "elev_m": float(rec["geo_meanElev"]) if pd.notna(rec.get("geo_meanElev")) else np.nan,
                "year_start": int(yr.min()),
                "year_end": int(yr.max()),
                "n_years": int(len(yr)),
                "span_years": int(yr.max() - yr.min() + 1),
                "database": str(rec.get("originalDatabase", "")),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    order = ["d18O", "d13C", "growth rate", "Mg/Ca"]
    out["proxy_type"] = pd.Categorical(out["proxy_type"], categories=order, ordered=True)
    return out.sort_values(["proxy_type", "dataset_id"]).reset_index(drop=True)


def speleothem_metadata_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Archive rows in DoD2k (matches official dod2k legend)."""
    sp = df.loc[df["archiveType"] == "Speleothem"]
    return (
        sp.groupby("paleoData_proxy")
        .size()
        .reset_index(name="n_metadata")
        .rename(columns={"paleoData_proxy": "proxy_type"})
    )


def speleothem_type_counts(records: pd.DataFrame) -> pd.DataFrame:
    """Series with valid observations in 0–2000 CE."""
    return (
        records.groupby("proxy_type", observed=False)
        .size()
        .reset_index(name="n_with_data")
        .sort_values("n_with_data", ascending=False)
    )


def speleothem_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    meta = speleothem_metadata_counts(df)
    rec = speleothem_type_counts(speleothem_records(df))
    labels = {
        "d18O": "δ18O",
        "d13C": "δ13C",
        "growth rate": "growth rate",
        "Mg/Ca": "Mg/Ca",
    }
    order = ["d18O", "d13C", "growth rate", "Mg/Ca"]
    out = meta.merge(rec, on="proxy_type", how="outer").fillna(0)
    out["proxy_type"] = pd.Categorical(out["proxy_type"], categories=order, ordered=True)
    out = out.sort_values("proxy_type").astype({"n_metadata": int, "n_with_data": int})
    out["proxy_label"] = out["proxy_type"].map(labels)
    return out[["proxy_label", "proxy_type", "n_metadata", "n_with_data"]]


def _style_for(proxy: str) -> dict:
    return SPELEO_PROXY_STYLE.get(proxy, {**DEFAULT_STYLE, "label": proxy})


def build_speleothem_map(records: pd.DataFrame) -> go.Figure:
    """Interactive world map — symbol & colour by paleoData_proxy."""
    if records.empty:
        raise ValueError("No speleothem records to plot.")

    records = records.copy()
    records["marker_symbol"] = records["proxy_type"].map(lambda p: _style_for(p)["symbol"])
    records["marker_color"] = records["proxy_type"].map(lambda p: _style_for(p)["color"])
    records["proxy_label"] = records["proxy_type"].map(lambda p: _style_for(p)["label"])
    records["hover"] = records.apply(
        lambda r: (
            f"<b>{r['site']}</b><br>"
            f"{r['dataset_id']}<br>"
            f"Proxy: {r['proxy_label']}<br>"
            f"Years: {r['year_start']}–{r['year_end']} (n={r['n_years']})<br>"
            f"Interpretation: {r['interpretation']}<br>"
            f"Source: {r['database']}"
        ),
        axis=1,
    )

    fig = go.Figure()

    for proxy, grp in records.groupby("proxy_type", observed=False):
        style = _style_for(str(proxy))
        fig.add_trace(
            go.Scattergeo(
                lon=grp["lon"],
                lat=grp["lat"],
                mode="markers",
                name=f"Speleothem: {style['label']} (n={len(grp)})",
                marker=dict(
                    size=9,
                    symbol=style["symbol"],
                    color=style["color"],
                    line=dict(width=0.6, color="#333333"),
                    opacity=0.92,
                ),
                text=grp["hover"],
                hoverinfo="text",
            )
        )

    fig.update_geos(
        projection_type="natural earth",
        showland=True,
        landcolor="#f0f0f0",
        showocean=True,
        oceancolor="#e6f2fa",
        showcountries=True,
        countrycolor="#bbbbbb",
        coastlinecolor="#888888",
        lataxis=dict(range=[-40, 75]),
    )

    n_types = records["proxy_type"].nunique()
    n_series = len(records)

    fig.update_layout(
        template="plotly_white",
        title=dict(
            text=(
                f"<b>DoD2k speleothem archive</b> — {n_series} series, "
                f"{n_types} proxy types"
            ),
            x=0.5,
            xanchor="center",
            font=dict(size=18),
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-0.08,
            xanchor="center",
            x=0.5,
            font=dict(size=12),
        ),
        margin=dict(l=0, r=0, t=70, b=90),
        height=560,
    )
    return fig


def build_speleothem_timeline(
    records: pd.DataFrame,
    max_rows: int = 80,
    min_span_years: int = 50,
) -> go.Figure:
    """Gantt-style timeline of temporal coverage (longest series first)."""
    if records.empty:
        raise ValueError("No speleothem records for timeline.")

    plot_df = records.loc[records["span_years"] >= min_span_years].copy()
    if plot_df.empty:
        plot_df = records.nlargest(min(max_rows, len(records)), "span_years").copy()
    elif len(plot_df) > max_rows:
        plot_df = plot_df.nlargest(max_rows, "span_years")

    plot_df = plot_df.sort_values(["proxy_type", "year_start"], ascending=[True, True])
    plot_df["label"] = plot_df.apply(
        lambda r: _series_label(r["site"], r["proxy_type"], r["dataset_id"]),
        axis=1,
    )

    fig = go.Figure()
    for _, row in plot_df.iterrows():
        style = _style_for(row["proxy_type"])
        fig.add_trace(
            go.Bar(
                x=[row["span_years"]],
                y=[row["label"]],
                base=[row["year_start"]],
                orientation="h",
                marker=dict(color=style["color"], line=dict(width=0.4, color="#333333")),
                hovertemplate=(
                    f"<b>{row['site']}</b><br>"
                    f"{row['dataset_id']}<br>"
                    f"{row['year_start']}–{row['year_end']} CE<br>"
                    f"n = {row['n_years']}<extra></extra>"
                ),
                showlegend=False,
            )
        )

    fig.update_layout(
        template="plotly_white",
        title=dict(
            text=(
                f"<b>Temporal coverage</b> — {len(plot_df)} longest series "
                f"(≥ {min_span_years} yr span)"
            ),
            x=0.5,
            font=dict(size=16),
        ),
        xaxis=dict(title="Year (CE)", range=[YEAR_MIN, YEAR_MAX], dtick=250),
        yaxis=dict(title="", automargin=True, tickfont=dict(size=10), categoryorder="array",
                   categoryarray=plot_df["label"].tolist()),
        barmode="overlay",
        height=max(480, 14 * len(plot_df)),
        margin=dict(l=10, r=20, t=50, b=40),
    )
    return fig


def build_speleothem_dashboard(df: pd.DataFrame) -> go.Figure:
    """Combined map + timeline in one interactive figure."""
    records = speleothem_records(df)
    map_fig = build_speleothem_map(records)
    time_fig = build_speleothem_timeline(records, max_rows=60)

    fig = make_subplots(
        rows=2,
        cols=1,
        row_heights=[0.55, 0.45],
        vertical_spacing=0.06,
        specs=[[{"type": "geo"}], [{"type": "xy"}]],
        subplot_titles=("Geographic distribution", "Time coverage (selected series)"),
    )

    for trace in map_fig.data:
        fig.add_trace(trace, row=1, col=1)

    for trace in time_fig.data:
        fig.add_trace(trace, row=2, col=1)

    fig.update_geos(projection_type="natural earth", row=1, col=1)
    fig.update_layout(
        template="plotly_white",
        title=dict(
            text="<b>DoD2k speleothems</b> — proxy types & time series",
            x=0.5,
            font=dict(size=20),
        ),
        height=980,
        legend=dict(orientation="h", y=1.02, x=0.5, xanchor="center"),
        margin=dict(t=90),
    )
    fig.update_xaxes(title_text="Year (CE)", row=2, col=1)
    return fig
