#!/usr/bin/env Rscript
# Offline ENSO method sweep — not part of analysis.qmd
# Tests lags, seasons, targets, filters on DoD2k speleothem × Niño3.4

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
})

root <- if (dir.exists("course/data")) "." else if (dir.exists("../course/data")) ".." else "../.."
proj <- if (file.exists("analysis.qmd")) "." else "proxy-cfr-comparison"
cache <- file.path(proj, "output", "cache")
cru_dir <- file.path(cache, "cru_ts")

y0 <- 1901L
y1 <- 2000L
min_n <- 10L

detrend_linear <- function(years, values) {
  ok <- is.finite(values)
  if (sum(ok) < 3L) return(rep(NA_real_, length(values)))
  fit <- lm(values[ok] ~ years[ok])
  pred <- predict(fit, newdata = data.frame(years = years))
  values - pred
}

assign_enso_region <- function(lon, lat) {
  case_when(
    lon >= 138 & lon <= 155 & lat >= -40 & lat <= -28 ~ "SE Australia",
    lon >= 130 & lon <= 180 & lat >= -45 & lat <= -10 ~ "Eastern Australasia",
    lon >= 110 & lon <= 180 & lat >= -30 & lat <= 5 ~ "SW Pacific / IASM",
    TRUE ~ NA_character_
  )
}

# --- load data ---
catalog <- readRDS(file.path(cache, "dod2k_timeseries.rds"))
cru_clim <- readRDS(file.path(cache, "cru_site_climate.rds"))
nino_monthly <- readRDS(file.path(cache, "nino34_index.rds"))

sites <- catalog |>
  filter(archive_type == "Speleothem") |>
  group_by(dataset_id, proxy_type, lon, lat) |>
  summarize(
    year_start = min(year), year_end = max(year),
    n_years = n(), span = max(year) - min(year) + 1L,
    pts_per_year = n() / pmax(span, 1L),
    .groups = "drop"
  ) |>
  mutate(
    enso_region = assign_enso_region(lon, lat),
    proxy_type = as.character(proxy_type)
  )

loc <- sites |> distinct(lon, lat) |> mutate(loc_id = row_number())

proxy_ann <- catalog |>
  filter(archive_type == "Speleothem", proxy_type %in% c("Mg/Ca", "d18O")) |>
  group_by(dataset_id, proxy_type, lon, lat, year) |>
  summarize(value = mean(value), .groups = "drop") |>
  inner_join(loc, by = c("lon", "lat")) |>
  inner_join(sites |> dplyr::select(dataset_id, enso_region, n_years, span, pts_per_year),
             by = "dataset_id")

# --- Niño3.4 indices ---
nino_idx <- function(monthly) {
  djf <- monthly |>
    filter(month_num %in% c(12L, 1L, 2L)) |>
    mutate(climate_year = if_else(month_num == 12L, year + 1L, year)) |>
    group_by(climate_year) |>
    summarize(nino_djf = mean(nino34, na.rm = TRUE), .groups = "drop") |>
    rename(year = climate_year)

  son <- monthly |>
    filter(month_num %in% c(9L, 10L, 11L)) |>
    group_by(year) |>
    summarize(nino_son = mean(nino34, na.rm = TRUE), .groups = "drop")

  jja <- monthly |>
    filter(month_num %in% c(6L, 7L, 8L)) |>
    group_by(year) |>
    summarize(nino_jja = mean(nino34, na.rm = TRUE), .groups = "drop")

  ann <- monthly |>
    group_by(year) |>
    summarize(nino_ann = mean(nino34, na.rm = TRUE), .groups = "drop")

  reduce(list(djf, son, jja, ann), full_join, by = "year")
}

nino <- nino_idx(nino_monthly)

# lag helper: positive lag = proxy leads nino (proxy year Y vs nino Y+lag)
lag_nino <- function(df, col, lag) {
  nino_lag <- nino |>
    dplyr::select(year, !!col) |>
    mutate(year = year - lag)
  inner_join(df, nino_lag, by = "year")
}

# --- seasonal CRU precipitation at sites (monthly nc) ---
seasonal_cru <- function(loc) {
  pre_nc <- file.path(cru_dir, "cru_ts4.07.1901.2022.pre.dat.nc")
  if (!file.exists(pre_nc)) {
    message("No CRU pre nc — skipping seasonal targets")
    return(NULL)
  }
  pts <- vect(loc, geom = c("lon", "lat"), crs = "EPSG:4326")
  mat <- as.matrix(
    terra::extract(rast(pre_nc, subds = "pre"), pts, method = "simple")[, -1]
  )
  years <- 1901:2022
  map_dfr(seq_len(nrow(loc)), function(i) {
    m <- mat[i, ]
    tibble(
      loc_id = loc$loc_id[i],
      year = years,
      pre_ann = vapply(years, function(yr) {
        j <- (yr - 1901L) * 12L + 1:12
        sum(m[j], na.rm = TRUE)
      }, numeric(1)),
      # Austral summer DJF labelled to year Y = Dec(Y-1)+Jan(Y)+Feb(Y)
      pre_djf = vapply(years, function(yr) {
        dec <- (yr - 1901L) * 12L + 12L
        if (yr == 1901L) dec <- NA_integer_
        jan <- (yr - 1901L) * 12L + 1L
        feb <- (yr - 1901L) * 12L + 2L
        sum(m[c(dec, jan, feb)], na.rm = TRUE)
      }, numeric(1)),
      pre_mam = vapply(years, function(yr) {
        j <- (yr - 1901L) * 12L + 3:5
        sum(m[j], na.rm = TRUE)
      }, numeric(1)),
      pre_jja = vapply(years, function(yr) {
        j <- (yr - 1901L) * 12L + 6:8
        sum(m[j], na.rm = TRUE)
      }, numeric(1))
    )
  })
}

clim_season <- seasonal_cru(loc)
clim_full <- cru_clim |>
  dplyr::select(loc_id, year, tmp_ann, tmp_jja, pre_ann) |>
  left_join(clim_season |> dplyr::select(loc_id, year, pre_djf, pre_mam, pre_jja),
            by = c("loc_id", "year"))

# --- correlation engine ---
corr_one <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < min_n) return(list(r = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- cor.test(x[ok], y[ok], method = "pearson")
  list(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

corr_spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < min_n) return(list(r = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  list(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

partial_r <- function(x, y, z) {
  ok <- complete.cases(x, y, z)
  if (sum(ok) < min_n + 2L) return(list(r = NA_real_, p = NA_real_, n = sum(ok)))
  rx <- residuals(lm(x[ok] ~ z[ok]))
  ry <- residuals(lm(y[ok] ~ z[ok]))
  corr_one(rx, ry)
}

first_diff <- function(x) c(NA_real_, diff(x))

screen_series <- function(df, nino_col, clim_cols = character(),
                          lags = -2:2, method = "pearson",
                          partial_col = NULL, transform = "detrend") {
  map_dfr(lags, function(lag) {
    d <- df |>
      filter(year >= y0, year <= y1)
    d <- lag_nino(d, nino_col, lag)

    if (transform == "detrend") {
      d <- d |> mutate(px = detrend_linear(year, value))
    } else if (transform == "diff") {
      d <- d |> arrange(year) |> mutate(px = first_diff(detrend_linear(year, value)))
    } else {
      d <- d |> mutate(px = value)
    }

    targets <- c(nino_col, clim_cols)
    map_dfr(targets, function(target) {
      yy <- d[[target]]
      xx <- d$px
      if (!is.null(partial_col) && target == nino_col) {
        res <- partial_r(xx, yy, d[[partial_col]])
        meth <- paste0("partial|", partial_col)
      } else if (method == "spearman") {
        res <- corr_spearman(xx, yy)
        meth <- "spearman"
      } else {
        res <- corr_one(xx, yy)
        meth <- "pearson"
      }
      tibble(
        lag = lag, target = target, method = meth,
        r = res$r, p = res$p, n = res$n
      )
    })
  })
}

# build per-series results for method grid
enso_pool <- proxy_ann |> filter(!is.na(enso_region))

nino_cols <- c("nino_djf", "nino_son", "nino_jja", "nino_ann")
clim_cols <- c("pre_ann", "pre_djf", "pre_mam", "pre_jja", "tmp_ann", "tmp_jja")
lags <- -2:2

cat("=== ENSO method sweep (offline) ===\n")
cat("Series in ENSO regions (Mg/Ca + d18O):", n_distinct(enso_pool$dataset_id), "\n\n")

results <- map_dfr(unique(enso_pool$dataset_id), function(id) {
  row <- enso_pool |> filter(dataset_id == id) |> slice(1)
  d0 <- enso_pool |>
    filter(dataset_id == id) |>
    left_join(clim_full, by = c("loc_id", "year"))

  base <- map_dfr(nino_cols, function(ncol) {
    screen_series(d0, ncol, clim_cols = character(), lags = lags)
  }) |> mutate(variant = "detrend_pearson")

  partial <- map_dfr(c("nino_djf", "nino_son"), function(ncol) {
    screen_series(d0, ncol, lags = -1:1, partial_col = "pre_ann")
  }) |> mutate(variant = "partial_pre_ann")

  spearman <- map_dfr(c("nino_djf", "nino_son"), function(ncol) {
    screen_series(d0, ncol, lags = -1:1, method = "spearman")
  }) |> mutate(variant = "spearman")

  diffed <- map_dfr(c("nino_djf", "nino_son"), function(ncol) {
    screen_series(d0, ncol, lags = -1:1, transform = "diff")
  }) |> mutate(variant = "first_diff")

  # climate targets at lag 0 and best nino lag per index
  clim_at0 <- screen_series(d0, "nino_djf", clim_cols = clim_cols, lags = 0) |>
    filter(target %in% clim_cols) |>
    mutate(variant = "vs_climate_lag0")

  bind_rows(
    base |> filter(target %in% nino_cols),
    partial |> filter(target %in% c("nino_djf", "nino_son")),
    spearman |> filter(target %in% c("nino_djf", "nino_son")),
    diffed |> filter(target %in% c("nino_djf", "nino_son")),
    clim_at0
  ) |>
    mutate(
      dataset_id = id,
      proxy_type = row$proxy_type,
      enso_region = row$enso_region,
      pts_per_year = row$pts_per_year
    )
})

# --- summaries ---
best_per_series <- results |>
  filter(grepl("^nino", target), variant %in% c("detrend_pearson", "spearman", "partial_pre_ann", "first_diff")) |>
  group_by(dataset_id, proxy_type, enso_region, variant) |>
  slice_max(abs(r), n = 1, with_ties = FALSE) |>
  ungroup()

cat("--- Current Step 3 baseline (lag 0, nino_djf, detrend pearson) ---\n")
baseline <- results |>
  filter(variant == "detrend_pearson", target == "nino_djf", lag == 0) |>
  arrange(desc(abs(r)))
print(baseline |> dplyr::select(dataset_id, enso_region, proxy_type, r, p, n))
cat("\nSignificant (p<0.05):", sum(baseline$p < 0.05, na.rm = TRUE), "\n\n")

cat("--- Best lag per series (detrend pearson, any nino season) ---\n")
best_lag <- results |>
  filter(variant == "detrend_pearson", grepl("^nino", target)) |>
  group_by(dataset_id) |>
  slice_max(abs(r), n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(abs(r)))
print(best_lag |> dplyr::select(dataset_id, enso_region, proxy_type, target, lag, r, p))
cat("\nSignificant (p<0.05):", sum(best_lag$p < 0.05, na.rm = TRUE), "\n")
cat("Improved |r| vs baseline (same series):",
    round(median(abs(best_lag$r) - abs(baseline$r[match(best_lag$dataset_id, baseline$dataset_id)]), na.rm = TRUE), 3), "\n\n")

cat("--- Method comparison: max |r| across all variants (nino targets only) ---\n")
method_best <- results |>
  filter(grepl("^nino", target)) |>
  mutate(method_label = paste(variant, target, "lag", lag)) |>
  group_by(dataset_id, proxy_type, enso_region) |>
  slice_max(abs(r), n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(abs(r)))
print(head(method_best, 15) |>
        dplyr::select(dataset_id, enso_region, proxy_type, method_label, r, p))
cat("\nSignificant (p<0.05) any method:", sum(method_best$p < 0.05, na.rm = TRUE), "\n\n")

cat("--- Proxy vs local climate (often stronger than Niño?) lag 0 ---\n")
clim_best <- results |>
  filter(variant == "vs_climate_lag0") |>
  group_by(dataset_id) |>
  slice_max(abs(r), n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(abs(r)))
print(head(clim_best, 12) |> dplyr::select(dataset_id, enso_region, proxy_type, target, r, p))

cat("\n--- High-resolution filter (pts_per_year >= 1) ---\n")
hi <- method_best |> filter(pts_per_year >= 1)
cat("n series:", nrow(hi), "| sig:", sum(hi$p < 0.05, na.rm = TRUE), "\n")
if (nrow(hi)) print(hi |> dplyr::select(dataset_id, proxy_type, method_label, r, p))

cat("\n--- Partial correlation: proxy vs Niño | controlling local pre_ann ---\n")
partial_sum <- results |>
  filter(variant == "partial_pre_ann", target %in% c("nino_djf", "nino_son")) |>
  group_by(dataset_id) |>
  slice_max(abs(r), n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(abs(r)))
print(head(partial_sum, 10) |> dplyr::select(dataset_id, enso_region, proxy_type, target, lag, r, p))
cat("Significant:", sum(partial_sum$p < 0.05, na.rm = TRUE), "\n")

cat("\n--- Climate targets with lag sweep (best per series) ---\n")
clim_lag <- results |>
  filter(variant == "vs_climate_lag0", target %in% clim_cols) |>
  group_by(dataset_id) |>
  slice_max(abs(r), n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(abs(r)))
print(head(clim_lag, 12) |> dplyr::select(dataset_id, enso_region, target, lag, r, p))
cat("Significant vs climate:", sum(clim_lag$p < 0.05, na.rm = TRUE), "\n")

cat("\n--- Composite-style: mean z-score by region×proxy, best nino setting ---\n")
composite_test <- function(region, ptype, nino_col = "nino_djf", lag = 0L) {
  ids <- enso_pool |>
    filter(enso_region == region, proxy_type == ptype) |>
    pull(dataset_id) |> unique()
  if (length(ids) < 2L) return(NULL)

  ts <- map_dfr(ids, function(id) {
    d <- enso_pool |>
      filter(dataset_id == id, year >= y0, year <= y1)
    d <- lag_nino(d, nino_col, lag) |>
      mutate(px = detrend_linear(year, value))
    tibble(year = d$year, px = d$px, id = id)
  })

  wide <- ts |>
    pivot_wider(names_from = id, values_from = px) |>
    left_join(nino |> dplyr::select(year, all_of(nino_col)), by = "year")

  ids_ok <- intersect(ids, setdiff(names(wide), c("year", nino_col)))
  if (length(ids_ok) < 2L) return(NULL)

  signs <- baseline |>
    filter(dataset_id %in% ids_ok) |>
    arrange(match(dataset_id, ids_ok)) |>
    mutate(s = sign(replace_na(r, 1))) |>
    pull(s)

  mat <- as.matrix(wide |> dplyr::select(all_of(ids_ok)))
  for (j in seq_along(ids_ok)) {
    if (signs[j] < 0) mat[, j] <- -mat[, j]
  }
  comp <- rowMeans(mat, na.rm = TRUE)
  ref_y0 <- 1901L; ref_y1 <- 1980L
  zscore <- function(x, yrs) {
    ref <- x[yrs >= ref_y0 & yrs <= ref_y1]
    (x - mean(ref, na.rm = TRUE)) / sd(ref, na.rm = TRUE)
  }
  yrs <- wide$year
  nino_v <- wide[[nino_col]]
  ok <- complete.cases(comp, nino_v)
  cal <- ok & yrs >= 1901 & yrs <= 1980
  val <- ok & yrs >= 1981 & yrs <= 2000
  tibble(
    region = region, proxy_type = ptype, nino_col = nino_col, lag = lag,
    n = length(ids_ok),
    r_cal = cor(zscore(comp, yrs)[cal], zscore(nino_v, yrs)[cal]),
    r_val = cor(zscore(comp, yrs)[val], zscore(nino_v, yrs)[val]),
    r_full = cor(zscore(comp[ok], yrs[ok]), zscore(nino_v[ok], yrs[ok]))
  )
}

comp_grid <- crossing(
  region = unique(enso_pool$enso_region),
  ptype = c("d18O", "Mg/Ca"),
  nino_col = c("nino_djf", "nino_son", "nino_ann"),
  lag = c(-1L, 0L, 1L)
) |>
  pmap_dfr(composite_test) |>
  arrange(desc(abs(r_val)))

print(comp_grid)

cat("\n--- SW Pacific d18O: per-series best vs composite (validation r) ---\n")
swp <- method_best |> filter(enso_region == "SW Pacific / IASM", proxy_type == "d18O")
print(swp)

out_rds <- file.path(cache, "enso_method_sweep_results.rds")
saveRDS(list(results = results, method_best = method_best, comp_grid = comp_grid), out_rds)
cat("\nSaved:", out_rds, "\n")
