# Data loading — CRUTEM, PAGES, N-TREND (DoD2k via Python cache; see dod2k_cache.R)

#' Load CRUTEM JJA anomaly grid (1901–2015)
load_crutem_jja_anomaly <- function(paths = project_paths()) {
  f <- file.path(paths$data_dir, "crutem_anom_jjamean_1901-2015.Rdata")
  env <- rlang::env()
  load(f, envir = env)
  rlang::env_get(env, "cru_anom_jjamean")
}

#' Crop a gridded list object to a lat/lon box
crop_grid_domain <- function(grid, domain) {
  lon_idx <- grid$lon >= domain$lon_w & grid$lon <= domain$lon_e
  lat_idx <- grid$lat >= domain$lat_s & grid$lat <= domain$lat_n

  grid$data <- grid$data[lon_idx, lat_idx, , drop = FALSE]
  grid$lon <- grid$lon[lon_idx]
  grid$lat <- grid$lat[lat_idx]
  grid
}

#' Build reconstruction target for Master_CFR.R
build_reconstruction_target <- function(domain = DOMAIN_NHET_JJA, paths = project_paths()) {
  load_crutem_jja_anomaly(paths) |>
    crop_grid_domain(domain)
}

#' Load raw PAGES 2.0 annual proxy table + metadata
load_pages2k <- function(paths = project_paths()) {
  proxy <- readr::read_table(
    file.path(paths$data_dir, "PAGES_proxy_ann_2.0.0.txt"),
    show_col_types = FALSE
  )
  meta <- readr::read_table(
    file.path(paths$data_dir, "PAGES_metadata_2.0.0.txt"),
    show_col_types = FALSE
  )

  list(
    years = proxy[[1L]],
    values = as.matrix(proxy[, -1L, drop = FALSE]),
    meta = meta,
    names = as.character(meta[1L, -1L]),
    archives = as.character(meta[4L, -1L]),
    lon = as.numeric(meta[3L, -1L]),
    lat = as.numeric(meta[2L, -1L])
  )
}

#' Load N-TREND 2015 tree-ring network
load_ntrend <- function(paths = project_paths()) {
  proxy <- readr::read_delim(
    file.path(paths$data_dir, "N-TREND2015_data.csv"),
    delim = ";",
    show_col_types = FALSE
  )
  meta <- readr::read_delim(
    file.path(paths$data_dir, "N-TREND2015_transformed_coordinates.csv"),
    delim = ";",
    show_col_types = FALSE
  )

  list(
    years = proxy[[1L]],
    values = as.matrix(proxy[, -1L, drop = FALSE]),
    lon = meta$Long_transformed,
    lat = meta$Lat_transformed
  )
}

#' Convert internal proxy table to Master_CFR `proxydata` list
as_proxydata <- function(years, values, lon, lat) {
  stopifnot(length(lon) == ncol(values), length(lat) == ncol(values))
  list(
    data = values,
    lon = lon,
    lat = lat,
    time = years
  )
}
