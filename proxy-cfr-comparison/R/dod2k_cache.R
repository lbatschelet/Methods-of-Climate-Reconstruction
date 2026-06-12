# Read DoD2k catalog exported by Quarto {python} chunks

#' Path to DoD2k Parquet cache (written by Python)
dod2k_cache_path <- function(paths = project_paths()) {
  file.path(paths$cache_dir, "dod2k_timeseries.parquet")
}

#' Read long-format DoD2k timeseries (tibble)
read_dod2k_catalog <- function(paths = project_paths()) {
  path <- dod2k_cache_path(paths)
  if (!file.exists(path)) {
    stop(
      "DoD2k cache not found at ", path, ". ",
      "Run the Python chunk in analysis.qmd first (exports Parquet)."
    )
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required to read the DoD2k cache.")
  }

  arrow::read_parquet(path) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      dataset_id = as.character(.data$dataset_id),
      archive_type = as.character(.data$archive_type),
      proxy_type = as.character(.data$proxy_type),
      interpretation = as.character(.data$interpretation),
      year = as.integer(.data$year),
      value = as.double(.data$value),
      lon = as.double(.data$lon),
      lat = as.double(.data$lat)
    )
}

#' Summarise archive × proxy counts from catalog
summarize_dod2k_catalog <- function(catalog = read_dod2k_catalog()) {
  catalog |>
    dplyr::distinct(.data$dataset_id, .data$archive_type, .data$proxy_type) |>
    dplyr::count(.data$archive_type, .data$proxy_type, name = "n", sort = TRUE)
}

#' Filter catalog and pivot to wide proxy matrix for Master_CFR
catalog_to_proxy_block <- function(catalog,
                                   archive_types = NULL,
                                   proxy_types = NULL,
                                   year_min = 0L,
                                   year_max = 2000L) {
  filtered <- catalog |>
    dplyr::filter(
      .data$year >= year_min,
      .data$year <= year_max
    )

  if (!is.null(archive_types)) {
    filtered <- dplyr::filter(filtered, .data$archive_type %in% archive_types)
  }
  if (!is.null(proxy_types)) {
    filtered <- dplyr::filter(filtered, .data$proxy_type %in% proxy_types)
  }

  if (nrow(filtered) == 0L) {
    stop("No DoD2k records match the requested filters.")
  }

  meta <- filtered |>
    dplyr::distinct(
      .data$dataset_id,
      .data$lon,
      .data$lat,
      .data$archive_type,
      .data$proxy_type,
      .data$interpretation
    )

  wide <- filtered |>
    tidyr::pivot_wider(
      id_cols = "year",
      names_from = "dataset_id",
      values_from = "value",
      values_fn = mean
    ) |>
    dplyr::arrange(.data$year)

  years <- wide$year
  values <- as.matrix(dplyr::select(wide, -"year"))

  list(
    years = years,
    values = values,
    lon = meta$lon[match(colnames(values), meta$dataset_id)],
    lat = meta$lat[match(colnames(values), meta$dataset_id)],
    archive_type = meta$archive_type[match(colnames(values), meta$dataset_id)],
    proxy_type = meta$proxy_type[match(colnames(values), meta$dataset_id)],
    interpretation = meta$interpretation[match(colnames(values), meta$dataset_id)],
    type = meta$proxy_type[match(colnames(values), meta$dataset_id)],
    names = colnames(values)
  )
}
