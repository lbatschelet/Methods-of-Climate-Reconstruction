# Data loading — CRUTEM target grids and proxy archives

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Load CRUTEM JJA anomaly grid (1901–2015)
load_crutem_jja_anomaly <- function(paths = project_paths()) {
  f <- file.path(paths$data_dir, "crutem_anom_jjamean_1901-2015.Rdata")
  env <- new.env(parent = emptyenv())
  load(f, envir = env)
  get("cru_anom_jjamean", envir = env)
}

#' Crop a gridded list object to a lat/lon box
crop_grid_domain <- function(grid, domain) {
  lon_idx <- which(grid$lon >= domain$lon_w & grid$lon <= domain$lon_e)
  lat_idx <- which(grid$lat >= domain$lat_s & grid$lat <= domain$lat_n)
  grid$data <- grid$data[lon_idx, lat_idx, , drop = FALSE]
  grid$lon <- grid$lon[lon_idx]
  grid$lat <- grid$lat[lat_idx]
  grid
}

#' Build reconstruction target for Master_CFR.R
build_reconstruction_target <- function(domain = DOMAIN_NHET_JJA, paths = project_paths()) {
  grid <- load_crutem_jja_anomaly(paths)
  crop_grid_domain(grid, domain)
}

#' Load raw PAGES 2.0 annual proxy table + metadata
load_pages2k <- function(paths = project_paths()) {
  proxy <- utils::read.table(
    file.path(paths$data_dir, "PAGES_proxy_ann_2.0.0.txt"),
    header = TRUE
  )
  meta <- utils::read.table(
    file.path(paths$data_dir, "PAGES_metadata_2.0.0.txt"),
    sep = "\t",
    stringsAsFactors = FALSE
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
  proxy <- utils::read.table(
    file.path(paths$data_dir, "N-TREND2015_data.csv"),
    sep = ";",
    header = TRUE
  )
  meta <- utils::read.csv(
    file.path(paths$data_dir, "N-TREND2015_transformed_coordinates.csv"),
    sep = ";"
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

#' Path to official DoD2k v2.0 compact CSV bundle
dod2k_v2_dir <- function(paths = project_paths()) {
  file.path(paths$data_dir, "dod2k_v2.0")
}

#' Load DoD2k v2.0 metadata table (LiPDverse `archiveType`, `paleoData_proxy`, …)
#'
#' See Evans et al. (2026) ESSD and https://lluecke.github.io/dod2k/
load_dod2k_metadata <- function(paths = project_paths()) {
  meta_path <- file.path(dod2k_v2_dir(paths), "dod2k_v2.0_compact_metadata.csv")
  if (!file.exists(meta_path)) {
    stop(
      "DoD2k metadata not found at ", meta_path, ". ",
      "Copy the compact CSV bundle from the dod2k repository (data/dod2k_v2.0/)."
    )
  }
  meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
  row.names(meta) <- meta$datasetId
  meta
}

#' Parse one row of a DoD2k compact values/year CSV into a numeric vector
parse_dod2k_compact_row <- function(fields) {
  vals <- suppressWarnings(as.numeric(fields))
  vals[!is.na(vals)]
}

#' Read a DoD2k compact column CSV (paleoData_values or year)
read_dod2k_compact_column <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2L) {
    stop("Empty compact CSV: ", path)
  }
  header <- strsplit(lines[[1L]], ",", fixed = TRUE)[[1L]]
  if (length(header) < 2L || header[[2L]] == "") {
    stop("Unexpected compact CSV header in ", path)
  }
  ids <- character(length(lines) - 1L)
  arrays <- vector("list", length(lines) - 1L)
  for (i in seq_along(arrays)) {
    parts <- strsplit(lines[[i + 1L]], ",", fixed = TRUE)[[1L]]
    ids[[i]] <- parts[[1L]]
    arrays[[i]] <- parse_dod2k_compact_row(parts[-1L])
  }
  names(arrays) <- ids
  arrays
}

#' Load official DoD2k v2.0 compact dataframe (4781 records, 1–2000 CE)
#'
#' Mirrors `dod2k_utilities.ut_functions.load_compact_dataframe_from_csv`.
load_dod2k_v2 <- function(paths = project_paths()) {
  base <- dod2k_v2_dir(paths)
  value_path <- file.path(base, "dod2k_v2.0_compact_paleoData_values.csv")
  year_path <- file.path(base, "dod2k_v2.0_compact_year.csv")
  if (!file.exists(value_path) || !file.exists(year_path)) {
    stop(
      "DoD2k v2.0 compact data missing under ", base, ". ",
      "See Data/dod2k_v2.0/ in the dod2k GitHub repository."
    )
  }

  meta <- load_dod2k_metadata(paths)
  values_by_id <- read_dod2k_compact_column(value_path)
  years_by_id <- read_dod2k_compact_column(year_path)
  ids <- intersect(intersect(meta$datasetId, names(values_by_id)), names(years_by_id))
  if (length(ids) < 2L) {
    stop("Fewer than 2 DoD2k v2.0 records could be joined from compact CSV files.")
  }

  # Align sparse compact series to the course/DoD2k CE grid (0–2000)
  years_ref <- seq.int(0L, 2000L)
  n_time <- length(years_ref)
  mat <- matrix(NA_real_, nrow = n_time, ncol = length(ids))
  colnames(mat) <- ids
  lon <- lat <- archive_type <- proxy_type <- interpretation <- character(length(ids))

  for (j in seq_along(ids)) {
    id <- ids[[j]]
    yr <- round(years_by_id[[id]])
    vals <- values_by_id[[id]]
    if (length(yr) != length(vals)) {
      stop("Year/value length mismatch for record ", id)
    }
    in_range <- yr >= years_ref[[1L]] & yr <= years_ref[[n_time]]
    if (!any(in_range)) {
      warning("No years in 0–2000 CE for record ", id, call. = FALSE)
      next
    }
    row_idx <- match(yr[in_range], years_ref)
    mat[row_idx, j] <- vals[in_range]
    m <- meta[id, ]
    lon[[j]] <- as.numeric(m$geo_meanLon)
    lat[[j]] <- as.numeric(m$geo_meanLat)
    archive_type[[j]] <- m$archiveType
    proxy_type[[j]] <- m$paleoData_proxy
    interpretation[[j]] <- m$interpretation_variable
  }

  list(
    variant = "v2",
    years = years_ref,
    values = mat,
    lon = lon,
    lat = lat,
    archive_type = archive_type,
    proxy_type = proxy_type,
    interpretation = interpretation,
    type = proxy_type,
    names = ids,
    metadata = meta[ids, , drop = FALSE]
  )
}

#' Attach DoD2k v2.0 metadata to a course RData subset (matched by `datasetId`)
attach_dod2k_metadata <- function(dod, paths = project_paths()) {
  meta <- load_dod2k_metadata(paths)
  idx <- match(dod$names, meta$datasetId)
  if (any(is.na(idx))) {
    missing <- dod$names[is.na(idx)]
    stop(
      length(missing), " record(s) in variant '", dod$variant,
      "' are absent from DoD2k v2.0 metadata, e.g. ",
      paste(head(missing, 3L), collapse = ", ")
    )
  }
  m <- meta[idx, ]
  dod$archive_type <- m$archiveType
  dod$proxy_type <- m$paleoData_proxy
  dod$interpretation <- m$interpretation_variable
  dod$metadata <- m
  dod
}

#' Load DoD2k proxy collection
#'
#' @param variant `v2` (official compact CSV, recommended), or course subsets
#'   `mean` / `tm` / `full` from `Data/DoD2k.RData` (pre-filtered; `mean` has no
#'   speleothems — see DoD2k documentation).
load_dod2k <- function(variant = c("v2", "mean", "tm", "full"), paths = project_paths()) {
  variant <- match.arg(variant)
  if (variant == "v2") {
    return(load_dod2k_v2(paths))
  }

  env <- new.env(parent = emptyenv())
  load(file.path(paths$data_dir, "DoD2k.RData"), envir = env)
  raw <- get(DOD2K_VARIANTS[[variant]], envir = env)

  dod <- list(
    variant = variant,
    years = raw$time,
    values = as.matrix(raw$data),
    lon = raw$lon,
    lat = raw$lat,
    type = raw$type,
    names = colnames(raw$data)
  )
  attach_dod2k_metadata(dod, paths)
}

#' Tabulate DoD2k records by archive and proxy measurement type
summarize_dod2k_archives <- function(variant = "v2", paths = project_paths()) {
  dod <- load_dod2k(variant, paths)
  tab <- as.data.frame(table(
    archive = dod$archive_type,
    proxy = dod$proxy_type,
    stringsAsFactors = FALSE
  ))
  tab <- tab[tab$Freq > 0L, ]
  tab[order(tab$Freq, decreasing = TRUE), c("archive", "proxy", "Freq")]
}

#' Backwards-compatible alias (proxy measurement types only — prefer `summarize_dod2k_archives`)
summarize_dod2k_types <- function(variant = "v2", paths = project_paths()) {
  dod <- load_dod2k(variant, paths)
  counts <- sort(table(dod$proxy_type), decreasing = TRUE)
  data.frame(
    proxy_type = names(counts),
    n = as.integer(counts),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Subset DoD2k using official metadata fields
#'
#' Filter on `archiveType` (e.g. Speleothem, Wood) and/or `paleoData_proxy`
#' (e.g. d18O, ring width). Do not equate d18O with speleothem — tree and coral
#' archives also carry d18O in DoD2k (Evans et al. 2026).
subset_dod2k <- function(dod,
                         archive_types = NULL,
                         proxy_types = NULL,
                         interpretation = NULL,
                         domain = NULL) {
  keep <- rep(TRUE, ncol(dod$values))

  if (!is.null(archive_types)) {
    keep <- keep & dod$archive_type %in% archive_types
  }
  if (!is.null(proxy_types)) {
    keep <- keep & dod$proxy_type %in% proxy_types
  }
  if (!is.null(interpretation)) {
    keep <- keep & dod$interpretation %in% interpretation
  }
  if (!is.null(domain)) {
    keep <- keep &
      dod$lat >= domain$lat_s & dod$lat <= domain$lat_n &
      dod$lon >= domain$lon_w & dod$lon <= domain$lon_e
  }

  if (!any(keep)) {
    stop(
      "No DoD2k records match filters. ",
      "archive_types=", paste(archive_types %||% "*", collapse = ", "),
      "; proxy_types=", paste(proxy_types %||% "*", collapse = ", ")
    )
  }

  list(
    years = dod$years,
    values = dod$values[, keep, drop = FALSE],
    lon = dod$lon[keep],
    lat = dod$lat[keep],
    archive_type = dod$archive_type[keep],
    proxy_type = dod$proxy_type[keep],
    interpretation = dod$interpretation[keep],
    type = dod$proxy_type[keep],
    names = dod$names[keep]
  )
}

#' @rdname subset_dod2k
subset_dod2k_types <- function(dod, proxy_types) {
  subset_dod2k(dod, proxy_types = proxy_types)
}
