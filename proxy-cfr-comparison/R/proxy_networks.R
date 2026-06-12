# Assemble and quality-control proxy networks (tidyverse)

#' Align proxy matrix rows to a master year sequence
align_proxy_years <- function(years, values, master_years) {
  idx <- match(master_years, years)
  if (any(is.na(idx))) {
    rlang::abort("Master years missing from proxy time axis.")
  }
  values[idx, , drop = FALSE]
}

#' Filter columns by calibration coverage; optionally keep top-N
filter_calibration_coverage <- function(values,
                                        years,
                                        calib_start,
                                        calib_end,
                                        min_coverage = 0.8,
                                        max_proxies = NULL) {
  calib_mask <- years >= calib_start & years <= calib_end
  if (!any(calib_mask)) {
    rlang::abort("No years inside the calibration window.")
  }

  coverage <- colMeans(!is.na(values[calib_mask, , drop = FALSE]))
  keep <- coverage >= min_coverage
  if (!any(keep)) {
    rlang::abort("No proxies pass the calibration coverage threshold.")
  }

  idx <- which(keep)
  if (!is.null(max_proxies) && length(idx) > max_proxies) {
    idx <- idx[order(coverage[idx], decreasing = TRUE)[seq_len(max_proxies)]]
    keep <- seq_along(keep) %in% idx
  }

  list(
    values = values[, keep, drop = FALSE],
    keep = keep,
    coverage = coverage[keep],
    n_dropped = sum(!keep)
  )
}

subset_pages_archives <- function(pages, archives) {
  archives <- unique(archives)
  keep <- pages$archives %in% archives
  if (!any(keep)) {
    rlang::abort("No PAGES records match the requested archives.")
  }

  list(
    years = pages$years,
    values = pages$values[, keep, drop = FALSE],
    lon = pages$lon[keep],
    lat = pages$lat[keep],
    names = pages$names[keep],
    archive_type = pages$archives[keep],
    type = pages$archives[keep]
  )
}

cbind_proxy_blocks <- function(block_a, block_b, master_years) {
  va <- align_proxy_years(block_a$years, block_a$values, master_years)
  vb <- align_proxy_years(block_b$years, block_b$values, master_years)

  list(
    years = master_years,
    values = cbind(va, vb),
    lon = c(block_a$lon, block_b$lon),
    lat = c(block_a$lat, block_b$lat),
    names = c(block_a$names, block_b$names),
    archive_type = c(block_a$archive_type, block_b$archive_type),
    type = c(block_a$type, block_b$type)
  )
}

dod2k_block_from_spec <- function(spec, catalog, domain = NULL) {
  block <- catalog_to_proxy_block(
    catalog,
    archive_types = spec$archive_types,
    proxy_types = spec$proxy_types
  )

  if (!is.null(domain)) {
    in_domain <- block$lat >= domain$lat_s &
      block$lat <= domain$lat_n &
      block$lon >= domain$lon_w &
      block$lon <= domain$lon_e
    if (!any(in_domain)) {
      rlang::abort("No DoD2k proxies inside the requested domain.")
    }
    block$values <- block$values[, in_domain, drop = FALSE]
    block$lon <- block$lon[in_domain]
    block$lat <- block$lat[in_domain]
    block$names <- block$names[in_domain]
    block$archive_type <- block$archive_type[in_domain]
    block$type <- block$type[in_domain]
    block$interpretation <- block$interpretation[in_domain]
  }

  block
}

assemble_proxy_block <- function(spec, master_years, paths, catalog = NULL, domain = NULL) {
  switch(
    spec$builder,
    dod2k = {
      if (is.null(catalog)) catalog <- read_dod2k_catalog(paths)
      dod2k_block_from_spec(spec, catalog, domain)
    },
    dod2k_combined = {
      if (is.null(catalog)) catalog <- read_dod2k_catalog(paths)
      purrr::imap(spec$proxy_groups, \(group, name) {
        block <- dod2k_block_from_spec(group, catalog, domain)
        cap <- spec$max_proxies_per_group[[name]]
        if (is.null(cap)) {
          return(block)
        }
        aligned <- align_proxy_years(block$years, block$values, master_years)
        filtered <- filter_calibration_coverage(
          aligned,
          master_years,
          spec$.calib_start,
          spec$.calib_end,
          min_coverage = spec$min_cal_coverage,
          max_proxies = cap
        )
        list(
          years = master_years,
          values = filtered$values,
          lon = block$lon[filtered$keep],
          lat = block$lat[filtered$keep],
          names = block$names[filtered$keep],
          archive_type = block$archive_type[filtered$keep],
          type = block$type[filtered$keep]
        )
      }) |>
        purrr::reduce(\(a, b) cbind_proxy_blocks(a, b, master_years))
    },
    ntrend = {
      nt <- load_ntrend(paths)
      list(
        years = nt$years,
        values = nt$values,
        lon = nt$lon,
        lat = nt$lat,
        names = colnames(nt$values),
        archive_type = rep("Wood", ncol(nt$values)),
        type = rep("ring width", ncol(nt$values))
      )
    },
    pages = subset_pages_archives(load_pages2k(paths), spec$archives),
    combined_speleothem_ntrend = {
      pages <- load_pages2k(paths)
      nt <- load_ntrend(paths)
      sp <- subset_pages_archives(pages, "speleothem")
      cbind_proxy_blocks(sp, nt, master_years)
    },
    rlang::abort(paste("Unknown builder:", spec$builder))
  )
}

build_proxy_network <- function(network_id,
                                recon_period = NULL,
                                calib_start,
                                calib_end,
                                catalog = NULL,
                                domain = NULL,
                                paths = project_paths()) {
  spec <- PROXY_NETWORKS[[network_id]]
  if (is.null(spec)) {
    rlang::abort(paste("Unknown network_id:", network_id))
  }

  recon_period <- spec$recon_period %||% recon_period %||% RECON_PERIOD
  master_years <- seq.int(recon_period$start, recon_period$end)

  spec_run <- spec
  if (spec$builder == "dod2k_combined") {
    spec_run$.calib_start <- calib_start
    spec_run$.calib_end <- calib_end
  }

  block <- assemble_proxy_block(
    spec_run, master_years, paths,
    catalog = catalog, domain = domain
  )
  aligned <- align_proxy_years(block$years, block$values, master_years)

  max_px <- if (spec$builder == "dod2k_combined") NULL else spec$max_proxies
  filtered <- filter_calibration_coverage(
    aligned,
    master_years,
    calib_start,
    calib_end,
    min_coverage = spec$min_cal_coverage,
    max_proxies = max_px
  )

  if (ncol(filtered$values) < 2L) {
    rlang::abort(paste0(
      "Fewer than 2 proxies remain for '", network_id, "' after filtering."
    ))
  }

  list(
    proxydata = as_proxydata(
      master_years,
      filtered$values,
      block$lon[filtered$keep],
      block$lat[filtered$keep]
    ),
    meta = list(
      network_id = network_id,
      label = spec$label,
      data_source = if (grepl("^dod2k", spec$builder)) "DoD2k v2.0" else spec$builder,
      n_proxies_raw = ncol(block$values),
      n_proxies_used = ncol(filtered$values),
      n_dropped = filtered$n_dropped,
      calib_window = c(calib_start, calib_end),
      recon_period = recon_period,
      proxy_names = block$names[filtered$keep],
      proxy_types = block$type[filtered$keep],
      archive_types = block$archive_type[filtered$keep]
    )
  )
}

proxy_network_map_df <- function(network) {
  tibble::tibble(
    lon = network$proxydata$lon,
    lat = network$proxydata$lat,
    name = network$meta$proxy_names,
    archive = network$meta$archive_types,
    proxy = network$meta$proxy_types,
    network = network$meta$network_id
  )
}
