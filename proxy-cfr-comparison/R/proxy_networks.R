# Assemble and quality-control proxy networks

#' Align a proxy matrix to a master year sequence
align_proxy_years <- function(years, values, master_years) {
  idx <- match(master_years, years)
  if (any(is.na(idx))) {
    missing <- master_years[is.na(idx)]
    stop("Master years not found in proxy time axis: ", paste(range(missing), collapse = "-"))
  }
  values[idx, , drop = FALSE]
}

#' Drop columns with insufficient calibration coverage; optionally keep top-N by coverage
filter_calibration_coverage <- function(values, years, calib_start, calib_end,
                                        min_coverage = 0.8,
                                        max_proxies = NULL) {
  calib_rows <- years >= calib_start & years <= calib_end
  if (!any(calib_rows)) {
    stop("No years fall inside calibration window.")
  }
  coverage <- colMeans(!is.na(values[calib_rows, , drop = FALSE]))
  keep <- coverage >= min_coverage
  if (!any(keep)) {
    stop("No proxies pass the calibration coverage threshold.")
  }

  idx <- which(keep)
  if (!is.null(max_proxies) && length(idx) > max_proxies) {
    ord <- order(coverage[idx], decreasing = TRUE)[seq_len(max_proxies)]
    idx <- idx[ord]
    keep <- rep(FALSE, length(keep))
    keep[idx] <- TRUE
  }

  list(
    values = values[, keep, drop = FALSE],
    keep = keep,
    coverage = coverage[keep],
    n_dropped = sum(!coverage >= min_coverage) + max(0L, sum(coverage >= min_coverage) - sum(keep))
  )
}

#' Subset PAGES 2.0 by archive type(s)
subset_pages_archives <- function(pages, archives) {
  archives <- unique(archives)
  keep <- pages$archives %in% archives
  if (!any(keep)) {
    stop("No PAGES records match archives: ", paste(archives, collapse = ", "))
  }
  list(
    years = pages$years,
    values = pages$values[, keep, drop = FALSE],
    lon = pages$lon[keep],
    lat = pages$lat[keep],
    names = pages$names[keep],
    archives = pages$archives[keep],
    type = pages$archives[keep]
  )
}

#' Horizontal bind of two proxy blocks on a shared timeline
cbind_proxy_blocks <- function(block_a, block_b, master_years) {
  va <- align_proxy_years(block_a$years, block_a$values, master_years)
  vb <- align_proxy_years(block_b$years, block_b$values, master_years)
  list(
    years = master_years,
    values = cbind(va, vb),
    lon = c(block_a$lon, block_b$lon),
    lat = c(block_a$lat, block_b$lat),
    names = c(block_a$names, block_b$names),
    archive_type = c(
      block_a$archive_type %||% block_a$type %||% rep(NA_character_, ncol(va)),
      block_b$archive_type %||% block_b$type %||% rep(NA_character_, ncol(vb))
    ),
    type = c(block_a$type, block_b$type)
  )
}

#' Apply DoD2k archive/proxy filters from a network spec
subset_dod2k_spec <- function(dod, spec, domain = NULL) {
  subset_dod2k(
    dod,
    archive_types = spec$archive_types,
    proxy_types = spec$proxy_types,
    interpretation = spec$interpretation,
    domain = domain
  )
}

#' Assemble proxy block from a network specification
assemble_proxy_block <- function(spec, master_years, paths, domain = NULL) {
  switch(
    spec$builder,
    dod2k = {
      dod <- load_dod2k(spec$dod2k_variant %||% "v2", paths)
      subset_dod2k_spec(dod, spec, domain)
    },
    dod2k_combined = {
      dod <- load_dod2k(spec$dod2k_variant %||% "v2", paths)
      groups <- spec$proxy_groups
      caps <- spec$max_proxies_per_group
      blocks <- lapply(names(groups), function(g) {
        block <- subset_dod2k_spec(dod, groups[[g]], domain)
        cap <- if (!is.null(caps)) caps[[g]] else NULL
        if (is.null(cap)) {
          return(block)
        }
        aligned <- align_proxy_years(block$years, block$values, master_years)
        filtered <- filter_calibration_coverage(
          aligned,
          master_years,
          calib_start = spec$.calib_start,
          calib_end = spec$.calib_end,
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
      })
      Reduce(function(a, b) cbind_proxy_blocks(a, b, master_years), blocks)
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
    pages = {
      pages <- load_pages2k(paths)
      subset_pages_archives(pages, spec$archives)
    },
    combined_speleothem_ntrend = {
      pages <- load_pages2k(paths)
      nt <- load_ntrend(paths)
      sp <- subset_pages_archives(pages, "speleothem")
      cbind_proxy_blocks(sp, nt, master_years)
    },
    stop("Unknown builder: ", spec$builder)
  )
}

#' Null-coalescing infix
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Build a proxy network from a preset definition
build_proxy_network <- function(network_id,
                                recon_period = NULL,
                                calib_start,
                                calib_end,
                                domain = NULL,
                                paths = project_paths()) {
  spec <- PROXY_NETWORKS[[network_id]]
  if (is.null(spec)) {
    stop("Unknown network_id: ", network_id)
  }

  recon_period <- recon_period %||% spec$recon_period %||% RECON_PERIOD
  master_years <- seq.int(recon_period$start, recon_period$end)

  spec_run <- spec
  if (spec$builder == "dod2k_combined") {
    spec_run$.calib_start <- calib_start
    spec_run$.calib_end <- calib_end
  }

  block <- assemble_proxy_block(spec_run, master_years, paths, domain = domain)
  aligned <- align_proxy_years(block$years, block$values, master_years)

  max_px <- spec$max_proxies
  if (spec$builder == "dod2k_combined") {
    max_px <- NULL
  }

  filtered <- filter_calibration_coverage(
    aligned,
    master_years,
    calib_start,
    calib_end,
    min_coverage = spec$min_cal_coverage,
    max_proxies = max_px
  )

  if (ncol(filtered$values) < 2L) {
    stop(
      "Fewer than 2 proxies remain after filtering for '", network_id,
      "'. Adjust calibration window or network definition."
    )
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
      dod2k_variant = spec$dod2k_variant,
      n_proxies_raw = ncol(block$values),
      n_proxies_used = ncol(filtered$values),
      n_dropped = filtered$n_dropped,
      calib_window = c(calib_start, calib_end),
      recon_period = recon_period,
      proxy_names = block$names[filtered$keep],
      proxy_types = block$type[filtered$keep],
      archive_types = if (!is.null(block$archive_type)) block$archive_type[filtered$keep] else NULL
    )
  )
}

#' Summarise proxy locations for mapping
proxy_network_map_df <- function(network) {
  data.frame(
    lon = network$proxydata$lon,
    lat = network$proxydata$lat,
    name = network$meta$proxy_names,
    archive = network$meta$archive_types,
    proxy = network$meta$proxy_types,
    network = network$meta$network_id,
    stringsAsFactors = FALSE
  )
}
