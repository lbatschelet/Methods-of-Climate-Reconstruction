# Experiment definitions for proxy CFR comparison
# Sourced by R/load_project.R — do not edit paths here; see R/paths.R

#' Target domain: NH extratropics JJA anomalies (same as course Exercise 14)
DOMAIN_NHET_JJA <- list(
  id = "nhet_jja",
  label = "NH extratropics (40–90°N), JJA anomalies",
  lon_w = -180,
  lon_e = 180,
  lat_s = 40,
  lat_n = 90,
  season = "jja",
  variable = "anomaly"
)

#' Reconstruction period — N-TREND / PAGES
RECON_PERIOD <- list(start = 750L, end = 2011L)

#' Reconstruction period — DoD2k (calendar years 0–2000 on file axis)
RECON_PERIOD_DOD2K <- list(start = 750L, end = 2000L)

#' DoD2k loaders — see https://lluecke.github.io/dod2k/
#' `v2` = official compact CSV (recommended); course RData subsets are legacy.
DOD2K_VARIANTS <- c(
  v2 = "v2",
  mean = "dod2k.m.data",
  tm = "dod2k.m.tm.data",
  full = "dod2kdata"
)

#' LiPDverse archive types (DoD2k `archiveType` field — NOT `paleoData_proxy`)
DOD2K_ARCHIVE_SPELEOTHEM <- "Speleothem"
DOD2K_ARCHIVE_WOOD <- "Wood"

#' Proxy measurement types (`paleoData_proxy` in DoD2k metadata)
DOD2K_PROXY_D18O <- "d18O"
DOD2K_PROXY_D13C <- "d13C"
DOD2K_PROXY_GROWTH_RATE <- "growth rate"
DOD2K_PROXY_MGCA <- "Mg/Ca"
DOD2K_PROXY_RING_WIDTH <- "ring width"

#' All speleothem proxy types in DoD2k v2.0 (archiveType Speleothem)
DOD2K_SPELEOTHEM_PROXIES <- c(
  DOD2K_PROXY_D18O,
  DOD2K_PROXY_D13C,
  DOD2K_PROXY_GROWTH_RATE,
  DOD2K_PROXY_MGCA
)

#' Proxy network presets
PROXY_NETWORKS <- list(
  # --- DoD2k v2.0 (archiveType + paleoData_proxy per ESSD / dod2k docs) -------
  dod2k_speleothems = list(
    id = "dod2k_speleothems",
    label = "DoD2k speleothem calcite δ18O (archiveType Speleothem)",
    builder = "dod2k",
    dod2k_variant = "v2",
    archive_types = DOD2K_ARCHIVE_SPELEOTHEM,
    proxy_types = DOD2K_PROXY_D18O,
    min_cal_coverage = 0.5,
    recon_period = RECON_PERIOD_DOD2K,
    max_proxies = NULL
  ),
  dod2k_trees = list(
    id = "dod2k_trees",
    label = "DoD2k tree-ring width (Wood, top 60 by cal. coverage)",
    builder = "dod2k",
    dod2k_variant = "v2",
    archive_types = DOD2K_ARCHIVE_WOOD,
    proxy_types = DOD2K_PROXY_RING_WIDTH,
    min_cal_coverage = 0.8,
    recon_period = RECON_PERIOD_DOD2K,
    max_proxies = 60L
  ),
  dod2k_speleo_trees = list(
    id = "dod2k_speleo_trees",
    label = "DoD2k speleothem δ18O + tree rings (combined)",
    builder = "dod2k_combined",
    dod2k_variant = "v2",
    proxy_groups = list(
      speleothem = list(
        archive_types = DOD2K_ARCHIVE_SPELEOTHEM,
        proxy_types = DOD2K_PROXY_D18O
      ),
      trees = list(
        archive_types = DOD2K_ARCHIVE_WOOD,
        proxy_types = DOD2K_PROXY_RING_WIDTH
      )
    ),
    min_cal_coverage = 0.5,
    recon_period = RECON_PERIOD_DOD2K,
    max_proxies_per_group = list(speleothem = NULL, trees = 60L)
  ),

  # --- Legacy course RData subsets (pre-filtered; limited speleothem coverage) --
  dod2k_course_mt = list(
    id = "dod2k_course_mt",
    label = "Course DoD2k MT subset (dod2k.m.data — no speleothems)",
    builder = "dod2k",
    dod2k_variant = "mean",
    archive_types = DOD2K_ARCHIVE_WOOD,
    proxy_types = DOD2K_PROXY_RING_WIDTH,
    min_cal_coverage = 0.8,
    recon_period = RECON_PERIOD_DOD2K,
    max_proxies = 60L
  ),

  # --- Other course datasets ---------------------------------------------------
  ntrend_trees = list(
    id = "ntrend_trees",
    label = "N-TREND tree rings (54 series)",
    builder = "ntrend",
    min_cal_coverage = 0.8,
    recon_period = RECON_PERIOD,
    max_proxies = NULL
  ),
  pages_speleothems = list(
    id = "pages_speleothems",
    label = "PAGES 2.0 speleothems (4 series only)",
    builder = "pages",
    archives = "speleothem",
    min_cal_coverage = 0.5,
    recon_period = RECON_PERIOD,
    max_proxies = NULL
  ),
  speleothem_ntrend = list(
    id = "speleothem_ntrend",
    label = "PAGES speleothems + N-TREND (calibration-filtered)",
    builder = "combined_speleothem_ntrend",
    min_cal_coverage = 0.5,
    recon_period = RECON_PERIOD,
    max_proxies = NULL
  )
)

#' Calibration / validation schemes
CALIBRATION_SCHEMES <- list(
  course_split = list(
    id = "course_split",
    label = "Course default: cal 1941–2000, val 1901–1940",
    calibration = c(1941L, 2000L),
    validation = c(1901L, 1940L)
  ),
  early_holdout = list(
    id = "early_holdout",
    label = "Early hold-out: cal 1951–2000, val 1901–1950",
    calibration = c(1951L, 2000L),
    validation = c(1901L, 1950L)
  ),
  recent_holdout = list(
    id = "recent_holdout",
    label = "Recent hold-out: cal 1941–1990, val 1991–2010",
    calibration = c(1941L, 1990L),
    validation = c(1991L, 2010L)
  )
)

#' Default experiment matrix
EXPERIMENT_PLAN <- expand.grid(
  network_id = c(
    "dod2k_speleothems",
    "dod2k_trees",
    "dod2k_speleo_trees",
    "ntrend_trees"
  ),
  calibration_id = names(CALIBRATION_SCHEMES),
  stringsAsFactors = FALSE
)
EXPERIMENT_PLAN$experiment_id <- paste(
  EXPERIMENT_PLAN$network_id,
  EXPERIMENT_PLAN$calibration_id,
  sep = "__"
)
EXPERIMENT_PLAN$run <- TRUE
