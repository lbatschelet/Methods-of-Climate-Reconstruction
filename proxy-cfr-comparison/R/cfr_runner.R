# Run Master_CFR.R experiments with consistent bookkeeping

#' List all planned experiment combinations
experiment_grid <- function() {
  EXPERIMENT_PLAN[EXPERIMENT_PLAN$run, c("experiment_id", "network_id", "calibration_id")]
}

#' Resolve configuration objects for one experiment
resolve_experiment <- function(network_id, calibration_id) {
  network <- PROXY_NETWORKS[[network_id]]
  calib <- CALIBRATION_SCHEMES[[calibration_id]]
  if (is.null(network) || is.null(calib)) {
    stop("Unknown network_id or calibration_id.")
  }
  recon_period <- network$recon_period %||% RECON_PERIOD
  list(
    experiment_id = paste(network_id, calibration_id, sep = "__"),
    network = network,
    calibration = calib,
    domain = DOMAIN_NHET_JJA,
    recon_period = recon_period
  )
}

#' Load CFR output written by Master_CFR.R
load_cfr_output <- function(out_suffix, startyear, paths = project_paths()) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required.")
  }
  nc_path <- file.path(paths$exercises_dir, out_suffix, paste0("output_", out_suffix, ".nc"))
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  field <- ncdf4::ncvar_get(nc, "tas")
  if (length(dim(field)) == 4L) {
    field <- field[, , 1L, ]
  }

  spat_path <- file.path(
    paths$exercises_dir, out_suffix,
    paste0("Recon_spatmean_ROSM_", out_suffix, ".txt")
  )
  spat <- utils::read.table(spat_path, sep = ";")

  list(
    field = list(
      data = field,
      lon = ncdf4::ncvar_get(nc, "lon"),
      lat = ncdf4::ncvar_get(nc, "lat"),
      time = seq.int(startyear, length.out = dim(field)[3L])
    ),
    spatmean = stats::ts(spat[, 2L], start = startyear),
    out_suffix = out_suffix
  )
}

#' Run a single CFR experiment
#'
#' @param network_id Key in PROXY_NETWORKS
#' @param calibration_id Key in CALIBRATION_SCHEMES
#' @param overwrite If FALSE, skip when output already recorded
#' @return List with reconstruction, target, proxy network metadata, and paths
run_cfr_experiment <- function(network_id,
                             calibration_id,
                             overwrite = FALSE,
                             paths = project_paths()) {
  ensure_output_dirs(paths)
  cfg <- resolve_experiment(network_id, calibration_id)
  manifest_path <- file.path(
    paths$cfr_output_dir,
    paste0(cfg$experiment_id, ".rds")
  )

  if (!overwrite && file.exists(manifest_path)) {
    message("Loading cached experiment: ", cfg$experiment_id)
    return(readRDS(manifest_path))
  }

  cal <- cfg$calibration$calibration
  network <- build_proxy_network(
    network_id = network_id,
    recon_period = cfg$recon_period,
    calib_start = cal[1L],
    calib_end = cal[2L],
    paths = paths
  )
  target <- build_reconstruction_target(cfg$domain, paths)

  message("Running CFR: ", cfg$experiment_id)
  message("  Network: ", network$meta$label, " (n = ", network$meta$n_proxies_used, ")")
  message("  Calibration: ", cal[1L], "-", cal[2L])

  recon_name <- paste0("PFR_", cfg$experiment_id)
  exp_dir <- paths$exercises_dir

  # Master_CFR.R sources Recon_workflow_clean.r with local = FALSE, so parameters
  # must live in the global environment during the run.
  cfr_globals <- list(
    reconstruction.target = target,
    proxydata = network$proxydata,
    recon.name = recon_name,
    do.cps = FALSE,
    do.field = TRUE,
    do.index = TRUE,
    calib.start = cal[1L],
    calib.end = cal[2L],
    startyear = cfg$recon_period$start,
    endyear = cfg$recon_period$end
  )
  old_globals <- mget(names(cfr_globals), envir = .GlobalEnv, ifnotfound = list(NULL))
  list2env(cfr_globals, envir = .GlobalEnv)

  old_wd <- getwd()
  on.exit({
    setwd(old_wd)
    for (nm in names(cfr_globals)) {
      if (is.null(old_globals[[nm]])) {
        if (exists(nm, envir = .GlobalEnv, inherits = FALSE)) {
          rm(list = nm, envir = .GlobalEnv)
        }
      } else {
        assign(nm, old_globals[[nm]], envir = .GlobalEnv)
      }
    }
  }, add = TRUE)
  setwd(exp_dir)

  source(file.path(exp_dir, "Master_CFR.R"), local = FALSE)
  out.suffix <- get("out.suffix", envir = .GlobalEnv)
  startyear <- get("startyear", envir = .GlobalEnv)

  cfr <- load_cfr_output(out.suffix, startyear, paths)

  result <- list(
    experiment_id = cfg$experiment_id,
    network_id = network_id,
    calibration_id = calibration_id,
    network = network,
    target = target,
    reconstruction = cfr$field,
    spatmean = cfr$spatmean,
    out_suffix = cfr$out_suffix,
    config = cfg,
    run_time = Sys.time()
  )
  attr(result, "manifest_path") <- manifest_path
  saveRDS(result, manifest_path)
  result
}

#' Run all experiments flagged in EXPERIMENT_PLAN
run_experiment_batch <- function(overwrite = FALSE, paths = project_paths()) {
  plan <- EXPERIMENT_PLAN[EXPERIMENT_PLAN$run, ]
  results <- vector("list", nrow(plan))
  for (i in seq_len(nrow(plan))) {
    results[[i]] <- tryCatch(
      run_cfr_experiment(
        plan$network_id[i],
        plan$calibration_id[i],
        overwrite = overwrite,
        paths = paths
      ),
      error = function(e) {
        warning("Experiment failed: ", plan$experiment_id[i], " — ", conditionMessage(e))
        list(experiment_id = plan$experiment_id[i], error = conditionMessage(e))
      }
    )
  }
  names(results) <- plan$experiment_id
  results
}
