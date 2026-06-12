# Validation metrics — mirrors Exercise 6 logic, extended to grids and multiple splits

#' Score one paired time series (observed vs reconstructed)
score_timeseries <- function(obs, pred, min_years = 30L) {
  ok <- stats::complete.cases(obs, pred)
  n <- sum(ok)
  if (n < min_years) {
    return(list(n = n, r = NA_real_, rmse = NA_real_, rmsess = NA_real_))
  }
  obs <- obs[ok]
  pred <- pred[ok]
  rmse <- sqrt(mean((obs - pred)^2))
  clim <- mean(obs)
  rmse_clim <- sqrt(mean((obs - clim)^2))
  list(
    n = n,
    r = stats::cor(obs, pred),
    rmse = rmse,
    rmsess = if (rmse_clim > 0) 1 - rmse / rmse_clim else NA_real_
  )
}

#' Latitude-weighted spatial mean (consistent with Master_CFR index)
lat_weighted_mean <- function(field_array, lat, time_idx) {
  w <- cos(lat * pi / 180)^0.5
  weight_mat <- outer(rep(1, nrow(field_array)), w)
  vapply(time_idx, function(t) {
    slice <- field_array[, , t]
    ok <- !is.na(slice)
    if (!any(ok)) return(NA_real_)
    ww <- weight_mat[ok]
    sum(slice[ok] * ww) / sum(ww)
  }, numeric(1))
}

#' Grid-cell skill maps for one evaluation period
validate_grid_period <- function(recon, target, years, min_years = 30L) {
  idx_r <- match(years, recon$time)
  idx_t <- match(years, target$time)
  ok_y <- stats::complete.cases(idx_r, idx_t)
  idx_r <- idx_r[ok_y]
  idx_t <- idx_t[ok_y]

  nx <- length(recon$lon)
  ny <- length(recon$lat)
  r <- rmse <- rmsess <- array(NA_real_, dim = c(nx, ny))

  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      s <- score_timeseries(
        target$data[i, j, idx_t],
        recon$data[i, j, idx_r],
        min_years = min_years
      )
      r[i, j] <- s$r
      rmse[i, j] <- s$rmse
      rmsess[i, j] <- s$rmsess
    }
  }

  list(years = years[ok_y], r = r, rmse = rmse, rmsess = rmsess, n_years = length(idx_t))
}

#' Full validation for one CFR result and calibration scheme
validate_cfr_result <- function(result,
                              calibration_id = result$calibration_id,
                              min_years = 30L) {
  scheme <- CALIBRATION_SCHEMES[[calibration_id]]
  if (is.null(scheme)) {
    stop("Unknown calibration_id: ", calibration_id)
  }

  recon <- result$reconstruction
  target <- result$target
  cal_years <- seq(scheme$calibration[1L], scheme$calibration[2L])
  val_years <- seq(scheme$validation[1L], scheme$validation[2L])

  grid_cal <- validate_grid_period(recon, target, cal_years, min_years)
  grid_val <- validate_grid_period(recon, target, val_years, min_years)

  idx_val <- match(val_years, recon$time)
  obs_ix <- lat_weighted_mean(target$data, target$lat, match(val_years, target$time))
  pred_ix <- as.numeric(result$spatmean[idx_val])
  index_val <- score_timeseries(obs_ix, pred_ix, min_years = min_years)

  list(
    experiment_id = result$experiment_id,
    network_id = result$network_id,
    calibration_id = calibration_id,
    scheme = scheme,
    grid_calibration = grid_cal,
    grid_validation = grid_val,
    index_validation = index_val,
    min_years = min_years
  )
}

#' Tidy summary table across experiments
print_skill_summary <- function(skill) {
  if (is.data.frame(skill)) {
    return(skill)
  }
  # Single result vs list of results
  if (!is.null(skill$experiment_id) && !is.null(skill$grid_validation)) {
    skill <- list(skill)
  }

  rows <- lapply(skill, function(s) {
    data.frame(
      experiment = s$experiment_id,
      network = s$network_id,
      calibration = s$calibration_id,
      val_median_r = median(s$grid_validation$r, na.rm = TRUE),
      val_median_rmse = median(s$grid_validation$rmse, na.rm = TRUE),
      val_median_rmsess = median(s$grid_validation$rmsess, na.rm = TRUE),
      cal_median_r = median(s$grid_calibration$r, na.rm = TRUE),
      index_val_r = s$index_validation$r,
      index_val_rmsess = s$index_validation$rmsess,
      n_cells_val = sum(!is.na(s$grid_validation$r)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], function(x) round(x, 3))
  out
}

#' Compare skill tables for multiple experiments
skill_comparison_table <- function(results, calibration_id) {
  skills <- lapply(results, function(res) {
    if (!is.null(res$error)) return(NULL)
    validate_cfr_result(res, calibration_id = calibration_id)
  })
  skills <- Filter(Negate(is.null), skills)
  print_skill_summary(skills)
}
