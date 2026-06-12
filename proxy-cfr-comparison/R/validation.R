# Validation metrics (tidyverse + Exercise 6 logic)

score_timeseries <- function(obs, pred, min_years = 30L) {
  tibble::tibble(obs = obs, pred = pred) |>
    tidyr::drop_na() |>
    (\(df) {
      n <- nrow(df)
      if (n < min_years) {
        return(tibble::tibble(n = n, r = NA_real_, rmse = NA_real_, rmsess = NA_real_))
      }
      rmse <- sqrt(mean((df$obs - df$pred)^2))
      rmse_clim <- sqrt(mean((df$obs - mean(df$obs))^2))
      tibble::tibble(
        n = n,
        r = cor(df$obs, df$pred),
        rmse = rmse,
        rmsess = if (rmse_clim > 0) 1 - rmse / rmse_clim else NA_real_
      )
    })()
}

lat_weighted_mean <- function(field_array, lat, time_idx) {
  w <- cos(lat * pi / 180)^0.5
  purrr::map_dbl(time_idx, \(t) {
    slice <- field_array[, , t]
    ok <- !is.na(slice)
    if (!any(ok)) return(NA_real_)
    ww <- outer(rep(1, nrow(slice)), w)[ok]
    sum(slice[ok] * ww) / sum(ww)
  })
}

validate_grid_period <- function(recon, target, years, min_years = 30L) {
  year_idx <- tibble::tibble(
    year = years,
    idx_r = match(years, recon$time),
    idx_t = match(years, target$time)
  ) |>
    tidyr::drop_na()

  idx_r <- year_idx$idx_r
  idx_t <- year_idx$idx_t
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

  list(
    years = year_idx$year,
    r = r,
    rmse = rmse,
    rmsess = rmsess,
    n_years = length(idx_t)
  )
}

validate_cfr_result <- function(result,
                              calibration_id = result$calibration_id,
                              min_years = 30L) {
  scheme <- CALIBRATION_SCHEMES[[calibration_id]]
  if (is.null(scheme)) {
    rlang::abort(paste("Unknown calibration_id:", calibration_id))
  }

  cal_years <- seq(scheme$calibration[1L], scheme$calibration[2L])
  val_years <- seq(scheme$validation[1L], scheme$validation[2L])

  grid_cal <- validate_grid_period(result$reconstruction, result$target, cal_years, min_years)
  grid_val <- validate_grid_period(result$reconstruction, result$target, val_years, min_years)

  idx_val <- match(val_years, result$reconstruction$time)
  obs_ix <- lat_weighted_mean(
    result$target$data,
    result$target$lat,
    match(val_years, result$target$time)
  )
  pred_ix <- as.numeric(result$spatmean[idx_val])

  list(
    experiment_id = result$experiment_id,
    network_id = result$network_id,
    calibration_id = calibration_id,
    scheme = scheme,
    grid_calibration = grid_cal,
    grid_validation = grid_val,
    index_validation = score_timeseries(obs_ix, pred_ix, min_years),
    min_years = min_years
  )
}

summarise_skill <- function(skill) {
  tibble::tibble(
    experiment = skill$experiment_id,
    network = skill$network_id,
    calibration = skill$calibration_id,
    val_median_r = median(skill$grid_validation$r, na.rm = TRUE),
    val_median_rmse = median(skill$grid_validation$rmse, na.rm = TRUE),
    val_median_rmsess = median(skill$grid_validation$rmsess, na.rm = TRUE),
    cal_median_r = median(skill$grid_calibration$r, na.rm = TRUE),
    index_val_r = dplyr::pull(skill$index_validation, r),
    index_val_rmsess = dplyr::pull(skill$index_validation, rmsess),
    n_cells_val = sum(!is.na(skill$grid_validation$r))
  ) |>
    dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 3)))
}

print_skill_summary <- function(skills) {
  if (!is.null(skills$experiment_id) && !is.null(skills$grid_validation)) {
    skills <- list(skills)
  }
  purrr::map_dfr(skills, summarise_skill)
}

skill_comparison_table <- function(results, calibration_id) {
  results |>
    purrr::keep(\(res) is.null(res$error)) |>
    purrr::map(\(res) validate_cfr_result(res, calibration_id)) |>
    print_skill_summary()
}
