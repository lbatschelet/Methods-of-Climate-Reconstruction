# Helper functions for Methods of Climate Reconstruction exercises

fill_na_linear <- function(x) {
  x <- as.numeric(x)
  if (!anyNA(x)) {
    return(x)
  }
  idx <- seq_along(x)
  ok <- !is.na(x)
  if (sum(ok) < 2) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  }
  x[!ok] <- stats::approx(idx[ok], x[ok], xout = idx[!ok], rule = 2)$y
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}

#' Butterworth low-pass filter for time series (single or multivariate)
tsfilt <- function(ts.data, period, filter.type = "bw") {
  if (!is.ts(ts.data)) {
    ts.data <- ts(ts.data)
  }
  if (filter.type != "bw") {
    stop("Only Butterworth low-pass filtering ('bw') is implemented.")
  }
  if (!requireNamespace("dplR", quietly = TRUE)) {
    stop("Package 'dplR' is required for tsfilt().")
  }
  filter_col <- function(x) {
    dplR::pass.filt(fill_na_linear(x), W = period, type = "low", method = "Butterworth")
  }
  if (NCOL(ts.data) == 1) {
    return(filter_col(ts.data))
  }
  out <- ts.data
  for (i in seq_len(NCOL(ts.data))) {
    out[, i] <- filter_col(ts.data[, i])
  }
  out
}

#' Find CRUTEM grid-box indices nearest to each station location
getgridboxnum <- function(station, grid) {
  xgridnum <- vapply(station$lon, function(lon) {
    which.min(abs(grid$lon - lon))
  }, integer(1))
  ygridnum <- vapply(station$lat, function(lat) {
    which.min(abs(grid$lat - lat))
  }, integer(1))
  list(xgridnum = xgridnum, ygridnum = ygridnum)
}

#' Score one station combination with a fixed train/test split.
#'
#' Procedure (same for every candidate set of stations):
#'   1. CALIBRATION: fit lm(Y ~ stations) on cal_idx only (1951-2015 in the course)
#'   2. VALIDATION:  apply those coefficients to val_idx only (1901-1950) — no re-fitting
#'   3. Compute r_cal, r_val, rmse_val, rmsess on the respective periods
#'
#' Returns NULL if too few complete years or the regression is not identifiable.
validate_station_regression <- function(Y, station_matrix, cal_idx, val_idx,
                                        station_cols = seq_len(ncol(station_matrix))) {
  station_cols <- unique(as.integer(station_cols))
  Ycal <- Y[cal_idx]
  Yval <- Y[val_idx]
  Xc <- as.data.frame(station_matrix[cal_idx, station_cols, drop = FALSE])
  Xv <- as.data.frame(station_matrix[val_idx, station_cols, drop = FALSE])
  ok_cal <- complete.cases(Ycal, Xc)
  n_obs <- sum(ok_cal)

  # Cap predictors so the regression is identifiable (more obs than coefficients)
  max_stations <- max(1L, n_obs - 2L)
  if (length(station_cols) > max_stations) {
    station_cols <- station_cols[seq_len(max_stations)]
    Xc <- as.data.frame(station_matrix[cal_idx, station_cols, drop = FALSE])
    Xv <- as.data.frame(station_matrix[val_idx, station_cols, drop = FALSE])
    ok_cal <- complete.cases(Ycal, Xc)
    n_obs <- sum(ok_cal)
  }

  if (n_obs < length(station_cols) + 3L) {
    return(NULL)
  }
  fit <- stats::lm(Ycal[ok_cal] ~ ., data = Xc[ok_cal, , drop = FALSE])
  pred_cal <- stats::fitted(fit)
  pred_val <- as.numeric(stats::predict(fit, newdata = Xv))
  if (length(pred_val) != length(Yval)) {
    return(NULL)
  }
  ok_v <- complete.cases(Yval, pred_val)
  if (sum(ok_v) < 10L) {
    return(NULL)
  }
  rmse_v <- sqrt(mean((Yval[ok_v] - pred_val[ok_v])^2))
  sd_v <- stats::sd(Yval[ok_v])
  list(
    idx = station_cols,
    n_used = length(station_cols),
    fit = fit,
    r_cal = stats::cor(Ycal[ok_cal], pred_cal),
    r_val = stats::cor(Yval[ok_v], pred_val[ok_v]),
    rmse_cal = sqrt(mean((Ycal[ok_cal] - pred_cal)^2)),
    rmse_val = rmse_v,
    rmsess = if (sd_v > 0) 1 - rmse_v / sd_v else NA_real_,
    n_cal = sum(ok_cal),
    n_val = sum(ok_v),
    pred_cal = pred_cal,
    pred_val = pred_val,
    ok_cal = ok_cal,
    ok_val = ok_v
  )
}

#' Exhaustive search: try every k-station combination and keep the one with highest `metric`.
#'
#' IMPORTANT: only ONE metric is optimized (no weighted combination).
#'   metric = "r_val"   -> maximize validation correlation (default in Ex. 7)
#'   metric = "rmsess"  -> maximize validation skill vs climatology (1 - RMSE/RMSE_clim)
#'
#' r_val and rmsess usually pick slightly different station triplets.
#' rmse_val is never used as the objective here (only reported).
find_best_stations <- function(Y, station_matrix, cal_idx, val_idx, k = 3L,
                               metric = c("r_val", "rmsess")) {
  metric <- match.arg(metric)
  n_st <- ncol(station_matrix)
  if (k < 1L || k > 3L) {
    stop("Exhaustive search is implemented for k = 1, 2, or 3 only.")
  }
  best <- list(score = -Inf)
  combos <- 0L

  eval_combo <- function(idx) {
    combos <<- combos + 1L
    res <- validate_station_regression(Y, station_matrix, cal_idx, val_idx, idx)
    if (is.null(res) || is.na(res[[metric]])) {
      return(invisible(NULL))
    }
    if (res[[metric]] > best$score) {
      best <<- c(list(score = res[[metric]]), res)
    }
    invisible(NULL)
  }

  if (k == 1L) {
    for (i in seq_len(n_st)) eval_combo(i)
  } else if (k == 2L) {
    for (i in seq_len(n_st - 1L)) {
      for (j in (i + 1L):n_st) eval_combo(c(i, j))
    }
  } else {
    for (i in seq_len(n_st - 2L)) {
      for (j in (i + 1L):(n_st - 1L)) {
        for (l in (j + 1L):n_st) eval_combo(c(i, j, l))
      }
    }
  }

  c(best, list(combos_tested = combos, metric = metric, k = k))
}

#' Generate an AR(1) pseudoproxy from an instrumental (or climate) signal
#'
#' P(t) = T(t) + N(t),  N(t) = rho * N(t-1) + eps(t)
#' SNR = sd(signal) / sd(noise)  =>  target var(noise) = (sd(signal)/SNR)^2
#' var(eps) = var(noise) * (1 - rho^2)
create.ar1.acf.pseudoproxy <- function(signal, rho = 0.6, snr = 0.4) {
  signal <- as.numeric(signal)
  l <- length(signal)
  sd_noise <- stats::sd(signal, na.rm = TRUE) / snr
  var_noise <- sd_noise^2
  var_eps <- var_noise * (1 - rho^2)
  eps <- stats::rnorm(l, sd = sqrt(var_eps))
  noise <- numeric(l)
  noise[1] <- eps[1]
  if (l > 1L) {
    for (t in 2:l) {
      noise[t] <- rho * noise[t - 1L] + eps[t]
    }
  }
  signal + noise
}

#' White-noise pseudoproxy (AR0): P = signal + N, N ~ N(0, (sd(signal)/SNR)^2)
create.white.pseudoproxy <- function(signal, snr = 0.4) {
  signal <- as.numeric(signal)
  sd_noise <- stats::sd(signal, na.rm = TRUE) / snr
  signal + stats::rnorm(length(signal), sd = sd_noise)
}

#' Correlation of each proxy with local CRUTEM grid box (1901-2012 overlap)
proxy_cors <- function(proxdata, gridnum, gridded) {
  n <- length(proxdata$lon)
  cor_v <- pval_v <- df_v <- rep(NA_real_, n)
  years_use <- intersect(proxdata$time, gridded$time)
  years_use <- years_use[years_use > 1900 & years_use < 2013]

  for (prox in seq_len(n)) {
    iy <- match(years_use, proxdata$time)
    ig <- match(years_use, gridded$time)
    px <- proxdata$data[iy, prox]
    gy <- gridded$data[gridnum$xgridnum[prox], gridnum$ygridnum[prox], ig]
    ok <- complete.cases(px, gy)
    if (sum(ok) > 5L) {
      ct <- tryCatch(
        stats::cor.test(px[ok], gy[ok]),
        error = function(e) NULL
      )
      if (!is.null(ct)) {
        cor_v[prox] <- unname(ct$estimate)
        df_v[prox] <- ct$parameter
        pval_v[prox] <- ct$p.value
      }
    }
  }
  list(r = cor_v, p = pval_v, df = df_v)
}

#' Build a year-aligned proxy matrix (years x proxies)
proxy_matrix <- function(proxdata, years) {
  iy <- match(years, proxdata$time)
  if (any(is.na(iy))) {
    stop("Some reconstruction years are missing from the proxy time axis.")
  }
  as.matrix(proxdata$data[iy, , drop = FALSE])
}

#' Simple PCR climate-field reconstruction (teaching version of Master_CFR)
#'
#' Steps (same logic as in the course PCR workflow, but without nests/ensemble):
#'   1. Standardise each proxy using the calibration period
#'   2. PCA on calibration proxies -> retain PCs explaining `var_explained`
#'   3. Regress each target grid cell on proxy PCs (calibration years)
#'   4. Project proxy PCs for all reconstruction years and predict the field
#'   5. Optionally rescale each grid cell to match calibration variance
run_pcr_cfr <- function(reconstruction.target,
                        proxydata,
                        calib.start = 1941L,
                        calib.end = 2000L,
                        startyear = 750L,
                        endyear = 2011L,
                        recon.name = "Test_reconstruction",
                        var_explained = 0.7,
                        min_cal_coverage = 0.8,
                        do.var.adj = TRUE,
                        write_nc = TRUE,
                        out_dir = getwd()) {
  years_recon <- seq.int(startyear, endyear)
  cal_years <- seq.int(calib.start, calib.end)
  cal_years <- cal_years[cal_years %in% years_recon]

  X <- proxy_matrix(proxydata, years_recon)
  n_prox <- ncol(X)
  cal_idx <- match(cal_years, years_recon)

  # Keep proxies with enough data in the calibration window
  cal_cov <- colMeans(!is.na(X[cal_idx, , drop = FALSE]))
  keep <- which(cal_cov >= min_cal_coverage)
  if (length(keep) < 2L) {
    stop("Fewer than 2 proxies pass the calibration coverage filter.")
  }
  X <- X[, keep, drop = FALSE]

  # Impute remaining gaps with calibration mean (simple infill for teaching)
  cal_means <- colMeans(X[cal_idx, , drop = FALSE], na.rm = TRUE)
  for (j in seq_len(ncol(X))) {
    na_j <- is.na(X[, j])
    X[na_j, j] <- cal_means[j]
  }

  mu <- colMeans(X[cal_idx, , drop = FALSE])
  sd_x <- apply(X[cal_idx, , drop = FALSE], 2L, stats::sd)
  sd_x[sd_x == 0 | !is.finite(sd_x)] <- 1
  Xs <- scale(X, center = mu, scale = sd_x)

  pca <- stats::prcomp(Xs[cal_idx, , drop = FALSE], center = FALSE, scale. = FALSE)
  cumvar <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  n_pc <- max(1L, which(cumvar >= var_explained)[1L])
  scores <- Xs %*% pca$rotation[, seq_len(n_pc), drop = FALSE]

  tgt <- reconstruction.target$data
  nx <- dim(tgt)[1L]
  ny <- dim(tgt)[2L]
  tgt_years <- reconstruction.target$time
  tgt_cal_idx <- match(cal_years, tgt_years)

  recon_arr <- array(NA_real_, dim = c(nx, ny, length(years_recon)))
  pc_names <- paste0("PC", seq_len(n_pc))
  df_all <- as.data.frame(scores)
  names(df_all) <- pc_names

  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      y <- tgt[i, j, tgt_cal_idx]
      df_cal <- df_all[cal_idx, , drop = FALSE]
      df_cal$y <- y
      ok <- complete.cases(df_cal)
      if (sum(ok) < n_pc + 3L) next
      fit <- stats::lm(y ~ ., data = df_cal[ok, , drop = FALSE])
      pred <- as.numeric(stats::predict(fit, newdata = df_all))
      if (do.var.adj) {
        sd_obs <- stats::sd(y[ok])
        sd_pred <- stats::sd(pred[cal_idx[ok]], na.rm = TRUE)
        if (is.finite(sd_obs) && is.finite(sd_pred) && sd_pred > 0) {
          pred <- pred * (sd_obs / sd_pred)
        }
      }
      recon_arr[i, j, ] <- pred
    }
  }

  recon <- list(
    data = recon_arr,
    lon = reconstruction.target$lon,
    lat = reconstruction.target$lat,
    time = years_recon
  )

  spatmean <- apply(recon_arr, 3L, mean, na.rm = TRUE)
  spatmean_recon <- stats::ts(spatmean, start = startyear)

  date_stamp <- format(Sys.time(), "%Y_%m_%d_%H%M%S")
  out.suffix <- paste0(recon.name, "_", date_stamp)
  out_folder <- file.path(out_dir, out.suffix)
  dir.create(out_folder, showWarnings = FALSE)

  spat_file <- file.path(out_folder, paste0("Recon_spatmean_ROSM_", out.suffix, ".txt"))
  spat_out <- cbind(years_recon, spatmean)
  write.table(spat_out, file = spat_file, sep = ";", row.names = FALSE, col.names = FALSE)

  if (write_nc && requireNamespace("ncdf4", quietly = TRUE)) {
    nc_file <- file.path(out_folder, paste0("output_", out.suffix, ".nc"))
    lon_dim <- ncdf4::ncdim_def("lon", "degrees_east", reconstruction.target$lon)
    lat_dim <- ncdf4::ncdim_def("lat", "degrees_north", reconstruction.target$lat)
    time_dim <- ncdf4::ncdim_def("time", "year", years_recon)
    var_def <- ncdf4::ncvar_def(
      "tas", "degC",
      dim = list(lon_dim, lat_dim, time_dim),
      missval = -9999
    )
    nc <- ncdf4::nc_create(nc_file, vars = var_def)
    ncdf4::ncvar_put(nc, var_def, recon_arr)
    ncdf4::nc_close(nc)
  }

  list(
    recon = recon,
    spatmean_recon = spatmean_recon,
    out.suffix = out.suffix,
    out_folder = out_folder,
    n_proxies_used = length(keep),
    n_pcs = n_pc
  )
}

#' Grid-wise validation metrics for a spatial reconstruction (Ex. 15)
validate_spatial_recon <- function(recon, target, cal_years = NULL, val_years = NULL) {
  nx <- length(recon$lon)
  ny <- length(recon$lat)

  period_scores <- function(years) {
    idx_recon <- match(years, recon$time)
    idx_tgt <- match(years, target$time)
    ok_y <- complete.cases(idx_recon, idx_tgt)
    idx_recon <- idx_recon[ok_y]
    idx_tgt <- idx_tgt[ok_y]
    r_out <- rmse_out <- rmsess_out <- array(NA_real_, dim = c(nx, ny))
    for (i in seq_len(nx)) {
      for (j in seq_len(ny)) {
        obs <- target$data[i, j, idx_tgt]
        pred <- recon$data[i, j, idx_recon]
        ok <- complete.cases(obs, pred)
        if (sum(ok) < 5L) next
        r_out[i, j] <- stats::cor(obs[ok], pred[ok])
        rmse_v <- sqrt(mean((obs[ok] - pred[ok])^2))
        rmse_out[i, j] <- rmse_v
        clim <- mean(obs[ok])
        rmse_clim <- sqrt(mean((obs[ok] - clim)^2))
        rmsess_out[i, j] <- if (rmse_clim > 0) 1 - rmse_v / rmse_clim else NA_real_
      }
    }
    list(r = r_out, rmse = rmse_out, rmsess = rmsess_out)
  }

  cal <- if (!is.null(cal_years)) period_scores(cal_years) else NULL
  val <- if (!is.null(val_years)) period_scores(val_years) else NULL

  list(
    r_cal = if (!is.null(cal)) cal$r else array(NA_real_, dim = c(nx, ny)),
    rmse_cal = if (!is.null(cal)) cal$rmse else array(NA_real_, dim = c(nx, ny)),
    rmsess_cal = if (!is.null(cal)) cal$rmsess else array(NA_real_, dim = c(nx, ny)),
    r_val = if (!is.null(val)) val$r else array(NA_real_, dim = c(nx, ny)),
    rmse_val = if (!is.null(val)) val$rmse else array(NA_real_, dim = c(nx, ny)),
    rmsess = if (!is.null(val)) val$rmsess else array(NA_real_, dim = c(nx, ny))
  )
}

#' Load world coastline segments for ggplot maps (NA rows = segment breaks)
load_world_coastline <- function(path = "../Data/world_coastline.dat") {
  raw <- utils::read.table(path, sep = "", header = FALSE)
  world <- data.frame(lon = raw[[1]], lat = raw[[2]])
  world$seg <- cumsum(is.na(world$lon) | is.na(world$lat))
  world[stats::complete.cases(world), , drop = FALSE]
}

#' Load field reconstruction output written by Master_CFR.R
load_cfr_output <- function(out.suffix, startyear = NULL) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required to load CFR NetCDF output.")
  }
  nc_path <- file.path(out.suffix, paste0("output_", out.suffix, ".nc"))
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  d <- ncdf4::ncvar_get(nc, "tas")
  if (length(dim(d)) == 4L) {
    d <- d[, , 1L, ]
  }
  lon <- ncdf4::ncvar_get(nc, "lon")
  lat <- ncdf4::ncvar_get(nc, "lat")

  spat_path <- file.path(out.suffix, paste0("Recon_spatmean_ROSM_", out.suffix, ".txt"))
  spat <- utils::read.table(spat_path, sep = ";")
  if (is.null(startyear)) {
    startyear <- as.integer(spat[1, 1])
  }
  recon <- list(
    data = d,
    lon = lon,
    lat = lat,
    time = seq.int(startyear, length.out = dim(d)[3])
  )
  spatmean_recon <- stats::ts(spat[, 2], start = startyear)
  list(recon = recon, spatmean_recon = spatmean_recon)
}
