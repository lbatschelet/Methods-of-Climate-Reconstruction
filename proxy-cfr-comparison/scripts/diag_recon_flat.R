library(tidyverse)
cache <- "output/cache"
catalog <- readRDS(file.path(cache, "dod2k_timeseries.rds"))
cru_clim <- readRDS(file.path(cache, "cru_site_climate.rds"))
speleo_corr <- readRDS(file.path(cache, "speleo_instrument_correlations.rds"))
cal_y0 <- 1901L; cal_y1 <- 1980L; val_y0 <- 1981L; val_y1 <- 2000L
min_overlap <- 10L; screen_p <- 0.05; min_screen_n <- 3L
climate_col <- c(tmp_ann = "tmp_ann", tmp_jja = "tmp_jja", pre_ann = "pre_ann")
r_col <- c(tmp_ann = "r_tmp_ann", tmp_jja = "r_tmp_jja", pre_ann = "r_pre_ann")

detrend_linear <- function(years, values) {
  ok <- is.finite(values)
  if (sum(ok) < 5L) return(rep(NA_real_, length(values)))
  fit <- lm(values[ok] ~ years[ok])
  values - predict(fit, newdata = data.frame(years = years)) + mean(values[ok])
}
zscore_ref <- function(x, years, ref_y0, ref_y1) {
  ref <- x[years >= ref_y0 & years <= ref_y1]
  ref <- ref[is.finite(ref)]
  if (length(ref) < 3L || sd(ref) == 0) return(rep(NA_real_, length(x)))
  (x - mean(ref)) / sd(ref)
}
assign_region <- function(lon, lat) {
  case_when(
    lat < -10 & lon > 100 ~ "Australasia",
    lat <= 15 & lon >= -85 & lon < -30 ~ "South America",
    lat > 15 & lon < -30 ~ "Americas (north)",
    lat >= 35 & lon >= -15 & lon < 45 ~ "Europe",
    lat >= -35 & lat < 40 & lon >= -20 & lon < 55 ~ "Africa / Med",
    lon >= 50 ~ "Asia",
    TRUE ~ "Other"
  )
}
speleothem_annual <- function(catalog) {
  catalog |>
    filter(archive_type == "Speleothem") |>
    group_by(dataset_id, proxy_type, lon, lat, year) |>
    summarize(value = mean(value), .groups = "drop")
}
speleo_sites <- catalog |>
  filter(archive_type == "Speleothem") |>
  group_by(dataset_id, proxy_type, lon, lat) |>
  summarize(.groups = "drop")

speleo_corr <- speleo_corr |>
  mutate(
    region = assign_region(lon, lat),
    p_best = case_when(
      best_var == "tmp_ann" ~ p_tmp_ann,
      best_var == "tmp_jja" ~ p_tmp_jja,
      best_var == "pre_ann" ~ p_pre_ann,
      TRUE ~ NA_real_
    ),
    screened = is.finite(best_r) & !is.na(p_best) & p_best < screen_p & n_overlap >= min_overlap
  )

loc <- speleo_sites |> distinct(lon, lat) |> mutate(loc_id = row_number())
proxy_ann <- speleothem_annual(catalog) |> inner_join(loc, by = c("lon", "lat"))

calibrate_one_site <- function(proxy_df, clim_df, climate_var) {
  col <- climate_col[[climate_var]]
  d <- proxy_df |>
    inner_join(clim_df |> transmute(year, climate = .data[[col]]), by = "year") |>
    filter(year >= cal_y0, year <= val_y1)
  if (nrow(d) < min_overlap) return(NULL)
  d <- d |>
    mutate(
      proxy_dt = detrend_linear(year, value),
      proxy_z = zscore_ref(proxy_dt, year, cal_y0, cal_y1),
      climate_z = zscore_ref(climate, year, cal_y0, cal_y1)
    )
  cal <- d |> filter(year >= cal_y0, year <= cal_y1)
  ok <- complete.cases(cal[, c("proxy_z", "climate_z")])
  if (sum(ok) < min_overlap) return(NULL)
  fit <- lm(climate_z ~ proxy_z, data = cal[ok, ])
  list(data = d |> mutate(recon_z = predict(fit, newdata = d)), coef = coef(fit))
}

screened <- speleo_corr |> filter(screened)
combos <- screened |> count(region, best_var, name = "n") |> filter(n >= min_screen_n)

cat("=== Composite validation diagnostics ===\n")
for (i in seq_len(nrow(combos))) {
  reg <- combos$region[i]
  cv <- combos$best_var[i]
  sub <- screened |> filter(region == reg, best_var == cv)

  site_recons <- map(sub$dataset_id, function(id) {
    row <- sub |> filter(dataset_id == id)
    px <- proxy_ann |> filter(dataset_id == id)
    cdf <- cru_clim |> filter(loc_id == px$loc_id[1])
    out <- calibrate_one_site(px, cdf, cv)
    if (is.null(out)) return(NULL)
    out$data |>
      mutate(
        dataset_id = id,
        sign_r = sign(row[[r_col[[cv]]]][1])
      )
  }) |> compact()

  all_years <- bind_rows(site_recons) |>
    mutate(recon_z_signed = if_else(sign_r < 0, -recon_z, recon_z))

  col <- climate_col[[cv]]
  regional_clim <- map_dfr(sub$dataset_id, function(id) {
    px <- proxy_ann |> filter(dataset_id == id) |> slice(1)
    cru_clim |>
      filter(loc_id == px$loc_id) |>
      transmute(year, climate = .data[[col]])
  }) |>
    group_by(year) |>
    summarize(climate = mean(climate, na.rm = TRUE), .groups = "drop") |>
    mutate(climate_z = zscore_ref(climate, year, cal_y0, cal_y1))

  composite <- all_years |>
    group_by(year) |>
    summarize(recon_z = mean(recon_z_signed, na.rm = TRUE), .groups = "drop") |>
    left_join(regional_clim |> select(year, climate_z), by = "year")

  val <- composite |> filter(year >= val_y0, year <= val_y1)
  cat(sprintf(
    "%s / %s | sd(CRU)=%.2f sd(recon)=%.4f r_val=%.3f | per-site sd(recon): %s\n",
    reg, cv,
    sd(val$climate_z, na.rm = TRUE),
    sd(val$recon_z, na.rm = TRUE),
    cor(val$climate_z, val$recon_z, use = "complete"),
    paste(round(map_dbl(site_recons, \(x) sd(x$recon_z[x$year >= val_y0], na.rm = TRUE)), 3), collapse = ", ")
  ))
}

cat("\n=== Test: composite without sign flip (proxy_z mean -> calibrate once) ===\n")
reg <- "Americas (north)"; cv <- "tmp_ann"
sub <- screened |> filter(region == reg, best_var == cv)
col <- climate_col[[cv]]
regional_clim <- map_dfr(sub$dataset_id, function(id) {
  px <- proxy_ann |> filter(dataset_id == id) |> slice(1)
  cru_clim |> filter(loc_id == px$loc_id) |> transmute(year, climate = .data[[col]])
}) |>
  group_by(year) |> summarize(climate = mean(climate), .groups = "drop")

proxy_comp <- map_dfr(sub$dataset_id, function(id) {
  px <- proxy_ann |> filter(dataset_id == id)
  px |>
    filter(year >= cal_y0, year <= val_y1) |>
    mutate(
      proxy_dt = detrend_linear(year, value),
      proxy_z = zscore_ref(proxy_dt, year, cal_y0, cal_y1)
    ) |>
    select(year, proxy_z)
}) |>
  group_by(year) |>
  summarize(proxy_z = mean(proxy_z, na.rm = TRUE), .groups = "drop") |>
  inner_join(regional_clim, by = "year") |>
  mutate(
    climate_z = zscore_ref(climate, year, cal_y0, cal_y1)
  )

cal <- proxy_comp |> filter(year >= cal_y0, year <= cal_y1)
fit <- lm(climate_z ~ proxy_z, data = cal)
proxy_comp <- proxy_comp |> mutate(recon_z = predict(fit, newdata = proxy_comp))
val <- proxy_comp |> filter(year >= val_y0)
cat(sprintf("Pooled calibration: sd(recon)=%.3f r_val=%.3f coef=%s\n",
  sd(val$recon_z), cor(val$climate_z, val$recon_z), paste(round(coef(fit),3), collapse=", ")))
