# ggplot2 visualisation helpers

theme_pfr <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(colour = "grey35", hjust = 0, size = rel(0.92)),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    )
}

load_coastline <- function(paths = project_paths()) {
  raw <- utils::read.table(file.path(paths$data_dir, "world_coastline.dat"), sep = "")
  world <- data.frame(lon = raw[[1L]], lat = raw[[2L]])
  world$seg <- cumsum(is.na(world$lon) | is.na(world$lat))
  world[stats::complete.cases(world), , drop = FALSE]
}

grid_to_df <- function(grid_array, lon, lat) {
  df <- expand.grid(lon = lon, lat = lat)
  df$value <- as.vector(grid_array)
  df
}

#' Map a skill grid (r, RMSE, RMSESS)
plot_skill_map <- function(grid_array, lon, lat, title,
                           limits = c(-1, 1), midpoint = 0,
                           domain = DOMAIN_NHET_JJA,
                           coastline = NULL) {
  if (is.null(coastline)) {
    coastline <- load_coastline()
  }
  df <- grid_to_df(grid_array, lon, lat)
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = df, ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data$value)) +
    ggplot2::geom_path(
      data = coastline,
      ggplot2::aes(x = .data$lon, y = .data$lat, group = .data$seg),
      inherit.aes = FALSE,
      colour = "grey35",
      linewidth = 0.15
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = midpoint, limits = limits,
      oob = scales::squish, na.value = "grey92",
      name = NULL
    ) +
    ggplot2::coord_fixed(
      ratio = 1.35,
      xlim = c(domain$lon_w, domain$lon_e),
      ylim = c(domain$lat_s, domain$lat_n)
    ) +
    ggplot2::labs(title = title, x = "Longitude", y = "Latitude") +
    theme_pfr()
}

#' Proxy network locations
plot_proxy_map <- function(network_df, title = "Proxy network") {
  ggplot2::ggplot(network_df, ggplot2::aes(x = .data$lon, y = .data$lat, colour = .data$network)) +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::borders("world", fill = NA, colour = "grey70") +
    ggplot2::coord_fixed(ratio = 1.3) +
    ggplot2::labs(title = title, x = "Longitude", y = "Latitude", colour = "Network") +
    theme_pfr()
}

#' Bar chart comparing index validation skill
plot_skill_comparison <- function(summary_df, metric = "index_val_r") {
  summary_df$experiment <- reorder(summary_df$experiment, summary_df[[metric]])
  ggplot2::ggplot(summary_df, ggplot2::aes(x = .data$experiment, y = .data[[metric]], fill = .data$network)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Comparison across experiments —", metric),
      x = NULL, y = metric
    ) +
    theme_pfr()
}
