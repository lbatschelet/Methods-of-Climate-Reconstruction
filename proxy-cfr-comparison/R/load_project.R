# Load project configuration and all R modules

load_project <- function(project_root = NULL) {
  if (is.null(project_root)) {
    candidates <- c(
      getwd(),
      normalizePath(file.path(getwd(), ".."), mustWork = FALSE),
      normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
    )
    project_root <- NULL
    for (path in unique(candidates[nzchar(candidates)])) {
      if (file.exists(file.path(path, "config", "experiments.R"))) {
        project_root <- normalizePath(path)
        break
      }
    }
    if (is.null(project_root)) {
      stop("Cannot find proxy-cfr-comparison project root.")
    }
  } else {
    project_root <- normalizePath(project_root, mustWork = TRUE)
  }

  assign("PROJECT_ROOT", project_root, envir = .GlobalEnv)

  r_dir <- file.path(project_root, "R")
  cfg <- file.path(project_root, "config", "experiments.R")

  source(file.path(r_dir, "paths.R"), local = .GlobalEnv)
  source(cfg, local = .GlobalEnv)
  source(file.path(r_dir, "io.R"), local = .GlobalEnv)
  source(file.path(r_dir, "proxy_networks.R"), local = .GlobalEnv)
  source(file.path(r_dir, "cfr_runner.R"), local = .GlobalEnv)
  source(file.path(r_dir, "validation.R"), local = .GlobalEnv)
  source(file.path(r_dir, "plotting.R"), local = .GlobalEnv)

  ensure_output_dirs()
  invisible(project_root)
}

if (interactive() && !exists("PROJECT_LOADED", envir = .GlobalEnv)) {
  load_project()
  assign("PROJECT_LOADED", TRUE, envir = .GlobalEnv)
}
