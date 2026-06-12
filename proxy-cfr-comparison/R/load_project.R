# Load project modules and dependencies

locate_project_root <- function() {
  starts <- unique(na.omit(c(
    Sys.getenv("PROXY_CFR_PROJECT_ROOT", unset = NA_character_),
    getwd(),
    {
      if (!requireNamespace("knitr", quietly = TRUE)) NULL
      else tryCatch({
        inp <- knitr::current_input()
        if (is.null(inp) || !nzchar(inp)) NULL else dirname(normalizePath(inp, mustWork = FALSE))
      }, error = function(e) NULL)
    }
  )))

  for (start in starts) {
    path <- start
    for (i in seq_len(10L)) {
      for (root in unique(c(path, file.path(path, "proxy-cfr-comparison")))) {
        if (file.exists(file.path(root, "config", "experiments.R"))) {
          return(normalizePath(root))
        }
      }
      parent <- dirname(path)
      if (identical(parent, path)) break
      path <- parent
    }
  }

  rlang::abort(
    "Cannot find proxy-cfr-comparison/. Setwd to the project or open analysis.qmd."
  )
}

load_project <- function(project_root = NULL) {
  project_root <- project_root %||% locate_project_root()
  project_root <- normalizePath(project_root, mustWork = TRUE)

  if (normalizePath(getwd(), mustWork = FALSE) != project_root) {
    setwd(project_root)
  }

  r_dir <- file.path(project_root, "R")

  source(file.path(r_dir, "paths.R"), local = .GlobalEnv)
  source(file.path(project_root, "config", "experiments.R"), local = .GlobalEnv)
  source(file.path(r_dir, "io.R"), local = .GlobalEnv)
  source(file.path(r_dir, "dod2k_cache.R"), local = .GlobalEnv)
  source(file.path(r_dir, "proxy_networks.R"), local = .GlobalEnv)
  source(file.path(r_dir, "cfr_runner.R"), local = .GlobalEnv)
  source(file.path(r_dir, "validation.R"), local = .GlobalEnv)
  source(file.path(r_dir, "plotting.R"), local = .GlobalEnv)

  ensure_output_dirs()
  invisible(project_root)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
