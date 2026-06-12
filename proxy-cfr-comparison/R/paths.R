# Path resolution for the Quarto project and course materials

#' Locate the proxy-cfr-comparison project root (`config/experiments.R`)
find_project_root <- function() {
  locate_project_root()
}

#' Locate repository root (contains `course/data` and `proxy-cfr-comparison/`)
find_repo_root <- function() {
  env_root <- Sys.getenv("CLIMATE_RECON_REPO_ROOT", unset = "")
  if (nzchar(env_root)) {
    root <- normalizePath(env_root, mustWork = FALSE)
    if (dir.exists(file.path(root, "course", "data"))) {
      return(root)
    }
  }

  proj <- find_project_root()
  repo <- normalizePath(file.path(proj, ".."), mustWork = FALSE)
  if (dir.exists(file.path(repo, "course", "data"))) {
    return(repo)
  }

  stop(
    "Cannot find repository root with course/data/. ",
    "Expected layout: <repo>/course/data and <repo>/proxy-cfr-comparison/."
  )
}

#' Standard paths used across the project
project_paths <- function() {
  repo <- find_repo_root()
  project_root <- find_project_root()

  list(
    repo_root = repo,
    project_root = project_root,
    course_root = file.path(repo, "course"),
    data_dir = file.path(repo, "course", "data"),
    exercises_dir = file.path(repo, "course", "exercises"),
    output_dir = file.path(project_root, "output"),
    cache_dir = file.path(project_root, "output", "cache"),
    cfr_output_dir = file.path(project_root, "output", "cfr"),
    tables_dir = file.path(project_root, "output", "tables"),
    figures_dir = file.path(project_root, "output", "figures")
  )
}

#' Create output directories if missing
ensure_output_dirs <- function(paths = project_paths()) {
  invisible(vapply(
    paths[c("output_dir", "cache_dir", "cfr_output_dir", "tables_dir", "figures_dir")],
    function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE),
    logical(1)
  ))
}
