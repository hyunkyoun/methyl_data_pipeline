# Keep each execution in its own directory, including two runs started in one second.
create_output_dirs <- function(results_root = "results",
                               run_id = paste0("analysis_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))) {
  if (length(run_id) != 1 || is.na(run_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id)) {
    stop("run_id must contain only letters, numbers, underscores, and hyphens.")
  }
  if (!dir.exists(results_root) && !dir.create(results_root, recursive = TRUE)) {
    stop("Cannot create results directory: ", results_root)
  }
  run_dir <- file.path(results_root, run_id)
  suffix <- 0L
  while (!dir.create(run_dir, showWarnings = FALSE)) {
    if (!file.exists(run_dir)) stop("Cannot create analysis directory: ", run_dir)
    suffix <- suffix + 1L
    run_dir <- file.path(results_root, paste0(run_id, sprintf("_%02d", suffix)))
  }
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  dirs <- list(
    root = run_dir,
    betas = file.path(run_dir, "tables", "betas"),
    summaries = file.path(run_dir, "tables", "summaries"),
    metadata = file.path(run_dir, "metadata")
  )
  for (path in unlist(dirs, use.names = FALSE)) {
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
      stop("Cannot create output directory: ", path)
    }
  }
  dirs
}

experiment_output_dirs <- function(output_dirs, experiment, same_parent = FALSE) {
  experiment <- as.character(experiment)
  if (length(experiment) != 1L || is.na(experiment) || !grepl("^[A-Za-z0-9_-]+$", experiment)) {
    stop("Experiment must be a non-empty, filesystem-safe identifier.")
  }
  root <- file.path(output_dirs$root, paste0("Experiment_", experiment, if (same_parent) "_sameparent" else ""))
  dirs <- list(root = root, metadata = file.path(root, "metadata"),
               dmp_tables = file.path(root, "tables", "dmp"),
               pca_tables = file.path(root, "tables", "pca"),
               fasta_details = file.path(root, "tables", "fasta_details"),
               pca_figures = file.path(root, "figures", "pca"),
               volcano_figures = file.path(root, "figures", "volcano"),
               heatmap_figures = file.path(root, "figures", "heatmap"),
               fasta = file.path(root, "fasta"))
  for (path in dirs) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dirs
}

# Relative paths keep the index useful when a results folder is moved or shared.
write_output_manifest <- function(output_dirs) {
  files <- sort(list.files(output_dirs$root, recursive = TRUE))
  files <- setdiff(files, "output_manifest.csv")
  manifest <- data.frame(
    Path = files,
    Bytes = unname(file.info(file.path(output_dirs$root, files))$size)
  )
  write.csv(manifest, file.path(output_dirs$root, "output_manifest.csv"), row.names = FALSE)
  invisible(manifest)
}
