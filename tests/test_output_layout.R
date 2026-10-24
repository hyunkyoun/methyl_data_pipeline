# Run from the repository root: Rscript --vanilla tests/test_output_layout.R
test_library <- Sys.getenv("METHYL_DMP_TEST_LIBRARY")
if (nzchar(test_library)) .libPaths(c(test_library, .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
})
source("R/output_paths.R")
source("R/heatmap_plots.R")

# Exercise the workflow's actual saving/plotting blocks with a tiny indexed genome.
workflow_lines <- readLines("workflow_new.R", warn = FALSE)
workflow_block <- function(start_marker, end_marker = NULL) {
  first <- match(start_marker, workflow_lines)
  last <- if (is.null(end_marker)) length(workflow_lines) + 1L else match(end_marker, workflow_lines)
  stopifnot(!is.na(first), !is.na(last), first < last)
  parse(text = workflow_lines[seq.int(first, last - 1L)])
}

run_tests <- function() {
  original_dir <- getwd()
  scratch <- tempfile("methyl-layout-test-")
  stopifnot(dir.create(scratch), dir.create(file.path(scratch, "R")),
            dir.create(file.path(scratch, "idats")))
  source_files <- c("workflow_new.R", "R/dmp_analysis.R", "R/output_paths.R",
                    "R/pca_analysis.R", "R/heatmap_plots.R")
  stopifnot(all(file.copy(source_files, file.path(scratch, source_files))))
  on.exit(setwd(original_dir))
  setwd(scratch)

  cohort <- function(experiment, tissue, genotypes) {
    data.frame(Experiment = experiment, Tissue.Cell.type = tissue, True.Genotype = genotypes)
  }
  samples <- rbind(
    cohort(2, "Tcells", c(rep("WT", 9), rep(c("DP10", "DP16", "DP10_DP16"), each = 3))),
    cohort(7, "Tcells", rep(c("WT", "DP10", "DP16", "DP10_DP16"), each = 3)),
    cohort(3, "thymus", rep(c("WT", "DP16", "RUNX1-DEL", "RUNX1-DEL_DP16"), each = 3)),
    cohort(4, "thymocytes", rep(c("WT", "DP16", "RUNX1-DEL", "RUNX1-DEL_DP16"), each = 3)),
    cohort(4, "thymus", "RUNX1-DEL_DP16"),
    cohort(5, "thymocytes", rep(c("WT", "Ts65Dn"), each = 3))
  )
  samples$age <- ifelse(samples$Experiment == 5, 2, 17.5)
  # An age-2 animal within Experiment 4 ensures age selection is not just an Exp4 shortcut.
  samples$age[which(samples$Experiment == 4 & samples$Tissue.Cell.type == "thymocytes")[1]] <- 2
  samples$Samples <- paste0("mouse", seq_len(nrow(samples)))
  # Use three of the actual excluded WT identifiers; all other WT animals must remain.
  samples$Samples[samples$Experiment == 5 & samples$True.Genotype == "WT"] <-
    c("11984_5", "11985_5", "11988_5")
  samples$group <- "repeated_legacy_group"
  samples$cross.type <- "Double"
  samples$cross.type[samples$Experiment == 2 & samples$True.Genotype == "WT"] <-
    rep(c("DP10", "DP16", "Double"), each = 3)
  samples$Genotype <- samples$True.Genotype
  samples$Genotype[samples$Experiment == 2 & samples$True.Genotype == "WT"] <-
    rep(c("WT_DP10", "WT_DP16", "WT_DP10_WT_DP16"), each = 3)
  samples$BATCH <- rep(c("BATCH1", "BATCH2"), length.out = nrow(samples))
  write.csv(samples, "idats/samplesheet_test.csv", row.names = FALSE)
  set.seed(20261014)
  betas <- matrix(runif(100 * nrow(samples), 0.45, 0.55), nrow = 100,
                  dimnames = list(paste0("cg", 1:100), samples$Samples))
  effect <- ifelse(samples$True.Genotype == "WT", 0, 0.3)
  betas[1, ] <- 0.5 + effect + rnorm(nrow(samples), sd = 0.002)
  betas[2, ] <- 0.5 - effect + rnorm(nrow(samples), sd = 0.002)
  # Opposite effects from the omitted RUNX1-DEL group must not label retained CpGs Mixed.
  mixed_effect <- ifelse(samples$True.Genotype == "RUNX1-DEL", -effect, effect)
  betas[3, ] <- 0.5 + mixed_effect + rnorm(nrow(samples), sd = 0.002)

  env <- new.env(parent = globalenv())
  env$idir <- "idats"
  env$samplesheet_path <- "idats/samplesheet_test.csv"
  invisible(capture.output(eval(workflow_block("# Saving files", 'cat("Reading IDAT files")'), env)))
  dirs <- env$output_dirs
  original_exp2 <- experiment_output_dirs(dirs, 2)
  writeLines("Preserve the earlier pooled-WT analysis.", file.path(original_exp2$root, "previous_run.txt"))
  env$betas_adjusted <- betas
  env$matching_samples <- samples
  env$skip_combat <- FALSE
  write.csv(betas, file.path(dirs$betas, "betas_combat_adjusted.csv"))
  write.csv(samples, file.path(dirs$metadata, "samples_after_qc.csv"), row.names = FALSE)

  env$heatmap_calls <- list()
  env$plot_dmp_heatmap <- function(probes, selected_samples, plot_title, output_file,
                                  directions, betas_data, targets_data, ...) {
    env$heatmap_calls[[output_file]] <- list(probes = probes, samples = selected_samples,
                                            directions = directions,
                                            genotypes = targets_data$True.Genotype[match(selected_samples, targets_data$Samples)])
    plot_dmp_heatmap(probes, selected_samples, plot_title, output_file,
                     directions, betas_data, targets_data, ...)
  }
  writeXStringSet(DNAStringSet(c(chr1 = paste(rep("ACGT", 1000), collapse = ""))),
                 "idats/reference.fasta")
  indexFa("idats/reference.fasta")
  env$BSgenome.Mmusculus.UCSC.mm10 <- FaFile("idats/reference.fasta")
  env$annoMouse <- data.frame(Name = paste0("cg", 1:3, "_BC21"), chr = "chr1",
                              pos = c(300L, 1000L, 2000L))
  messages <- capture.output(suppressWarnings(eval(workflow_block("# DMP analysis"), env)))
  stopifnot(length(env$results_list) == 13L, length(env$heatmap_calls) == 17L)
  stopifnot(basename(env$experiment_dirs[["2"]]$root) == "Experiment_2_sameparent",
            readLines(file.path(original_exp2$root, "previous_run.txt")) == "Preserve the earlier pooled-WT analysis.")

  pooled_root <- file.path(dirs$root, "All_Experiments")
  pooled_scores <- read.csv(file.path(pooled_root, "tables", "pca", "scores",
                                      "Genotypes_AllExperiments_all_tissues.csv"))
  pooled_expected <- samples[samples$Experiment != 5, ]
  stopifnot(setequal(pooled_scores$Samples, pooled_expected$Samples),
            !any(pooled_scores$True.Genotype == "Ts65Dn"),
            sum(pooled_scores$True.Genotype == "WT") == sum(pooled_expected$True.Genotype == "WT"),
            length(unique(pooled_scores$Experiment)) == 4L,
            all(file.exists(file.path(pooled_root, "figures", "pca", c("genotype", "batch", "tissue"),
                                      "Genotypes_AllExperiments_all_tissues.pdf"))))
  for (tissue in unique(pooled_expected$Tissue.Cell.type)) {
    prefix <- paste0("Genotypes_AllExperiments_", tissue)
    scores <- read.csv(file.path(pooled_root, "tables", "pca", "scores", paste0(prefix, ".csv")))
    stopifnot(setequal(scores$Samples, pooled_expected$Samples[pooled_expected$Tissue.Cell.type == tissue]),
              all(scores$Tissue.Cell.type == tissue),
              all(file.exists(file.path(pooled_root, "figures", "pca", c("genotype", "batch"),
                                        paste0(prefix, ".pdf")))))
  }

  for (id in names(env$results_list)) {
    result <- env$results_list[[id]]
    exp_dirs <- env$experiment_dirs[[result$experiment]]
    expected <- c(file.path(exp_dirs$dmp_tables, paste0(id, ".csv")),
                  file.path(exp_dirs$volcano_figures, paste0(id, ".pdf")),
                  file.path(exp_dirs$heatmap_figures, paste0(id, ".pdf")),
                  file.path(exp_dirs$fasta, c("GOM", "LOM"), paste0(id, ".fasta")),
                  file.path(exp_dirs$fasta_details, c("GOM", "LOM"), paste0(id, ".csv")))
    stopifnot(all(file.exists(expected)),
              isTRUE(all.equal(read.csv(expected[1], row.names = 1, colClasses = c(Experiment = "character")), result$table,
                               check.attributes = FALSE)))
    call <- env$heatmap_calls[[expected[3]]]
    stopifnot(setequal(call$samples, c(result$control_samples, result$experimental_samples)),
              identical(call$probes, result$dmps$Probe))
    members <- samples[match(call$samples, samples$Samples), ]
    stopifnot(all(members$Experiment == result$experiment),
              all(members$Tissue.Cell.type == result$tissue))
    if (result$experiment == "2") {
      controls <- members[members$True.Genotype == "WT", ]
      cases <- members[members$True.Genotype != "WT", ]
      stopifnot(nrow(controls) == 3L, all(controls$cross.type == result$wt_cross),
                nrow(cases) == 3L, all(cases$cross.type == "Double"))
    }
  }

  age_dirs <- env$experiment_dirs[["4"]]
  age_scores <- read.csv(file.path(age_dirs$pca_tables, "Age17.5", "scores", "Genotypes_Exp4_thymocytes.csv"))
  age_expected <- samples$Samples[samples$Experiment == 4 & samples$Tissue.Cell.type == "thymocytes" & samples$age == 17.5]
  stopifnot(setequal(age_scores$Samples, age_expected), nrow(age_scores) == 11L,
            all(age_scores$age == 17.5), all(age_scores$Experiment == 4),
            all(file.exists(file.path(age_dirs$pca_figures, "Age17.5", c("genotype", "batch"),
                                      "Genotypes_Exp4_thymocytes.pdf"))),
            !dir.exists(file.path(env$experiment_dirs[["5"]]$pca_figures, "Age17.5")))

  for (experiment in names(env$experiment_dirs)) {
    exp_dirs <- env$experiment_dirs[[experiment]]
    metadata <- read.csv(file.path(exp_dirs$metadata, "samples.csv"))
    stopifnot(all(metadata$Experiment == experiment))
    for (tissue in unique(metadata$Tissue.Cell.type)) {
      prefix <- paste0("Genotypes_Exp", experiment, "_", tissue)
      score_file <- file.path(exp_dirs$pca_tables, "scores", paste0(prefix, ".csv"))
      if (experiment == "4" && tissue == "thymus") {
        stopifnot(!file.exists(score_file))
        next
      }
      scores <- read.csv(score_file)
      stopifnot(all(scores$Experiment == experiment), all(scores$Tissue.Cell.type == tissue),
                setequal(scores$Samples, metadata$Samples[metadata$Tissue.Cell.type == tissue]),
                all(file.exists(file.path(exp_dirs$pca_figures, c("genotype", "batch"),
                                          paste0(prefix, ".pdf")))))
    }
  }

  for (experiment in c("3", "4")) {
    tissue <- if (experiment == "3") "thymus" else "thymocytes"
    exp_dirs <- env$experiment_dirs[[experiment]]
    filename <- paste0("RUNX_vs_WT_Exp", experiment, "_", tissue, ".pdf")
    union_call <- env$heatmap_calls[[file.path(exp_dirs$heatmap_figures, "combined", filename)]]
    focused <- env$heatmap_calls[[file.path(exp_dirs$heatmap_figures, "top1000", filename)]]
    expected_samples <- samples$Samples[samples$Experiment == experiment & samples$Tissue.Cell.type == tissue &
                                        samples$True.Genotype %in% c("WT", "DP16", "RUNX1-DEL_DP16")]
    stopifnot(setequal(union_call$samples, expected_samples),
              setequal(focused$samples, expected_samples),
              identical(unique(union_call$genotypes), c("WT", "Dp16", "Dp16Runx1")),
              identical(unique(focused$genotypes), c("WT", "Dp16", "Dp16Runx1")),
              union_call$directions["cg3"] == "GOM", focused$directions["cg3"] == "GOM",
              setequal(union_call$probes, paste0("cg", 1:3)))
  }
  manifest <- read.csv(file.path(dirs$root, "output_manifest.csv"))
  stopifnot(!anyDuplicated(manifest$Path),
            setequal(manifest$Path, setdiff(list.files(dirs$root, recursive = TRUE), "output_manifest.csv")),
            all(manifest$Bytes == file.info(file.path(dirs$root, manifest$Path))$size),
            !any(grepl("equal_mouse|by_tissue|heatmap_mean_centered|__|experiments_3_4", manifest$Path)),
            sum(endsWith(manifest$Path, ".pdf")) == 51L)
  fasta_summary <- read.csv(file.path(dirs$summaries, "FASTA_summary.csv"))
  stopifnot(nrow(fasta_summary) == 13L, all(fasta_summary$Sequences.Written == 3L))

  # Another run with the same name cannot overwrite any of these outputs.
  saved_files <- file.path(dirs$root, c(manifest$Path, "output_manifest.csv"))
  before <- tools::md5sum(saved_files)
  next_dirs <- create_output_dirs("results", basename(dirs$root))
  stopifnot(identical(next_dirs$root, paste0(dirs$root, "_01")),
            identical(before, tools::md5sum(saved_files)))

  cat("Output layout tests passed: 13 comparisons, 51 PDFs, exact pooled/experiment/tissue and age membership,",
      "RUNX direction calls, FASTA exports, short filenames, manifest, and preserved earlier runs.\n")
}
run_tests()
