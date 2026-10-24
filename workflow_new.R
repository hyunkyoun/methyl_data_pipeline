# Plan: IDAT -> QC -> ComBAT -> DMP analysis -> PCA/Plotss -> FASTA extraction

# Setup
package_dir <- file.path(getwd(), "packages")
if (!dir.exists(package_dir)) {
  stop("Run source(\"install_packages.R\") from methyl_data_pipeline before starting the workflow.")
}
.libPaths(c(package_dir, .libPaths()))

required_libs <- c("sesame", "minfi", "IlluminaMouseMethylationmanifest", 
                   "IlluminaMouseMethylationanno.12.v1.mm10", "limma", "sva", 
                   "ggsci", "wateRmelon", "snow", "BiocParallel", "parallel", 
                   "RColorBrewer", "readxl", "openxlsx", "ggplot2", "gridExtra", "pheatmap",
                   "GenomicRanges", "GenomicFeatures", "BSgenome.Mmusculus.UCSC.mm10", 
                   "Biostrings")
lapply(required_libs, library, character.only = TRUE)
annoMouse <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)

# IDAT
idir <- "idats"
samplesheet_path <- file.path(idir, "samplesheet_169_samples.csv")
cores <- 4

# Load analysis samples before choosing IDAT pairs; brain is excluded throughout.
targets_raw <- read.csv(samplesheet_path, header = TRUE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
names(targets_raw) <- make.names(names(targets_raw))
if (!all(c("Basename", "Tissue.Cell.type") %in% names(targets_raw))) {
  stop("Samplesheet requires Basename and Tissue.Cell.type columns.")
}
targets_raw <- targets_raw[!tolower(trimws(targets_raw$Tissue.Cell.type)) %in% "brain", , drop = FALSE]
targets_raw$Basename <- as.character(targets_raw$Basename)
if (anyNA(targets_raw$Basename) || any(!nzchar(trimws(targets_raw$Basename)))) {
  stop("Samplesheet contains missing or blank Basename identifiers.")
}
if (anyDuplicated(targets_raw$Basename)) {
  stop("Samplesheet contains duplicate Basename identifiers; each IDAT pair must be unique.")
}
targets_raw$Samples <- paste0(targets_raw$Sample.ID, "_", targets_raw$Experiment)
targets_raw$EPIC_RUN_CLEAN <- targets_raw$EPIC.Run
idat_prefixes <- searchIDATprefixes(idir)
if (any(!targets_raw$Basename %in% names(idat_prefixes))) {
  stop("Samplesheet contains Basenames absent from complete IDAT pairs.")
}
idat_prefixes <- idat_prefixes[names(idat_prefixes) %in% targets_raw$Basename]

# Saving files
source("R/output_paths.R")
results_root <- "results"
output_dirs <- create_output_dirs(results_root)
cat("Saving this analysis under:", output_dirs$root, "\n")
for (source_path in c(samplesheet_path, "workflow_new.R", "R/dmp_analysis.R", "R/output_paths.R",
                      "R/pca_analysis.R", "R/heatmap_plots.R")) {
  if (!file.copy(source_path, file.path(output_dirs$metadata, basename(source_path)))) {
    stop("Could not save run metadata copy: ", source_path)
  }
}

cat("Reading IDAT files")
betas <- openSesame(idat_prefixes, BPPARAM = MulticoreParam(workers = cores), collapseToPfx = TRUE)
detP <- openSesame(idat_prefixes, func = pOOBAH, return.pval = TRUE, BPPARAM = MulticoreParam(workers = cores))
detP <- betasCollapseToPfx(detP)

# Alignment
common_samples <- intersect(colnames(betas), colnames(detP))
betas <- betas[, common_samples, drop = FALSE]
detP <- detP[, common_samples, drop = FALSE]

write.csv(betas, file.path(output_dirs$betas, "betas_raw_idat.csv"), row.names = TRUE)
cat("Raw data saved:", nrow(betas), "probes x", ncol(betas), "samples\n")

# QC
cat("Applying QC")

qc_filter_run <- function(betas_run, detP_run, targets_run, run_id) {
  good_samples <- which(colMeans(detP_run, na.rm = TRUE) < 0.05)
  if (length(good_samples) == 0) stop(paste("No good samples in RUN", run_id))

  betas_run <- betas_run[, good_samples, drop = FALSE]
  detP_run <- detP_run[, good_samples, drop = FALSE]
  targets_run <- targets_run[good_samples, , drop = FALSE]

  keep_probes <- rowSums(detP_run < 0.05, na.rm = TRUE) == ncol(betas_run)
  betas_run <- betas_run[keep_probes, , drop = FALSE]
  detP_run <- detP_run[keep_probes, , drop = FALSE]

  complete_idx <- complete.cases(betas_run)
  betas_run <- betas_run[complete_idx, , drop = FALSE]
  detP_run <- detP_run[complete_idx, , drop = FALSE]

  cpg_idx <- grep("^cg", rownames(betas_run))
  betas_run <- betas_run[cpg_idx, , drop = FALSE]
  sex_cpgs <- annoMouse$Name[annoMouse$chr %in% c("chrX", "chrY")]
  sex_cpgs_base <- sub("_.*$", "", sex_cpgs)
  betas_run <- betas_run[!rownames(betas_run) %in% sex_cpgs_base, , drop = FALSE]
  detP_run <- detP_run[rownames(betas_run), , drop = FALSE]

  return(list(betas = betas_run, detP = detP_run, targets = targets_run))
}

common_sample_ids <- intersect(colnames(betas), targets_raw$Basename)
betas <- betas[, common_sample_ids, drop = FALSE]
detP <- detP[, common_sample_ids, drop = FALSE]
targets_matched <- targets_raw[match(common_sample_ids, targets_raw$Basename), , drop = FALSE]

# By run QC
filtered_data <- list()
unique_runs <- sort(unique(targets_matched$EPIC_RUN_CLEAN))

for (run_id in unique_runs) {
  run_idx <- targets_matched$EPIC_RUN_CLEAN == run_id
  result <- qc_filter_run(betas[, run_idx, drop = FALSE],
                          detP[, run_idx, drop = FALSE],
                          targets_matched[run_idx, , drop = FALSE],
                          run_id)
  filtered_data[[paste0("RUN_", gsub("_", "", run_id))]] <- result
}

common_cpgs_qc <- Reduce(intersect, lapply(filtered_data, function(x) rownames(x$betas)))
betas_qc_combined <- do.call(cbind, lapply(filtered_data, function(x) x$betas[common_cpgs_qc, , drop = FALSE]))
detP_qc_combined <- do.call(cbind, lapply(filtered_data, function(x) x$detP[common_cpgs_qc, , drop = FALSE]))
targets_qc_combined <- do.call(rbind, lapply(filtered_data, function(x) x$targets))

# Save filtered data
write.csv(betas_qc_combined, file.path(output_dirs$betas, "betas_qc_filtered.csv"), row.names = TRUE)
cat("QC filtered data:", nrow(betas_qc_combined), "probes x", ncol(betas_qc_combined), "samples\n")

betas <- betas_qc_combined
detP <- detP_qc_combined
targets <- targets_qc_combined

# Samplesheet matching
cat("Matching to samplesheet")

targets <- targets_raw
current_samples <- colnames(betas)
matched_indices <- match(current_samples, targets$Basename)
if (anyNA(matched_indices)) {
  stop("QC retained an IDAT sample without a matching samplesheet Basename.")
}
targets <- targets[matched_indices, , drop = FALSE]

colnames(betas) <- targets$Samples

# Technology between old vs new
technology_by_run <- c("1_2" = "BATCH1", "3_4" = "BATCH2", "5" = "BATCH2", "6" = "BATCH2")
targets$BATCH <- unname(technology_by_run[as.character(targets$EPIC_RUN_CLEAN)])

betas_combined <- betas
matching_samples <- targets
cat("Matched dataset:", nrow(betas_combined), "probes x", ncol(betas_combined), "samples\n")

# ComBAT
biological_columns <- c("Experiment", "Tissue.Cell.type", "True.Genotype")
if (!all(biological_columns %in% names(matching_samples)) ||
    any(vapply(matching_samples[biological_columns], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1)))) {
  stop("ComBat requires non-missing Experiment, tissue, and True.Genotype metadata.")
}
sample_ids <- as.character(matching_samples$Samples)
batch <- factor(matching_samples$BATCH)
group <- interaction(matching_samples$Experiment, matching_samples$Tissue.Cell.type,
                     matching_samples$True.Genotype, drop = TRUE)
skip_combat <- nlevels(batch) < 2

if (skip_combat) {
  cat("Skipping ComBat (single technology)\n")
  betas_adjusted <- betas_combined
} else {
  cat("Running ComBat while preserving biological group effects...\n")

  # Protect the mean of each experiment/tissue/genotype combination.
  if (nlevels(group) > 1) {
    mod <- model.matrix(~ group)
  } else {
    mod <- model.matrix(~ 1, data = matching_samples)
  }

  # Match ComBat's joint design: batch indicators plus non-intercept group terms.
  batch_design <- model.matrix(~ 0 + batch)
  group_design <- mod[, colnames(mod) != "(Intercept)", drop = FALSE]
  combat_design <- cbind(batch_design, group_design)
  if (qr(combat_design)$rank < ncol(combat_design)) {
    stop("Group and technology are fully confounded after QC; ComBat cannot separate them.")
  }
  if (nrow(combat_design) <= ncol(combat_design)) {
    stop("Too few samples remain to estimate residual variation for ComBat.")
  }

  # Keep the existing policy: drop CpGs with non-finite M-values, including beta 0/1.
  mvals_combined <- log2(betas_combined / (1 - betas_combined))
  finite_probes <- rowSums(!is.finite(mvals_combined)) == 0
  mvals_combined <- mvals_combined[finite_probes, , drop = FALSE]
  if (nrow(mvals_combined) < 2) {
    stop("ComBat requires at least two CpGs with finite M-values in every sample.")
  }

  mvals_combat <- ComBat(
    dat = mvals_combined,
    batch = batch,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  if (any(!is.finite(mvals_combat))) {
    stop("ComBat returned non-finite M-values; inspect the retained probes and design.")
  }

  # Numerically stable inverse of log2(beta / (1 - beta)); values stay within 0-1.
  betas_adjusted <- plogis(mvals_combat * log(2))
}

write.csv(betas_adjusted, file.path(output_dirs$betas, "betas_combat_adjusted.csv"), row.names = TRUE)
write.csv(matching_samples, file.path(output_dirs$metadata, "samples_after_qc.csv"), row.names = FALSE)
cat("Final dataset:", nrow(betas_adjusted), "probes x", ncol(betas_adjusted), "samples\n")

# DMP analysis
source("R/dmp_analysis.R")

# Every experiment has its own output folder and genotype-versus-WT comparisons.
experiments <- sort(unique(as.character(matching_samples$Experiment)))
same_parent_experiments <- c("1", "2")
experiment_dirs <- setNames(lapply(experiments, function(experiment) {
  dirs <- experiment_output_dirs(output_dirs, experiment, same_parent = experiment %in% same_parent_experiments)
  write.csv(matching_samples[matching_samples$Experiment == experiment, , drop = FALSE],
            file.path(dirs$metadata, "samples.csv"), row.names = FALSE)
  dirs
}), experiments)

# Shared significance rules for result tables, volcanoes, heatmaps, and FASTA.
dmp_fdr_cutoff <- 0.05
dmp_delta_beta_cutoff <- 0.1

dmp_analysis <- run_dmp_analysis(
  betas_adjusted, matching_samples,
  fdr_cutoff = dmp_fdr_cutoff, delta_beta_cutoff = dmp_delta_beta_cutoff,
  same_parent_experiments = same_parent_experiments
)
results_list <- dmp_analysis$results
write.csv(dmp_analysis$contrast_weights, file.path(output_dirs$summaries, "DMP_contrast_weights.csv"), row.names = FALSE)

dmp_summary_rows <- list()
for (analysis_id in names(results_list)) {
  result <- results_list[[analysis_id]]
  dirs <- experiment_dirs[[result$experiment]]
  write.csv(result$table, file.path(dirs$dmp_tables, paste0(analysis_id, ".csv")), row.names = TRUE)
  dmp_summary_rows[[analysis_id]] <- data.frame(
    Contrast = result$contrast_name,
    Comparison = result$contrast_expression,
    Experiment = result$experiment,
    Tissue = result$tissue,
    WT.Cross = result$wt_cross,
    Weighting = result$weighting,
    Weighting.Label = result$weighting_label,
    Experimental.Samples = length(result$experimental_samples),
    Control.Samples = length(result$control_samples),
    Tested.CpGs = nrow(result$table),
    GOM = sum(result$dmps$diff > 0),
    LOM = sum(result$dmps$diff < 0),
    FDR.Cutoff = dmp_fdr_cutoff,
    Delta.Beta.Cutoff = dmp_delta_beta_cutoff
  )
}
if (length(dmp_summary_rows) > 0) {
  dmp_summary <- do.call(rbind, dmp_summary_rows)
  write.csv(dmp_summary, file.path(output_dirs$summaries, "DMP_summary.csv"), row.names = FALSE)
}

# PCA from top 10000 CpGs
source("R/pca_analysis.R")
pca_ntop <- 10000L
pca_palettes <- pca_color_palettes(matching_samples)

# Restore pooled overviews, excluding Ts65Dn and its four Experiment 5 WT controls.
pooled_pca_excluded_wt <- c("11984_5", "11985_5", "11988_5", "11989_5")
pooled_pca_samples <- matching_samples[
  !tolower(trimws(matching_samples$Tissue.Cell.type)) %in% "brain" &
    !matching_samples$True.Genotype %in% "Ts65Dn" &
    !matching_samples$Samples %in% pooled_pca_excluded_wt, , drop = FALSE
]
pooled_pca_root <- file.path(output_dirs$root, "All_Experiments")
pooled_pca_result <- create_simple_pca(
  betas_adjusted, pooled_pca_samples, pca_ntop,
  output_prefix = "Genotypes_AllExperiments_all_tissues",
  figure_dir = file.path(pooled_pca_root, "figures", "pca"),
  table_dir = file.path(pooled_pca_root, "tables", "pca"),
  scope_label = "All experiments | all tissues", color_palettes = pca_palettes
)
pooled_tissues <- unique(as.character(pooled_pca_samples$Tissue.Cell.type))
pooled_tissue_pca_results <- setNames(lapply(pooled_tissues, function(tissue) {
  create_simple_pca(
    betas_adjusted, pooled_pca_samples[pooled_pca_samples$Tissue.Cell.type == tissue, , drop = FALSE],
    pca_ntop, output_prefix = paste0("Genotypes_AllExperiments_", tissue),
    figure_dir = file.path(pooled_pca_root, "figures", "pca"),
    table_dir = file.path(pooled_pca_root, "tables", "pca"),
    scope_label = paste0("All experiments | ", tissue),
    plot_tissue = FALSE, color_palettes = pca_palettes
  )
}), pooled_tissues)

# Refit PCA and select variable CpGs using only one experiment/tissue at a time.
tissue_pca_results <- setNames(lapply(experiments, function(experiment) {
  dirs <- experiment_dirs[[experiment]]
  experiment_samples <- matching_samples[matching_samples$Experiment == experiment, , drop = FALSE]
  tissue_results <- create_tissue_pcas(
    betas_adjusted, experiment_samples,
    pca_ntop, output_prefix = paste0("Genotypes_Exp", experiment),
    figure_dir = dirs$pca_figures, table_dir = dirs$pca_tables, color_palettes = pca_palettes
  )
  # Additional age-matched thymocyte PCA, keeping each experiment separate.
  if ("age" %in% names(experiment_samples)) {
    age_samples <- experiment_samples[experiment_samples$Tissue.Cell.type == "thymocytes" &
                                        experiment_samples$age %in% 17.5, , drop = FALSE]
    if (nrow(age_samples) > 0L) {
      create_simple_pca(
        betas_adjusted, age_samples, pca_ntop,
        output_prefix = paste0("Genotypes_Exp", experiment, "_thymocytes"),
        figure_dir = file.path(dirs$pca_figures, "Age17.5"),
        table_dir = file.path(dirs$pca_tables, "Age17.5"),
        scope_label = paste0("Experiment ", experiment, " | thymocytes"),
        plot_title = "Age17.5", plot_tissue = FALSE, color_palettes = pca_palettes
      )
    }
  }
  tissue_results
}), experiments)

# Stored DMP results plotted
source("R/heatmap_plots.R")
plot_mode <- "both" # both includes volcano and heatmaps
runx_heatmap_top_n <- 1000L

# Pairwise plots reuse the stored DMP calls and comparison samples.
for (analysis_id in names(results_list)) {
  result <- results_list[[analysis_id]]
  contrast_name <- result$contrast_name
  results <- result$table
  DMPs <- result$dmps
  dirs <- experiment_dirs[[result$experiment]]
  plot_title <- sprintf("%s vs WT | Experiment %s | %s\nSamples: %d vs %d",
                        result$genotype, result$experiment, result$tissue,
                        length(result$experimental_samples), length(result$control_samples))
  if (result$wt_cross != "All") {
    plot_title <- paste0(plot_title, "\nWT parental cross: ", result$wt_cross)
  }

  cat("\n--- Plotting", contrast_name, "-", result$weighting_label, "---\n")

  # Volcano colors and the FDR axis use the same calls as heatmaps and FASTA.
  if (plot_mode %in% c("both", "volcano")) {
    xx <- results
    xx$col <- "ivory3"
    xx$col[xx$Significant & xx$diff > 0] <- "indianred"
    xx$col[xx$Significant & xx$diff < 0] <- "deepskyblue3"
    xx <- xx[order(xx$adj.P.Val, decreasing = TRUE), , drop = FALSE]

    # Keep zero centered and expand to include every effect and both cutoff lines.
    x_limit <- max(0.5, abs(xx$diff), dmp_delta_beta_cutoff, na.rm = TRUE)
    pdf(file.path(dirs$volcano_figures, paste0(analysis_id, ".pdf")))
    plot(xx$diff, -log10(pmax(xx$adj.P.Val, .Machine$double.xmin)),
         xlim = c(-x_limit, x_limit),
         col = xx$col, pch = 16,
         ylab = "-log10(FDR)",
         xlab = "Beta difference (experimental - control)",
         main = plot_title)
    abline(h = -log10(dmp_fdr_cutoff), lty = 3)
    abline(v = c(-dmp_delta_beta_cutoff, dmp_delta_beta_cutoff),
           lty = 2, col = "gray50")
    dev.off()
  }

  if (plot_mode %in% c("both", "heatmap") && nrow(DMPs) > 0) {
    plot_dmp_heatmap(
      DMPs$Probe, c(result$control_samples, result$experimental_samples), plot_title,
      file.path(dirs$heatmap_figures, paste0(analysis_id, ".pdf")),
      directions = dmp_heatmap_directions(DMPs$Probe, list(DMPs)),
      betas_data = betas_adjusted, targets_data = matching_samples
    )
  }

  cat("Significant GOM:", sum(DMPs$diff > 0), "\n")
  cat("Significant LOM:", sum(DMPs$diff < 0), "\n")
  cat("DMP cutoff: FDR <", dmp_fdr_cutoff,
      "and |delta beta| >", dmp_delta_beta_cutoff, "\n")
}

# Combined RUNX heatmaps
if (plot_mode %in% c("both", "heatmap")) {
  runx_genotypes <- c("WT", "DP16", "RUNX1-DEL_DP16")
  runx_labels <- c(WT = "WT", DP16 = "Dp16", "RUNX1-DEL_DP16" = "Dp16Runx1")
  runx_cohorts <- unique(matching_samples[matching_samples$True.Genotype %in% "RUNX1-DEL_DP16",
                                          c("Experiment", "Tissue.Cell.type")])
  for (i in seq_len(nrow(runx_cohorts))) {
    experiment <- as.character(runx_cohorts$Experiment[i])
    tissue <- runx_cohorts$Tissue.Cell.type[i]
    runx_results <- Filter(function(x) x$experiment == experiment && x$tissue == tissue &&
                            x$genotype %in% runx_genotypes, results_list)
    probes <- unique(unlist(lapply(runx_results, function(x) x$dmps$Probe), use.names = FALSE))
    if (length(probes) == 0) {
      cat("Skipping combined RUNX", tissue, "heatmap: no significant DMPs.\n")
      next
    }
    dirs <- experiment_dirs[[experiment]]
    runx_heatmap_dir <- file.path(dirs$heatmap_figures, "combined")
    focused_heatmap_dir <- file.path(dirs$heatmap_figures, paste0("top", runx_heatmap_top_n))
    dir.create(runx_heatmap_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(focused_heatmap_dir, recursive = TRUE, showWarnings = FALSE)
    runx_filename <- paste0("RUNX_vs_WT_Exp", experiment, "_", tissue, ".pdf")
    panel <- matching_samples[matching_samples$Experiment == experiment &
                                matching_samples$Tissue.Cell.type == tissue &
                                matching_samples$True.Genotype %in% runx_genotypes, , drop = FALSE]
    panel <- panel[order(match(panel$True.Genotype, runx_genotypes)), , drop = FALSE]
    counts <- table(factor(panel$True.Genotype, levels = runx_genotypes))
    names(counts) <- unname(runx_labels[names(counts)])
    panel$True.Genotype <- unname(runx_labels[as.character(panel$True.Genotype)])
    plot_title <- paste0(
      "Dp16 and Dp16Runx1 vs WT | Experiment ", experiment, " | ", tissue, "\n",
      paste0(names(counts), " (n=", counts, ")", collapse = " | "), "\n",
      length(probes), " unique DMPs from Dp16 and Dp16Runx1 comparisons vs WT"
    )
    cat("Plotting combined RUNX", tissue, ":", nrow(panel), "samples and", length(probes), "DMPs\n")
    plot_dmp_heatmap(
      probes, panel$Samples, plot_title,
      file.path(runx_heatmap_dir, runx_filename),
      directions = dmp_heatmap_directions(probes, lapply(runx_results, `[[`, "dmps")),
      betas_data = betas_adjusted, targets_data = panel,
      direction_note = "GOM/LOM: Dp16 and Dp16Runx1 vs WT; Mixed: opposite significant directions"
    )

    # Select rows using only DP16 vs shared WT, then show the three retained genotypes.
    dp16_result <- Filter(function(x) x$genotype == "DP16", runx_results)
    if (length(dp16_result) == 0L || nrow(dp16_result[[1]]$dmps) == 0L) {
      cat("Skipping focused RUNX", tissue, "heatmap: no significant DP16-vs-WT DMPs.\n")
      next
    }
    ranked_dmps <- dp16_result[[1]]$dmps
    ranked_dmps <- ranked_dmps[order(ranked_dmps$adj.P.Val, ranked_dmps$Probe), , drop = FALSE]
    top_probes <- head(unique(ranked_dmps$Probe), runx_heatmap_top_n)
    plot_dmp_heatmap(
      top_probes, panel$Samples,
      paste0("Dp16 and Dp16Runx1 vs WT | Experiment ", experiment, " | ", tissue, "\n",
             paste0(names(counts), " (n=", counts, ")", collapse = " | "), "\n",
             "Top ", length(top_probes), " of ", nrow(ranked_dmps),
             " significant DP16 vs WT CpGs, ranked by FDR"),
      file.path(focused_heatmap_dir, runx_filename),
      directions = dmp_heatmap_directions(top_probes, list(ranked_dmps)),
      betas_data = betas_adjusted, targets_data = panel,
      direction_note = "GOM: gain of methylation; LOM: loss of methylation (DP16 vs WT)"
    )
  }
}

# FASTA sequences from DMPS
extract_and_write_FASTA <- function(DMPs, contrast_name, annoMouse, genome, fasta_dir, details_dir) {
  cat("\n--- Extracting FASTA sequences for", contrast_name, "---\n")

  # Clean annotation IDs: remove suffix (e.g. "_BC21", "_TC21")
  annoMouse$CleanName <- sub("_.*", "", annoMouse$Name)

  # Match probes to annotation
  anno_subset <- annoMouse[match(DMPs$Probe, annoMouse$CleanName), c("chr","pos")]
  DMPs$chr <- anno_subset$chr
  DMPs$pos <- anno_subset$pos

  # Keep only valid CpGs
  valid <- !is.na(DMPs$chr) & !is.na(DMPs$pos)
  if (sum(valid) == 0) {
    cat("No valid genomic coordinates for", contrast_name, "\n")
    return(list(GOM = 0, LOM = 0, written = 0))
  }
  DMPs <- DMPs[valid, ]

  # Split into hyper/hypo
  GOM <- DMPs[DMPs$diff > 0, ]
  LOM <- DMPs[DMPs$diff < 0, ]

  # Build ±100 bp windows
  GOM.gr <- with(GOM, GRanges(chr, IRanges(pos - 99, pos + 100)))
  LOM.gr <- with(LOM, GRanges(chr, IRanges(pos - 99, pos + 100)))

  # Reduce overlaps
  GOM.gr <- GenomicRanges::reduce(GOM.gr)
  LOM.gr <- GenomicRanges::reduce(LOM.gr)

  # Extract and name sequences only when that direction has genomic regions.
  GOM.fasta <- Biostrings::DNAStringSet()
  LOM.fasta <- Biostrings::DNAStringSet()
  if (length(GOM.gr) > 0) {
    GOM.fasta <- getSeq(genome, GOM.gr)
    names(GOM.fasta) <- paste0(seqnames(GOM.gr), ":", start(GOM.gr), "-", end(GOM.gr))
  } else {
    cat("No methylation gains (GOM) found for", contrast_name, "\n")
  }
  if (length(LOM.gr) > 0) {
    LOM.fasta <- getSeq(genome, LOM.gr)
    names(LOM.fasta) <- paste0(seqnames(LOM.gr), ":", start(LOM.gr), "-", end(LOM.gr))
  } else {
    cat("No methylation losses (LOM) found for", contrast_name, "\n")
  }

  # Empty categories write empty outputs, clearing any sequences from an earlier run.
  for (direction in c("GOM", "LOM")) {
    dir.create(file.path(fasta_dir, direction), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(details_dir, direction), recursive = TRUE, showWarnings = FALSE)
  }
  writeXStringSet(GOM.fasta, file.path(fasta_dir, "GOM", paste0(contrast_name, ".fasta")))
  writeXStringSet(LOM.fasta, file.path(fasta_dir, "LOM", paste0(contrast_name, ".fasta")))
  write.csv(GOM, file.path(details_dir, "GOM", paste0(contrast_name, ".csv")), row.names = FALSE)
  write.csv(LOM, file.path(details_dir, "LOM", paste0(contrast_name, ".csv")), row.names = FALSE)

  cat("Hypermethylated (GOM):", nrow(GOM), "\n")
  cat("Hypomethylated (LOM):", nrow(LOM), "\n")
  cat("Total FASTA sequences written:", length(GOM.fasta) + length(LOM.fasta), "\n")

  return(list(GOM = nrow(GOM), LOM = nrow(LOM),
              written = length(GOM.fasta) + length(LOM.fasta)))
}

# Extract FASTA from the stored DMP calls
genome <- BSgenome.Mmusculus.UCSC.mm10
summary_list <- list()

for (analysis_id in names(results_list)) {
  result <- results_list[[analysis_id]]
  dirs <- experiment_dirs[[result$experiment]]
  cat("\n=== Extracting FASTA:", result$contrast_name, "-",
      result$weighting_label, "===\n")

  significant_tt <- result$dmps
  if (nrow(significant_tt) > 0) {
    summary_list[[analysis_id]] <- extract_and_write_FASTA(
      significant_tt, analysis_id, annoMouse, genome,
      fasta_dir = dirs$fasta,
      details_dir = dirs$fasta_details
    )
  } else {
    cat("No significant DMPs found for", analysis_id, "\n")
    summary_list[[analysis_id]] <- list(GOM = 0, LOM = 0, written = 0)
  }
}

# Summary
fasta_summary <- data.frame(
  Analysis = as.character(names(summary_list)),
  GOM.Annotated.CpGs = vapply(summary_list, `[[`, numeric(1), "GOM"),
  LOM.Annotated.CpGs = vapply(summary_list, `[[`, numeric(1), "LOM"),
  Sequences.Written = vapply(summary_list, `[[`, numeric(1), "written")
)
write.csv(fasta_summary, file.path(output_dirs$summaries, "FASTA_summary.csv"), row.names = FALSE)
cat("\n=== FASTA EXTRACTION SUMMARY ===\n")
total_written <- 0
for (contrast_name in names(summary_list)) {
  res <- summary_list[[contrast_name]]
  cat(contrast_name, ":",
      res$GOM, "GOM,", res$LOM, "LOM,",
      res$written, "sequences written\n")
  total_written <- total_written + res$written
}
if (total_written == 0) {
  cat("No FASTA files were created.\n")
  cat("Possible issues:\n")
  cat("1. No significant DMPs found (try relaxing thresholds)\n")
  cat("2. Tissue/genotype matching issues (check the samplesheet)\n")
  cat("3. Annotation matching problems (check annoMouse object)\n")
} else {
  cat("FASTA files have been created in the experiment folders under:", output_dirs$root, "\n")
}

# ---- Save run metadata and output index ----
run_settings <- data.frame(
  Setting = c("IDAT.Directory", "Samplesheet", "DMP.FDR.Cutoff", "DMP.Delta.Beta.Cutoff",
              "PCA.Top.CpGs.Requested", "RUNX.Heatmap.Top.DMPs.Requested",
              "Plot.Mode", "ComBat.Skipped", "Completed.UTC"),
  Value = c(normalizePath(idir, winslash = "/"), normalizePath(samplesheet_path, winslash = "/"),
            dmp_fdr_cutoff, dmp_delta_beta_cutoff, pca_ntop, runx_heatmap_top_n, plot_mode, skip_combat,
            format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
)
write.csv(run_settings, file.path(output_dirs$metadata, "run_settings.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dirs$metadata, "sessionInfo.txt"))
write_output_manifest(output_dirs)
cat("All results saved under:", output_dirs$root, "\n")
