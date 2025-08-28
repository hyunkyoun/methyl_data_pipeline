# Methylation Data Analysis Pipeline (Clean)
# IDAT -> QC -> ComBat -> PCA

# ---- Setup 
.libPaths("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/packages")
setwd("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline")

required_libs <- c("sesame", "minfi", "IlluminaMouseMethylationmanifest", 
                   "IlluminaMouseMethylationanno.12.v1.mm10", "limma", "sva", 
                   "ggsci", "wateRmelon", "snow", "BiocParallel", "parallel", 
                   "RColorBrewer", "readxl", "openxlsx", "ggplot2", "gridExtra", "pheatmap",
                   "GenomicRanges", "GenomicFeatures", "BSgenome.Mmusculus.UCSC.mm10")
lapply(required_libs, library, character.only = TRUE)

source("bmiq/DoBMIQ.R")
annoMouse <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)

# ---- Step 1: Read IDAT files
idir <- "idat"
cores <- 4

cat("Reading IDAT files...\n")
betas <- openSesame(idir, BPPARAM = MulticoreParam(workers = cores), collapseToPfx = TRUE)
detP <- openSesame(idir, func = pOOBAH, return.pval = TRUE, BPPARAM = MulticoreParam(workers = cores))
detP <- betasCollapseToPfx(detP)

# Align datasets
common_samples <- intersect(colnames(betas), colnames(detP))
betas <- betas[, common_samples]
detP <- detP[, common_samples]

write.csv(betas, "betas_raw_idat.csv", row.names = TRUE)
cat("Raw data saved:", nrow(betas), "probes x", ncol(betas), "samples\n")

# ---- Step 2: Per-run QC filtering ----
cat("Applying per-run QC filtering...\n")

# Define QC filtering function per run
qc_filter_run <- function(betas_run, detP_run, targets_run, run_id) {
  good_samples <- which(colMeans(detP_run, na.rm = TRUE) < 0.05)
  if (length(good_samples) == 0) stop(paste("No good samples in RUN", run_id))

  betas_run <- betas_run[, good_samples, drop = FALSE]
  detP_run <- detP_run[, good_samples, drop = FALSE]
  targets_run <- targets_run[good_samples, , drop = FALSE]

  keep_probes <- rowSums(detP_run < 0.05, na.rm = TRUE) == ncol(betas_run)
  betas_run <- betas_run[keep_probes, ]
  detP_run <- detP_run[keep_probes, ]

  complete_idx <- complete.cases(betas_run)
  betas_run <- betas_run[complete_idx, ]
  detP_run <- detP_run[complete_idx, ]

  cpg_idx <- grep("^cg", rownames(betas_run))
  betas_run <- betas_run[cpg_idx, ]
  sex_cpgs <- annoMouse$Name[annoMouse$chr %in% c("chrX", "chrY")]
  sex_cpgs_base <- sub("_.*$", "", sex_cpgs)
  betas_run <- betas_run[!rownames(betas_run) %in% sex_cpgs_base, ]
  detP_run <- detP_run[rownames(betas_run), ]

  return(list(betas = betas_run, detP = detP_run, targets = targets_run))
}

# Load targets and prepare for filtering
targets_raw <- read.csv("data/samplesheet.csv", header = TRUE)
names(targets_raw) <- make.names(names(targets_raw))
targets_raw$Beadchip <- as.character(targets_raw$Beadchip)
targets_raw$Basename <- paste0(targets_raw$Beadchip, "_", targets_raw$Position)
targets_raw$Samples <- paste0(targets_raw$Sample.ID, "_", targets_raw$Experiment)
targets_raw$EPIC_RUN_CLEAN <- targets_raw$EPIC.Run
targets_raw <- targets_raw[!duplicated(targets_raw$Basename), ]

# Match samples
common_sample_ids <- intersect(colnames(betas), targets_raw$Basename)
betas <- betas[, common_sample_ids]
detP <- detP[, common_sample_ids]
targets_matched <- targets_raw[match(common_sample_ids, targets_raw$Basename), ]

# Apply QC filtering by run
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

# Find common CpGs across runs and merge
common_cpgs_qc <- Reduce(intersect, lapply(filtered_data, function(x) rownames(x$betas)))
betas_qc_combined <- do.call(cbind, lapply(filtered_data, function(x) x$betas[common_cpgs_qc, , drop = FALSE]))
detP_qc_combined <- do.call(cbind, lapply(filtered_data, function(x) x$detP[common_cpgs_qc, , drop = FALSE]))
targets_qc_combined <- do.call(rbind, lapply(filtered_data, function(x) x$targets))

# Save filtered data
write.csv(betas_qc_combined, "betas_qc_filtered.csv", row.names = TRUE)
cat("QC filtered data:", nrow(betas_qc_combined), "probes x", ncol(betas_qc_combined), "samples\n")

# Update global variables for consistency with later steps
betas <- betas_qc_combined
detP <- detP_qc_combined
targets <- targets_qc_combined

# ---- Step 3: Match to samplesheet
cat("Matching to samplesheet...\n")

targets <- read.csv("data/samplesheet.csv", header = TRUE)
names(targets) <- make.names(names(targets))
targets$Beadchip <- as.character(targets$Beadchip)
targets$Basename <- paste0(targets$Beadchip, "_", targets$Position)
targets$Samples <- paste0(targets$Sample.ID, "_", targets$Experiment)
targets$EPIC_RUN_CLEAN <- targets$EPIC.Run
targets <- targets[!duplicated(targets$Basename), ]

# Match samples
current_samples <- colnames(betas)
matched_indices <- match(current_samples, targets$Basename)
valid_matches <- !is.na(matched_indices)

matched_samples <- current_samples[valid_matches]
betas <- betas[, matched_samples]
detP <- detP[, matched_samples]
targets <- targets[matched_indices[valid_matches], ]
targets <- targets[match(colnames(betas), targets$Basename), ]

# Rename to Sample IDs
colnames(betas) <- targets$Samples
targets$BATCH <- ifelse(targets$EPIC_RUN_CLEAN == "1_2", "BATCH1", "BATCH2")

betas_combined <- betas
matching_samples <- targets

cat("Matched dataset:", nrow(betas_combined), "probes x", ncol(betas_combined), "samples\n")

# ---- Step 4: ComBat normalization
use_tissue_design <- length(unique(targets$Tissue.Cell.type)) >= 2
skip_combat <- length(unique(targets$BATCH)) < 2

if (skip_combat) {
  cat("Skipping ComBat (single batch)\n")
  betas_adjusted <- betas_combined
} else {
  cat("Running ComBat normalization...\n")
  
  # Convert to M-values
  mvals_combined <- log2(betas_combined / (1 - betas_combined))
  mvals_combined[is.infinite(mvals_combined)] <- NA
  mvals_combined <- mvals_combined[complete.cases(mvals_combined), ]
  
  # Design matrix
  if (use_tissue_design) {
    mod <- model.matrix(~as.factor(Tissue.Cell.type), data=matching_samples)
  } else {
    mod <- model.matrix(~1, data=matching_samples)
  }
  
  # Run ComBat
  mvals_combat <- ComBat(
    dat = mvals_combined,
    batch = as.numeric(factor(matching_samples$BATCH)),
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  # Convert back to beta values
  betas_adjusted <- 2^mvals_combat / (1 + 2^mvals_combat)
  betas_adjusted[betas_adjusted < 0] <- 0
  betas_adjusted[betas_adjusted > 1] <- 1
}

write.csv(betas_adjusted, "betas_combat_adjusted.csv", row.names = TRUE)
cat("Final dataset:", nrow(betas_adjusted), "probes x", ncol(betas_adjusted), "samples\n")

# ---- Step 5 + Step 6: Extract top variable CpGs & Generate PCA plots 
create_simple_pca <- function(betas_data, targets_data, ntop = 100, output_prefix) {
  cat("Processing", output_prefix, "data...\n")
  
  # Transpose to match Python: rows = samples, columns = CpGs
  df_T <- t(betas_data)
  
  # Select top N most variable CpGs 
  cpg_variances <- apply(betas_data, 1, var, na.rm = TRUE)
  top_cpgs <- names(sort(cpg_variances, decreasing = TRUE)[1:ntop])
  df_T_reduced <- df_T[, top_cpgs]
  
  cat("Selected top", ntop, "CpGs based on variance\n")
  
  # Ensure sample alignment with targets
  df_T_reduced <- df_T_reduced[match(targets_data$Samples, rownames(df_T_reduced)), ]
  rownames(df_T_reduced) <- targets_data$Samples
  
  # Double-check alignment
  stopifnot(all(rownames(df_T_reduced) == targets_data$Samples))
  
  # Standardize (center and scale)
  scaled_data <- scale(df_T_reduced, center = TRUE, scale = TRUE)
  
  # Handle any NaN/Inf values
  scaled_data[is.na(scaled_data)] <- 0
  scaled_data[is.infinite(scaled_data)] <- 0
  
  # Run PCA
  pca <- prcomp(scaled_data, center = FALSE, scale. = FALSE)  
  
  # Create PCA dataframe with metadata
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Sample = rownames(pca$x)
  )
  pca_df <- cbind(pca_df, targets_data)
  
  # Calculate variance explained percentages
  var_explained <- summary(pca)$importance[2, ] * 100
  pc1_var <- round(var_explained[1], 2)
  pc2_var <- round(var_explained[2], 2)
  
  # Helper function to create individual plots
  create_plot <- function(color_var, color_palette, label) {
    ggplot(pca_df, aes_string(x = "PC1", y = "PC2", color = color_var)) +
      geom_point(size = 2.5, alpha = 0.7) +
      scale_color_manual(values = color_palette) +
      theme_minimal() +
      theme(
        panel.grid = element_line(color = "grey90"),
        plot.title = element_text(size = 14, hjust = 0.5)
      ) +
      labs(
        title = paste("PCA of Top", ntop, "CpGs - PC1 vs PC2 by", label),
        x = paste0("PC1 (", pc1_var, "%)"),
        y = paste0("PC2 (", pc2_var, "%)")
      )
  }
  
  # Batch colors
  batch_colors <- c("red", "blue")[1:length(unique(pca_df$BATCH))]
  names(batch_colors) <- unique(pca_df$BATCH)
  
  p_batch <- create_plot("BATCH", batch_colors, "Batch")
  ggsave(paste0(output_prefix, "_PCA_Batch_PC1_PC2.pdf"), p_batch, width = 8, height = 6)
  
  # Tissue colors
  tissue_levels <- unique(pca_df$Tissue.Cell.type)
  tissue_colors <- rainbow(length(tissue_levels))
  names(tissue_colors) <- tissue_levels
  
  p_tissue <- create_plot("Tissue.Cell.type", tissue_colors, "Tissue")
  ggsave(paste0(output_prefix, "_PCA_Tissue_PC1_PC2.pdf"), p_tissue, width = 8, height = 6)
  
  # Genotype colors 
  genotype_levels <- unique(pca_df$True.Genotype)
  genotype_colors <- c("green", "orange", "purple", "red", "blue", "brown", 
                      "pink", "gray", "cyan", "magenta", "yellow", "black")[1:length(genotype_levels)]
  names(genotype_colors) <- genotype_levels
  
  p_genotype <- create_plot("True.Genotype", genotype_colors, "Genotype")
  ggsave(paste0(output_prefix, "_PCA_Genotype_PC1_PC2.pdf"), p_genotype, width = 8, height = 6)
  
  # Print summary info
  cat("PC1 explains", pc1_var, "% of variance\n")
  cat("PC2 explains", pc2_var, "% of variance\n")
  cat("Total variance explained:", round(sum(var_explained[1:2]), 2), "%\n")
  
  # Save the filtered beta values for the top CpGs
  betas_top <- betas_data[top_cpgs, ]
  write.csv(betas_top, paste0(output_prefix, "_betas_top", ntop, "_cpgs.csv"), row.names = TRUE)
  cat("Saved filtered beta values:", nrow(betas_top), "CpGs x", ncol(betas_top), "samples\n")
  
  write.csv(data.frame(CpG = top_cpgs, Variance = cpg_variances[top_cpgs]), 
            paste0(output_prefix, "_top", ntop, "_cpgs_simple.csv"), row.names = FALSE)
  
  return(list(
    pca_df = pca_df, 
    pca = pca, 
    top_cpgs = top_cpgs,
    betas_top = betas_top  
  ))
}

# Now run the analysis and capture betas_top
cat("Running PCA analysis and extracting top CpGs...\n")

# For ComBat adjusted data  
combat_result <- create_simple_pca(betas_adjusted, matching_samples, 10000, "COMBAT_ADJUSTED")

# Extract betas_top from the results
betas_top_combat <- combat_result$betas_top
cat("Extracted betas_top with dimensions:", nrow(betas_top_combat), "x", ncol(betas_top_combat), "\n")

# ---- Step 7: Analysis Heatmaps + Volcano with Mode and Group Control ----

plot_mode <- "both"  # Options: "both", "volcano", "heatmap"

# Define sample groups for each contrast (for heatmaps)
heatmap_group_map <- list(
  dp10_vs_WT_thymus1 = c("groupdp10_thymus_dp10", "groupdouble_thymus_dp10", "groupdp10_thymus_wt", "groupdouble_thymus_wt", "groupdp16_thymus_wt"),
  dp16_vs_WT_thymus1 = c("groupdp16_thymus_dp16", "groupdouble_thymus_dp16", "groupdp10_thymus_wt", "groupdouble_thymus_wt", "groupdp16_thymus_wt"),
  dp16_dp10_vs_WT_thymus1 = c("groupdouble_thymus_dp10_dp16", "groupdp10_thymus_wt", "groupdouble_thymus_wt", "groupdp16_thymus_wt"),
  dp10_vs_WT_tcells2 = c("groupdp10_tcells_dp10", "groupdouble_tcells_dp10", "groupdp10_tcells_wt", "groupdouble_tcells_wt", "groupdp16_tcells_wt"),
  dp16_vs_WT_tcells2 = c("groupdp16_tcells_dp16", "groupdouble_tcells_dp16", "groupdp10_tcells_wt", "groupdouble_tcells_wt", "groupdp16_tcells_wt"),
  dp16_dp10_vs_WT_tcells2 = c("groupdouble_tcells_dp10_dp16", "groupdp10_tcells_wt", "groupdouble_tcells_wt", "groupdp16_tcells_wt"),
  dp16_vs_WT_thymus3 = c("grouprunx_thymus_dp16", "grouprunx_thymus_wt"),
  runx1_dp16_vs_WT_thymus3 = c("grouprunx_thymus_2copy", "grouprunx_thymus_wt"),
  runx1_vs_WT_thymus3 = c("grouprunx_thymus_1copy", "grouprunx_thymus_wt"),
  runx1_dp16_vs_WT_thymocytes4 = c("grouprunx_thymocytes_2copy", "grouprunx_thymocytes_wt"),
  runx1_vs_WT_thymocytes4 = c("grouprunx_thymocytes_1copy", "grouprunx_thymocytes_wt"),
  dp16_vs_WT_thymocytes4 = c("grouprunx_thymocytes_dp16", "grouprunx_thymocytes_wt"),
  Ts65Dn_vs_WT_thymocytes5 = c("groupTs65Dn", "groupwt_for_Ts65Dn"),
  dp10_vs_WT_thymus7 = c("groupdp10_tcells_dp10_5", "groupdouble_tcells_dp10_5", "groupdp10_tcells_wt_5", "groupdouble_tcells_wt_5", "groupdp16_tcells_wt_5"),
  dp16_vs_WT_thymus7 = c("groupdp16_tcells_dp16_5", "groupdouble_tcells_dp16_5", "groupdp10_tcells_wt_5", "groupdouble_tcells_wt_5", "groupdp16_tcells_wt_5"),
  dp16_dp10_vs_WT_thymus7 = c("groupdouble_tcells_dp10_dp16_5", "groupdp10_tcells_wt_5", "groupdouble_tcells_wt_5", "groupdp16_tcells_wt_5")
)

# Define group mapping for each contrast (experimental vs control groups)
group_labels_map <- list(
  dp10_vs_WT_thymus1 = list(
    experimental = c("groupdp10_thymus_dp10", "groupdouble_thymus_dp10"),
    control = c("groupdp10_thymus_wt", "groupdouble_thymus_wt", "groupdp16_thymus_wt"),
    exp_label = "dp10", ctrl_label = "WT"
  ),
  dp16_vs_WT_thymus1 = list(
    experimental = c("groupdp16_thymus_dp16", "groupdouble_thymus_dp16"),
    control = c("groupdp10_thymus_wt", "groupdouble_thymus_wt", "groupdp16_thymus_wt"),
    exp_label = "dp16", ctrl_label = "WT"
  ),
  dp16_dp10_vs_WT_thymus1 = list(
    experimental = c("groupdouble_thymus_dp10_dp16"),
    control = c("groupdp10_thymus_wt", "groupdouble_thymus_wt", "groupdp16_thymus_wt"),
    exp_label = "dp16_dp10", ctrl_label = "WT"
  ),
  dp10_vs_WT_tcells2 = list(
    experimental = c("groupdp10_tcells_dp10", "groupdouble_tcells_dp10"),
    control = c("groupdp10_tcells_wt", "groupdouble_tcells_wt", "groupdp16_tcells_wt"),
    exp_label = "dp10", ctrl_label = "WT"
  ),
  dp16_vs_WT_tcells2 = list(
    experimental = c("groupdp16_tcells_dp16", "groupdouble_tcells_dp16"),
    control = c("groupdp10_tcells_wt", "groupdouble_tcells_wt", "groupdp16_tcells_wt"),
    exp_label = "dp16", ctrl_label = "WT"
  ),
  dp16_dp10_vs_WT_tcells2 = list(
    experimental = c("groupdouble_tcells_dp10_dp16"),
    control = c("groupdp10_tcells_wt", "groupdouble_tcells_wt", "groupdp16_tcells_wt"),
    exp_label = "dp16_dp10", ctrl_label = "WT"
  ),
  dp16_vs_WT_thymus3 = list(
    experimental = c("grouprunx_thymus_dp16"),
    control = c("grouprunx_thymus_wt"),
    exp_label = "dp16", ctrl_label = "WT"
  ),
  runx1_dp16_vs_WT_thymus3 = list(
    experimental = c("grouprunx_thymus_2copy"),
    control = c("grouprunx_thymus_wt"),
    exp_label = "runx1_dp16", ctrl_label = "WT"
  ),
  runx1_vs_WT_thymus3 = list(
    experimental = c("grouprunx_thymus_1copy"),
    control = c("grouprunx_thymus_wt"),
    exp_label = "runx1", ctrl_label = "WT"
  ),
  runx1_dp16_vs_WT_thymocytes4 = list(
    experimental = c("grouprunx_thymocytes_2copy"),
    control = c("grouprunx_thymocytes_wt"),
    exp_label = "runx1_dp16", ctrl_label = "WT"
  ),
  runx1_vs_WT_thymocytes4 = list(
    experimental = c("grouprunx_thymocytes_1copy"),
    control = c("grouprunx_thymocytes_wt"),
    exp_label = "runx1", ctrl_label = "WT"
  ),
  dp16_vs_WT_thymocytes4 = list(
    experimental = c("grouprunx_thymocytes_dp16"),
    control = c("grouprunx_thymocytes_wt"),
    exp_label = "dp16", ctrl_label = "WT"
  ),
  Ts65Dn_vs_WT_thymocytes5 = list(
    experimental = c("groupTs65Dn"),
    control = c("groupwt_for_Ts65Dn"),
    exp_label = "Ts65Dn", ctrl_label = "WT"
  ),
  dp10_vs_WT_thymus7 = list(
    experimental = c("groupdp10_tcells_dp10_5", "groupdouble_tcells_dp10_5"),
    control = c("groupdp10_tcells_wt_5", "groupdouble_tcells_wt_5", "groupdp16_tcells_wt_5"),
    exp_label = "dp10", ctrl_label = "WT"
  ),
  dp16_vs_WT_thymus7 = list(
    experimental = c("groupdp16_tcells_dp16_5", "groupdouble_tcells_dp16_5"),
    control = c("groupdp10_tcells_wt_5", "groupdouble_tcells_wt_5", "groupdp16_tcells_wt_5"),
    exp_label = "dp16", ctrl_label = "WT"
  ),
  dp16_dp10_vs_WT_thymus7 = list(
    experimental = c("groupdouble_tcells_dp10_dp16_5"),
    control = c("groupdp10_tcells_wt_5", "groupdouble_tcells_wt_5", "groupdp16_tcells_wt_5"),
    exp_label = "dp16_dp10", ctrl_label = "WT"
  )
)

betas <- betas_top_combat
targets <- matching_samples

contrast_names <- c(
  "dp10_vs_WT_thymus1", "dp16_vs_WT_thymus1", "dp16_dp10_vs_WT_thymus1",
  "dp10_vs_WT_tcells2", "dp16_vs_WT_tcells2", "dp16_dp10_vs_WT_tcells2",
  "dp16_vs_WT_thymus3", "runx1_dp16_vs_WT_thymus3", "runx1_vs_WT_thymus3",
  "runx1_dp16_vs_WT_thymocytes4", "runx1_vs_WT_thymocytes4", "dp16_vs_WT_thymocytes4",
  "Ts65Dn_vs_WT_thymocytes5", "dp10_vs_WT_thymus7", "dp16_vs_WT_thymus7", "dp16_dp10_vs_WT_thymus7"
)

contrasts_list <- c(
  "(groupdp10_thymus_dp10 + groupdouble_thymus_dp10)/2 - (groupdp10_thymus_wt + groupdouble_thymus_wt + groupdp16_thymus_wt)/3",
  "(groupdp16_thymus_dp16 + groupdouble_thymus_dp16)/2 - (groupdp10_thymus_wt + groupdouble_thymus_wt + groupdp16_thymus_wt)/3",
  "groupdouble_thymus_dp10_dp16 - (groupdp10_thymus_wt + groupdouble_thymus_wt + groupdp16_thymus_wt)/3",
  "(groupdp10_tcells_dp10 + groupdouble_tcells_dp10)/2 - (groupdp10_tcells_wt + groupdouble_tcells_wt + groupdp16_tcells_wt)/3",
  "(groupdp16_tcells_dp16 + groupdouble_tcells_dp16)/2 - (groupdp10_tcells_wt + groupdouble_tcells_wt + groupdp16_tcells_wt)/3",
  "groupdouble_tcells_dp10_dp16 - (groupdp10_tcells_wt + groupdouble_tcells_wt + groupdp16_tcells_wt)/3",
  "grouprunx_thymus_dp16 - grouprunx_thymus_wt",
  "grouprunx_thymus_2copy - grouprunx_thymus_wt",
  "grouprunx_thymus_1copy - grouprunx_thymus_wt",
  "grouprunx_thymocytes_2copy - grouprunx_thymocytes_wt",
  "grouprunx_thymocytes_1copy - grouprunx_thymocytes_wt",
  "grouprunx_thymocytes_dp16 - grouprunx_thymocytes_wt",
  "groupTs65Dn - groupwt_for_Ts65Dn",
  "(groupdp10_tcells_dp10_5 + groupdouble_tcells_dp10_5)/2 - (groupdp10_tcells_wt_5 + groupdouble_tcells_wt_5 + groupdp16_tcells_wt_5)/3",
  "(groupdp16_tcells_dp16_5 + groupdouble_tcells_dp16_5)/2 - (groupdp10_tcells_wt_5 + groupdouble_tcells_wt_5 + groupdp16_tcells_wt_5)/3",
  "groupdouble_tcells_dp10_dp16_5 - (groupdp10_tcells_wt_5 + groupdouble_tcells_wt_5 + groupdp16_tcells_wt_5)/3"
)

mod <- model.matrix(~ 0 + group, data = targets)
fit <- lmFit(betas, mod)

for (i in seq_along(contrasts_list)) {
  contrast_name <- contrast_names[i]
  contrast_formula <- contrasts_list[i]
  
  cat("\n--- Running contrast", i, ":", contrast_name, "---\n")
  
  cont_matrix <- makeContrasts(contrasts = contrast_formula, levels = mod)
  fit_cont <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit_cont)
  results <- topTable(fit2, num = Inf, coef = 1, sort.by = 'none')
  write.csv(results, paste0("DM_", contrast_name, ".csv"), row.names = TRUE)
  
  sigUp <- results[which(results$P.Value < 1e-4 & results$logFC > 0.1), ]
  sigDown <- results[which(results$P.Value < 1e-4 & results$logFC < -0.1), ]
  notSig <- results[which(!(results$P.Value < 1e-4 & abs(results$logFC) > 0.1)), ]
  DMPs <- results[which(results$adj.P.Val < 1e-4 & abs(results$logFC) > 0.1), ]
  
  # Volcano plot
  if (plot_mode %in% c("both", "volcano")) {
    pdf(paste0("volcano_", contrast_name, ".pdf"))
    if (nrow(sigUp) > 0) sigUp$col <- "indianred"
    if (nrow(sigDown) > 0) sigDown$col <- "deepskyblue3"
    if (nrow(notSig) > 0) notSig$col <- "ivory3"
    
    xx <- rbind(sigUp, sigDown, notSig)
    xx <- xx[order(xx$P.Value), ]
    
    plot(xx$logFC, -log10(xx$P.Value),
         xlim = c(-0.5, 0.5),
         col = xx$col, pch = 16,
         ylab = "-log10(pval)",
         xlab = paste("Effect Size (", contrast_name, ")"),
         main = paste("Volcano Plot -", contrast_name))
    
    abline(h = -log10(1e-4), lty = 3)
    abline(v = c(-0.1, 0.1), lty = 2, col = "gray50")
    dev.off()
  }
  
  # --- Z-score Heatmap Block ---
  if (plot_mode %in% c("both", "heatmap") && nrow(DMPs) > 0) {
    common_cpgs <- intersect(rownames(DMPs), rownames(betas))
    if (length(common_cpgs) == 0) next
    
    DMPs.beta <- betas[common_cpgs, , drop = FALSE]
    DMPs.beta <- DMPs.beta[complete.cases(DMPs.beta), , drop = FALSE]
    
    if (contrast_name %in% names(group_labels_map)) {
      group_info <- group_labels_map[[contrast_name]]
      
      exp_groups <- gsub("^group", "", group_info$experimental)
      ctrl_groups <- gsub("^group", "", group_info$control)
      
      exp_samples <- targets$Samples[targets$group %in% exp_groups]
      ctrl_samples <- targets$Samples[targets$group %in% ctrl_groups]
      
      selected_samples <- c(exp_samples, ctrl_samples)
      selected_samples <- intersect(selected_samples, colnames(DMPs.beta))
      
      if (length(selected_samples) > 0) {
        DMPs.beta <- DMPs.beta[, selected_samples, drop = FALSE]
        
        col_annotation <- data.frame(
          Group = c(rep(group_info$exp_label, length(exp_samples)), 
                    rep(group_info$ctrl_label, length(ctrl_samples))),
          row.names = selected_samples
        )
        
        group_colors <- c("red", "blue")
        names(group_colors) <- c(group_info$exp_label, group_info$ctrl_label)
        annotation_colors <- list(Group = group_colors)
        
        if (nrow(DMPs.beta) > 0 && ncol(DMPs.beta) > 0) {
          zscore_matrix <- t(scale(t(DMPs.beta)))
          
          row_height <- 0.02         # inches per CpG row (try 0.25 for better vertical spacing)
          extra_padding <- 2         # inches for annotations, legends, labels
          n_cpgs <- nrow(DMPs.beta)
          pdf_height <- n_cpgs * row_height + extra_padding
          
          pdf(paste0("heatmap_zscore_", contrast_name, ".pdf"), width = 12, height = pdf_height)
          pheatmap(zscore_matrix,
                  cluster_rows = TRUE,
                  cluster_cols = FALSE,
                  show_rownames = FALSE,
                  show_colnames = FALSE,
                  scale = "none",
                  treeheight_col = 0,
                  treeheight_row = 0,
                  annotation_col = col_annotation,
                  annotation_colors = annotation_colors,
                  annotation_names_col = TRUE,
                  main = paste("Z-Score Heatmap of Significant DMPs –", contrast_name),
                  fontsize = 10,
                  border_color = NA)  # <--- this line removes the cell borders
                    dev.off()
        }
      }
    } else {
      if (contrast_name %in% names(heatmap_group_map)) {
        selected_groups <- heatmap_group_map[[contrast_name]]
        group_names <- gsub("^group", "", selected_groups)
        sample_mask <- targets$group %in% group_names
        DMPs.beta <- DMPs.beta[, sample_mask, drop = FALSE]
      }
      
      if (nrow(DMPs.beta) > 0 && ncol(DMPs.beta) > 0) {
        zscore_matrix <- t(scale(t(DMPs.beta)))
        
        # Fixed row height
        row_height <- 0.05  # inches per CpG row
        n_cpgs <- nrow(zscore_matrix)
        pdf_height <- n_cpgs * row_height
        
        pdf(paste0("heatmap_zscore_", contrast_name, ".pdf"), width = 12, height = pdf_height)
        pheatmap(zscore_matrix,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                show_rownames = FALSE,
                scale = "none",
                treeheight_col = 1,
                treeheight_row = 0,
                main = paste("Z-Score Heatmap of Significant DMPs –", contrast_name))
        dev.off()
      }
    }
  }

  cat("Hypermethylated (p < 1e-4):", nrow(sigUp), "\n")
  cat("Hypomethylated (p < 1e-4):", nrow(sigDown), "\n")
  cat("Not significant:", nrow(notSig), "\n")
  cat("Significant DMPs (FDR < 1e-4 & |logFC| > 0.1):", nrow(DMPs), "\n")
}