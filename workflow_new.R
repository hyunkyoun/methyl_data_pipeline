# Methylation Data Analysis Pipeline (Clean)
# IDAT -> QC -> ComBat -> PCA -> DMP Analysis -> FASTA Extraction for AME

# ---- Setup 
.libPaths("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/packages")
setwd("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline")

required_libs <- c("sesame", "minfi", "IlluminaMouseMethylationmanifest", 
                   "IlluminaMouseMethylationanno.12.v1.mm10", "limma", "sva", 
                   "ggsci", "wateRmelon", "snow", "BiocParallel", "parallel", 
                   "RColorBrewer", "readxl", "openxlsx", "ggplot2", "gridExtra", "pheatmap",
                   "GenomicRanges", "GenomicFeatures", "BSgenome.Mmusculus.UCSC.mm10", 
                   "Biostrings")  # Added Biostrings for writeXStringSet
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

# ---- Function for extracting FASTA sequences from DMPs ----

extract_and_write_FASTA <- function(DMPs, contrast_name, annoMouse, genome) {
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

  # Extract sequences
  GOM.fasta <- getSeq(genome, GOM.gr)
  LOM.fasta <- getSeq(genome, LOM.gr)

  names(GOM.fasta) <- paste0(seqnames(GOM.gr), ":", start(GOM.gr), "-", end(GOM.gr))
  names(LOM.fasta) <- paste0(seqnames(LOM.gr), ":", start(LOM.gr), "-", end(LOM.gr))

  # Save outputs
  dir.create("FASTA_sequences", showWarnings = FALSE)
  writeXStringSet(GOM.fasta, file.path("FASTA_sequences", paste0(contrast_name, "_GOM.fasta")))
  writeXStringSet(LOM.fasta, file.path("FASTA_sequences", paste0(contrast_name, "_LOM.fasta")))
  write.csv(GOM, file.path("FASTA_sequences", paste0(contrast_name, "_GOM_details.csv")), row.names = FALSE)
  write.csv(LOM, file.path("FASTA_sequences", paste0(contrast_name, "_LOM_details.csv")), row.names = FALSE)

  cat("Hypermethylated (GOM):", nrow(GOM), "\n")
  cat("Hypomethylated (LOM):", nrow(LOM), "\n")
  cat("Total FASTA sequences written:", length(GOM.fasta) + length(LOM.fasta), "\n")

  return(list(GOM = nrow(GOM), LOM = nrow(LOM),
              written = length(GOM.fasta) + length(LOM.fasta)))
}

# ---- Run DMP analysis and FASTA extraction (CORRECTED VERSION) ----

library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10

results_list <- list()
summary_list <- list()

for (contrast_name in names(group_labels_map)) {
  cat("\n=== Processing contrast:", contrast_name, "===\n")

  groups <- group_labels_map[[contrast_name]]
  
  # CORRECTED: Match samples based on group column, not sample names
  # Find samples that belong to experimental groups
  exp_samples <- matching_samples$Samples[matching_samples$group %in% gsub("^group", "", groups$experimental)]
  ctrl_samples <- matching_samples$Samples[matching_samples$group %in% gsub("^group", "", groups$control)]
  
  # Debug output to see what's happening
  cat("Experimental groups:", paste(gsub("^group", "", groups$experimental), collapse = ", "), "\n")
  cat("Control groups:", paste(gsub("^group", "", groups$control), collapse = ", "), "\n")
  cat("Available groups in data:", paste(unique(matching_samples$group), collapse = ", "), "\n")
  cat("Found experimental samples:", length(exp_samples), "\n")
  cat("Found control samples:", length(ctrl_samples), "\n")
  
  if (length(exp_samples) < 2 | length(ctrl_samples) < 2) {
    cat("Skipping", contrast_name, "(not enough samples: exp=", length(exp_samples), ", ctrl=", length(ctrl_samples), ")\n")
    
    # Additional debugging - show first few sample IDs and groups
    cat("Sample group examples:\n")
    print(head(matching_samples[, c("Samples", "group")], 10))
    next
  }

  # Subset betas
  betas_sub <- betas_adjusted[, c(exp_samples, ctrl_samples), drop = FALSE]

  # Design matrix
  design <- model.matrix(~ 0 + factor(c(rep("exp", length(exp_samples)),
                                        rep("ctrl", length(ctrl_samples)))))
  colnames(design) <- c("ctrl", "exp")

  fit <- lmFit(betas_sub, design)
  contrast.matrix <- makeContrasts(exp - ctrl, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  tt <- topTable(fit2, number = Inf, adjust = "BH")
  tt$Probe <- rownames(tt)
  tt$diff <- tt$logFC

  # Save results
  results_list[[contrast_name]] <- tt
  write.csv(tt, paste0("DM_", contrast_name, "_custom.csv"), row.names = TRUE)

  # Extract FASTA - only for significant DMPs
  significant_tt <- tt[tt$adj.P.Val < 0.05 & abs(tt$diff) > 0.1, ]
  if (nrow(significant_tt) > 0) {
    res <- extract_and_write_FASTA(significant_tt, contrast_name, annoMouse, genome)
    summary_list[[contrast_name]] <- res
  } else {
    cat("No significant DMPs found for", contrast_name, "\n")
    summary_list[[contrast_name]] <- list(GOM = 0, LOM = 0, written = 0)
  }
}

# ---- Print summary ----
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
  cat("2. Group matching issues (check sample group assignments)\n")
  cat("3. Annotation matching problems (check annoMouse object)\n")
} else {
  cat("FASTA files have been created in 'FASTA_sequences/'.\n")
}

# ---- Additional debugging function ----
debug_sample_matching <- function() {
  cat("\n=== SAMPLE MATCHING DEBUG INFO ===\n")
  cat("Total samples in matching_samples:", nrow(matching_samples), "\n")
  cat("Unique groups in data:", length(unique(matching_samples$group)), "\n")
  cat("Group distribution:\n")
  print(table(matching_samples$group))
  
  cat("\nFirst 10 samples and their groups:\n")
  print(head(matching_samples[, c("Samples", "group")], 10))
  
  cat("\nExpected groups from group_labels_map (first contrast):\n")
  first_contrast <- group_labels_map[[1]]
  cat("Experimental:", paste(gsub("^group", "", first_contrast$experimental), collapse = ", "), "\n")
  cat("Control:", paste(gsub("^group", "", first_contrast$control), collapse = ", "), "\n")
}

# Run debugging
debug_sample_matching()
