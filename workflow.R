# # Methylation Data Analysis Pipeline
# # IDAT -> QC -> BMIQ (per run) -> ComBat-MET (across all runs) -> PCA

# # ==============================================================================
# # SETUP AND INITIALIZATION
# # ==============================================================================

# # Set library paths and working directory
# .libPaths("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/packages")
# setwd("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline")

# # Load required libraries
# print("Loading libraries...")
# required_libs <- c("sesame", "minfi", "IlluminaMouseMethylationmanifest", 
#                    "IlluminaMouseMethylationanno.12.v1.mm10", "limma", "sva", 
#                    "ggsci", "wateRmelon", "snow", "BiocParallel", "parallel", 
#                    "RColorBrewer", "readxl", "openxlsx", "ggplot2")

# lapply(required_libs, library, character.only = TRUE)

# # Install and load ComBat-MET if not available
# if (!requireNamespace("ComBatMet", quietly = TRUE)) {
#   if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#   devtools::install_github("JmWangBio/ComBatMet")
# }
# library(ComBatMet)

# # Source BMIQ function
# print("Sourcing DoBMIQ function...")
# source("bmiq/DoBMIQ.R")

# # Load annotation data
# print("Loading annotation data...")
# annoMouse <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)

# # ==============================================================================
# # STEP 1: READ IDAT FILES AND GENERATE DETECTION P-VALUES
# # ==============================================================================

# print("=== STEP 1: Reading IDAT files ===")
# idir <- "idat"
# cores <- 4

# # Read raw data and calculate detection p-values
# betas <- openSesame(idir, BPPARAM = MulticoreParam(workers = cores), collapseToPfx = TRUE)
# detP <- openSesame(idir, func = pOOBAH, return.pval = TRUE, BPPARAM = MulticoreParam(workers = cores))
# detP <- betasCollapseToPfx(detP)

# # Load and prepare sample metadata
# targets <- read.csv("data/master_samplesheet.csv", header = TRUE)
# targets$Basename <- paste0(targets$BEADCHIP, "_", targets$POSITION)
# targets$Samples <- paste0(targets$tb, "_", targets$RUN)
# targets$BEADCHIP <- as.character(as.numeric(targets$BEADCHIP))
# targets <- targets[!duplicated(targets$Basename), ]

# # Align metadata with beta matrix
# common_samples <- intersect(targets$Basename, colnames(betas))
# targets <- targets[targets$Basename %in% common_samples, ]
# betas <- betas[, match(targets$Basename, colnames(betas))]
# detP <- detP[, match(targets$Basename, colnames(detP))]

# # Verify sample matching
# if (!identical(targets$Basename, colnames(betas))) {
#   stop("Sample matching failed!")
# }

# # Generate initial QC plots
# print("Generating initial QC plots...")
# pal <- pal_simpsons()(12)
# pdf("meanDP_initial.pdf")
# barplot(colMeans(detP), col = pal[factor(targets$Samples)], las = 2, 
#         cex.names = 0.4, ylab = "Mean detection p-values")
# abline(h = 0.05, col = "red")
# legend("topleft", legend = levels(factor(targets$Samples)), fill = pal, bg = "white")
# dev.off()

# pdf("density_initial.pdf")
# densityPlot(betas, sampGroups = targets$Samples, main = "Beta Density (Initial)", 
#             legend = FALSE, pal = pal)
# legend("top", legend = levels(factor(targets$Samples)), text.col = pal)
# dev.off()

# # ==============================================================================
# # STEP 2: QC FILTERING BY RUN
# # ==============================================================================

# print("=== STEP 2: QC filtering by run ===")

# qc_filter_run <- function(betas_run, detP_run, targets_run, run_id) {
#   print(paste("Processing RUN", run_id))
  
#   # Filter samples with poor detection p-values
#   good_samples <- which(colMeans(detP_run) < 0.05)
#   if (length(good_samples) == 0) {
#     warning(paste("No good samples in RUN", run_id))
#     return(NULL)
#   }
  
#   betas_run <- betas_run[, good_samples, drop = FALSE]
#   detP_run <- detP_run[, good_samples, drop = FALSE]
#   targets_run <- targets_run[good_samples, ]
  
#   # Keep probes detected in all remaining samples
#   keep_probes <- rowSums(detP_run < 0.05) == ncol(betas_run)
#   betas_run <- betas_run[keep_probes, ]
#   detP_run <- detP_run[keep_probes, ]
  
#   # Remove probes with NA/NaN values
#   complete_idx <- complete.cases(betas_run)
#   betas_run <- betas_run[complete_idx, ]
#   detP_run <- detP_run[complete_idx, ]
  
#   # Keep only CpG probes
#   cpg_idx <- grep("^cg", rownames(betas_run))
#   betas_run <- betas_run[cpg_idx, ]
#   detP_run <- detP_run[rownames(betas_run), ]
  
#   # Remove sex chromosome probes
#   sex_cpgs <- annoMouse$Name[annoMouse$chr %in% c("chrX", "chrY")]
#   keep_sex <- !rownames(betas_run) %in% sex_cpgs
#   betas_run <- betas_run[keep_sex, ]
#   detP_run <- detP_run[rownames(betas_run), ]
  
#   print(paste("RUN", run_id, "after QC:", nrow(betas_run), "probes x", ncol(betas_run), "samples"))
  
#   return(list(betas = betas_run, detP = detP_run, targets = targets_run))
# }

# # Apply QC filtering to each run
# unique_runs <- sort(unique(targets$RUN))
# filtered_data <- list()

# for (run_id in unique_runs) {
#   idx <- targets$RUN == run_id
#   result <- qc_filter_run(betas[, idx], detP[, idx], targets[idx, ], run_id)
#   if (!is.null(result)) {
#     filtered_data[[paste0("RUN", run_id)]] <- result
#   }
# }

# # ==============================================================================
# # STEP 3: BMIQ NORMALIZATION (PER RUN)
# # ==============================================================================

# print("=== STEP 3: BMIQ normalization per run ===")

# run_bmiq_normalization <- function(betas_run, run_name) {
#   print(paste("Running BMIQ on", run_name))
  
#   temp_file <- paste0("temp_betas_", run_name, ".csv")
#   write.csv(data.frame(CpG_ID = rownames(betas_run), betas_run), 
#             file = temp_file, row.names = FALSE)
  
#   tryCatch({
#     DoBMIQ(probes_path = "data/probesample.xlsx", beta_path = temp_file)
#     load("bmiq/results/bmiq.Rd")
    
#     # Save BMIQ results
#     write.csv(bmiq.m, paste0("betas_bmiq_", run_name, ".csv"))
#     print(paste("BMIQ completed for", run_name, "- Probes:", nrow(bmiq.m), "Samples:", ncol(bmiq.m)))
    
#     return(bmiq.m)
#   }, error = function(e) {
#     print(paste("BMIQ failed for", run_name, ":", e$message))
#     print("Using raw data for this run...")
#     return(betas_run)
#   }, finally = {
#     if (file.exists(temp_file)) file.remove(temp_file)
#   })
# }

# # Apply BMIQ to each run
# betas_bmiq_by_run <- list()
# for (run_name in names(filtered_data)) {
#   betas_bmiq_by_run[[run_name]] <- run_bmiq_normalization(
#     filtered_data[[run_name]]$betas, run_name
#   )
# }

# # ==============================================================================
# # STEP 4: COMBAT NORMALIZATION USING COMBAT-MET
# # ==============================================================================

# print("=== STEP 4: Logit-M-value ComBat Normalization ===")

# # --- Step 4.1: Combine runs as before ---
# common_cpgs_12 <- Reduce(intersect, lapply(betas_bmiq_by_run[c("RUN1", "RUN2")], rownames))
# betas_RUN1_2 <- cbind(
#   betas_bmiq_by_run$RUN1[common_cpgs_12, , drop = FALSE],
#   betas_bmiq_by_run$RUN2[common_cpgs_12, , drop = FALSE]
# )

# common_cpgs_345 <- Reduce(intersect, lapply(betas_bmiq_by_run[c("RUN3", "RUN4", "RUN5")], rownames))
# betas_RUN3_4_5 <- cbind(
#   betas_bmiq_by_run$RUN3[common_cpgs_345, , drop = FALSE],
#   betas_bmiq_by_run$RUN4[common_cpgs_345, , drop = FALSE],
#   betas_bmiq_by_run$RUN5[common_cpgs_345, , drop = FALSE]
# )

# common_cpgs_all <- intersect(rownames(betas_RUN1_2), rownames(betas_RUN3_4_5))
# beta_combined <- cbind(
#   betas_RUN1_2[common_cpgs_all, , drop = FALSE],
#   betas_RUN3_4_5[common_cpgs_all, , drop = FALSE]
# )

# # --- Step 4.2: Match sample metadata ---
# samplesheet <- read.csv("data/master_samplesheet.csv", stringsAsFactors = FALSE)
# samplesheet$Basename <- paste0("X", samplesheet$BEADCHIP, "_", samplesheet$POSITION)
# samplesheet$Samples <- paste0(samplesheet$tb, "_", samplesheet$RUN)

# sample_names <- colnames(beta_combined)
# matching_samples <- samplesheet[samplesheet$Basename %in% sample_names, ]
# matching_samples <- matching_samples[match(sample_names, matching_samples$Basename), ]

# if (!identical(sample_names, matching_samples$Basename)) {
#   stop("Sample matching failed in batch correction step!")
# }

# # --- Step 4.3: Define batch and group ---
# batch <- factor(ifelse(matching_samples$RUN %in% c(1, 2), "batch1", "batch2"))
# group <- factor(matching_samples$cell.type)

# # --- Step 4.4: Mild filtering (ensure safe logit transform) ---
# epsilon <- 1e-5
# beta_combined <- beta_combined[
#   apply(beta_combined, 1, function(x) all(x > epsilon & x < (1 - epsilon))),
# ]

# cat("Remaining probes after logit-safe filtering:", nrow(beta_combined), "\n")

# # --- Step 4.5: Convert to logit M-values ---
# M <- log2(beta_combined / (1 - beta_combined))

# # --- Step 4.6: Build design matrix ---
# mod <- model.matrix(~ group)

# # --- Step 4.7: Apply ComBat ---
# library(sva)

# M_corrected <- ComBat(
#   dat = beta_combined,
#   batch = batch,
#   mod = mod,
#   ref.batch = "batch1",   # <- anchor to batch1
#   par.prior = TRUE,       # enable empirical Bayes
#   prior.plots = FALSE
# )

# # --- Step 4.8: Back-transform to β-values ---
# beta_corrected <- 2^M_corrected / (1 + 2^M_corrected)

# # --- Optional: Save corrected β-values ---
# dir.create("output", showWarnings = FALSE)
# write.csv(beta_corrected, "output/betas_corrected_logit_combat.csv")

# print("✅ Logit-M-value ComBat normalization complete.")


# # ==============================================================================
# # STEP 5: PCA ANALYSIS
# # ==============================================================================

# print("=== STEP 5: PCA analysis and visualization ===")

# # Load sample sheet
# samplesheet <- read.csv("data/master_samplesheet.csv", stringsAsFactors = FALSE)

# # Construct Basename (must match colnames in beta matrix)
# samplesheet$Basename <- paste0("X", samplesheet$BEADCHIP, "_", samplesheet$POSITION)

# # Ensure your beta matrix uses Basename-style columns
# # colnames(beta_corrected) should be like "X205243950045_R06C02"
# if (!all(colnames(beta_corrected) %in% samplesheet$Basename)) {
#   stop("Some column names in beta_corrected not found in sample sheet Basename.")
# }

# # Match sample names using Basename → tb
# new_names <- samplesheet$tb[match(colnames(beta_corrected), samplesheet$Basename)]

# # Apply new names
# colnames(beta_corrected) <- new_names

# # Use betas_combat as the final adjusted matrix
# betas_final <- beta_corrected

# # Clamp values to avoid Inf/NaN in log2 transformation (important for MDS)
# betas_final[betas_final <= 0.001] <- 0.001
# betas_final[betas_final >= 0.999] <- 0.999

# # Subset for top 10K most variable probes
# probe_vars <- apply(betas_final, 1, var)
# top_probes <- names(sort(probe_vars, decreasing = TRUE))[1:10000]
# betas_final_top <- betas_final[top_probes, ]

# # Get sample metadata again and align
# samplesheet <- read.csv("data/master_samplesheet.csv", header = TRUE)
# samplesheet$Basename <- paste0(samplesheet$BEADCHIP, "_", samplesheet$POSITION)
# samplesheet$Samples <- paste0(samplesheet$tb, "_", samplesheet$RUN)
# samplesheet$Samples_unique <- make.unique(as.character(samplesheet$Samples))
# rownames(samplesheet) <- samplesheet$Samples_unique

# # Match metadata to columns
# sample_match <- samplesheet[colnames(betas_final), ]

# # Define color palettes
# celltype_colors <- c("t-cells-DNA/RNA" = "#E31A1C", "THYMOCYTES" = "#1F78B4", "UNFRACTIONATED" = "#33A02C")
# batch_colors <- c("Batch1" = "red", "Batch2" = "blue")

# # Build color vectors
# sample_celltype_colors <- celltype_colors[as.character(sample_match$cell.type)]
# sample_match$Batch <- ifelse(sample_match$RUN <= 2, "Batch1", "Batch2")
# sample_batch_colors <- batch_colors[sample_match$Batch]

# # === Begin PCA Visualizations ===
# print("Generating PCA plots...")

# pdf("pca_simple_combat.pdf", width = 16, height = 8)
# par(mfrow = c(2, 3))

# # PCA by Cell Type (Top 10K probes)
# plotMDS(betas_final_top, top = 10000, gene.selection = "common", col = sample_celltype_colors,
#         main = "Top 10K CpGs - Cell Type", pch = 19, cex = 1.2)
# legend("topright", legend = names(celltype_colors), col = celltype_colors, pch = 19, cex = 0.8)

# # PCA by Cell Type (All probes)
# plotMDS(betas_final, gene.selection = "common", col = sample_celltype_colors,
#         main = "All CpGs - Cell Type", pch = 19, cex = 1.2)
# legend("topright", legend = names(celltype_colors), col = celltype_colors, pch = 19, cex = 0.8)

# # PCA by Batch (Top 10K probes)
# plotMDS(betas_final_top, top = 10000, gene.selection = "common", col = sample_batch_colors,
#         main = "Top 10K CpGs - Batch (should be mixed)", pch = 19, cex = 1.2)
# legend("topright", legend = names(batch_colors), col = batch_colors, pch = 19)

# # PCA by Batch (All CpGs)
# plotMDS(betas_final, gene.selection = "common", col = sample_batch_colors,
#         main = "All CpGs - Batch (should be mixed)", pch = 19, cex = 1.2)
# legend("topright", legend = names(batch_colors), col = batch_colors, pch = 19)

# # PCA by Experiment (if present)
# if ("Experiment" %in% colnames(sample_match)) {
#   exp_colors <- rainbow(length(unique(sample_match$Experiment)))
#   names(exp_colors) <- unique(sample_match$Experiment)
#   sample_exp_colors <- exp_colors[as.character(sample_match$Experiment)]
  
#   plotMDS(betas_final, gene.selection = "common", col = sample_exp_colors,
#           main = "All CpGs - Experiment", pch = 19, cex = 1.2)
#   legend("topright", legend = names(exp_colors), col = exp_colors, pch = 19, cex = 0.7)
# }

# # Sample Label View
# mds <- plotMDS(betas_final, gene.selection = "common", dim = c(1, 2), plot = FALSE)
# plot(mds$x, mds$y,
#      xlab = paste("Dimension 1 (", round(mds$var.explained[1], 1), "%)", sep = ""),
#      ylab = paste("Dimension 2 (", round(mds$var.explained[2], 1), "%)", sep = ""),
#      main = "Sample Labels", pch = 19, col = "black", cex = 0.8)
# text(mds$x, mds$y, labels = colnames(betas_final), pos = 3, cex = 0.5)

# dev.off()

# print("✅ PCA plots saved to: pca_simple_combat.pdf")

# # ==============================================================================
# # STEP 5 ALT: PCA ANALYSIS BEFORE COMBAT NORMALIZATION
# # ==============================================================================

# print("=== STEP 5 ALT: PCA analysis BEFORE batch normalization ===")

# # Load sample sheet
# samplesheet <- read.csv("data/master_samplesheet.csv", stringsAsFactors = FALSE)

# # Construct Basename (should match pre-ComBat colnames)
# samplesheet$Basename <- paste0("X", samplesheet$BEADCHIP, "_", samplesheet$POSITION)
# samplesheet$Samples <- samplesheet$tb  # readable names

# # Ensure no duplicated Basename rows
# samplesheet <- samplesheet[!duplicated(samplesheet$Basename), ]

# # Rename beta matrix columns using tb names (from Basename)
# sample_names <- colnames(beta_combined)
# rename_map <- samplesheet$Samples[match(sample_names, samplesheet$Basename)]

# # Check that renaming works
# if (any(is.na(rename_map))) {
#   stop("❌ Sample renaming failed: some Basename entries in beta_combined not found in sample sheet.")
# }

# # Apply renaming
# colnames(beta_combined) <- rename_map

# # Assign as final matrix for PCA
# betas_uncorrected <- beta_combined

# # Clamp extreme beta values to avoid issues with log2 distance
# betas_uncorrected[betas_uncorrected <= 0.001] <- 0.001
# betas_uncorrected[betas_uncorrected >= 0.999] <- 0.999

# # Subset to top 10,000 most variable probes
# probe_vars <- apply(betas_uncorrected, 1, var)
# top_probes <- names(sort(probe_vars, decreasing = TRUE))[1:10000]
# betas_uncorrected_top <- betas_uncorrected[top_probes, ]

# # Match metadata to renamed columns
# sample_match <- samplesheet[match(colnames(betas_uncorrected), samplesheet$Samples), ]

# # Build color vectors for cell type and batch
# celltype_colors <- c("t-cells-DNA/RNA" = "#E31A1C", "THYMOCYTES" = "#1F78B4", "UNFRACTIONATED" = "#33A02C")
# sample_celltype_colors <- celltype_colors[as.character(sample_match$cell.type)]

# sample_match$Batch <- ifelse(sample_match$RUN <= 2, "Batch1", "Batch2")
# batch_colors <- c("Batch1" = "red", "Batch2" = "blue")
# sample_batch_colors <- batch_colors[sample_match$Batch]

# # === PCA PLOTTING ===
# print("Generating PCA plots BEFORE batch correction...")

# pdf("pca_before_combat.pdf", width = 16, height = 8)
# par(mfrow = c(2, 3))

# # PCA by Cell Type (Top 10K CpGs)
# plotMDS(betas_uncorrected_top, top = 10000, gene.selection = "common", col = sample_celltype_colors,
#         main = "Top 10K CpGs (Pre-ComBat) - Cell Type", pch = 19, cex = 1.2)
# legend("topright", legend = names(celltype_colors), col = celltype_colors, pch = 19, cex = 0.8)

# # PCA by Cell Type (All CpGs)
# plotMDS(betas_uncorrected, gene.selection = "common", col = sample_celltype_colors,
#         main = "All CpGs (Pre-ComBat) - Cell Type", pch = 19, cex = 1.2)
# legend("topright", legend = names(celltype_colors), col = celltype_colors, pch = 19, cex = 0.8)

# # PCA by Batch (Top 10K CpGs)
# plotMDS(betas_uncorrected_top, top = 10000, gene.selection = "common", col = sample_batch_colors,
#         main = "Top 10K CpGs (Pre-ComBat) - Batch", pch = 19, cex = 1.2)
# legend("topright", legend = names(batch_colors), col = batch_colors, pch = 19)

# # PCA by Batch (All CpGs)
# plotMDS(betas_uncorrected, gene.selection = "common", col = sample_batch_colors,
#         main = "All CpGs (Pre-ComBat) - Batch", pch = 19, cex = 1.2)
# legend("topright", legend = names(batch_colors), col = batch_colors, pch = 19)

# # PCA by Experiment (if present)
# if ("Experiment" %in% colnames(sample_match)) {
#   exp_colors <- rainbow(length(unique(sample_match$Experiment)))
#   names(exp_colors) <- unique(sample_match$Experiment)
#   sample_exp_colors <- exp_colors[as.character(sample_match$Experiment)]
  
#   plotMDS(betas_uncorrected, gene.selection = "common", col = sample_exp_colors,
#           main = "All CpGs (Pre-ComBat) - Experiment", pch = 19, cex = 1.2)
#   legend("topright", legend = names(exp_colors), col = exp_colors, pch = 19, cex = 0.7)
# }

# # PCA by Sample Label
# mds <- plotMDS(betas_uncorrected, gene.selection = "common", dim = c(1, 2), plot = FALSE)
# plot(mds$x, mds$y,
#      xlab = paste("Dimension 1 (", round(mds$var.explained[1], 1), "%)", sep = ""),
#      ylab = paste("Dimension 2 (", round(mds$var.explained[2], 1), "%)", sep = ""),
#      main = "Sample Labels (Pre-ComBat)", pch = 19, col = "black", cex = 0.8)
# text(mds$x, mds$y, labels = colnames(betas_uncorrected), pos = 3, cex = 0.5)

# dev.off()

# print("✅ PCA BEFORE ComBat saved to: pca_before_combat.pdf")

# ==============================================================================
# STEP 6: SAVE RESULTS
# ==============================================================================

print("=== STEP 6: Saving results ===")

# Save beta matrix
write.csv(betas_final_top, "betas_combat_adjusted_final.csv")

# Rebuild metadata to save
meta_final <- sample_match[match(colnames(betas_final_top), rownames(sample_match)), ]
write.csv(meta_final, "metadata_final.csv", row.names = FALSE)

# Create summary statistics
summary_stats <- data.frame(
  Metric = c("Initial samples", "Initial probes", "Final samples", "Final probes",
             "Runs processed", "Common CpGs"),
  Value = c(ncol(betas), nrow(betas), ncol(betas_final_top), nrow(betas_final_top),
            length(unique_runs), length(common_cpgs_all))
)
write.csv(summary_stats, "processing_summary.csv", row.names = FALSE)

print("=== PIPELINE COMPLETED SUCCESSFULLY ===")
print("Final dataset:")
print(paste("Samples:", ncol(betas_adjusted)))
print(paste("Probes:", nrow(betas_adjusted)))
print(paste("Runs processed:", length(unique_runs)))
