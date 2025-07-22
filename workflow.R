# Methylation Data Analysis Pipeline
# IDAT -> QC -> BMIQ (per run) -> ComBat (across all runs) -> PCA

# ==============================================================================
# SETUP AND INITIALIZATION
# ==============================================================================

# Set library paths and working directory
.libPaths("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/packages")
setwd("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline")

# Load required libraries
print("Loading libraries...")
required_libs <- c("sesame", "minfi", "IlluminaMouseMethylationmanifest", 
                   "IlluminaMouseMethylationanno.12.v1.mm10", "limma", "sva", 
                   "ggsci", "wateRmelon", "snow", "BiocParallel", "parallel", 
                   "RColorBrewer", "readxl", "openxlsx", "ggplot2")

lapply(required_libs, library, character.only = TRUE)

# Source BMIQ function
print("Sourcing DoBMIQ function...")
source("bmiq/DoBMIQ.R")

# Load annotation data
print("Loading annotation data...")
annoMouse <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)

# ==============================================================================
# STEP 1: READ IDAT FILES AND GENERATE DETECTION P-VALUES
# ==============================================================================

print("=== STEP 1: Reading IDAT files ===")
idir <- "idat"
cores <- 4

# Read raw data and calculate detection p-values
betas <- openSesame(idir, BPPARAM = MulticoreParam(workers = cores), collapseToPfx = TRUE)
detP <- openSesame(idir, func = pOOBAH, return.pval = TRUE, BPPARAM = MulticoreParam(workers = cores))
detP <- betasCollapseToPfx(detP)

# Load and prepare sample metadata
targets <- read.csv("data/samplesheet.csv", header = TRUE)
targets$Beadchip <- as.character(targets$Beadchip)
targets$Basename <- paste0(targets$Beadchip, "_", targets$Position)

# Create sample names using Sample.ID instead of tb
targets$Samples <- paste0("TB", targets$Sample.ID, "_", targets$Experiment)

# Keep EPIC Run as is for batch grouping
targets$EPIC_RUN_CLEAN <- targets$EPIC.Run

# Remove duplicated samples
targets <- targets[!duplicated(targets$Basename), ]

print("Sample names in beta matrix (first 5):")
print(head(colnames(betas), 5))
print("Sample names in metadata (first 5):")
print(head(targets$Basename, 5))

# Align metadata with beta matrix
common_samples <- intersect(targets$Basename, colnames(betas))
print(paste("Common samples found:", length(common_samples)))
print(paste("Samples in beta matrix:", ncol(betas)))
print(paste("Samples in metadata:", nrow(targets)))

targets <- targets[targets$Basename %in% common_samples, ]
betas <- betas[, match(targets$Basename, colnames(betas))]
detP <- detP[, match(targets$Basename, colnames(detP))]

# Verify sample matching
if (!identical(targets$Basename, colnames(betas))) {
  stop("Sample matching failed!")
}

# Generate initial QC plots
print("Generating initial QC plots...")
pal <- pal_simpsons()(12)
pdf("meanDP_initial.pdf")
barplot(colMeans(detP), col = pal[factor(targets$Samples)], las = 2, 
        cex.names = 0.4, ylab = "Mean detection p-values")
abline(h = 0.05, col = "red")
legend("topleft", legend = levels(factor(targets$Samples)), fill = pal, bg = "white")
dev.off()

pdf("density_initial.pdf")
densityPlot(betas, sampGroups = targets$Samples, main = "Beta Density (Initial)", 
            legend = FALSE, pal = pal)
legend("top", legend = levels(factor(targets$Samples)), text.col = pal)
dev.off()

# ==============================================================================
# STEP 2: QC FILTERING BY RUN
# ==============================================================================

print("=== STEP 2: QC filtering by run ===")

qc_filter_run <- function(betas_run, detP_run, targets_run, run_id) {
  print(paste("Processing RUN", run_id))
  
  # Filter samples with poor detection p-values
  good_samples <- which(colMeans(detP_run) < 0.05)
  if (length(good_samples) == 0) {
    warning(paste("No good samples in RUN", run_id))
    return(NULL)
  }
  
  betas_run <- betas_run[, good_samples, drop = FALSE]
  detP_run <- detP_run[, good_samples, drop = FALSE]
  targets_run <- targets_run[good_samples, ]
  
  # Keep probes detected in all remaining samples
  keep_probes <- rowSums(detP_run < 0.05) == ncol(betas_run)
  betas_run <- betas_run[keep_probes, ]
  detP_run <- detP_run[keep_probes, ]
  
  # Remove probes with NA/NaN values
  complete_idx <- complete.cases(betas_run)
  betas_run <- betas_run[complete_idx, ]
  detP_run <- detP_run[complete_idx, ]
  
  # Keep only CpG probes
  cpg_idx <- grep("^cg", rownames(betas_run))
  betas_run <- betas_run[cpg_idx, ]
  detP_run <- detP_run[rownames(betas_run), ]
  
  # Remove sex chromosome probes
  sex_cpgs <- annoMouse$Name[annoMouse$chr %in% c("chrX", "chrY")]
  keep_sex <- !rownames(betas_run) %in% sex_cpgs
  betas_run <- betas_run[keep_sex, ]
  detP_run <- detP_run[rownames(betas_run), ]
  
  print(paste("RUN", run_id, "after QC:", nrow(betas_run), "probes x", ncol(betas_run), "samples"))
  
  return(list(betas = betas_run, detP = detP_run, targets = targets_run))
}

# Apply QC filtering to each run
unique_runs <- sort(unique(targets$EPIC_RUN_CLEAN))
filtered_data <- list()

for (run_id in unique_runs) {
  idx <- targets$EPIC_RUN_CLEAN == run_id
  result <- qc_filter_run(betas[, idx], detP[, idx], targets[idx, ], run_id)
  if (!is.null(result)) {
    filtered_data[[paste0("RUN_", gsub("_", "", run_id))]] <- result
  }
}

# ==============================================================================
# STEP 3: BMIQ NORMALIZATION (PER RUN)
# ==============================================================================

print("=== STEP 3: BMIQ normalization per run ===")

run_bmiq_normalization <- function(betas_run, run_name) {
  print(paste("Running BMIQ on", run_name))
  
  temp_file <- paste0("temp_betas_", run_name, ".csv")
  write.csv(data.frame(CpG_ID = rownames(betas_run), betas_run), 
            file = temp_file, row.names = FALSE)
  
  tryCatch({
    DoBMIQ(probes_path = "data/probesample.xlsx", beta_path = temp_file)
    load("bmiq/results/bmiq.Rd")
    
    # Save BMIQ results
    write.csv(bmiq.m, paste0("betas_bmiq_", run_name, ".csv"))
    print(paste("BMIQ completed for", run_name, "- Probes:", nrow(bmiq.m), "Samples:", ncol(bmiq.m)))
    
    return(bmiq.m)
  }, error = function(e) {
    print(paste("BMIQ failed for", run_name, ":", e$message))
    print("Using raw data for this run...")
    return(betas_run)
  }, finally = {
    if (file.exists(temp_file)) file.remove(temp_file)
  })
}

# Apply BMIQ to each run
betas_bmiq_by_run <- list()
for (run_name in names(filtered_data)) {
  betas_bmiq_by_run[[run_name]] <- run_bmiq_normalization(
    filtered_data[[run_name]]$betas, run_name
  )
}

# ==============================================================================
# STEP 4: COMBINE RUNS AND COMBAT BATCH CORRECTION
# ==============================================================================

print("=== STEP 4: ComBat batch correction ===")

# Combine all BMIQ-normalized runs
common_cpgs <- Reduce(intersect, lapply(betas_bmiq_by_run, rownames))
betas_combined <- do.call(cbind, lapply(betas_bmiq_by_run, function(mat) {
  mat[common_cpgs, , drop = FALSE]
}))

# Load samplesheet and match samples
samplesheet <- read.csv("data/samplesheet.csv", header = TRUE)
samplesheet$Beadchip <- as.character(samplesheet$Beadchip)
samplesheet$Basename <- paste0(samplesheet$Beadchip, "_", samplesheet$Position)
samplesheet$Samples <- paste0("TB", samplesheet$Sample.ID, "_", samplesheet$Experiment)

# Keep EPIC Run for batch grouping
samplesheet$EPIC_RUN_CLEAN <- samplesheet$EPIC.Run

# Match samples (remove X prefix from beta matrix column names)
betas_colnames_clean <- gsub("^X", "", colnames(betas_combined))
matching_samples <- samplesheet[match(betas_colnames_clean, samplesheet$Basename), ]
matching_samples <- matching_samples[!is.na(matching_samples$Basename), ]
betas_combined <- betas_combined[, !is.na(samplesheet[match(betas_colnames_clean, samplesheet$Basename), ]$Basename)]
colnames(betas_combined) <- matching_samples$Samples

# Create batch groups: BATCH1 (1_2), BATCH2 (3_4 + 5)
matching_samples$BATCH <- ifelse(matching_samples$EPIC_RUN_CLEAN == "1_2", "BATCH1", "BATCH2")
batch_vector <- factor(matching_samples$BATCH)

# Use Tissue Cell Type as biological groups
group_vector <- factor(matching_samples$Tissue.Cell.type)

print("EPIC Run distribution:")
print(table(matching_samples$EPIC_RUN_CLEAN))
print("Tissue Cell Type distribution:")
print(table(matching_samples$Tissue.Cell.type))
print("Batch x Tissue Cell Type design:")
print(table(batch_vector, group_vector))

# Ensure ComBat normalization is done properly
print("Preparing data for ComBat normalization...")

# Remove any samples/probes with missing data
print("Checking for missing data...")
missing_samples <- colSums(is.na(betas_combined)) > 0
missing_probes <- rowSums(is.na(betas_combined)) > 0

if(any(missing_samples)) {
  print(paste("Removing", sum(missing_samples), "samples with missing data"))
  betas_combined <- betas_combined[, !missing_samples]
  batch_vector <- batch_vector[!missing_samples]
  group_vector <- group_vector[!missing_samples]
  matching_samples <- matching_samples[!missing_samples, ]
}

if(any(missing_probes)) {
  print(paste("Removing", sum(missing_probes), "probes with missing data"))
  betas_combined <- betas_combined[!missing_probes, ]
}

# Clamp beta values to avoid infinite M-values
print("Clamping beta values...")
betas_combined[betas_combined <= 0.001] <- 0.001
betas_combined[betas_combined >= 0.999] <- 0.999

# Convert to M-values for ComBat
print("Converting to M-values...")
mvals_combined <- log2(betas_combined / (1 - betas_combined))

# Check for any remaining infinite or NaN values
inf_check <- sum(is.infinite(mvals_combined))
nan_check <- sum(is.nan(mvals_combined))
print(paste("Infinite values:", inf_check))
print(paste("NaN values:", nan_check))

if(inf_check > 0 || nan_check > 0) {
  stop("Still have infinite or NaN values after clamping!")
}

# Run ComBat normalization properly
print("Running ComBat batch correction...")

# Create design matrix with tissue cell types as biological groups
mod <- model.matrix(~ group_vector)
print("Design matrix for tissue cell types:")
print(colnames(mod))
print(paste("Design matrix dimensions:", nrow(mod), "x", ncol(mod)))
print(paste("Batch vector length:", length(batch_vector)))
print(paste("M-values matrix dimensions:", nrow(mvals_combined), "x", ncol(mvals_combined)))

# Show the biological groups being preserved
print("Biological groups (Tissue Cell Types) being preserved:")
print(levels(group_vector))

# Verify everything matches
if(length(batch_vector) != ncol(mvals_combined)) {
  stop("Batch vector length doesn't match number of samples!")
}
if(length(group_vector) != ncol(mvals_combined)) {
  stop("Group vector length doesn't match number of samples!")
}

# Run ComBat
tryCatch({
  mvals_adjusted <- ComBat(
    dat = mvals_combined,
    batch = batch_vector,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE,
  )
  
  # Convert back to beta values
  betas_adjusted <- 2^mvals_adjusted / (1 + 2^mvals_adjusted)
  
  print("ComBat normalization completed successfully!")
  print(paste("Final dimensions:", nrow(betas_adjusted), "probes x", ncol(betas_adjusted), "samples"))
  
}, error = function(e) {
  print(paste("ComBat failed with error:", e$message))
  print("This might be due to:")
  print("1. Insufficient samples per batch")
  print("2. Confounded batch-group design")
  print("3. Singular design matrix")
  
  # Show batch/group breakdown
  print("Batch x Tissue Cell Type table:")
  print(table(batch_vector, group_vector))
  
  print("Using uncorrected data...")
  betas_adjusted <<- betas_combined
})

# ==============================================================================
# STEP 4.5: SELECT TOP 10K MOST VARIABLE CPGS (DESeq2 METHOD)
# ==============================================================================

print("=== STEP 4.5: Selecting top 10k most variable CpGs (DESeq2 method) ===")

# Transform to M-values first (equivalent to DESeq2's variance stabilizing transformation)
print("Converting to M-values for variance calculation...")
mvals_for_selection <- log2(betas_adjusted / (1 - betas_adjusted))

# Calculate row variance on transformed data (M-values) - DESeq2 approach
print("Calculating variance on M-values (transformed space)...")
cpg_variances <- apply(mvals_for_selection, 1, var, na.rm = TRUE)

# Select ntop (10,000) genes by variance in transformed space - exactly like DESeq2
print("Selecting top 10k CpGs by variance in transformed space...")
ntop <- 10000
select_var <- order(cpg_variances, decreasing = TRUE)[seq_len(min(ntop, length(cpg_variances)))]
top_10k_cpgs <- rownames(mvals_for_selection)[select_var]

# Subset both beta and M-value matrices to top 10k CpGs
betas_top10k <- betas_adjusted[top_10k_cpgs, ]
mvals_top10k <- mvals_for_selection[top_10k_cpgs, ]

print(paste("Selected", nrow(betas_top10k), "CpGs out of", nrow(betas_adjusted), "total CpGs"))
print(paste("Variance range of selected CpGs (in M-value space):", 
            round(min(cpg_variances[select_var]), 4), 
            "to", 
            round(max(cpg_variances[select_var]), 4)))

# Save the top 10k CpGs dataset
print("Saving top 10k CpGs dataset...")
write.csv(betas_top10k, "betas_combat_adjusted_top10k_cpgs_deseq2method.csv")
write.csv(mvals_top10k, "mvals_combat_adjusted_top10k_cpgs_deseq2method.csv")

# Save the list of selected CpGs and their variances (in M-value space)
top10k_info <- data.frame(
  CpG_ID = top_10k_cpgs,
  Variance_Mvals = cpg_variances[select_var],
  Rank = seq_len(length(top_10k_cpgs))
)
write.csv(top10k_info, "top10k_cpgs_info_deseq2method.csv", row.names = FALSE)

print("Top 10k CpGs selection completed using DESeq2 methodology!")
print("Note: Variance calculated on M-values (transformed space) as per DESeq2 best practices")

# ==============================================================================
# STEP 5A: PCA ANALYSIS (USING TOP 10K CPGS - DESeq2 METHOD)
# ==============================================================================

print("=== STEP 5A: PCA analysis (using top 10k CpGs - DESeq2 method) ===")

# Use M-values for PCA (already calculated above)
betas_for_pca <- betas_top10k
mvals_for_pca <- mvals_top10k
group_for_pca <- group_vector
batch_for_pca <- batch_vector

print("Debugging PCA setup...")
print(paste("Samples in data:", ncol(betas_for_pca)))
print(paste("CpGs in data:", nrow(betas_for_pca)))
print(paste("Length of group_for_pca:", length(group_for_pca)))
print(paste("Length of batch_for_pca:", length(batch_for_pca)))

print("Sample names in data (first 5):")
print(head(colnames(betas_for_pca), 5))
print("Group vector (first 5):")
print(head(as.character(group_for_pca), 5))

# Perform PCA on M-values (centered, not scaled - DESeq2 default)
pca_result <- prcomp(t(mvals_for_pca), center = TRUE, scale. = FALSE)

print("PCA completed successfully!")
print(paste("PCA sample names (first 5):"))
print(head(rownames(pca_result$x), 5))

# Create additional grouping vectors for PCA
genotype_for_pca <- factor(matching_samples$Genotype)
cross_type_for_pca <- factor(matching_samples$cross.type)

# Create PCA dataframe with proper alignment
pca_df_top10k <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = rownames(pca_result$x),
  Tissue_Cell_Type = group_for_pca[match(rownames(pca_result$x), colnames(betas_for_pca))],
  Batch = batch_for_pca[match(rownames(pca_result$x), colnames(betas_for_pca))],
  Genotype = genotype_for_pca[match(rownames(pca_result$x), colnames(betas_for_pca))],
  Cross_Type = cross_type_for_pca[match(rownames(pca_result$x), colnames(betas_for_pca))]
)

# Check for NAs in the PCA dataframe
print("Checking PCA dataframe:")
print(paste("Rows in pca_df:", nrow(pca_df_top10k)))
print(paste("NAs in Tissue_Cell_Type:", sum(is.na(pca_df_top10k$Tissue_Cell_Type))))
print(paste("NAs in Batch:", sum(is.na(pca_df_top10k$Batch))))
print(paste("NAs in Genotype:", sum(is.na(pca_df_top10k$Genotype))))
print(paste("NAs in Cross_Type:", sum(is.na(pca_df_top10k$Cross_Type))))

# Plot by tissue cell type
pdf("PCA_Combat_adjusted_tissue_type_top10k_deseq2.pdf", width = 10, height = 8)
p1 <- ggplot(pca_df_top10k, aes(x = PC1, y = PC2, color = Tissue_Cell_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on ComBat Adjusted M-values (by Tissue Cell Type) - Top 10k CpGs (DESeq2 method)",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p1)
dev.off()

# Plot by batch
pdf("PCA_Combat_adjusted_batch_top10k_deseq2.pdf", width = 10, height = 8)
p2 <- ggplot(pca_df_top10k, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on ComBat Adjusted M-values (by Batch) - Top 10k CpGs (DESeq2 method)",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p2)
dev.off()

# Plot by genotype
pdf("PCA_Combat_adjusted_genotype_top10k_deseq2.pdf", width = 10, height = 8)
p3 <- ggplot(pca_df_top10k, aes(x = PC1, y = PC2, color = Genotype)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on ComBat Adjusted M-values (by Genotype) - Top 10k CpGs (DESeq2 method)",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p3)
dev.off()

# Plot by cross type
pdf("PCA_Combat_adjusted_cross_type_top10k_deseq2.pdf", width = 10, height = 8)
p4 <- ggplot(pca_df_top10k, aes(x = PC1, y = PC2, color = Cross_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "PCA on ComBat Adjusted M-values (by Cross Type) - Top 10k CpGs (DESeq2 method)",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p4)
dev.off()

print("PCA plots (Top 10k CpGs - DESeq2 method) saved successfully!")

# ==============================================================================
# STEP 5B: PCA ANALYSIS (USING ALL CPGS)
# ==============================================================================

print("=== STEP 5B: PCA analysis (using ALL CpGs) ===")

# Use ALL CpGs for PCA comparison
betas_all_cpgs <- betas_adjusted
mvals_all_cpgs <- log2(betas_all_cpgs / (1 - betas_all_cpgs))

print("Debugging PCA setup for ALL CpGs...")
print(paste("Samples in ALL CpG data:", ncol(betas_all_cpgs)))
print(paste("CpGs in ALL CpG data:", nrow(betas_all_cpgs)))

# Perform PCA on M-values using ALL CpGs
print("Performing PCA on ALL CpGs...")
pca_result_all <- prcomp(t(mvals_all_cpgs), center = TRUE, scale. = FALSE)

print("PCA on ALL CpGs completed successfully!")

# Create PCA dataframe for ALL CpGs
pca_df_all_cpgs <- data.frame(
  PC1 = pca_result_all$x[, 1],
  PC2 = pca_result_all$x[, 2],
  Sample = rownames(pca_result_all$x),
  Tissue_Cell_Type = group_for_pca[match(rownames(pca_result_all$x), colnames(betas_all_cpgs))],
  Batch = batch_for_pca[match(rownames(pca_result_all$x), colnames(betas_all_cpgs))],
  Genotype = genotype_for_pca[match(rownames(pca_result_all$x), colnames(betas_all_cpgs))],
  Cross_Type = cross_type_for_pca[match(rownames(pca_result_all$x), colnames(betas_all_cpgs))]
)

# Plot by tissue cell type (ALL CpGs)
pdf("PCA_Combat_adjusted_tissue_type_ALL_cpgs.pdf", width = 10, height = 8)
p5 <- ggplot(pca_df_all_cpgs, aes(x = PC1, y = PC2, color = Tissue_Cell_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on ComBat Adjusted M-values (by Tissue Cell Type) - ALL CpGs",
    x = paste0("PC1 (", round(summary(pca_result_all)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_all)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p5)
dev.off()

# Plot by batch (ALL CpGs)
pdf("PCA_Combat_adjusted_batch_ALL_cpgs.pdf", width = 10, height = 8)
p6 <- ggplot(pca_df_all_cpgs, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on ComBat Adjusted M-values (by Batch) - ALL CpGs",
    x = paste0("PC1 (", round(summary(pca_result_all)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_all)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p6)
dev.off()

# Plot by genotype (ALL CpGs)
pdf("PCA_Combat_adjusted_genotype_ALL_cpgs.pdf", width = 10, height = 8)
p7 <- ggplot(pca_df_all_cpgs, aes(x = PC1, y = PC2, color = Genotype)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on ComBat Adjusted M-values (by Genotype) - ALL CpGs",
    x = paste0("PC1 (", round(summary(pca_result_all)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_all)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p7)
dev.off()

# Plot by cross type (ALL CpGs)
pdf("PCA_Combat_adjusted_cross_type_ALL_cpgs.pdf", width = 10, height = 8)
p8 <- ggplot(pca_df_all_cpgs, aes(x = PC1, y = PC2, color = Cross_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "PCA on ComBat Adjusted M-values (by Cross Type) - ALL CpGs",
    x = paste0("PC1 (", round(summary(pca_result_all)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_all)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p8)
dev.off()

print("PCA plots (ALL CpGs) saved successfully!")

# ==============================================================================
# COMPARISON SUMMARY
# ==============================================================================

print("=== PCA COMPARISON SUMMARY ===")
print("Generated 8 PCA plots total:")
print("Top 10k CpGs (DESeq2 method):")
print("  - PCA_Combat_adjusted_tissue_type_top10k_deseq2.pdf")
print("  - PCA_Combat_adjusted_batch_top10k_deseq2.pdf") 
print("  - PCA_Combat_adjusted_genotype_top10k_deseq2.pdf")
print("  - PCA_Combat_adjusted_cross_type_top10k_deseq2.pdf")
print("")
print("ALL CpGs:")
print("  - PCA_Combat_adjusted_tissue_type_ALL_cpgs.pdf")
print("  - PCA_Combat_adjusted_batch_ALL_cpgs.pdf")
print("  - PCA_Combat_adjusted_genotype_ALL_cpgs.pdf") 
print("  - PCA_Combat_adjusted_cross_type_ALL_cpgs.pdf")
print("")
print("Variance explained comparison:")
print(paste("Top 10k CpGs - PC1:", round(summary(pca_result)$importance[2, 1] * 100, 1), "%, PC2:", round(summary(pca_result)$importance[2, 2] * 100, 1), "%"))
print(paste("ALL CpGs - PC1:", round(summary(pca_result_all)$importance[2, 1] * 100, 1), "%, PC2:", round(summary(pca_result_all)$importance[2, 2] * 100, 1), "%"))

# Save both PCA coordinate sets for comparison
write.csv(pca_df_top10k, "pca_coordinates_top10k_deseq2.csv", row.names = FALSE)
write.csv(pca_df_all_cpgs, "pca_coordinates_ALL_cpgs.csv", row.names = FALSE)