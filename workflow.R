# Methylation Data Analysis Pipeline (Improved)
# IDAT -> QC -> BMIQ (per run) -> ComBat (across all runs) -> PCA

# ---- Setup 

.libPaths("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/packages")
setwd("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline")

required_libs <- c("sesame", "minfi", "IlluminaMouseMethylationmanifest", 
                   "IlluminaMouseMethylationanno.12.v1.mm10", "limma", "sva", 
                   "ggsci", "wateRmelon", "snow", "BiocParallel", "parallel", 
                   "RColorBrewer", "readxl", "openxlsx", "ggplot2", "ComBatMet")
lapply(required_libs, library, character.only = TRUE)

source("bmiq/DoBMIQ.R")
annoMouse <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)

# ---- Step 1: Read idat files

idir <- "idat"
cores <- 4

betas <- openSesame(idir, BPPARAM = MulticoreParam(workers = cores), collapseToPfx = TRUE)
detP <- openSesame(idir, func = pOOBAH, return.pval = TRUE, BPPARAM = MulticoreParam(workers = cores))
detP <- betasCollapseToPfx(detP)

targets <- read.csv("data/samplesheet.csv", header = TRUE)
targets$Beadchip <- as.character(targets$Beadchip)
targets$Basename <- paste0(targets$Beadchip, "_", targets$Position)
targets$Samples <- paste0(targets$Sample.ID, "_", targets$Experiment)
targets$EPIC_RUN_CLEAN <- targets$EPIC.Run
targets <- targets[!duplicated(targets$Basename), ]

common_samples <- intersect(targets$Basename, colnames(betas))
targets <- targets[targets$Basename %in% common_samples, ]
betas <- betas[, match(targets$Basename, colnames(betas))]
detP <- detP[, match(targets$Basename, colnames(detP))]

stopifnot(identical(targets$Basename, colnames(betas)))

# ---- Step 2: Run by run Qc checks

qc_filter_run <- function(betas_run, detP_run, targets_run, run_id) {
  good_samples <- which(colMeans(detP_run) < 0.05)
  if (length(good_samples) == 0) stop(paste("No good samples in RUN", run_id))

  betas_run <- betas_run[, good_samples, drop = FALSE]
  detP_run <- detP_run[, good_samples, drop = FALSE]
  targets_run <- targets_run[good_samples, ]

  keep_probes <- rowSums(detP_run < 0.05) == ncol(betas_run)
  betas_run <- betas_run[keep_probes, ]
  detP_run <- detP_run[keep_probes, ]

  complete_idx <- complete.cases(betas_run)
  betas_run <- betas_run[complete_idx, ]
  detP_run <- detP_run[complete_idx, ]

  cpg_idx <- grep("^cg", rownames(betas_run))
  betas_run <- betas_run[cpg_idx, ]
  sex_cpgs <- annoMouse$Name[annoMouse$chr %in% c("chrX", "chrY")]
  betas_run <- betas_run[!rownames(betas_run) %in% sex_cpgs, ]
  detP_run <- detP_run[rownames(betas_run), ]

  return(list(betas = betas_run, detP = detP_run, targets = targets_run))
}

unique_runs <- sort(unique(targets$EPIC_RUN_CLEAN))
filtered_data <- list()

for (run_id in unique_runs) {
  idx <- targets$EPIC_RUN_CLEAN == run_id
  result <- qc_filter_run(betas[, idx], detP[, idx], targets[idx, ], run_id)
  filtered_data[[paste0("RUN_", gsub("_", "", run_id))]] <- result
}

# Combine betas after QC filtering (pre-BMIQ)
common_cpgs_qc <- Reduce(intersect, lapply(filtered_data, function(x) rownames(x$betas)))
betas_qc_combined <- do.call(cbind, lapply(filtered_data, function(x) x$betas[common_cpgs_qc, , drop = FALSE]))
write.csv(betas_qc_combined, "betas_qc_combined.csv", row.names = TRUE)

# ---- Step 3: Run by Run BMIQ checks

run_bmiq_normalization <- function(betas_run, run_name) {
  temp_file <- paste0("temp_betas_", run_name, ".csv")
  write.csv(data.frame(CpG_ID = rownames(betas_run), betas_run), file = temp_file, row.names = FALSE)

  tryCatch({
    DoBMIQ(probes_path = "data/probesample.xlsx", beta_path = temp_file)
    load("bmiq/results/bmiq.Rd")
    return(bmiq.m)
  }, error = function(e) {
    stop(paste("BMIQ failed for", run_name, ":", e$message))
  }, finally = {
    if (file.exists(temp_file)) file.remove(temp_file)
  })
}

# combining BMIQ runs

betas_bmiq_by_run <- lapply(names(filtered_data), function(run_name) {
  run_bmiq_normalization(filtered_data[[run_name]]$betas, run_name)
})
names(betas_bmiq_by_run) <- names(filtered_data)

common_cpgs <- Reduce(intersect, lapply(betas_bmiq_by_run, rownames))
betas_combined <- do.call(cbind, lapply(betas_bmiq_by_run, function(mat) mat[common_cpgs, , drop = FALSE]))
write.csv(betas_combined, "betas_bmiq.csv", row.names = TRUE)

# ---- Step 4: ComBat normalization

betas_colnames_clean <- gsub("^X", "", colnames(betas_combined))
matching_samples <- targets[match(betas_colnames_clean, targets$Basename), ]
matching_samples <- matching_samples[!is.na(matching_samples$Basename), ]
betas_combined <- betas_combined[, !is.na(targets[match(betas_colnames_clean, targets$Basename), ]$Basename)]
colnames(betas_combined) <- matching_samples$Samples

missing_samples <- colSums(is.na(betas_combined)) > 0
missing_probes <- rowSums(is.na(betas_combined)) > 0
if (any(missing_samples)) {
  betas_combined <- betas_combined[, !missing_samples]
  matching_samples <- matching_samples[!missing_samples, ]
}
if (any(missing_probes)) {
  betas_combined <- betas_combined[!missing_probes, ]
}

matching_samples$BATCH <- ifelse(matching_samples$EPIC_RUN_CLEAN == "1_2", "BATCH1", "BATCH2")
batch_vector <- factor(matching_samples$BATCH)
group_vector <- factor(matching_samples$Tissue.Cell.type)

mod <- model.matrix(~as.factor(Tissue.Cell.type), data=matching_samples)
stopifnot(identical(colnames(betas_combined), matching_samples$Samples))

mvals_adjusted <- log2(betas_combined / (1 - betas_combined))
mvals_adjusted[is.infinite(mvals_adjusted)] <- NA
mvals_adjusted <- mvals_adjusted[complete.cases(mvals_adjusted), ]

mvals_combat <- ComBat(
  dat = mvals_adjusted,
  batch = as.numeric(batch_vector),
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

betas_adjusted <- 2^mvals_combat / (1 + 2^mvals_combat)
betas_adjusted[betas_adjusted < 0] <- 0
betas_adjusted[betas_adjusted > 1] <- 1
write.csv(betas_adjusted, "betas_combat_adjusted.csv", row.names = TRUE)

# ---- Step 5: Choosing top 10k CpGs

# Combat
cpg_variances <- apply(mvals_combat, 1, var, na.rm = TRUE)
ntop <- 10000
select_var <- order(cpg_variances, decreasing = TRUE)[seq_len(min(ntop, length(cpg_variances)))]
top_10k_cpgs <- rownames(mvals_combat)[select_var]

betas_top10k <- betas_adjusted[top_10k_cpgs, ]
mvals_top10k <- mvals_combat[top_10k_cpgs, ]

write.csv(betas_top10k, "betas_combat_adjusted_top10k_cpgs_deseq2method.csv")
write.csv(mvals_top10k, "mvals_combat_adjusted_top10k_cpgs_deseq2method.csv")

# BMIQ 
mvals_bmiq <- log2(betas_combined / (1 - betas_combined))
mvals_bmiq[is.infinite(mvals_bmiq)] <- NA

cpg_variances_bmiq <- apply(mvals_bmiq, 1, var, na.rm = TRUE)
ntop <- 10000
select_var_bmiq <- order(cpg_variances_bmiq, decreasing = TRUE)[seq_len(min(ntop, length(cpg_variances_bmiq)))]
top_10k_cpgs_bmiq <- rownames(betas_combined)[select_var_bmiq]

betas_top10k_bmiq <- betas_combined[top_10k_cpgs_bmiq, ]
mvals_top10k_bmiq <- mvals_bmiq[top_10k_cpgs_bmiq, ]

write.csv(betas_top10k_bmiq, "betas_BMIQ_adjusted_top10k_cpgs_deseq2method.csv")
write.csv(mvals_top10k_bmiq, "mvals_BMIQ_adjusted_top10k_cpgs_deseq2method.csv")

# Raw
mvals_qc <- log2(betas_qc_combined / (1 - betas_qc_combined))
mvals_qc[is.infinite(mvals_qc)] <- NA

cpg_variances_qc <- apply(mvals_qc, 1, var, na.rm = TRUE)
ntop <- 10000
select_var_qc <- order(cpg_variances_qc, decreasing = TRUE)[seq_len(min(ntop, length(cpg_variances_qc)))]
top_10k_cpgs_qc <- rownames(betas_qc_combined)[select_var_qc]

betas_top10k_qc <- betas_qc_combined[top_10k_cpgs_qc, ]
mvals_top10k_qc <- mvals_qc[top_10k_cpgs_qc, ]

write.csv(betas_top10k_qc, "betas_QC_adjusted_top10k_cpgs_deseq2method.csv")
write.csv(mvals_top10k_qc, "mvals_QC_adjusted_top10k_cpgs_deseq2method.csv")

# ---- Step 6: PCA analysis (top 10k CpGs)

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
genotype_for_pca <- factor(matching_samples$True.Genotype)
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

# ---- Step 5B: PCA analysis (All CpGs)

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

# ---- Step 5C: PCA analysis (top 10k for BMIQ + Qc (Raw data))

print("=== PCA analysis (using top 10k CpGs - BMIQ data) ===")

# Use M-values for PCA (already calculated above)
betas_for_pca_bmiq <- betas_top10k_bmiq
mvals_for_pca_bmiq <- mvals_top10k_bmiq
group_for_pca_bmiq <- group_vector
batch_for_pca_bmiq <- batch_vector

print("Debugging PCA setup for BMIQ...")
print(paste("Samples in data:", ncol(betas_for_pca_bmiq)))
print(paste("CpGs in data:", nrow(betas_for_pca_bmiq)))
print(paste("Length of group_for_pca:", length(group_for_pca_bmiq)))
print(paste("Length of batch_for_pca:", length(batch_for_pca_bmiq)))

print("Sample names in data (first 5):")
print(head(colnames(betas_for_pca_bmiq), 5))
print("Group vector (first 5):")
print(head(as.character(group_for_pca_bmiq), 5))

# Perform PCA on M-values (centered, not scaled - DESeq2 default)
pca_result_bmiq <- prcomp(t(mvals_for_pca_bmiq), center = TRUE, scale. = FALSE)

print("PCA completed successfully for BMIQ!")
print(paste("PCA sample names (first 5):"))
print(head(rownames(pca_result_bmiq$x), 5))

# Create additional grouping vectors for PCA
genotype_for_pca_bmiq <- factor(matching_samples$True.Genotype)
cross_type_for_pca_bmiq <- factor(matching_samples$cross.type)

# Create PCA dataframe with proper alignment
pca_df_top10k_bmiq <- data.frame(
  PC1 = pca_result_bmiq$x[, 1],
  PC2 = pca_result_bmiq$x[, 2],
  Sample = rownames(pca_result_bmiq$x),
  Tissue_Cell_Type = group_for_pca_bmiq[match(rownames(pca_result_bmiq$x), colnames(betas_for_pca_bmiq))],
  Batch = batch_for_pca_bmiq[match(rownames(pca_result_bmiq$x), colnames(betas_for_pca_bmiq))],
  Genotype = genotype_for_pca_bmiq[match(rownames(pca_result_bmiq$x), colnames(betas_for_pca_bmiq))],
  Cross_Type = cross_type_for_pca_bmiq[match(rownames(pca_result_bmiq$x), colnames(betas_for_pca_bmiq))]
)

# Check for NAs in the PCA dataframe
print("Checking PCA dataframe for BMIQ:")
print(paste("Rows in pca_df:", nrow(pca_df_top10k_bmiq)))
print(paste("NAs in Tissue_Cell_Type:", sum(is.na(pca_df_top10k_bmiq$Tissue_Cell_Type))))
print(paste("NAs in Batch:", sum(is.na(pca_df_top10k_bmiq$Batch))))
print(paste("NAs in Genotype:", sum(is.na(pca_df_top10k_bmiq$Genotype))))
print(paste("NAs in Cross_Type:", sum(is.na(pca_df_top10k_bmiq$Cross_Type))))

# Plot by tissue cell type
pdf("PCA_BMIQ_tissue_type_top10k.pdf", width = 10, height = 8)
p1_bmiq <- ggplot(pca_df_top10k_bmiq, aes(x = PC1, y = PC2, color = Tissue_Cell_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on BMIQ Normalized M-values (by Tissue Cell Type) - Top 10k CpGs",
    x = paste0("PC1 (", round(summary(pca_result_bmiq)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_bmiq)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p1_bmiq)
dev.off()

# Plot by batch
pdf("PCA_BMIQ_batch_top10k.pdf", width = 10, height = 8)
p2_bmiq <- ggplot(pca_df_top10k_bmiq, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on BMIQ Normalized M-values (by Batch) - Top 10k CpGs",
    x = paste0("PC1 (", round(summary(pca_result_bmiq)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_bmiq)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p2_bmiq)
dev.off()

print("PCA plots (Top 10k CpGs - BMIQ data) saved successfully!")

print("=== PCA analysis (using top 10k CpGs - QC data) ===")

# Use M-values for PCA (already calculated above)
betas_for_pca_qc <- betas_top10k_qc
mvals_for_pca_qc <- mvals_top10k_qc
group_for_pca_qc <- group_vector
batch_for_pca_qc <- batch_vector

print("Debugging PCA setup for QC...")
print(paste("Samples in data:", ncol(betas_for_pca_qc)))
print(paste("CpGs in data:", nrow(betas_for_pca_qc)))
print(paste("Length of group_for_pca:", length(group_for_pca_qc)))
print(paste("Length of batch_for_pca:", length(batch_for_pca_qc)))

print("Sample names in data (first 5):")
print(head(colnames(betas_for_pca_qc), 5))
print("Group vector (first 5):")
print(head(as.character(group_for_pca_qc), 5))

# Perform PCA on M-values (centered, not scaled - DESeq2 default)
pca_result_qc <- prcomp(t(mvals_for_pca_qc), center = TRUE, scale. = FALSE)

print("PCA completed successfully for QC!")
print(paste("PCA sample names (first 5):"))
print(head(rownames(pca_result_qc$x), 5))

# Create additional grouping vectors for PCA
genotype_for_pca_qc <- factor(matching_samples$True.Genotype)
cross_type_for_pca_qc <- factor(matching_samples$cross.type)

# Create PCA dataframe with proper alignment
pca_df_top10k_qc <- data.frame(
  PC1 = pca_result_qc$x[, 1],
  PC2 = pca_result_qc$x[, 2],
  Sample = rownames(pca_result_qc$x),
  Tissue_Cell_Type = group_for_pca_qc[match(rownames(pca_result_qc$x), colnames(betas_for_pca_qc))],
  Batch = batch_for_pca_qc[match(rownames(pca_result_qc$x), colnames(betas_for_pca_qc))],
  Genotype = genotype_for_pca_qc[match(rownames(pca_result_qc$x), colnames(betas_for_pca_qc))],
  Cross_Type = cross_type_for_pca_qc[match(rownames(pca_result_qc$x), colnames(betas_for_pca_qc))]
)

# Check for NAs in the PCA dataframe
print("Checking PCA dataframe for QC:")
print(paste("Rows in pca_df:", nrow(pca_df_top10k_qc)))
print(paste("NAs in Tissue_Cell_Type:", sum(is.na(pca_df_top10k_qc$Tissue_Cell_Type))))
print(paste("NAs in Batch:", sum(is.na(pca_df_top10k_qc$Batch))))
print(paste("NAs in Genotype:", sum(is.na(pca_df_top10k_qc$Genotype))))
print(paste("NAs in Cross_Type:", sum(is.na(pca_df_top10k_qc$Cross_Type))))

# Plot by tissue cell type
pdf("PCA_QC_tissue_type_top10k.pdf", width = 10, height = 8)
p1_qc <- ggplot(pca_df_top10k_qc, aes(x = PC1, y = PC2, color = Tissue_Cell_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on QC-filtered M-values (by Tissue Cell Type) - Top 10k CpGs",
    x = paste0("PC1 (", round(summary(pca_result_qc)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_qc)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p1_qc)
dev.off()

# Plot by batch
pdf("PCA_QC_batch_top10k.pdf", width = 10, height = 8)
p2_qc <- ggplot(pca_df_top10k_qc, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on QC-filtered M-values (by Batch) - Top 10k CpGs",
    x = paste0("PC1 (", round(summary(pca_result_qc)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_qc)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p2_qc)
dev.off()

print("PCA plots (Top 10k CpGs - QC data) saved successfully!")

# ---- Step 5C: PCA analysis (ALL CpGs for BMIQ + Qc (Raw data))

print("=== PCA analysis (using ALL CpGs - BMIQ data) ===")

# Use M-values for PCA (all CpGs)
betas_for_pca_bmiq_all <- betas_combined
mvals_for_pca_bmiq_all <- mvals_bmiq
group_for_pca_bmiq_all <- group_vector
batch_for_pca_bmiq_all <- batch_vector

print("Debugging PCA setup for BMIQ (ALL CpGs)...")
print(paste("Samples in data:", ncol(betas_for_pca_bmiq_all)))
print(paste("CpGs in data:", nrow(betas_for_pca_bmiq_all)))
print(paste("Length of group_for_pca:", length(group_for_pca_bmiq_all)))
print(paste("Length of batch_for_pca:", length(batch_for_pca_bmiq_all)))

print("Sample names in data (first 5):")
print(head(colnames(betas_for_pca_bmiq_all), 5))
print("Group vector (first 5):")
print(head(as.character(group_for_pca_bmiq_all), 5))

# Perform PCA on M-values (centered, not scaled - DESeq2 default)
pca_result_bmiq_all <- prcomp(t(mvals_for_pca_bmiq_all), center = TRUE, scale. = FALSE)

print("PCA completed successfully for BMIQ (ALL CpGs)!")
print(paste("PCA sample names (first 5):"))
print(head(rownames(pca_result_bmiq_all$x), 5))

# Create additional grouping vectors for PCA
genotype_for_pca_bmiq_all <- factor(matching_samples$True.Genotype)
cross_type_for_pca_bmiq_all <- factor(matching_samples$cross.type)

# Create PCA dataframe with proper alignment
pca_df_all_bmiq <- data.frame(
  PC1 = pca_result_bmiq_all$x[, 1],
  PC2 = pca_result_bmiq_all$x[, 2],
  Sample = rownames(pca_result_bmiq_all$x),
  Tissue_Cell_Type = group_for_pca_bmiq_all[match(rownames(pca_result_bmiq_all$x), colnames(betas_for_pca_bmiq_all))],
  Batch = batch_for_pca_bmiq_all[match(rownames(pca_result_bmiq_all$x), colnames(betas_for_pca_bmiq_all))],
  Genotype = genotype_for_pca_bmiq_all[match(rownames(pca_result_bmiq_all$x), colnames(betas_for_pca_bmiq_all))],
  Cross_Type = cross_type_for_pca_bmiq_all[match(rownames(pca_result_bmiq_all$x), colnames(betas_for_pca_bmiq_all))]
)

# Check for NAs in the PCA dataframe
print("Checking PCA dataframe for BMIQ (ALL CpGs):")
print(paste("Rows in pca_df:", nrow(pca_df_all_bmiq)))
print(paste("NAs in Tissue_Cell_Type:", sum(is.na(pca_df_all_bmiq$Tissue_Cell_Type))))
print(paste("NAs in Batch:", sum(is.na(pca_df_all_bmiq$Batch))))
print(paste("NAs in Genotype:", sum(is.na(pca_df_all_bmiq$Genotype))))
print(paste("NAs in Cross_Type:", sum(is.na(pca_df_all_bmiq$Cross_Type))))

# Plot by tissue cell type
pdf("PCA_BMIQ_tissue_type_ALL_CpGs.pdf", width = 10, height = 8)
p1_bmiq_all <- ggplot(pca_df_all_bmiq, aes(x = PC1, y = PC2, color = Tissue_Cell_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on BMIQ Normalized M-values (by Tissue Cell Type) - ALL CpGs",
    x = paste0("PC1 (", round(summary(pca_result_bmiq_all)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_bmiq_all)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p1_bmiq_all)
dev.off()

# Plot by batch
pdf("PCA_BMIQ_batch_ALL_CpGs.pdf", width = 10, height = 8)
p2_bmiq_all <- ggplot(pca_df_all_bmiq, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on BMIQ Normalized M-values (by Batch) - ALL CpGs",
    x = paste0("PC1 (", round(summary(pca_result_bmiq_all)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_bmiq_all)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p2_bmiq_all)
dev.off()

print("PCA plots (ALL CpGs - BMIQ data) saved successfully!")

print("=== PCA analysis (using ALL CpGs - QC data) ===")

# Use M-values for PCA (all CpGs)
betas_for_pca_qc_all <- betas_qc_combined
mvals_for_pca_qc_all <- mvals_qc
group_for_pca_qc_all <- group_vector
batch_for_pca_qc_all <- batch_vector

print("Debugging PCA setup for QC (ALL CpGs)...")
print(paste("Samples in data:", ncol(betas_for_pca_qc_all)))
print(paste("CpGs in data:", nrow(betas_for_pca_qc_all)))
print(paste("Length of group_for_pca:", length(group_for_pca_qc_all)))
print(paste("Length of batch_for_pca:", length(batch_for_pca_qc_all)))

print("Sample names in data (first 5):")
print(head(colnames(betas_for_pca_qc_all), 5))
print("Group vector (first 5):")
print(head(as.character(group_for_pca_qc_all), 5))

# Perform PCA on M-values (centered, not scaled - DESeq2 default)
pca_result_qc_all <- prcomp(t(mvals_for_pca_qc_all), center = TRUE, scale. = FALSE)

print("PCA completed successfully for QC (ALL CpGs)!")
print(paste("PCA sample names (first 5):"))
print(head(rownames(pca_result_qc_all$x), 5))

# Create additional grouping vectors for PCA
genotype_for_pca_qc_all <- factor(matching_samples$True.Genotype)
cross_type_for_pca_qc_all <- factor(matching_samples$cross.type)

# Create PCA dataframe with proper alignment
pca_df_all_qc <- data.frame(
  PC1 = pca_result_qc_all$x[, 1],
  PC2 = pca_result_qc_all$x[, 2],
  Sample = rownames(pca_result_qc_all$x),
  Tissue_Cell_Type = group_for_pca_qc_all[match(rownames(pca_result_qc_all$x), colnames(betas_for_pca_qc_all))],
  Batch = batch_for_pca_qc_all[match(rownames(pca_result_qc_all$x), colnames(betas_for_pca_qc_all))],
  Genotype = genotype_for_pca_qc_all[match(rownames(pca_result_qc_all$x), colnames(betas_for_pca_qc_all))],
  Cross_Type = cross_type_for_pca_qc_all[match(rownames(pca_result_qc_all$x), colnames(betas_for_pca_qc_all))]
)

# Check for NAs in the PCA dataframe
print("Checking PCA dataframe for QC (ALL CpGs):")
print(paste("Rows in pca_df:", nrow(pca_df_all_qc)))
print(paste("NAs in Tissue_Cell_Type:", sum(is.na(pca_df_all_qc$Tissue_Cell_Type))))
print(paste("NAs in Batch:", sum(is.na(pca_df_all_qc$Batch))))
print(paste("NAs in Genotype:", sum(is.na(pca_df_all_qc$Genotype))))
print(paste("NAs in Cross_Type:", sum(is.na(pca_df_all_qc$Cross_Type))))

# Plot by tissue cell type
pdf("PCA_QC_tissue_type_ALL_CpGs.pdf", width = 10, height = 8)
p1_qc_all <- ggplot(pca_df_all_qc, aes(x = PC1, y = PC2, color = Tissue_Cell_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on QC-filtered M-values (by Tissue Cell Type) - ALL CpGs",
    x = paste0("PC1 (", round(summary(pca_result_qc_all)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_qc_all)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p1_qc_all)
dev.off()

# Plot by batch
pdf("PCA_QC_batch_ALL_CpGs.pdf", width = 10, height = 8)
p2_qc_all <- ggplot(pca_df_all_qc, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA on QC-filtered M-values (by Batch) - ALL CpGs",
    x = paste0("PC1 (", round(summary(pca_result_qc_all)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_qc_all)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p2_qc_all)
dev.off()

print("PCA plots (ALL CpGs - QC data) saved successfully!")