# .libPaths("C:/Users/Elliott/Documents/github/methyl_data_pipeline/packages") # Elliott PC

.libPaths("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/packages") # Elliott Mac
setwd("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline") # Elliott Mac

# load required packages
library(sesame)
library(minfi)
library(IlluminaMouseMethylationmanifest)
library(IlluminaMouseMethylationanno.12.v1.mm10)
library(limma)
library(sva)
library(ggsci)
library(wateRmelon)
library(snow)
library(BiocParallel)
library(parallel)
cat("Library Paths:\n")
cat(.libPaths(), sep = "\n")

# annotation
annoMouse <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)

# set up directory where all idat files are located
print("Reading in data...")
idir <- "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/idat"

# read in beta values and detection p-values
print("Processing idat files...")
cores <- 4
betas <- openSesame(idir, BPPARAM = BiocParallel::MulticoreParam(workers=cores), collapseToPfx=TRUE) # using SnowParam instead
detP <- openSesame(idir, func = pOOBAH, return.pval=TRUE, BPPARAM = BiocParallel::MulticoreParam(workers=cores))
detP <- betasCollapseToPfx(detP)
detP2 <- detP

# Save combined data to CSV files before reading phenotypic data
print("Saving combined data to CSV files...")

# Save beta values to CSV
write.csv(betas, "combined_betas.csv", row.names = TRUE)
print("Beta values saved to combined_betas.csv")

# Save detection p-values to CSV
write.csv(detP, "combined_detP.csv", row.names = TRUE)
print("Detection p-values saved to combined_detP.csv")

# read in phenotypic data
print("Reading in sample sheet...")
targets <- read.csv("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/master_samplesheet.csv",header=T)
targets$Basename <- paste0(targets$BEADCHIP, "_", targets$POSITION)
targets$Samples <- paste0(targets$tb, "_", targets$RUN)

write.csv(targets, "targets_pre.csv", row.names = FALSE)
print("Targets saved to targets.csv")

# targets$Basename <- paste0(targets$BEADCHIP, "_", targets$POSITION)
# targets$Samples <- paste0(targets$tb, "_", targets$Experiment)

targets2 <- targets

# Fix samplesheet and align with data
targets$BEADCHIP <- as.character(as.numeric(targets$BEADCHIP))  # Remove scientific notation
targets <- targets[!duplicated(paste0(targets$BEADCHIP, "_", targets$POSITION)), ]  # Remove duplicates
targets$Basename <- paste0(targets$BEADCHIP, "_", targets$POSITION)  # Add X prefix

# Keep only common samples and get everything in order
common_samples <- intersect(targets$Basename, colnames(betas))
targets <- targets[targets$Basename %in% common_samples, ]
betas <- betas[, match(targets$Basename, colnames(betas))]
detP <- detP[, match(targets$Basename, colnames(detP))]
betas2 <- betas

# Save combined data to CSV files before reading phenotypic data
print("Saving combined data to CSV files...")

# Save beta values to CSV
write.csv(betas, "betas.csv", row.names = TRUE)
print("Beta values saved to betas.csv")

# Save detection p-values to CSV
write.csv(detP, "detP.csv", row.names = TRUE)
print("Detection p-values saved to detP.csv")

write.csv(targets, "targets.csv", row.names = FALSE)
print("Targets saved to targets.csv")

### Diagnostic plots

# Barplots of mean detection p-value
print("plotting mean detection p-values...")
pal <- pal_simpsons()(12)
pdf("meanDP.pdf")
barplot(colMeans(detP), col=pal[factor(targets$Samples)], las=2,cex.names=0.4, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Samples)), fill=pal,bg="white")
dev.off()

# Density plot of beta distribution
print("plotting beta distribution...")
pdf("density.pdf")
densityPlot(betas, sampGroups=targets$Samples,main="Beta Density", legend=FALSE, pal = pal)
legend("top", legend = levels(factor(targets$Samples)),text.col=pal)
dev.off()

print("=== PROCESSING BY RUN ===")

# Split and filter each RUN separately
unique_runs <- sort(unique(targets$RUN))
betas_by_run <- list()
detP_by_run <- list()
targets_by_run <- list()

for(run_id in unique_runs) {
  # Split by RUN
  run_idx <- targets$RUN == run_id
  targets_run <- targets[run_idx, ]
  betas_run <- betas[, run_idx]
  detP_run <- detP[, run_idx]
  
  # Apply your original filtering to each RUN
  # Filter samples if their mean detection p-value is > 0.05
  good_samples <- which(colMeans(detP_run) < 0.05)
  betas_run <- betas_run[, good_samples]
  targets_run <- targets_run[good_samples, ]
  detP_run <- detP_run[, good_samples]
  
  # Filter CpGs if at most 1 sample has detection p-value > 0.05
  keep <- rowSums(detP_run < 0.05) == ncol(betas_run)
  betas_run <- betas_run[keep, ]
  detP_run <- detP_run[keep, ]
  
  # Filter CpGs for SNPs and other nonsense
  keep <- grep("cg", rownames(betas_run))
  betas_run <- betas_run[keep, ]
  detP_run <- detP_run[keep, ]
  
  # Remove SNPs and CpH probes (additional filtering)
  betas_run <- na.exclude(betas_run)
  betas_run <- betas_run[grep("cg", rownames(betas_run)), ]
  detP_run <- detP_run[rownames(betas_run), ]
  
  # Store filtered data
  betas_by_run[[paste0("RUN", run_id)]] <- betas_run
  detP_by_run[[paste0("RUN", run_id)]] <- detP_run
  targets_by_run[[paste0("RUN", run_id)]] <- targets_run
}

# Combine RUNs 1 & 2
if("RUN1" %in% names(betas_by_run) && "RUN2" %in% names(betas_by_run)) {
  common_cpgs_12 <- intersect(rownames(betas_by_run$RUN1), rownames(betas_by_run$RUN2))
  betas_RUN1_2 <- cbind(betas_by_run$RUN1[common_cpgs_12, ], betas_by_run$RUN2[common_cpgs_12, ])
  detP_RUN1_2 <- cbind(detP_by_run$RUN1[common_cpgs_12, ], detP_by_run$RUN2[common_cpgs_12, ])
  write.csv(betas_RUN1_2, "betas_RUN1_2.csv")
}

# Combine RUNs 3, 4 & 5
run_names_345 <- paste0("RUN", 3:5)
available_345 <- run_names_345[run_names_345 %in% names(betas_by_run)]
if(length(available_345) > 1) {
  common_cpgs_345 <- Reduce(intersect, lapply(available_345, function(x) rownames(betas_by_run[[x]])))
  betas_list_345 <- lapply(available_345, function(x) betas_by_run[[x]][common_cpgs_345, ])
  betas_RUN3_4_5 <- do.call(cbind, betas_list_345)
  write.csv(betas_RUN3_4_5, "betas_RUN3_4_5.csv")
}

# PCA
pdf("pca.pdf")
betas3 <- betas
colnames(betas3) <- targets$Samples
plotMDS(betas3, top=1000, gene.selection="common",col=pal[factor(targets$Samples)], dim=c(1,2))
legend("right", legend=levels(factor(targets$Samples)), text.col=pal,cex=0.7, bg="white")
plotMDS(betas3, top=1000, gene.selection="common",col=pal[factor(targets$Samples)], dim=c(1,3))
legend("right", legend=levels(factor(targets$Samples)), text.col=pal,cex=0.7, bg="white")
plotMDS(betas3, top=1000, gene.selection="common",col=pal[factor(targets$Samples)], dim=c(2,3))
legend("right", legend=levels(factor(targets$Samples)), text.col=pal,cex=0.7, bg="white")
dev.off()

save(betas, betas2, detP, detP2, targets2, targets, annoMouse, file = "processedBetas.rdata")

# ComBat normalization to combine RUNs 1,2 + RUNs 3,4,5
print("=== COMBAT NORMALIZATION ===")

if(exists("betas_RUN1_2") && exists("betas_RUN3_4_5")) {
  # Find common CpGs between the two combined datasets
  common_cpgs_all <- intersect(rownames(betas_RUN1_2), rownames(betas_RUN3_4_5))
  print(paste("Common CpGs across all RUNs:", length(common_cpgs_all)))
  
  # Combine the datasets
  betas_combined_all <- cbind(betas_RUN1_2[common_cpgs_all, ], betas_RUN3_4_5[common_cpgs_all, ])
  
  # Create batch variable (1 for RUNs 1&2, 2 for RUNs 3,4&5)
  batch <- c(rep(1, ncol(betas_RUN1_2)), rep(2, ncol(betas_RUN3_4_5)))
  
  # Apply ComBat normalization
  print("Applying ComBat normalization...")
  betas_combat <- ComBat(dat = betas_combined_all, batch = batch, mod = NULL)
  
  # Save ComBat normalized data
  write.csv(betas_combat, "betas_combat_normalized.csv")
  print("ComBat normalized data saved to betas_combat_normalized.csv")
  
  print(paste("Final normalized dataset:", nrow(betas_bmiq$nbeta), "CpGs x", ncol(betas_bmiq$nbeta), "samples"))
} else {
  print("Cannot perform ComBat normalization - missing combined RUN datasets")
}


