.libPaths("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/packages")
setwd("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline")

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

# Annotation
annoMouse <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)

# Read raw data
print("Reading in IDATs...")
idir <- "idat"
cores <- 4
betas <- openSesame(idir, BPPARAM = BiocParallel::MulticoreParam(workers=cores), collapseToPfx=TRUE)
detP <- openSesame(idir, func = pOOBAH, return.pval=TRUE, BPPARAM = BiocParallel::MulticoreParam(workers=cores))
detP <- betasCollapseToPfx(detP)

# write.csv(betas, "combined_betas.csv", row.names = TRUE)
# write.csv(detP, "combined_detP.csv", row.names = TRUE)

# Read metadata
targets <- read.csv("master_samplesheet.csv", header=TRUE)
targets$Basename <- paste0(targets$BEADCHIP, "_", targets$POSITION)
targets$Samples <- paste0(targets$tb, "_", targets$RUN)
write.csv(targets, "targets_pre.csv", row.names = FALSE)

# Align
targets$BEADCHIP <- as.character(as.numeric(targets$BEADCHIP))
targets <- targets[!duplicated(paste0(targets$BEADCHIP, "_", targets$POSITION)), ]
targets$Basename <- paste0(targets$BEADCHIP, "_", targets$POSITION)
common_samples <- intersect(targets$Basename, colnames(betas))
targets <- targets[targets$Basename %in% common_samples, ]
betas <- betas[, match(targets$Basename, colnames(betas))]
detP <- detP[, match(targets$Basename, colnames(detP))]

# write.csv(betas, "betas.csv")
# write.csv(detP, "detP.csv")
# write.csv(targets, "targets.csv")

# QC Plots
pal <- pal_simpsons()(12)
pdf("meanDP.pdf")
barplot(colMeans(detP), col=pal[factor(targets$Samples)], las=2, cex.names=0.4, ylab="Mean detection p-values")
abline(h=0.05, col="red")
legend("topleft", legend=levels(factor(targets$Samples)), fill=pal, bg="white")
dev.off()

pdf("density.pdf")
densityPlot(betas, sampGroups=targets$Samples, main="Beta Density", legend=FALSE, pal=pal)
legend("top", legend=levels(factor(targets$Samples)), text.col=pal)
dev.off()

# Filter by RUN
print("=== FILTERING BY RUN ===")
unique_runs <- sort(unique(targets$RUN))
betas_by_run <- list()
detP_by_run <- list()
targets_by_run <- list()

for(run_id in unique_runs) {
  idx <- targets$RUN == run_id
  targets_run <- targets[idx, ]
  betas_run <- betas[, idx]
  detP_run <- detP[, idx]
  
  # Sample filter
  good_samples <- which(colMeans(detP_run) < 0.05)
  betas_run <- betas_run[, good_samples]
  detP_run <- detP_run[, good_samples]
  targets_run <- targets_run[good_samples, ]
  
  # CpG filter
  keep <- rowSums(detP_run < 0.05) == ncol(betas_run)
  betas_run <- betas_run[keep, ]
  detP_run <- detP_run[keep, ]
  
  # Remove NA and non-CpGs
  betas_run <- na.exclude(betas_run)
  betas_run <- betas_run[grep("cg", rownames(betas_run)), ]
  detP_run <- detP_run[rownames(betas_run), ]
  
  # Remove sex CpGs
  sex_cpgs <- annoMouse$Name[annoMouse$chr %in% c("chrX", "chrY")]
  betas_run <- betas_run[!rownames(betas_run) %in% sex_cpgs, ]
  detP_run <- detP_run[rownames(betas_run), ]
  
  # Store
  betas_by_run[[paste0("RUN", run_id)]] <- betas_run
  detP_by_run[[paste0("RUN", run_id)]] <- detP_run
  targets_by_run[[paste0("RUN", run_id)]] <- targets_run
}

# Combine RUNs 1 & 2
if ("RUN1" %in% names(betas_by_run) && "RUN2" %in% names(betas_by_run)) {
  common_cpgs_12 <- intersect(rownames(betas_by_run$RUN1), rownames(betas_by_run$RUN2))
  betas_RUN1_2 <- cbind(betas_by_run$RUN1[common_cpgs_12, ], betas_by_run$RUN2[common_cpgs_12, ])
  detP_RUN1_2 <- cbind(detP_by_run$RUN1[common_cpgs_12, ], detP_by_run$RUN2[common_cpgs_12, ])
  write.csv(betas_RUN1_2, "betas_RUN1_2.csv")
}

# Combine RUNs 3, 4 & 5
run_names_345 <- paste0("RUN", 3:5)
available_345 <- run_names_345[run_names_345 %in% names(betas_by_run)]
if (length(available_345) > 1) {
  common_cpgs_345 <- Reduce(intersect, lapply(available_345, function(x) rownames(betas_by_run[[x]])))
  betas_list_345 <- lapply(available_345, function(x) betas_by_run[[x]][common_cpgs_345, ])
  betas_RUN3_4_5 <- do.call(cbind, betas_list_345)
  write.csv(betas_RUN3_4_5, "betas_RUN3_4_5.csv")
}

# === COMBAT NORMALIZATION ===
print("=== COMBAT NORMALIZATION ===")
if (exists("betas_RUN1_2") && exists("betas_RUN3_4_5")) {
  common_cpgs_all <- intersect(rownames(betas_RUN1_2), rownames(betas_RUN3_4_5))
  betas_combined <- cbind(
    betas_RUN1_2[common_cpgs_all, ],
    betas_RUN3_4_5[common_cpgs_all, ]
  )
  
  batch <- c(rep(1, ncol(betas_RUN1_2)), rep(2, ncol(betas_RUN3_4_5)))
  print("Applying ComBat...")
  betas_combat <- ComBat(dat = betas_combined, batch = batch, mod = NULL)
  write.csv(betas_combat, "betas_combat_normalized.csv")
  
  # === SELECT & SAVE TOP 10,000 CpGs ===
  print("Selecting top 10,000 most variable CpGs...")
  top_cpgs <- names(sort(apply(betas_combat, 1, var), decreasing = TRUE))[1:10000]
  betas_combat_top <- betas_combat[top_cpgs, ]
  write.csv(betas_combat_top, "betas_combat_top10000.csv")
  
  print("Saved: betas_combat_top10000.csv")
  
} else {
  print("Cannot apply ComBat — missing combined datasets.")
}

# PCA using top 10k
pdf("pca.pdf")
top_cpgs <- names(sort(apply(betas, 1, var), decreasing = TRUE))[1:10000]
betas3 <- betas[top_cpgs, ]
colnames(betas3) <- targets$Samples
plotMDS(betas3, top=10000, gene.selection="common", col=pal[factor(targets$Samples)], dim=c(1,2))
legend("right", legend=levels(factor(targets$Samples)), text.col=pal, cex=0.7, bg="white")
dev.off()
