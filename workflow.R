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
library(RColorBrewer)
library(readxl)
library(openxlsx)

cat("Library Paths:\n")
cat(.libPaths(), sep = "\n")

# Source the DoBMIQ function
source("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/bmiq/DoBMIQ.R")

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
targets <- read.csv("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/data/master_samplesheet.csv", header=TRUE)
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
  # Combine RUNs 1 & 2 and 3–5 (this is above here already)
  betas_combined <- cbind(
    betas_RUN1_2[common_cpgs_all, ],
    betas_RUN3_4_5[common_cpgs_all, ]
  )

  # === Rename columns of betas_combined to sample names from master_samplesheet ===
  print("Renaming columns...")
  samplesheet <- read.csv("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/data/master_samplesheet.csv", header=TRUE)
  samplesheet$Basename <- paste0(samplesheet$BEADCHIP, "_", samplesheet$POSITION)
  samplesheet$Samples <- paste0(samplesheet$tb, "_", samplesheet$RUN)

  # Match sample names based on column names in betas_combined
  matching_samples <- samplesheet[samplesheet$Basename %in% colnames(betas_combined), ]
  matching_samples <- matching_samples[match(colnames(betas_combined), matching_samples$Basename), ]

  # Check for mismatches
  if (!all(matching_samples$Basename == colnames(betas_combined))) {
    stop("Mismatch in sample order between betas_combined and master sample sheet.")
  }

  # Rename columns
  colnames(betas_combined) <- matching_samples$Samples

  # Continue with batch definition for ComBat
  batch <- c(rep(1, ncol(betas_RUN1_2)), rep(2, ncol(betas_RUN3_4_5)))
  print("Applying ComBat...")
  betas_combat <- ComBat(dat = betas_combined, batch = batch, mod = NULL)
  write.csv(betas_combat, "betas_combat_normalized.csv")

  
  # === BMIQ NORMALIZATION ===
  print("=== BMIQ NORMALIZATION ===")
  
  # Create a temporary CSV file with the ComBat-normalized data
  temp_beta_file <- "temp_betas_combat.csv"
  betas_combat_df <- data.frame(CpG_ID = rownames(betas_combat), betas_combat, stringsAsFactors = FALSE)
  write.csv(betas_combat_df, temp_beta_file, row.names = FALSE)
  
  # You'll need to specify the path to your probe annotation file
  probes_file <- "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/data/probesample.xlsx"  # UPDATE THIS PATH
  
 # Check if probe file exists
  if (file.exists(probes_file)) {
    print("Applying BMIQ normalization...")
    
    # Run BMIQ normalization with error handling
    tryCatch({
      DoBMIQ(probes_path = probes_file, beta_path = temp_beta_file)
      
      # Check if bmiq.Rd was created successfully
      if (file.exists("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/bmiq/results/bmiq.Rd")) {
        # Load the BMIQ-normalized data
        load("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/bmiq/results/bmiq.Rd")  # This loads bmiq.m
        print("Loaded bmiq.Rd")
        
        # Convert back to the same format as betas_combat
        betas_combat_bmiq <- bmiq.m
        
        # Save the BMIQ-normalized data
        write.csv(betas_combat_bmiq, "betas_combat_bmiq_normalized.csv")
        
        # Use BMIQ-normalized data for downstream analysis
        final_betas <- betas_combat_bmiq
        
        print("BMIQ normalization completed successfully!")
        
      } else {
        print("BMIQ normalization failed - bmiq.Rd file not created.")
        print("Using ComBat-normalized data for downstream analysis.")
        final_betas <- betas_combat
      }
    }, error = function(e) {
      print(paste("Error during BMIQ normalization:", e$message))
      print("Using ComBat-normalized data for downstream analysis.")
      final_betas <<- betas_combat
    })
    
    # Clean up temporary file
    if (file.exists(temp_beta_file)) {
      file.remove(temp_beta_file)
    }
    
  } else {
    print(paste("Probe annotation file not found:", probes_file))
    print("Using ComBat-normalized data for downstream analysis.")
    final_betas <- betas_combat
  }

  # === SELECT & SAVE TOP 10,000 CpGs ===
  print("Selecting top 10,000 most variable CpGs...")
  top_cpgs <- names(sort(apply(final_betas, 1, var), decreasing = TRUE))[1:10000]
  betas_final_top <- final_betas[top_cpgs, ]
  
  # Save with appropriate filename
  if (exists("betas_combat_bmiq")) {
    write.csv(betas_final_top, "betas_combat_bmiq_top10000.csv")
    print("Saved: betas_combat_bmiq_top10000.csv")
  } else {
    write.csv(betas_final_top, "betas_combat_top10000.csv")
    print("Saved: betas_combat_top10000.csv")
  }
  
} else {
  print("Cannot apply ComBat — missing combined datasets.")
}

# === PCA PLOTS ===

# Prepare color definitions
celltype_colors <- c("t-cells-DNA/RNA" = "#E31A1C", "THYMOCYTES" = "#1F78B4", "UNFRACTIONATED" = "#33A02C")
sample_celltype_colors <- celltype_colors[targets$cell.type]

targets$strain <- tolower(targets$strain)
targets$strain <- gsub("[()_]", "", targets$strain)
strain_names <- unique(targets$strain)
n_strains <- length(strain_names)
strain_palette <- c(brewer.pal(8, "Set1"), brewer.pal(3, "Set2"))[1:n_strains]
names(strain_palette) <- strain_names
sample_strain_colors <- strain_palette[targets$strain]

# Top 10,000 most variable CpGs
top_cpgs <- names(sort(apply(betas, 1, var), decreasing = TRUE))[1:10000]
betas_top <- betas[top_cpgs, ]

# Use all CpGs
betas_all <- final_betas

# === 1. PCA - Top 10,000 CpGs - Color by Cell Type ===
pdf("pca_top10000_by_celltype.pdf", width=10, height=8)
plotMDS(betas_top, top=10000, gene.selection="common", col=sample_celltype_colors, dim=c(1,2), pch=19, cex=1.2)
legend("topright", legend=names(celltype_colors), col=celltype_colors, pch=19, cex=0.9, bg="white", title="Cell Type")
dev.off()

# === 2. PCA - All CpGs - Color by Cell Type ===
pdf("pca_all_by_celltype.pdf", width=10, height=8)
plotMDS(betas_all, gene.selection="common", col=sample_celltype_colors, dim=c(1,2), pch=19, cex=1.2)
legend("topright", legend=names(celltype_colors), col=celltype_colors, pch=19, cex=0.9, bg="white", title="Cell Type")
dev.off()

# === 3. PCA - Top 10,000 CpGs - Color by Strain ===
pdf("pca_top10000_by_strain.pdf", width=12, height=8)
plotMDS(betas_top, top=10000, gene.selection="common", col=sample_strain_colors, dim=c(1,2), pch=19, cex=1.2)
if (n_strains <= 6) {
  legend("topright", legend=strain_names, col=strain_palette, pch=19, cex=0.8, bg="white", title="Strain")
} else {
  mid_point <- ceiling(n_strains/2)
  legend("topright", legend=strain_names[1:mid_point], col=strain_palette[1:mid_point], pch=19, cex=0.7, bg="white", title="Strain")
  legend("bottomright", legend=strain_names[(mid_point+1):n_strains], col=strain_palette[(mid_point+1):n_strains], pch=19, cex=0.7, bg="white")
}
dev.off()

# === 4. PCA - All CpGs - Color by Strain ===
pdf("pca_all_by_strain.pdf", width=12, height=8)
plotMDS(betas_all, gene.selection="common", col=sample_strain_colors, dim=c(1,2), pch=19, cex=1.2)
if (n_strains <= 6) {
  legend("topright", legend=strain_names, col=strain_palette, pch=19, cex=0.8, bg="white", title="Strain")
} else {
  mid_point <- ceiling(n_strains/2)
  legend("topright", legend=strain_names[1:mid_point], col=strain_palette[1:mid_point], pch=19, cex=0.7, bg="white", title="Strain")
  legend("bottomright", legend=strain_names[(mid_point+1):n_strains], col=strain_palette[(mid_point+1):n_strains], pch=19, cex=0.7, bg="white")
}

# === 5. PCA - All CpGs - Labels Only (No Color) ===
pdf("pca_all_labeled.pdf", width=12, height=10)
mds <- plotMDS(betas_all, gene.selection="common", dim=c(1,2), plot=FALSE)

plot(mds$x, mds$y,
     xlab=paste("Dimension 1 (", round(mds$var.explained[1], 1), "%)", sep=""),
     ylab=paste("Dimension 2 (", round(mds$var.explained[2], 1), "%)", sep=""),
     main="PCA: All CpGs with Sample Labels",
     pch=19, col="black", cex=0.8)

text(mds$x, mds$y, labels=targets$Samples, pos=3, cex=0.7)
dev.off()
