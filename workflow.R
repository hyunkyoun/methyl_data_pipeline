.libPaths("C:/Users/Elliott/Documents/github/methyl_data_pipeline/packages")

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
idir <- "C:/Users/Elliott/Documents/github/methyl_data_pipeline/idat"

# read in beta values and detection p-values
print("Processing idat files...")
cores <- 1
betas <- openSesame(idir, BPPARAM = SnowParam(workers=cores), collapseToPfx=TRUE) # using SnowParam instead
detP <- openSesame(idir, func = pOOBAH, return.pval=TRUE, BPPARAM = SnowParam(workers=cores))
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
targets <- read.csv("C:/Users/Elliott/Documents/github/methyl_data_pipeline/samplesheet.csv",header=T)
targets$Basename <- paste0(targets$BEADCHIP, "_", targets$POSITION)
targets$Samples <- paste0(targets$tb, "_", targets$RUN)

# targets$Basename <- paste0(targets$BEADCHIP, "_", targets$POSITION)
# targets$Samples <- paste0(targets$tb, "_", targets$Experiment)

targets2 <- targets

# get everything in order
betas <- betas[,match(targets$Basename, colnames(betas))]
betas2 <- betas

### Diagnostic plots

# Barplots of mean detection p-value
pal <- pal_simpsons()(12)
pdf("meanDP.pdf")
barplot(colMeans(detP), col=pal[factor(targets$Samples)], las=2,cex.names=0.4, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Samples)), fill=pal,bg="white")
dev.off()

# Density plot of beta distribution
pdf("density.pdf")
densityPlot(betas, sampGroups=targets$Samples,main="Beta Density", legend=FALSE, pal = pal)
legend("top", legend = levels(factor(targets$Samples)),text.col=pal)
dev.off()

# Filter samples if their mean detection p-value is > 0.05
betas <- betas[,which(colMeans(detP)<0.05)]
targets <- targets[which(colMeans(detP)<0.05),]
detP <- detP[,which(colMeans(detP)<0.05)]

# Filter CpGs if at most 1 sample has detection p-value > 0.05
keep <- rowSums(detP < 0.05) == ncol(betas)
#keep <- rowSums(detP < 0.01) == ncol(betas)
betas <- betas[keep,]

# Filter CpGs for SNPs and other nonsense
keep <- grep("cg", rownames(betas))
betas <- betas[keep,]

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

# remove SNPs and CpH probes
betas <- na.exclude(betas)
betas <- betas[grep("cg", rownames(betas)),]

save(betas, betas2, detP, detP2, targets2, targets, annoMouse, file = "processedBetas.rdata")
