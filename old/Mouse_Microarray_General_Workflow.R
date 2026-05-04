# load required packages
library(sesame)
library(minfi)
library(IlluminaMouseMethylationmanifest)
library(IlluminaMouseMethylationanno.12.v1.mm10)
library(limma)
library(sva)
library(ggsci)
library(wateRmelon)

# annotation
annoMouse <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)

# set up directory where all idat files are located
idir <- "/media/Data/Microarray/Mouse/Released_Data_Mouse/Data/"

# read in beta values and detection p-values
cores <- 50
betas <- openSesame(idir, BPPARAM = BiocParallel::MulticoreParam(workers=cores), collapseToPfx=TRUE)
detP <- openSesame(idir, func = pOOBAH, return.pval=TRUE, BPPARAM = BiocParallel::MulticoreParam(workers=cores))
detP <- betasCollapseToPfx(detP)
detP2 <- detP

# read in phenotypic data
targets <- read.table("/media/Data/Microarray/Mouse/Released_Data_Mouse/Data/samplesheet.csv",header=T)
targets$Basename <- paste0(targets$Plate,"_",targets$Position)
targets$Group2 <- paste0(targets$Tissue,"_",targets$Group)
targets2 <- targets

# get everything in order
betas <- betas[,match(targets$Basename, colnames(betas))]
betas2 <- betas

### Diagnostic plots

# Barplots of mean detection p-value
pal <- pal_simpsons()(12)
pdf("/media/Data/meanDP.pdf")
barplot(colMeans(detP), col=pal[factor(targets$Group2)], las=2,cex.names=0.4, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Group2)), fill=pal,bg="white")
dev.off()

# Density plot of beta distribution
pdf("/media/Data/density.pdf")
densityPlot(betas, sampGroups=targets$Group2,main="Beta Density", legend=FALSE, pal = pal)
legend("top", legend = levels(factor(targets$Group2)),text.col=pal)
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
pdf("/media/Data/pca.pdf")
betas3 <- betas
colnames(betas3) <- targets$Group2
plotMDS(betas3, top=1000, gene.selection="common",col=pal[factor(targets$Group2)], dim=c(1,2))
legend("right", legend=levels(factor(targets$Group2)), text.col=pal,cex=0.7, bg="white")
plotMDS(betas3, top=1000, gene.selection="common",col=pal[factor(targets$Group2)], dim=c(1,3))
legend("right", legend=levels(factor(targets$Group2)), text.col=pal,cex=0.7, bg="white")
plotMDS(betas3, top=1000, gene.selection="common",col=pal[factor(targets$Group2)], dim=c(2,3))
legend("right", legend=levels(factor(targets$Group2)), text.col=pal,cex=0.7, bg="white")
dev.off()

#betas.mac <- betas[,which(targets$Tissue == "Macrophage")]
#betas.mic <- betas[,which(targets$Tissue == "Microglia")]
#targets.mac <- targets[which(targets$Tissue == "Macrophage"), ]
#targets.mic <- targets[which(targets$Tissue == "Microglia"), ]
#colnames(betas.mac) <- targets.mac$Group2
#colnames(betas.mic) <- targets.mic$Group2
#pdf("pcaMac.pdf")
#plotMDS(betas.mac, top=1000, gene.selection="common",col=pal[factor(targets.mac$Group2)], dim=c(1,2))
#legend("right", legend=levels(factor(targets.mac$Group2)), text.col=pal,cex=0.7, bg="white")
#dev.off()
#pdf("pcaMic.pdf")
#plotMDS(betas.mic, top=1000, gene.selection="common",col=pal[factor(targets.mic$Group2)], dim=c(1,2))
#legend("right", legend=levels(factor(targets.mic$Group2)), text.col=pal,cex=0.7, bg="white")
#dev.off()

# remove SNPs and CpH probes
betas <- na.exclude(betas)
betas <- betas[grep("cg", rownames(betas)),]

#means <- colMeans(betas)
#plot(factor(targets$Group2), means, col = pal, xlab = "", ylab = "Mean beta", xaxt = "n")
#legend("bottomright", legend=levels(factor(targets$Group2)), text.col=pal,cex=0.7, bg="white")

save(betas, betas2, detP, detP2, targets2, targets, annoMouse, file = "/media/Data/Microarray/Mouse/Released_Data_Mouse/Data/processedBetas.rdata")

# General research questions:
# 1) Does malnutrition affect DNAm in macrophages and/or microglia? (LP vs Ctrl in each tissue)
# 2) Do recovery diets reverse these changes? (LP vs milk. LP vs wheat, LP vs peanut in each tissue)
# 3) Do microglia and macrophages differ? (Ctrl Mac vs Ctrl Micro, Ctrl2 Mac vs Ctrl2 Micro)

# very basic analysis to compare two main groups
targets$Group <- as.factor(targets$Group)
targets$Group2 <- as.factor(targets$Group2)
targets$Age <- gsub("wk","",targets$Age)
targets$Age <- as.numeric(targets$Age)
targets$Tissue <- as.factor(targets$Tissue)
#mod <- model.matrix(~ 0 + Group + Age + Tissue + Group:Tissue, data = targets) 
mod <- model.matrix(~ 0 + Group2, data = targets) 

# SVA
#mod0 <- model.matrix(~ 1, data = targets) 
#nCols <- ncol(mod)
#n.sv <- num.sv(betas,mod) # 16
#svobj <- sva(betas,mod,mod0,n.sv=n.sv)
#mod <- cbind(mod,svobj$sv)
#colnames(mod) <- c(colnames(mod)[1:nCols], paste0("sv", 1:n.sv))

fit <- lmFit(betas,mod)
conts <- c("Group2Macrophage_LP - Group2Macrophage_Ctrl", 
  "Group2Macrophage_LP - Group2Macrophage_Wheat", 
  "Group2Macrophage_LP - Group2Macrophage_Peanut", 
  "Group2Macrophage_LP - Group2Macrophage_Milk", 
  "Group2Microglia_LP - Group2Microglia_Ctrl", 
  "Group2Microglia_LP - Group2Microglia_Wheat", 
  "Group2Microglia_LP - Group2Microglia_Peanut", 
  "Group2Microglia_LP - Group2Microglia_Milk", 
  "Group2Macrophage_Ctrl - Group2Microglia_Ctrl", 
  "Group2Macrophage_Ctrl2 - Group2Microglia_Ctrl2",
  "Group2Macrophage_Ctrl2 - Group2Macrophage_Wheat",
  "Group2Macrophage_Ctrl2 - Group2Macrophage_Peanut",
  "Group2Macrophage_Ctrl2 - Group2Macrophage_Milk",
  "Group2Microglia_Ctrl2 - Group2Microglia_Wheat", 
  "Group2Microglia_Ctrl2 - Group2Microglia_Peanut", 
  "Group2Microglia_Ctrl2 - Group2Microglia_Milk") 

#contMatrix <- makeContrasts(Group2Macrophage_LP - Group2Macrophage_Ctrl, levels = mod)
contMatrix <- makeContrasts(contrasts = conts, levels = mod)
fit.cont <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit.cont)
#which(p.adjust(fit2$p.value, method = "fdr")<0.01)
#DMPs <- topTable(fit2, num = Inf, sort.by = 'none') # this gives all results, but pvalues only from first contrast
