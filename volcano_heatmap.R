################ 
### Volcano Plot
################ 

# clarifiers
# results = data.frame with significance (FDR) and effect size (diff)methyl) columns

# split into hyper/hypo and non-significant sites
sigUp <- results[which(results$FDR < 0.05 & results$diff_methyl > 0),]
sigDown <- results[which(results$FDR < 0.05 & results$diff_methyl < 0),]
notSig <- results[which(results$FDR > 0.05),]

# add colors
sigDown$col <- "deepskyblue3"
sigUp$col <- "indianred"
notSig$col <- "ivory3"

# bring it all back home
toPlot <- rbind(sigUp,sigDown,notSig)

# plot it out
plot(toPlot$diff_methyl, -log10(toPlot$FDR),
	ylim = c(0,10),
	xlim = c(-1,1),
	col = toPlot$col,
	pch=16,
	ylab = "-log10(padj)",
	xlab = "Beta (KO - WT)")
abline(h = -log10(0.05), lty = 3)

################ 
### Heatmap
################

# # clarifiers
# # bValues = data.frame/matrix of beta-values
# # results = data.frame with test results and CpG IDs (e.g., cg123456) as rownames

# # load package
# library(pheatmap)

# # filter for DMPs from results (FDR < 0.05)
# DMPs <- results[which(results$FDR < 0.05), ]
  
# # filter for beta-values of DMPs
# DMPs.beta <- bValues[match(rownames(DMPs), rownames(bValues)),]
# pheatmap(as.matrix(DMPs.beta), cluster_rows = TRUE, show_rownames = FALSE, cluster_cols = TRUE, treeheight_col = 1, scale="row", treeheight_row = 0)
