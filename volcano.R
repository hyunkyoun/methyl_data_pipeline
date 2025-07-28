# ---- User input for groups

# Prompt user for file paths
cat("Separate multiple file paths with commas\n")
group1_input <- readline(prompt = "Enter path(s) for Group 1 CSV(s): ")
group2_input <- readline(prompt = "Enter path(s) for Group 2 CSV(s): ")
save_dir <- readline(prompt = "Enter directory to save plots (leave blank to skip saving): ")

# Parse file paths
group1_files <- trimws(unlist(strsplit(group1_input, ",")))
group2_files <- trimws(unlist(strsplit(group2_input, ",")))
save_dir <- trimws(save_dir)

# Remove empty strings
group1_files <- group1_files[group1_files != ""]
group2_files <- group2_files[group2_files != ""]

# Validate inputs
if(length(group1_files) == 0 | length(group2_files) == 0) {
  stop("No valid files provided for one or both groups")
}

cat("Group 1 files:", paste(group1_files, collapse = ", "), "\n")
cat("Group 2 files:", paste(group2_files, collapse = ", "), "\n")

# ---- Loading and combining data

# Function to load and merge files
load_and_merge <- function(files) {
  dfs <- list()
  for(i in seq_along(files)) {
    file_path <- files[i]
    if(file.exists(file_path)) {
      df <- read.csv(file_path, row.names = 1, check.names = FALSE)
      dfs[[i]] <- df
      cat("Loaded:", file_path, "- Shape:", nrow(df), "x", ncol(df), "\n")
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  
  if(length(dfs) == 0) {
    stop("No valid files could be loaded")
  }
  
  # Combine all dataframes by columns
  combined_df <- do.call(cbind, dfs)
  return(combined_df)
}

cat("Loading data...\n")
group1_data <- load_and_merge(group1_files)
group2_data <- load_and_merge(group2_files)

# ---- Data alignment

# Find common CpGs between groups
common_cpgs <- intersect(rownames(group1_data), rownames(group2_data))
group1_data <- group1_data[common_cpgs, , drop = FALSE]
group2_data <- group2_data[common_cpgs, , drop = FALSE]

cat("Group 1:", ncol(group1_data), "samples\n")
cat("Group 2:", ncol(group2_data), "samples\n") 
cat("Common CpGs:", length(common_cpgs), "\n")

# ---- Differential Methylation Analysis with limma

cat("Performing differential methylation analysis...\n")

# Load required packages
if(!require(limma, quietly = TRUE)) {
  install.packages("limma")
  library(limma)
}

# Combine data for analysis
combined_data <- cbind(group1_data, group2_data)
cat("Combined dataset:", nrow(combined_data), "CpGs x", ncol(combined_data), "samples\n")

# Create group labels
group_labels <- c(rep("Group1", ncol(group1_data)), rep("Group2", ncol(group2_data)))
cat("Group1 samples:", sum(group_labels == "Group1"), "\n")
cat("Group2 samples:", sum(group_labels == "Group2"), "\n")

# Clamp beta values to avoid infinite M-values
beta_clamped <- pmax(pmin(combined_data, 0.999), 0.001)

# Convert beta values to M-values for proper statistical analysis
mvals <- log2(beta_clamped / (1 - beta_clamped))
cat("Converted to M-values for statistical analysis\n")

# Create design matrix for limma
group_factor <- factor(group_labels, levels = c("Group1", "Group2"))
design <- model.matrix(~ group_factor)
colnames(design) <- c("Intercept", "Group2_vs_Group1")

cat("Design matrix created:\n")
print(head(design))

# Fit linear model
cat("Fitting linear models...\n")
fit <- lmFit(mvals, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract results for Group2 vs Group1 comparison
limma_results <- topTable(fit, coef = "Group2_vs_Group1", number = Inf, sort.by = "none")

# Calculate actual beta value differences
cat("Calculating beta value differences...\n")
group1_means <- rowMeans(group1_data[common_cpgs, ], na.rm = TRUE)
group2_means <- rowMeans(group2_data[common_cpgs, ], na.rm = TRUE)
beta_diff <- group2_means - group1_means  # Group2 (KO) - Group1 (WT)

# Create results dataframe with beta differences and p-values
results <- data.frame(
  diff_methyl = beta_diff,  # Actual beta differences
  FDR = limma_results$P.Value,
  row.names = rownames(limma_results)
)

# Summary statistics
cat("\n=== DIFFERENTIAL METHYLATION RESULTS ===\n")
cat("Total CpGs analyzed:", nrow(results), "\n")
cat("Significant CpGs (p-value < 0.05):", sum(results$FDR < 0.05), "\n")
cat("CpGs with diff_methyl > 0 and p-value < 0.05:", 
    sum(results$FDR < 0.05 & results$diff_methyl > 0), "\n")
cat("CpGs with diff_methyl < 0 and p-value < 0.05:", 
    sum(results$FDR < 0.05 & results$diff_methyl < 0), "\n")

# Beta difference range
cat("Beta difference range:", round(min(results$diff_methyl, na.rm = TRUE), 3), 
    "to", round(max(results$diff_methyl, na.rm = TRUE), 3), "\n")

# Show top differentially methylated CpGs
cat("\nTop 10 differentially methylated CpGs:\n")
top_dmps <- results[order(results$FDR)[1:10], c("diff_methyl", "FDR")]
print(round(top_dmps, 4))

cat("\nDifferential analysis complete! Results stored in 'results' dataframe.\n")

# Create combined beta values for heatmap
bValues <- cbind(group1_data, group2_data)

# ---- Volcano plots

# clarifiers
# results = data.frame with significance (FDR) and effect size (diff_methyl) columns
# diff_methyl now contains actual beta value differences

# split into hyper/hypo and non-significant sites
sigUp <- results[which(results$FDR < 0.05 & results$diff_methyl > 0),]
sigDown <- results[which(results$FDR < 0.05 & results$diff_methyl < 0),]
notSig <- results[which(results$FDR > 0.05),]

# add colors
if(nrow(sigDown) > 0) sigDown$col <- "deepskyblue3"
if(nrow(sigUp) > 0) sigUp$col <- "indianred"
if(nrow(notSig) > 0) notSig$col <- "ivory3"

# bring it all back home
toPlot <- rbind(sigUp,sigDown,notSig)

# plot it out with corrected axis limits for beta differences
plot(toPlot$diff_methyl, -log10(toPlot$FDR),
	ylim = c(0,10),
	xlim = c(-1,1),  # Appropriate range for beta differences
	col = toPlot$col,
	pch=16,
	ylab = "-log10(p-value)",
	xlab = "Beta Difference (Group2 - Group1)")
abline(h = -log10(0.05), lty = 3)
abline(v = 0, lty = 2, col = "gray")  # Add reference line at 0

# Add summary text to plot
legend("topright", 
       legend = c(paste("Hypermethylated:", nrow(sigUp)),
                  paste("Hypomethylated:", nrow(sigDown)),
                  paste("Not significant:", nrow(notSig))),
       col = c("indianred", "deepskyblue3", "ivory3"),
       pch = 16,
       cex = 0.8)

# ---- Heat maps

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

# ---- Summary

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results and plots have been generated!\n")
cat("- Results dataframe: 'results' (with FDR and beta difference columns)\n")
cat("- Combined beta values: 'bValues'\n")
cat("- Beta differences range from", round(min(results$diff_methyl, na.rm = TRUE), 3), 
    "to", round(max(results$diff_methyl, na.rm = TRUE), 3), "\n")