# Load & preprocess methylation IDAT files using SeSAMe
# Filters probes with detection p-value > 0.05 using pOOBAH

# 1. Install sesame if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sesame", ask = FALSE, update = FALSE)

# 2. Load library
library(sesame)

# 3. Cache required annotation data (e.g., for pOOBAH)
sesameDataCache()

# 4. Set paths
idat_dir <- "new_idat"
output_file <- "filtered_beta_matrix.csv"

# 5. Load & preprocess samples
#    Returns a list of SigDFs (one per sample)
sdfs <- openSesame(idat_dir, func = NULL)  # No masking yet

# 6. Apply pOOBAH detection mask and extract β-values
beta_list <- lapply(sdfs, function(sdf) {
    masked_sdf <- pOOBAH(sdf, return.pval = FALSE, pval.threshold = 0.05)
    getBetas(masked_sdf)
})

# 7. Combine into matrix
beta_matrix <- do.call(cbind, beta_list)
colnames(beta_matrix) <- names(beta_list)

# 8. Remove rows with all NA (filtered out)
beta_matrix <- beta_matrix[rowSums(is.na(beta_matrix)) < ncol(beta_matrix), ]

# 9. Save result
write.csv(beta_matrix, file = output_file, quote = FALSE)

cat("Done. Saved to:", output_file, "\n")
