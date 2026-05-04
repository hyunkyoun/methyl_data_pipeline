# Load & preprocess methylation IDAT files using SeSAMe
# Filters probes with detection p-value > 0.05 using pOOBAH
# Processes samples by runs, then combines results

# 1. Install sesame if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sesame", ask = FALSE, update = FALSE)

# 2. Load library
library(sesame)

# 3. Cache required annotation data (e.g., for pOOBAH)
sesameDataCache()

# 4. Set paths
idat_dir <- "idat"
output_file <- "filtered_beta_matrix.csv"

# 5. Define run patterns based on IDAT file prefixes
run_patterns <- list(
    RUN_1 = c("205243950018", "205243950019", "205243950041", "205243950045"),
    RUN_2 = c("205730000016", "205730000017", "205730000018", "205730000029"),
    RUN_3 = c("20578610008"), # Add RUN_3 prefixes here
    RUN_4 = c("208319490008", "208319490009"), # Add RUN_4 prefixes here
    RUN_5 = c("209160440001", "209171910075")  # Add RUN_5 prefixes here
)

# 6. Function to get sample names for a specific run
get_run_samples <- function(run_prefixes, idat_directory) {
    if (length(run_prefixes) == 0) {
        return(character(0))
    }
    
    # Get all IDAT files in directory
    all_files <- list.files(idat_directory, pattern = "\\.idat$", full.names = FALSE)
    
    # Extract sample names (remove _Red.idat or _Grn.idat suffixes)
    sample_names <- unique(gsub("_(Red|Grn)\\.idat$", "", all_files))
    
    # Filter samples that match run prefixes
    run_samples <- sample_names[sapply(sample_names, function(x) {
        any(sapply(run_prefixes, function(prefix) startsWith(x, prefix)))
    })]
    
    return(run_samples)
}

# 7. Process each run separately
all_beta_matrices <- list()

for (run_name in names(run_patterns)) {
    cat("Processing", run_name, "...\n")
    
    # Get samples for this run
    run_samples <- get_run_samples(run_patterns[[run_name]], idat_dir)
    
    if (length(run_samples) == 0) {
        cat("No samples found for", run_name, ". Skipping...\n")
        next
    }
    
    cat("Found", length(run_samples), "samples for", run_name, "\n")
    
    # Load samples for this run
    run_sdfs <- openSesame(idat_dir, func = NULL, samples = run_samples)
    
    # Apply pOOBAH detection mask and extract β-values
    beta_list <- lapply(run_sdfs, function(sdf) {
        masked_sdf <- pOOBAH(sdf, return.pval = FALSE, pval.threshold = 0.05)
        getBetas(masked_sdf)
    })
    
    # Combine into matrix for this run
    run_beta_matrix <- do.call(cbind, beta_list)
    colnames(run_beta_matrix) <- names(beta_list)
    
    # Remove rows with all NA (filtered out) for this run
    run_beta_matrix <- run_beta_matrix[rowSums(is.na(run_beta_matrix)) < ncol(run_beta_matrix), ]
    
    # Save individual run file
    run_output_file <- paste0(run_name, ".csv")
    write.csv(run_beta_matrix, file = run_output_file, quote = FALSE)
    cat("Saved", run_name, "to:", run_output_file, "\n")
    
    # Store for final combination
    all_beta_matrices[[run_name]] <- run_beta_matrix
}

# 8. Combine all runs into one matrix
if (length(all_beta_matrices) > 0) {
    cat("Combining all runs...\n")
    
    # Get all unique probe names across runs
    all_probes <- unique(unlist(lapply(all_beta_matrices, rownames)))
    
    # Create a master matrix with all probes
    combined_matrix <- matrix(NA, nrow = length(all_probes), ncol = 0)
    rownames(combined_matrix) <- all_probes
    
    # Add each run's data to the master matrix
    for (run_name in names(all_beta_matrices)) {
        run_matrix <- all_beta_matrices[[run_name]]
        
        # Create a temporary matrix with all probes for this run
        temp_matrix <- matrix(NA, nrow = length(all_probes), ncol = ncol(run_matrix))
        rownames(temp_matrix) <- all_probes
        colnames(temp_matrix) <- colnames(run_matrix)
        
        # Fill in the values where probes match
        common_probes <- intersect(all_probes, rownames(run_matrix))
        temp_matrix[common_probes, ] <- run_matrix[common_probes, ]
        
        # Combine with master matrix
        combined_matrix <- cbind(combined_matrix, temp_matrix)
    }
    
    # Remove rows that are all NA (no data across any run)
    combined_matrix <- combined_matrix[rowSums(!is.na(combined_matrix)) > 0, ]
    
    # Save combined result
    write.csv(combined_matrix, file = output_file, quote = FALSE)
    cat("Done. Combined matrix saved to:", output_file, "\n")
    cat("Total samples:", ncol(combined_matrix), "\n")
    cat("Total probes:", nrow(combined_matrix), "\n")
} else {
    cat("No data processed. Please check run patterns and IDAT files.\n")
}