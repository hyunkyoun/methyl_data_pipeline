# CSV HEALTH CHECK FUNCTION
# Usage: csv_health_check("your_file.csv")

csv_health_check <- function(file_path, sample_rows = 5) {
  cat("=== CSV HEALTH CHECK ===\n")
  cat("File:", file_path, "\n\n")
  
  # Check if file exists
  if (!file.exists(file_path)) {
    cat("ERROR: File does not exist!\n")
    return(NULL)
  }
  
  # Read the CSV
  tryCatch({
    data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE)
    
    # Basic dimensions
    cat("BASIC INFO:\n")
    cat("- Dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
    cat("- File size:", round(file.size(file_path) / 1024^2, 2), "MB\n")
    cat("- Data type:", class(data), "\n\n")
    
    # Missing values
    cat("MISSING VALUES:\n")
    total_cells <- nrow(data) * ncol(data)
    na_count <- sum(is.na(data))
    na_percentage <- round((na_count / total_cells) * 100, 2)
    
    cat("- Total NA values:", na_count, "out of", total_cells, "cells\n")
    cat("- Percentage NA:", na_percentage, "%\n")
    
    # NA values per column
    col_na <- colSums(is.na(data))
    if (any(col_na > 0)) {
      cat("- Columns with NA values:\n")
      na_cols <- col_na[col_na > 0]
      for (i in 1:length(na_cols)) {
        col_name <- names(na_cols)[i]
        na_pct <- round((na_cols[i] / nrow(data)) * 100, 2)
        cat("  ", col_name, ":", na_cols[i], "(", na_pct, "%)\n")
      }
    } else {
      cat("- No NA values found in any column\n")
    }
    
    # NA values per row
    row_na <- rowSums(is.na(data))
    completely_missing_rows <- sum(row_na == ncol(data))
    partially_missing_rows <- sum(row_na > 0 & row_na < ncol(data))
    
    cat("- Completely missing rows:", completely_missing_rows, "\n")
    cat("- Partially missing rows:", partially_missing_rows, "\n\n")
    
    # Data type analysis
    cat("DATA TYPES:\n")
    if (ncol(data) > 0) {
      # Check if data is numeric (typical for methylation data)
      numeric_cols <- sapply(data, is.numeric)
      cat("- Numeric columns:", sum(numeric_cols), "out of", ncol(data), "\n")
      
      if (any(numeric_cols)) {
        numeric_data <- data[, numeric_cols, drop = FALSE]
        
        # Value range for numeric data
        cat("- Value range: [", round(min(numeric_data, na.rm = TRUE), 4), 
            " to ", round(max(numeric_data, na.rm = TRUE), 4), "]\n")
        
        # Check for infinite values
        inf_count <- sum(is.infinite(as.matrix(numeric_data)))
        cat("- Infinite values:", inf_count, "\n")
        
        # Values outside 0-1 range (important for methylation beta values)
        outside_01 <- sum(numeric_data < 0 | numeric_data > 1, na.rm = TRUE)
        cat("- Values outside [0,1] range:", outside_01, "\n")
        
        # Summary statistics
        cat("- Mean:", round(mean(as.matrix(numeric_data), na.rm = TRUE), 4), "\n")
        cat("- Median:", round(median(as.matrix(numeric_data), na.rm = TRUE), 4), "\n")
        cat("- Standard deviation:", round(sd(as.matrix(numeric_data), na.rm = TRUE), 4), "\n")
      }
    }
    
    cat("\n")
    
    # Sample preview
    cat("SAMPLE DATA PREVIEW:\n")
    cat("First", sample_rows, "rows and columns:\n")
    preview_data <- data[1:min(sample_rows, nrow(data)), 1:min(sample_rows, ncol(data)), drop = FALSE]
    print(preview_data)
    
    cat("\n")
    
    # Row and column names check
    cat("NAMING CHECK:\n")
    cat("- Row names start with:", paste(head(rownames(data), 3), collapse = ", "), "...\n")
    cat("- Column names start with:", paste(head(colnames(data), 3), collapse = ", "), "...\n")
    
    # Check for duplicate row/column names
    dup_rows <- sum(duplicated(rownames(data)))
    dup_cols <- sum(duplicated(colnames(data)))
    cat("- Duplicate row names:", dup_rows, "\n")
    cat("- Duplicate column names:", dup_cols, "\n\n")
    
    # Quality assessment
    cat("QUALITY ASSESSMENT:\n")
    if (na_percentage < 5) {
      cat("✓ Low missing data (< 5%)\n")
    } else if (na_percentage < 20) {
      cat("⚠ Moderate missing data (5-20%)\n")
    } else {
      cat("✗ High missing data (> 20%)\n")
    }
    
    if (inf_count == 0) {
      cat("✓ No infinite values\n")
    } else {
      cat("✗ Contains infinite values\n")
    }
    
    if (exists("outside_01") && outside_01 == 0) {
      cat("✓ All values in [0,1] range (good for beta values)\n")
    } else if (exists("outside_01") && outside_01 > 0) {
      cat("⚠ Some values outside [0,1] range\n")
    }
    
    if (dup_rows == 0 && dup_cols == 0) {
      cat("✓ No duplicate names\n")
    } else {
      cat("⚠ Contains duplicate names\n")
    }
    
    cat("\n=== HEALTH CHECK COMPLETE ===\n")
    return(data)
    
  }, error = function(e) {
    cat("ERROR reading file:", e$message, "\n")
    return(NULL)
  })
}

# QUICK HEALTH CHECK FUNCTION (less detailed)
quick_health_check <- function(file_path) {
  if (!file.exists(file_path)) {
    cat("File does not exist:", file_path, "\n")
    return()
  }
  
  data <- read.csv(file_path, header = TRUE, row.names = 1)
  na_pct <- round((sum(is.na(data)) / (nrow(data) * ncol(data))) * 100, 2)
  
  cat("File:", basename(file_path), "\n")
  cat("Size:", nrow(data), "x", ncol(data), "\n")
  cat("Missing:", na_pct, "%\n")
  cat("Range: [", round(min(data, na.rm = TRUE), 3), 
      " to ", round(max(data, na.rm = TRUE), 3), "]\n\n")
}

# USAGE EXAMPLES:
# Full health check:
result <- csv_health_check("betas_combat_normalized.csv")

# Quick check of multiple files:
# quick_health_check("betas_RUN1_2.csv")
# quick_health_check("betas_RUN3_4_5.csv")
# quick_health_check("betas_combat_normalized.csv")
# quick_health_check("betas_final_combat_bmiq.csv")