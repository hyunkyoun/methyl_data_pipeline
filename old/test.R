# ============================================================================
# DEEP DIAGNOSIS OF VOLCANO PLOT ISSUES
# ============================================================================

library(limma)
library(ggplot2)

# Load both the corrected and uncorrected data for comparison
cat("=== LOADING DATA FOR DIAGNOSIS ===\n")

# Load ComBat-corrected data
betas_corrected <- read.csv("betas_final_combat_simple.csv", row.names=1)
metadata <- read.csv("metadata_final_simple.csv")

# Try to load uncorrected data for comparison
# We'll reconstruct it from the BMIQ results if needed
if(file.exists("betas_bmiq_RUN1_2.csv") && file.exists("betas_bmiq_RUN3_4_5.csv")) {
  bmiq_12 <- read.csv("betas_bmiq_RUN1_2.csv", row.names=1)
  bmiq_345 <- read.csv("betas_bmiq_RUN3_4_5.csv", row.names=1)
  
  # Remove the first column if it's CpG_ID
  if(colnames(bmiq_12)[1] == "CpG_ID") bmiq_12 <- bmiq_12[, -1]
  if(colnames(bmiq_345)[1] == "CpG_ID") bmiq_345 <- bmiq_345[, -1]
  
  common_probes <- intersect(rownames(bmiq_12), rownames(bmiq_345))
  betas_uncorrected <- cbind(bmiq_12[common_probes, ], bmiq_345[common_probes, ])
  
  cat("Loaded uncorrected BMIQ data for comparison\n")
} else {
  betas_uncorrected <- NULL
  cat("Could not load uncorrected data\n")
}

cat("Corrected data:", nrow(betas_corrected), "probes x", ncol(betas_corrected), "samples\n")

# === STEP 1: VERIFY WHAT COMPARISON IS ACTUALLY BEING MADE ===
cat("\n=== STEP 1: VERIFY THE COMPARISON ===\n")

# Show exactly what's in your volcano plot data
cat("What comparison produced your volcano plot?\n")
cat("Please tell me:\n")
cat("1. Which two groups are you comparing?\n")
cat("2. Are you using betas_final_combat_simple.csv?\n")
cat("3. What differential analysis method are you using?\n")

# Let's look at the available groups
cat("\nAvailable genotype combinations:\n")
metadata$genotype_combo <- paste(metadata$genotype_dp10, "|", metadata$genotype_dp16, "|", metadata$Runx1)
combo_table <- table(metadata$genotype_combo)
print(combo_table[combo_table > 0])

# Check batch distribution of each combination
cat("\nBatch distribution by genotype combination:\n")
metadata$Batch <- ifelse(metadata$RUN <= 2, "Batch1", "Batch2")

for(combo in names(combo_table)[combo_table > 2]) {  # Only show combos with >2 samples
  samples_in_combo <- metadata$genotype_combo == combo
  batch_dist <- table(metadata$Batch[samples_in_combo])
  
  if(length(batch_dist) > 1) {
    cat(sprintf("'%s': Batch1=%d, Batch2=%d\n", combo, batch_dist["Batch1"], batch_dist["Batch2"]))
  } else {
    batch_name <- names(batch_dist)[1]
    cat(sprintf("'%s': %s=%d (ONLY ONE BATCH!)\n", combo, batch_name, batch_dist[1]))
  }
}

# === STEP 2: TEST A DEFINITELY SAFE COMPARISON ===
cat("\n=== STEP 2: TESTING A GUARANTEED SAFE COMPARISON ===\n")

# Let's compare samples WITHIN the same batch to eliminate batch effects entirely
batch1_samples <- metadata$Batch == "Batch1"
batch1_metadata <- metadata[batch1_samples, ]

cat("Batch 1 only - available combinations:\n")
batch1_combos <- table(batch1_metadata$genotype_combo)
print(batch1_combos[batch1_combos > 2])

# Find the two most common genotype combinations in Batch 1
if(length(batch1_combos[batch1_combos > 2]) >= 2) {
  
  top_combos <- names(sort(batch1_combos[batch1_combos > 2], decreasing = TRUE))[1:2]
  
  cat(sprintf("\nTesting comparison within Batch 1: '%s' vs '%s'\n", top_combos[1], top_combos[2]))
  
  # Get samples
  group1_samples <- batch1_metadata$genotype_combo == top_combos[1]
  group2_samples <- batch1_metadata$genotype_combo == top_combos[2]
  
  cat(sprintf("Group 1: %d samples\n", sum(group1_samples)))
  cat(sprintf("Group 2: %d samples\n", sum(group2_samples)))
  
  # Get corresponding beta values (using corrected data)
  batch1_betas <- betas_corrected[, batch1_samples]
  analysis_betas <- batch1_betas[, group1_samples | group2_samples]
  
  # Create group labels
  group_labels <- c(rep("Group1", sum(group1_samples)), rep("Group2", sum(group2_samples)))
  
  # Quick differential analysis
  analysis_betas[analysis_betas <= 0.001] <- 0.001
  analysis_betas[analysis_betas >= 0.999] <- 0.999
  
  M_values <- log2(analysis_betas / (1 - analysis_betas))
  M_values <- M_values[complete.cases(M_values), ]
  
  if(nrow(M_values) > 1000) {
    design <- model.matrix(~ 0 + factor(group_labels))
    colnames(design) <- c("Group1", "Group2")
    
    fit <- lmFit(M_values, design)
    contrast_matrix <- makeContrasts(Group1 - Group2, levels = design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    results <- topTable(fit2, number = Inf, sort.by = "none")
    results$adj_pval <- p.adjust(results$P.Value, method = "BH")
    
    # Count significant results
    sig_up <- sum(results$adj_pval < 0.05 & results$logFC > 0.5, na.rm = TRUE)
    sig_down <- sum(results$adj_pval < 0.05 & results$logFC < -0.5, na.rm = TRUE)
    
    cat(sprintf("\nWITHIN-BATCH COMPARISON RESULTS:\n"))
    cat(sprintf("Significant up: %d\n", sig_up))
    cat(sprintf("Significant down: %d\n", sig_down))
    
    if(sig_up > 0 && sig_down > 0) {
      ratio <- max(sig_up, sig_down) / min(sig_up, sig_down)
      cat(sprintf("Ratio: %.2f:1\n", ratio))
      
      if(ratio < 3) {
        cat("✅ WITHIN-BATCH comparison looks balanced!\n")
        cat("This confirms the issue is batch-related, not biological.\n")
      } else {
        cat("⚠️ Even within-batch comparison shows bias.\n")
        cat("This suggests a deeper issue with the comparison groups.\n")
      }
    } else {
      cat("No significant differences found within batch.\n")
    }
  }
} else {
  cat("Not enough different genotype combinations in Batch 1 for testing.\n")
}

# === STEP 3: CHECK IF COMBAT ACTUALLY WORKED ===
cat("\n=== STEP 3: VERIFY COMBAT EFFECTIVENESS ===\n")

if(!is.null(betas_uncorrected)) {
  # Compare batch differences before and after ComBat
  batch1_indices <- metadata$Batch == "Batch1"
  batch2_indices <- metadata$Batch == "Batch2"
  
  # Before ComBat
  batch1_means_before <- rowMeans(betas_uncorrected[, batch1_indices], na.rm = TRUE)
  batch2_means_before <- rowMeans(betas_uncorrected[, batch2_indices], na.rm = TRUE)
  batch_diff_before <- batch1_means_before - batch2_means_before
  
  # After ComBat
  batch1_means_after <- rowMeans(betas_corrected[, batch1_indices], na.rm = TRUE)
  batch2_means_after <- rowMeans(betas_corrected[, batch2_indices], na.rm = TRUE)
  batch_diff_after <- batch1_means_after - batch2_means_after
  
  cat("Batch differences (Batch1 - Batch2 mean methylation):\n")
  cat("Before ComBat - Range:", round(range(batch_diff_before, na.rm = TRUE), 4), "\n")
  cat("Before ComBat - SD:", round(sd(batch_diff_before, na.rm = TRUE), 4), "\n")
  cat("After ComBat - Range:", round(range(batch_diff_after, na.rm = TRUE), 4), "\n")
  cat("After ComBat - SD:", round(sd(batch_diff_after, na.rm = TRUE), 4), "\n")
  
  # ComBat should reduce the SD dramatically
  reduction <- sd(batch_diff_before, na.rm = TRUE) / sd(batch_diff_after, na.rm = TRUE)
  cat("Batch effect reduction factor:", round(reduction, 2), "x\n")
  
  if(reduction > 5) {
    cat("✅ ComBat worked well (>5x reduction)\n")
  } else if(reduction > 2) {
    cat("⚠️ ComBat had moderate effect (2-5x reduction)\n")
  } else {
    cat("❌ ComBat had minimal effect (<2x reduction)\n")
  }
}

# === STEP 4: RECOMMENDATIONS ===
cat("\n=== RECOMMENDATIONS ===\n")
cat("Based on the persistent volcano plot pattern:\n\n")

cat("1. 🔍 VERIFY YOUR COMPARISON:\n")
cat("   - Make sure you're comparing biologically relevant groups\n")
cat("   - Check that both groups have samples in both batches\n\n")

cat("2. 🧪 TRY WITHIN-BATCH ANALYSIS:\n")
cat("   - Analyze Batch1 and Batch2 separately\n")
cat("   - See if the same biological pattern appears in both\n\n")

cat("3. 📊 INCLUDE BATCH AS COVARIATE:\n")
cat("   - Instead of removing batch effects, model them:\n")
cat("   - design <- model.matrix(~ genotype + batch + cell_type)\n\n")

cat("4. 🔄 TRY DIFFERENT NORMALIZATION:\n")
cat("   - Maybe BMIQ + batch covariate instead of ComBat\n")
cat("   - Or more aggressive ComBat on M-values\n\n")

cat("Please run this diagnosis and tell me:\n")
cat("- What specific comparison is giving you the volcano plot?\n")
cat("- What do the within-batch results show?\n")
cat("- How much did ComBat reduce batch effects?\n")