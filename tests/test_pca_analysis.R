# Run from the repository root: Rscript --vanilla tests/test_pca_analysis.R
test_library <- Sys.getenv("METHYL_DMP_TEST_LIBRARY")
if (nzchar(test_library)) .libPaths(c(test_library, .libPaths()))
source("R/pca_analysis.R")

samples <- data.frame(
  Samples = paste0("sample", 1:8),
  Experiment = 1,
  Tissue.Cell.type = rep(c("tissue_A", "tissue_B"), each = 4),
  BATCH = rep(c("BATCH1", "BATCH2"), each = 4),
  True.Genotype = c("WT", "DP16", "WT", "DP16", "DP16", "WT", "DP16", "WT")
)
betas <- rbind(
  between_tissues = c(rep(0.1, 4), rep(0.9, 4)),
  A_signal = c(0.1, 0.2, 0.8, 0.9, rep(0.5, 4)),
  A_second = c(0.4, 0.1, 0.6, 0.5, rep(0.5, 4)),
  B_signal = c(rep(0.5, 4), 0.1, 0.2, 0.8, 0.9),
  B_second = c(rep(0.5, 4), 0.4, 0.1, 0.6, 0.5),
  constant = rep(0.5, 8)
)
colnames(betas) <- samples$Samples
scratch <- tempfile("methyl-tissue-pca-test-")
palettes <- pca_color_palettes(samples)
overall <- create_simple_pca(betas, samples, 2L, "overall",
                             file.path(scratch, "figures"), file.path(scratch, "tables"))
stopifnot("between_tissues" %in% overall$top_cpgs)

# Shuffled input columns must still produce correctly aligned sample scores and exports.
shuffled <- betas[, c(8, 2, 6, 4, 3, 7, 1, 5)]
tissues <- create_tissue_pcas(shuffled, samples, 10000L,
                             figure_dir = file.path(scratch, "figures"),
                             table_dir = file.path(scratch, "tables"))
stopifnot(identical(tissues$tissue_A$top_cpgs, c("A_signal", "A_second")),
          identical(tissues$tissue_B$top_cpgs, c("B_signal", "B_second")))
for (tissue in names(tissues)) {
  result <- tissues[[tissue]]
  idx <- which(samples$Tissue.Cell.type == tissue)
  stopifnot(identical(result$pca_df$Sample, samples$Samples[idx]),
            identical(result$pca_df$True.Genotype, samples$True.Genotype[idx]),
            identical(result$color_palettes, palettes),
            identical(colnames(result$betas_top), samples$Samples[idx]),
            all(is.finite(result$pca$x)))
  expected <- prcomp(t(betas[result$top_cpgs, idx]), center = TRUE, scale. = TRUE)
  stopifnot(isTRUE(all.equal(result$pca$sdev, expected$sdev)),
            isTRUE(all.equal(unname(as.matrix(dist(result$pca$x))),
                             unname(as.matrix(dist(expected$x))))))
  prefix <- paste0("COMBAT_ADJUSTED_", tissue)
  saved <- read.csv(file.path(scratch, "tables", "betas", paste0(prefix, ".csv")), row.names = 1)
  stopifnot(isTRUE(all.equal(as.matrix(saved), result$betas_top)),
            all(file.exists(file.path(scratch, "figures", c("genotype", "batch"), paste0(prefix, ".pdf")))))
}

# A tissue PCA must not silently combine animals from different experiments.
mixed_experiments <- samples
mixed_experiments$Experiment[5:8] <- 2
error <- tryCatch(create_tissue_pcas(betas, mixed_experiments,
                                    figure_dir = scratch, table_dir = scratch), error = identity)
stopifnot(inherits(error, "error"), grepl("exactly one Experiment", conditionMessage(error)))

# A small or invariant tissue should not prevent other tissues from being plotted.
warnings <- character()
expect_skip <- function(matrix, metadata) {
  withCallingHandlers(create_simple_pca(matrix, metadata, 2L, "skipped",
                                        file.path(scratch, "skipped_figures"),
                                        file.path(scratch, "skipped_tables")),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
}
stopifnot(is.null(expect_skip(betas[, 1:2, drop = FALSE], samples[1:2, ])),
          is.null(expect_skip(betas["constant", , drop = FALSE], samples)),
          length(warnings) == 2L,
          !dir.exists(file.path(scratch, "skipped_figures")))
bad_samples <- samples
bad_samples$Samples[1] <- "not_in_matrix"
error <- tryCatch(create_simple_pca(betas, bad_samples, 2L, "invalid",
                                    file.path(scratch, "invalid"), file.path(scratch, "invalid")),
                   error = identity)
stopifnot(inherits(error, "error"), grepl("absent", conditionMessage(error)))
cat("PCA tests passed: within-tissue feature selection, sample alignment, matching colors,",
    "independent PCA fits, limited variable CpGs, and skipped uninformative subsets.\n")
