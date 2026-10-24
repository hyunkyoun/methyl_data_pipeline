# Run from the repository root: Rscript --vanilla tests/test_sameparent_dmp.R
test_library <- Sys.getenv("METHYL_DMP_TEST_LIBRARY")
if (nzchar(test_library)) .libPaths(c(test_library, .libPaths()))
source("R/dmp_analysis.R")

cohort <- function(experiment) {
  data.frame(Experiment = experiment, Tissue.Cell.type = "Tcells",
             True.Genotype = c(rep("WT", 9), rep(c("DP10", "DP16", "DP10_DP16"), each = 4)),
             Genotype = c(rep(c("WT_DP10", "WT_DP16", "WT_DP10_WT_DP16"), each = 3),
                          "DP10", "DP10_WT_DP16", "DP10", "DP10_WT_DP16",
                          "DP16", "DP16_WT_DP10", "DP16", "DP16_WT_DP10", rep("DP10_DP16", 4)),
             cross.type = c(rep(c("DP10", "DP16", "Double"), each = 3),
                            "DP10", "Double", "DP10", "Double",
                            "DP16", "Double", "DP16", "Double", rep("Double", 4)))
}
samples <- rbind(cohort(1), cohort(2), cohort(7))
samples$Samples <- paste0("mouse", seq_len(nrow(samples)))
set.seed(20261016)
betas <- matrix(runif(200 * nrow(samples), 0.2, 0.8), nrow = 200,
                dimnames = list(paste0("cg", 1:200), samples$Samples))
# Distinct WT-cross means expose any accidental pooling of controls.
wt <- samples$True.Genotype == "WT"
betas[1, wt] <- rep(rep(c(0.2, 0.5, 0.8), each = 3), 3)
analysis <- run_dmp_analysis(betas, samples[nrow(samples):1, ])
baseline <- run_dmp_analysis(betas, samples, same_parent_experiments = character())
stopifnot(length(analysis$results) == 9L)
wt_genotypes <- c(DP10 = "WT_DP10", DP16 = "WT_DP16", DP10_DP16 = "WT_DP10_WT_DP16")
for (result in analysis$results) {
  expected_cases <- samples$Samples[samples$Experiment == result$experiment & samples$True.Genotype == result$genotype]
  expected_wt <- samples$Samples[samples$Experiment == result$experiment & wt &
                                  (result$experiment == "7" | samples$Genotype == wt_genotypes[[result$genotype]])]
  stopifnot(identical(result$experimental_samples, expected_cases),
            identical(result$control_samples, expected_wt))
  difference <- rowMeans(betas[, expected_cases]) - rowMeans(betas[, expected_wt])
  stopifnot(isTRUE(all.equal(unname(result$table$diff), unname(difference), tolerance = 1e-10)))
  if (result$experiment == "7") {
    stopifnot(identical(result, baseline$results[[result$contrast_name]]))
    next
  }
  # Independently fit the exact two arms and verify the moderated p-values too.
  selected <- c(expected_wt, expected_cases)
  design <- cbind(Intercept = 1, Dp = c(rep(0, length(expected_wt)), rep(1, length(expected_cases))))
  fit <- limma::eBayes(limma::lmFit(betas[, selected], design))
  expected <- limma::topTable(fit, coef = "Dp", number = Inf, sort.by = "none", adjust.method = "BH")
  stopifnot(isTRUE(all.equal(result$table$P.Value, expected$P.Value, tolerance = 1e-10)),
            isTRUE(all.equal(result$table$adj.P.Val, expected$adj.P.Val, tolerance = 1e-10)))
  # Samples outside this comparison cannot influence its effect or p-values.
  changed <- betas
  outside <- !colnames(changed) %in% selected
  changed[, outside] <- matrix(runif(nrow(changed) * sum(outside), 0.01, 0.99), nrow(changed))
  rerun <- run_dmp_analysis(changed, samples)
  stopifnot(identical(result$table, rerun$results[[result$contrast_name]]$table))
}

# cross.type cannot override the explicit WT Genotype rule.
changed_cross <- samples
changed_cross$cross.type <- "unrelated_cross_label"
stopifnot(identical(analysis$results, run_dmp_analysis(betas, changed_cross)$results),
          identical(analysis$results, run_dmp_analysis(betas, samples[, names(samples) != "cross.type"])$results))
# Only True.Genotype defines the Dp arm, even if its Genotype label is missing or misleading.
changed_case_labels <- samples
changed_case_labels$Genotype[!wt] <- NA_character_
stopifnot(identical(analysis$results, run_dmp_analysis(betas, changed_case_labels)$results))
changed_case_labels$Genotype[!wt] <- "WT_DP10"
stopifnot(identical(analysis$results, run_dmp_analysis(betas, changed_case_labels)$results))
expect_error <- function(metadata) {
  error <- tryCatch(run_dmp_analysis(betas, metadata), error = identity)
  stopifnot(inherits(error, "error"), grepl("requires Genotype", conditionMessage(error), fixed = TRUE))
}
expect_error(samples[, names(samples) != "Genotype"])
invalid <- samples
invalid$Genotype[1] <- NA_character_
expect_error(invalid)
invalid$Genotype[1] <- "unknown"
expect_error(invalid)

# Having WT animals from other crosses must not rescue an under-replicated arm.
keep <- !samples$Samples %in% samples$Samples[2:3]
warnings <- character()
small <- withCallingHandlers(run_dmp_analysis(betas[, keep], samples[keep, ]),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
stopifnot(!"DP10_vs_WT_Exp1_Tcells" %in% names(small$results),
          length(small$results) == 8L, length(warnings) == 1L,
          grepl("fewer than two", warnings[1], fixed = TRUE))
cat("Same-parent DMP tests passed: explicit WT Genotype selection, all True.Genotype members retained,",
    "independent pairwise fits, unchanged other experiments, metadata validation, and no control fallback.\n")
