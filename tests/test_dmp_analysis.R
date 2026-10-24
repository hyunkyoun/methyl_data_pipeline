# Run from the repository root: Rscript --vanilla tests/test_dmp_analysis.R
test_library <- Sys.getenv("METHYL_DMP_TEST_LIBRARY")
if (nzchar(test_library)) .libPaths(c(test_library, .libPaths()))
stopifnot(requireNamespace("limma", quietly = TRUE))
source("R/dmp_analysis.R")

cohort <- function(experiment, tissue, genotypes) {
  data.frame(Experiment = experiment, Tissue.Cell.type = tissue, True.Genotype = genotypes)
}
targets <- rbind(
  cohort(7, "Tcells", rep(c("WT", "DP10", "DP16"), each = 3)),
  cohort(8, "Tcells", rep(c("WT", "DP10"), each = 3)),
  cohort(7, "thymus", rep(c("WT", "DP10"), each = 3)),
  cohort(3, "thymus", "RUNX1-DEL_DP16"),
  cohort(4, "Tcells", c(rep("WT", 3), "DP10"))
)
targets$Samples <- paste0("mouse", seq_len(nrow(targets)))
# Legacy groups intentionally overlap and disagree with the actual genotype arms.
targets$group <- rep(c("legacy_A", "legacy_B"), length.out = nrow(targets))
targets$Genotype <- "misleading_legacy_label"
set.seed(20261014)
betas <- matrix(runif(100 * nrow(targets), 0.45, 0.55), nrow = 100,
                dimnames = list(paste0("cg", 1:100), targets$Samples))
effect <- ifelse(targets$True.Genotype == "WT", 0,
                 ifelse(targets$Experiment == 8, -0.35, 0.35))
betas[1, ] <- 0.5 + effect + rnorm(nrow(targets), sd = 0.002)
betas[2, ] <- 0.5 - effect + rnorm(nrow(targets), sd = 0.002)

warnings <- character()
analysis <- withCallingHandlers(run_dmp_analysis(betas, targets[nrow(targets):1, ]),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
expected_ids <- c("DP10_vs_WT_Exp7_Tcells", "DP16_vs_WT_Exp7_Tcells",
                  "DP10_vs_WT_Exp8_Tcells", "DP10_vs_WT_Exp7_thymus")
stopifnot(setequal(names(analysis$results), expected_ids), length(warnings) == 2L,
          any(grepl("Experiment 3 / thymus", warnings, fixed = TRUE)),
          nrow(analysis$contrast_weights) == 8L,
          all(analysis$contrast_weights$Sample.Count == 3L),
          all(analysis$contrast_weights$Coefficient %in% c(-1, 1)))
for (result in analysis$results) {
  cohort_samples <- targets[targets$Experiment == result$experiment &
                              targets$Tissue.Cell.type == result$tissue, ]
  controls <- cohort_samples$Samples[cohort_samples$True.Genotype == "WT"]
  treatment <- cohort_samples$Samples[cohort_samples$True.Genotype == result$genotype]
  stopifnot(identical(result$control_samples, controls),
            identical(result$experimental_samples, treatment),
            all(result$table$Experiment == result$experiment),
            all(result$table$Tissue == result$tissue),
            identical(result$dmps$Probe, c("cg1", "cg2")))
  expected_diff <- rowMeans(betas[, treatment]) - rowMeans(betas[, controls])
  stopifnot(isTRUE(all.equal(unname(result$table$diff), unname(expected_diff), tolerance = 1e-10)))
}
stopifnot(analysis$results$DP10_vs_WT_Exp7_Tcells$table$diff[1] > 0,
          analysis$results$DP10_vs_WT_Exp8_Tcells$table$diff[1] < 0)

# Other cohorts must not affect either effect sizes or moderated p-values.
changed <- betas
outside <- targets$Experiment != 7 | targets$Tissue.Cell.type != "Tcells"
changed[, outside] <- matrix(runif(nrow(betas) * sum(outside), 0.01, 0.99), nrow(betas))
changed_analysis <- suppressWarnings(run_dmp_analysis(changed, targets))
for (id in expected_ids[1:2]) {
  stopifnot(identical(analysis$results[[id]]$table, changed_analysis$results[[id]]$table))
}

# Metadata omissions and unmatched samples must fail before any comparisons are fitted.
expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
}
expect_error(run_dmp_analysis(betas, targets[, names(targets) != "Experiment"]), "Experiment")
incomplete <- targets
incomplete$Experiment[1] <- NA
expect_error(run_dmp_analysis(betas, incomplete), "non-missing")
expect_error(run_dmp_analysis(betas, targets[-1, ]), "invalid sample metadata")
duplicated <- targets
duplicated$Samples[2] <- duplicated$Samples[1]
expect_error(run_dmp_analysis(betas, duplicated), "invalid sample metadata")

# Zero significant hits still retain full results with the correct cohort labels.
empty <- suppressWarnings(run_dmp_analysis(betas, targets, delta_beta_cutoff = 1))
stopifnot(all(vapply(empty$results, function(x) nrow(x$dmps) == 0L, logical(1))))

# Eligible local metadata should produce 19 comparisons, including 9 T-cell contrasts.
sheet_path <- "idats/samplesheet_169_samples.csv"
if (file.exists(sheet_path)) {
  sheet <- read.csv(sheet_path, fileEncoding = "UTF-8-BOM")
  sheet <- sheet[!tolower(trimws(sheet$Tissue.Cell.type)) %in% "brain", , drop = FALSE]
  smoke_betas <- matrix(runif(100 * nrow(sheet), 0.2, 0.8), nrow = 100,
                        dimnames = list(paste0("cg", 1:100), sheet$Samples))
  smoke <- suppressWarnings(run_dmp_analysis(smoke_betas, sheet))
  stopifnot(length(smoke$results) == 19L,
            sum(vapply(smoke$results, function(x) x$tissue == "Tcells", logical(1))) == 9L)
  dp10 <- smoke$results$DP10_vs_WT_Exp8_Tcells
  double <- smoke$results$DP10_DP16_vs_WT_Exp8_Tcells
  stopifnot(length(dp10$experimental_samples) == 7L,
            length(double$experimental_samples) == 7L,
            "12926_8" %in% double$experimental_samples,
            !"12926_8" %in% dp10$experimental_samples,
            !any(vapply(smoke$results, function(x) x$tissue == "brain", logical(1))))
  for (result in smoke$results) {
    members <- sheet[match(c(result$experimental_samples, result$control_samples), sheet$Samples), ]
    stopifnot(all(members$Experiment == result$experiment),
              all(members$Tissue.Cell.type == result$tissue))
  }
  expected_counts <- list("1" = c(DP10 = 9L, DP16 = 10L, DP10_DP16 = 4L),
                          "2" = c(DP10 = 6L, DP16 = 7L, DP10_DP16 = 3L))
  expected_wt <- list("1" = c(DP10 = 6L, DP16 = 4L, DP10_DP16 = 5L),
                      "2" = c(DP10 = 4L, DP16 = 4L, DP10_DP16 = 3L))
  expected_wt_genotype <- c(DP10 = "WT_DP10", DP16 = "WT_DP16", DP10_DP16 = "WT_DP10_WT_DP16")
  for (result in Filter(function(x) x$experiment %in% c("1", "2"), smoke$results)) {
    controls <- sheet[match(result$control_samples, sheet$Samples), ]
    cases <- sheet[match(result$experimental_samples, sheet$Samples), ]
    stopifnot(length(result$experimental_samples) == expected_counts[[result$experiment]][[result$genotype]],
              length(result$control_samples) == expected_wt[[result$experiment]][[result$genotype]],
              all(controls$Genotype == expected_wt_genotype[[result$genotype]]),
              all(controls$True.Genotype == "WT"), all(cases$True.Genotype == result$genotype))
    if (result$genotype %in% c("DP10", "DP16")) stopifnot(any(cases$cross.type == "Double"))
  }
  stopifnot("11157_1" %in% smoke$results$DP10_vs_WT_Exp1_thymus$experimental_samples)
}
cat("DMP tests passed: same-experiment/tissue controls, independent fits, short names,",
    "sample alignment, minimum replication, corrected 12926_8, and 19 non-brain comparisons.\n")
