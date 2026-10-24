# Run from the repository root: Rscript --vanilla tests/test_combat.R
test_library <- Sys.getenv("METHYL_DMP_TEST_LIBRARY")
if (nzchar(test_library)) .libPaths(c(test_library, .libPaths()))
stopifnot(requireNamespace("sva", quietly = TRUE))

# Exercise the actual workflow, replacing only CSV output with a capture.
workflow_lines <- readLines("workflow_new.R", warn = FALSE)
first <- match("# Technology between old vs new", workflow_lines)
last <- match("# DMP analysis", workflow_lines)
stopifnot(length(first) == 1, !is.na(first), length(last) == 1, first < last)
combat_block <- parse(text = workflow_lines[seq.int(first, last - 1)])

run_combat_block <- function(betas, samples) {
  env <- new.env(parent = globalenv())
  env$betas <- betas
  env$targets <- samples
  env$output_dirs <- list(betas = file.path("test-output", "tables", "betas"),
                          metadata = file.path("test-output", "metadata"))
  env$combat_calls <- 0L
  env$ComBat <- function(...) {
    env$combat_calls <- env$combat_calls + 1L
    sva::ComBat(...)
  }
  env$saved_outputs <- list()
  env$write.csv <- function(x, file, row.names) env$saved_outputs[[file]] <- x
  invisible(capture.output(suppressMessages(eval(combat_block, env))))
  env
}

expect_error <- function(expr, text) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}

# Known biological and technology effects, with deliberately unequal group counts.
set.seed(20261001)
samples <- data.frame(
  Samples = paste0("mouse", seq_len(40)),
  Experiment = 1,
  True.Genotype = rep(c("WT", "DP16"), each = 20),
  group = "legacy_group_is_not_biology",
  EPIC_RUN_CLEAN = c(rep("1_2", 16), rep("3_4", 4), rep("1_2", 4),
                     rep("5", 8), rep("6", 8)),
  Tissue.Cell.type = "Tcells"
)
is_new <- samples$EPIC_RUN_CLEAN != "1_2"
is_dp16 <- samples$True.Genotype == "DP16"
n_probes <- 1000L
true_group_effect <- runif(n_probes, 0.7, 1.1)
true_technology_effect <- runif(n_probes, 0.6, 1.0)
noise <- matrix(rnorm(n_probes * nrow(samples), sd = 0.12), nrow = n_probes)
noise[, is_new] <- noise[, is_new] * 1.6
# Center within each group/technology cell so the known mean effects are exact.
for (cell in unique(paste(is_new, is_dp16))) {
  idx <- which(paste(is_new, is_dp16) == cell)
  noise[, idx] <- noise[, idx, drop = FALSE] - rowMeans(noise[, idx, drop = FALSE])
}
mvals <- rnorm(n_probes, -0.8, 0.2) +
  true_group_effect %o% as.numeric(is_dp16) +
  true_technology_effect %o% as.numeric(is_new) + noise
betas <- plogis(mvals * log(2))
dimnames(betas) <- list(paste0("cg", seq_len(n_probes)), samples$Samples)

result <- run_combat_block(betas, samples)
stopifnot(result$combat_calls == 1L,
          identical(dimnames(result$betas_adjusted), dimnames(betas)),
          all(is.finite(result$betas_adjusted)),
          all(result$betas_adjusted >= 0 & result$betas_adjusted <= 1),
          identical(as.character(result$batch), ifelse(is_new, "BATCH2", "BATCH1")),
          identical(result$saved_outputs[[file.path(result$output_dirs$betas,
                                                    "betas_combat_adjusted.csv")]],
                    result$betas_adjusted),
          identical(result$saved_outputs[[file.path(result$output_dirs$metadata,
                                                    "samples_after_qc.csv")]],
                    result$matching_samples))

adjusted_mvals <- log2(result$betas_adjusted / (1 - result$betas_adjusted))
for (technology in c(FALSE, TRUE)) {
  adjusted_group_effect <- rowMeans(adjusted_mvals[, is_dp16 & is_new == technology]) -
    rowMeans(adjusted_mvals[, !is_dp16 & is_new == technology])
  stopifnot(max(abs(adjusted_group_effect - true_group_effect)) < 1e-8)
}
adjusted_technology_effect <- rowMeans(adjusted_mvals[, !is_dp16 & is_new]) -
  rowMeans(adjusted_mvals[, !is_dp16 & !is_new])
stopifnot(median(abs(adjusted_technology_effect)) < 0.05)
cat(sprintf("Synthetic effects: median technology shift %.3f -> %.3f M units; group effects retained.\n",
            median(true_technology_effect), median(abs(adjusted_technology_effect))))

# Complete confounding must stop rather than silently removing group protection.
confounded <- samples
confounded$EPIC_RUN_CLEAN <- ifelse(is_dp16, "5", "1_2")
expect_error(run_combat_block(betas, confounded), "fully confounded")

# A full-rank model still needs residual degrees of freedom.
saturated_idx <- c(1, 17, 25)
expect_error(run_combat_block(betas[, saturated_idx], samples[saturated_idx, ]),
             "residual variation")

# Fail clearly on incomplete biological metadata.
for (column in c("Experiment", "Tissue.Cell.type", "True.Genotype")) {
  for (bad_value in c(NA_character_, "", "  ")) {
    incomplete <- samples
    incomplete[[column]][1] <- bad_value
    expect_error(run_combat_block(betas, incomplete), "non-missing Experiment, tissue, and True.Genotype")
  }
}
expect_error(run_combat_block(betas, samples[, names(samples) != "True.Genotype"]),
             "non-missing Experiment, tissue, and True.Genotype")

# One technology passes through unchanged, even if it contains several run labels.
single_technology <- run_combat_block(betas[, is_new], samples[is_new, ])
stopifnot(single_technology$combat_calls == 0L,
          identical(single_technology$betas_adjusted, betas[, is_new]))

# One biological group across two technologies uses an intercept-only model.
one_group <- samples
one_group$True.Genotype <- "WT"
single_group <- run_combat_block(betas, one_group)
stopifnot(single_group$combat_calls == 1L, ncol(single_group$mod) == 1L,
          identical(colnames(single_group$mod), "(Intercept)"),
          all(is.finite(single_group$betas_adjusted)))

# The same genotype in another tissue gets its own protected mean.
two_tissues <- samples
two_tissues$Tissue.Cell.type[seq(1, nrow(samples), 2)] <- "thymus"
tissue_result <- run_combat_block(betas, two_tissues)
stopifnot(nlevels(tissue_result$group) == 4L, ncol(tissue_result$mod) == 4L)

# The same tissue/genotype in another experiment also gets a protected mean.
two_experiments <- samples
two_experiments$Experiment <- rep(1:2, length.out = nrow(samples))
experiment_result <- run_combat_block(betas, two_experiments)
stopifnot(nlevels(experiment_result$group) == 4L, ncol(experiment_result$mod) == 4L,
          qr(experiment_result$combat_design)$rank == ncol(experiment_result$combat_design))

# Boundary/missing values keep the existing whole-CpG exclusion policy.
invalid_betas <- betas
invalid_betas[1, 1] <- 0
invalid_betas[2, 1] <- 1
invalid_betas[3, 1] <- NA_real_
filtered <- run_combat_block(invalid_betas, samples)
stopifnot(identical(rownames(filtered$betas_adjusted), rownames(betas)[-(1:3)]),
          identical(colnames(filtered$betas_adjusted), colnames(betas)))
expect_error(run_combat_block(invalid_betas[1:3, ], samples), "at least two CpGs")
expect_error(run_combat_block(invalid_betas[1:4, ], samples), "at least two CpGs")

# Optional smoke check with the user's local metadata and synthetic values only.
sheet_path <- "idats/samplesheet_169_samples.csv"
if (file.exists(sheet_path)) {
  sheet <- read.csv(sheet_path, fileEncoding = "UTF-8-BOM")
  sheet <- sheet[!tolower(trimws(sheet$Tissue.Cell.type)) %in% "brain", , drop = FALSE]
  smoke_betas <- matrix(runif(200 * nrow(sheet), 0.2, 0.8), nrow = 200,
                       dimnames = list(paste0("cg", seq_len(200)), sheet$Samples))
  smoke <- run_combat_block(smoke_betas, sheet)
  stopifnot(smoke$combat_calls == 1L,
            qr(smoke$combat_design)$rank == ncol(smoke$combat_design),
            all(is.finite(smoke$betas_adjusted)))
  cat(sprintf("Samplesheet smoke check: %d samples, %d groups, %d technologies, design rank %d/%d (before QC).\n",
              nrow(sheet), nlevels(smoke$group), nlevels(smoke$batch),
              qr(smoke$combat_design)$rank, ncol(smoke$combat_design)))
}

cat("ComBat tests passed with sva", as.character(utils::packageVersion("sva")), "\n")
