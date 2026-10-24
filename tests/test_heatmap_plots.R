# Run from the repository root: Rscript --vanilla tests/test_heatmap_plots.R
test_library <- Sys.getenv("METHYL_DMP_TEST_LIBRARY")
if (nzchar(test_library)) .libPaths(c(test_library, .libPaths()))
source("R/heatmap_plots.R")

samples <- data.frame(Samples = paste0("mouse", 1:6),
                      True.Genotype = rep(c("DP16", "WT", "RUNX1-DEL"), each = 2))
betas <- rbind(
  gain = c(0.8, 0.7, 0.2, 0.3, 0.7, 0.8),
  loss = c(0.1, 0.2, 0.8, 0.7, 0.3, 0.2),
  mixed = c(0.8, 0.7, 0.4, 0.5, 0.1, 0.2),
  gain2 = c(0.7, 0.8, 0.2, 0.1, 0.6, 0.7)
)
colnames(betas) <- samples$Samples
calls <- list(data.frame(Probe = rownames(betas), diff = c(0.5, -0.6, 0.3, 0.6)),
              data.frame(Probe = c("gain", "mixed"), diff = c(0.5, -0.3)))
directions <- dmp_heatmap_directions(rownames(betas), calls)
stopifnot(identical(unname(directions), c("GOM", "LOM", "Mixed", "GOM")))
# Reordering metadata, calls, or direction entries cannot mislabel CpGs or samples.
prepared <- prepare_dmp_heatmap(rownames(betas), samples$Samples, betas,
                                samples[6:1, ], directions[4:1])
stopifnot(identical(colnames(prepared$matrix), samples$Samples[c(3, 4, 1, 2, 5, 6)]),
          identical(unname(prepared$directions), c("GOM", "GOM", "Mixed", "LOM")),
          identical(as.integer(prepared$counts), c(2L, 1L, 1L)),
          identical(prepared$gaps_row, c(2L, 3L)),
          identical(prepared$gaps_col, c(2L, 4L)),
          sum(prepared$counts) == nrow(betas), prepared$height == 4)
original <- betas[rownames(prepared$matrix), colnames(prepared$matrix), drop = FALSE]
stopifnot(isTRUE(all.equal(prepared$matrix, sweep(original, 1, rowMeans(original), "-"))),
          max(abs(rowMeans(prepared$matrix))) < 1e-12)
stopifnot(identical(directions, dmp_heatmap_directions(rownames(betas), rev(calls))))
# A focused selection follows only its own comparison, even if another reverses a CpG.
focused <- dmp_heatmap_directions(c("mixed", "loss"), calls[1])
stopifnot(identical(unname(focused), c("GOM", "LOM")))

scratch <- tempfile("methyl-heatmap-test-")
dir.create(scratch)
fixtures <- list(mixed = rownames(betas), gain_only = c("gain", "gain2"),
                 loss_only = "loss", single_gain = "gain")
for (fixture in names(fixtures)) {
  probes <- fixtures[[fixture]]
  path <- file.path(scratch, paste0(fixture, ".pdf"))
  result <- plot_dmp_heatmap(probes, samples$Samples, paste("Bracket layout:", fixture),
                             path, directions, betas, samples)
  stopifnot(sum(result$counts) == length(probes),
            identical(readChar(path, 4, useBytes = TRUE), "%PDF"), file.info(path)$size > 1000)
}
# Identical displayed values still get a valid symmetric color scale.
constant <- betas["gain", , drop = FALSE]
constant[,] <- 0.5
zero <- plot_dmp_heatmap("gain", samples$Samples, "Zero centered values", file.path(scratch, "zero.pdf"),
                         directions, constant, samples)
stopifnot(zero$limit > 0, all(zero$matrix == 0))
# Missing or zero directions must not silently turn into a loss label.
stopifnot(inherits(try(dmp_heatmap_directions("absent", calls), silent = TRUE), "try-error"),
          inherits(try(dmp_heatmap_directions("gain", list(data.frame(Probe = "gain", diff = 0))),
                       silent = TRUE), "try-error"))
cat("Heatmap tests passed: CpG direction, mixed calls, focused selection, WT-first columns,",
    "alignment, centering, block counts, one-direction/single-row PDFs, and zero-value scale.\n")
