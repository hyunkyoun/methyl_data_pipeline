# Direction comes from significant experimental-minus-WT effects, not tile colors.
# A union of comparisons can contain CpGs with opposite significant directions.
dmp_heatmap_directions <- function(probes, dmp_tables) {
  calls <- do.call(rbind, lapply(dmp_tables, function(x) x[, c("Probe", "diff"), drop = FALSE]))
  if (is.null(calls) || anyNA(calls$diff) || any(!is.finite(calls$diff))) {
    stop("Heatmap directions require finite DMP effects.")
  }
  gains <- probes %in% calls$Probe[calls$diff > 0]
  losses <- probes %in% calls$Probe[calls$diff < 0]
  if (any(!gains & !losses)) stop("Every heatmap CpG must have a nonzero significant DMP effect.")
  setNames(ifelse(gains & losses, "Mixed", ifelse(gains, "GOM", "LOM")), probes)
}

prepare_dmp_heatmap <- function(probes, selected_samples, betas_data, targets_data, directions) {
  if (!length(probes) || anyDuplicated(probes) || anyDuplicated(selected_samples) ||
      !all(probes %in% rownames(betas_data)) || !all(selected_samples %in% colnames(betas_data))) {
    stop("Heatmap requires unique CpGs and samples present in the beta matrix.")
  }
  sample_indices <- match(selected_samples, targets_data$Samples)
  if (!length(selected_samples) || anyNA(sample_indices)) stop("Heatmap samples need matching metadata.")
  genotypes <- as.character(targets_data$True.Genotype[sample_indices])
  if (anyNA(genotypes) || any(!nzchar(genotypes))) stop("Heatmap samples need genotype labels.")
  genotype_order <- unique(genotypes)
  genotype_order <- c(intersect("WT", genotype_order), setdiff(genotype_order, "WT"))
  column_order <- order(match(genotypes, genotype_order))
  selected_samples <- selected_samples[column_order]
  genotypes <- genotypes[column_order]

  direction_order <- c("GOM", "Mixed", "LOM")
  if (is.null(names(directions)) || anyDuplicated(names(directions))) {
    stop("Heatmap directions must be named by unique CpG IDs.")
  }
  directions <- directions[probes]
  if (anyNA(directions) || !all(directions %in% direction_order)) {
    stop("Every heatmap CpG needs a GOM, LOM, or Mixed label.")
  }
  beta_matrix <- betas_data[probes, selected_samples, drop = FALSE]
  if (any(!is.finite(beta_matrix))) stop("Heatmap beta values must be finite.")
  centered <- sweep(beta_matrix, 1L, rowMeans(beta_matrix), FUN = "-")
  # Cluster within each block so GOM stays at the top and LOM at the bottom.
  row_order <- unlist(lapply(direction_order, function(direction) {
    indices <- which(directions == direction)
    if (length(indices) > 1L) {
      indices <- indices[stats::hclust(stats::dist(centered[indices, , drop = FALSE]))$order]
    }
    indices
  }), use.names = FALSE)
  centered <- centered[row_order, , drop = FALSE]
  directions <- directions[row_order]
  counts <- table(factor(directions, levels = direction_order))
  genotype_counts <- table(factor(genotypes, levels = genotype_order))
  list(
    matrix = centered, directions = directions, counts = counts,
    genotypes = genotypes, genotype_counts = genotype_counts,
    gaps_row = as.integer(head(cumsum(counts[counts > 0L]), -1L)),
    gaps_col = as.integer(head(cumsum(genotype_counts), -1L)),
    limit = max(0.01, abs(centered)),
    height = min(8.5, max(4, nrow(centered) * 0.02 + 3))
  )
}

draw_dmp_heatmap <- function(prepared, plot_title, direction_note) {
  colors <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  direction_colors <- c(GOM = "#B2182B", Mixed = "#666666", LOM = "#2166AC")
  heatmap <- pheatmap::pheatmap(
    prepared$matrix, cluster_rows = FALSE, cluster_cols = FALSE,
    show_rownames = FALSE, show_colnames = FALSE, scale = "none",
    color = colors, breaks = seq(-prepared$limit, prepared$limit, length.out = 101),
    gaps_col = prepared$gaps_col, gaps_row = prepared$gaps_row,
    border_color = NA, legend = FALSE, silent = TRUE
  )
  # Draw the matrix in its own viewport. Brackets share pheatmap's four-point gaps.
  matrix_grob <- heatmap$gtable$grobs[[which(heatmap$gtable$layout$name == "matrix")]]
  block_edges <- function(counts) {
    counts <- counts[counts > 0L]
    size <- (grid::unit(1, "npc") - grid::unit(4 * (length(counts) - 1L), "bigpts")) / sum(counts)
    ends <- cumsum(counts)
    list(start = (ends - counts) * size + grid::unit(4 * (seq_along(counts) - 1L), "bigpts"),
         end = ends * size + grid::unit(4 * (seq_along(counts) - 1L), "bigpts"))
  }
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    7, 3,
    heights = grid::unit.c(grid::unit(0.76, "inches"), grid::unit(0.52, "inches"),
                          grid::unit(1, "null"), grid::unit(0.38, "inches"),
                          grid::unit(0.58, "inches"), grid::unit(0.30, "inches"),
                          grid::unit(0.12, "inches")),
    widths = grid::unit.c(grid::unit(0.48, "inches"), grid::unit(1, "null"),
                         grid::unit(0.48, "inches"))
  )))
  on.exit(grid::popViewport())
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.text(plot_title, gp = grid::gpar(fontsize = 10, fontface = "bold", lineheight = 1.25))
  grid::popViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  column_edges <- block_edges(prepared$genotype_counts)
  for (i in seq_along(prepared$genotype_counts)) {
    left <- column_edges$start[i] + grid::unit(1, "bigpts")
    right <- column_edges$end[i] - grid::unit(1, "bigpts")
    grid::grid.lines(x = grid::unit.c(left, left, right, right), y = c(0.10, 0.26, 0.26, 0.10),
                     gp = grid::gpar(col = "#333333", lwd = 0.8))
    grid::grid.text(paste0(names(prepared$genotype_counts)[i], " (n=", prepared$genotype_counts[i], ")"),
                    x = (left + right) / 2, y = 0.60, gp = grid::gpar(fontsize = 10))
  }
  grid::popViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
  grid::grid.draw(matrix_grob)
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 3))
  row_edges <- block_edges(prepared$counts)
  present <- names(prepared$counts)[prepared$counts > 0L]
  for (i in seq_along(present)) {
    top <- grid::unit(1, "npc") - row_edges$start[i]
    bottom <- grid::unit(1, "npc") - row_edges$end[i]
    grid::grid.lines(x = grid::unit(c(3, 8, 8, 3), "bigpts"),
                     y = grid::unit.c(top, top, bottom, bottom),
                     gp = grid::gpar(col = direction_colors[present[i]], lwd = 1.3))
  }
  grid::popViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = 2))
  labels <- c("GOM", if (prepared$counts["Mixed"] > 0L) "Mixed", "LOM")
  label_counts <- as.integer(prepared$counts[labels])
  label_text <- paste(format(label_counts, big.mark = ",", trim = TRUE),
                       ifelse(label_counts == 1L, "CpG", "CpGs"), labels)
  label_grobs <- lapply(seq_along(labels), function(i) {
    grid::textGrob(label_text[i], gp = grid::gpar(col = direction_colors[labels[i]], fontsize = 11))
  })
  widths <- lapply(label_grobs, grid::grobWidth)
  gap <- grid::unit(16, "bigpts")
  total <- Reduce(`+`, widths) + (length(widths) - 1L) * gap
  x <- (grid::unit(1, "npc") - total) / 2
  for (i in seq_along(labels)) {
    grid::grid.draw(grid::editGrob(label_grobs[[i]], x = x + widths[[i]] / 2))
    x <- x + widths[[i]] + gap
  }
  grid::popViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 5, layout.pos.col = 2))
  grid::pushViewport(grid::viewport(width = grid::unit(3, "inches")))
  grid::grid.raster(as.raster(matrix(colors, nrow = 1)), x = 0.5, y = 0.78,
                    width = 1, height = grid::unit(6, "bigpts"), interpolate = FALSE)
  ticks <- seq(-prepared$limit, prepared$limit, length.out = 5)
  grid::grid.text(format(round(ticks, 3), trim = TRUE), x = seq(0, 1, length.out = 5),
                  y = 0.44, gp = grid::gpar(fontsize = 8))
  grid::grid.text("Mean-centered beta", y = 0.06, gp = grid::gpar(fontsize = 9))
  grid::popViewport(2)
  grid::pushViewport(grid::viewport(layout.pos.row = 6, layout.pos.col = 2))
  grid::grid.text(direction_note, gp = grid::gpar(fontsize = 8, col = "#555555"))
  grid::popViewport()
}

plot_dmp_heatmap <- function(probes, selected_samples, plot_title, output_file,
                             directions, betas_data, targets_data,
                             direction_note = "GOM: gain of methylation; LOM: loss of methylation (experimental vs WT)") {
  prepared <- prepare_dmp_heatmap(probes, selected_samples, betas_data, targets_data, directions)
  grDevices::pdf(output_file, width = 11, height = prepared$height)
  on.exit(grDevices::dev.off())
  draw_dmp_heatmap(prepared, plot_title, direction_note)
  invisible(prepared)
}
