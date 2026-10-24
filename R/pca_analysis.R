# Use one palette across experiment/tissue subsets.
pca_color_palettes <- function(targets_data) {
  batch_levels <- unique(as.character(targets_data$BATCH))
  tissue_levels <- unique(as.character(targets_data$Tissue.Cell.type))
  genotype_levels <- unique(as.character(targets_data$True.Genotype))
  list(
    batch = setNames(c("red", "blue")[seq_along(batch_levels)], batch_levels),
    tissue = setNames(grDevices::rainbow(length(tissue_levels)), tissue_levels),
    genotype = setNames(c("green", "orange", "purple", "red", "blue", "brown",
                         "pink", "gray", "cyan", "magenta", "yellow", "black")
                       [seq_along(genotype_levels)], genotype_levels)
  )
}

create_simple_pca <- function(betas_data, targets_data, ntop = 100, output_prefix,
                              figure_dir, table_dir, scope_label = NULL,
                              plot_tissue = TRUE, color_palettes = NULL, plot_title = NULL) {
  required_metadata <- c("Samples", "BATCH", "Tissue.Cell.type", "True.Genotype")
  if (!all(required_metadata %in% names(targets_data))) {
    stop("PCA metadata must contain: ", paste(required_metadata, collapse = ", "))
  }
  if (anyNA(targets_data$Samples) || anyDuplicated(targets_data$Samples) ||
      is.null(colnames(betas_data)) || anyDuplicated(colnames(betas_data))) {
    stop("PCA requires unique, non-missing sample identifiers.")
  }
  sample_index <- match(targets_data$Samples, colnames(betas_data))
  if (anyNA(sample_index)) stop("PCA metadata contains samples absent from the beta matrix.")
  betas_data <- betas_data[, sample_index, drop = FALSE]
  if (ncol(betas_data) < 3L) {
    warning("Skipping PCA for ", output_prefix, ": fewer than three samples.")
    return(NULL)
  }
  if (length(ntop) != 1L || !is.finite(ntop) || ntop < 2L || ntop != as.integer(ntop)) {
    stop("ntop must be an integer of at least two.")
  }
  if (is.null(rownames(betas_data)) || anyDuplicated(rownames(betas_data))) {
    stop("PCA requires unique CpG identifiers as beta-matrix row names.")
  }

  cat("Processing", output_prefix, "data...\n")
  # Rank CpGs using only the samples supplied to this PCA.
  cpg_variances <- apply(betas_data, 1L, var, na.rm = TRUE)
  variable_cpgs <- cpg_variances[is.finite(cpg_variances) & cpg_variances > 0]
  top_cpgs <- names(head(sort(variable_cpgs, decreasing = TRUE), ntop))
  n_selected <- length(top_cpgs)
  if (n_selected < 2L) {
    warning("Skipping PCA for ", output_prefix, ": fewer than two variable CpGs.")
    return(NULL)
  }
  betas_top <- betas_data[top_cpgs, , drop = FALSE]
  scaled_data <- scale(t(betas_top), center = TRUE, scale = TRUE)
  scaled_data[!is.finite(scaled_data)] <- 0
  pca <- prcomp(scaled_data, center = FALSE, scale. = FALSE)
  pca_df <- cbind(data.frame(PC1 = pca$x[, 1L], PC2 = pca$x[, 2L],
                            Sample = rownames(pca$x)), targets_data)
  var_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  if (is.null(color_palettes)) color_palettes <- pca_color_palettes(targets_data)

  for (path in c(figure_dir, table_dir)) {
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
      stop("Could not create PCA output directory: ", path)
    }
  }
  create_plot <- function(color_var, color_palette, label) {
    title <- if (!is.null(plot_title)) {
      plot_title
    } else if (is.null(scope_label)) {
      paste("PCA of Top", n_selected, "CpGs - PC1 vs PC2 by", label)
    } else {
      paste0(scope_label, ": PCA by ", label)
    }
    subtitle <- if (is.null(scope_label)) NULL else {
      paste0(if (!is.null(plot_title)) paste0(scope_label, " | ") else "",
             paste(nrow(targets_data), "samples | Top", format(n_selected, big.mark = ","),
                   "variable CpGs | ComBat adjusted"))
    }
    ggplot2::ggplot(pca_df, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]],
                                       color = .data[[color_var]])) +
      ggplot2::geom_point(size = 2.5, alpha = 0.7) +
      ggplot2::scale_color_manual(values = color_palette) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid = ggplot2::element_line(color = "grey90"),
                     plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = title, subtitle = subtitle,
                    x = paste0("PC1 (", round(var_explained[1L], 2), "%)"),
                    y = paste0("PC2 (", round(var_explained[2L], 2), "%)"))
  }
  plots <- list(
    Batch = create_plot("BATCH", color_palettes$batch, "Batch"),
    Genotype = create_plot("True.Genotype", color_palettes$genotype, "Genotype")
  )
  if (plot_tissue) {
    plots$Tissue <- create_plot("Tissue.Cell.type", color_palettes$tissue, "Tissue")
  }
  for (label in names(plots)) {
    plot_dir <- file.path(figure_dir, tolower(label))
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(plot_dir, paste0(output_prefix, ".pdf")),
                    plots[[label]], width = 8, height = 6)
  }

  for (kind in c("betas", "cpgs", "scores", "variance")) {
    dir.create(file.path(table_dir, kind), recursive = TRUE, showWarnings = FALSE)
  }
  write.csv(betas_top, file.path(table_dir, "betas", paste0(output_prefix, ".csv")),
            row.names = TRUE)
  write.csv(data.frame(CpG = top_cpgs, Variance = cpg_variances[top_cpgs]),
            file.path(table_dir, "cpgs", paste0(output_prefix, ".csv")),
            row.names = FALSE)
  write.csv(pca_df, file.path(table_dir, "scores", paste0(output_prefix, ".csv")), row.names = FALSE)
  write.csv(data.frame(PC = paste0("PC", seq_along(var_explained)),
                       Variance.Percent = var_explained,
                       Cumulative.Percent = cumsum(var_explained)),
            file.path(table_dir, "variance", paste0(output_prefix, ".csv")), row.names = FALSE)
  cat("Selected", n_selected, "CpGs for", nrow(targets_data), "samples; PC1 + PC2 explain",
      round(sum(var_explained[1:2]), 2), "% of variance.\n")
  list(pca_df = pca_df, pca = pca, top_cpgs = top_cpgs, betas_top = betas_top,
       color_palettes = color_palettes)
}

create_tissue_pcas <- function(betas_data, targets_data, ntop = 10000,
                               output_prefix = "COMBAT_ADJUSTED", figure_dir, table_dir,
                               color_palettes = pca_color_palettes(targets_data)) {
  experiments <- unique(as.character(targets_data$Experiment))
  if (length(experiments) != 1L || anyNA(experiments) || !nzchar(trimws(experiments))) {
    stop("Tissue-specific PCA requires exactly one Experiment.")
  }
  tissues <- unique(as.character(targets_data$Tissue.Cell.type))
  if (length(tissues) == 0L || anyNA(tissues) || any(!nzchar(trimws(tissues)))) {
    stop("Tissue-specific PCA requires non-empty Tissue.Cell.type labels.")
  }
  folder_names <- make.unique(gsub("[^A-Za-z0-9_-]", "_", tissues), sep = "_")
  results <- setNames(vector("list", length(tissues)), tissues)
  for (i in seq_along(tissues)) {
    tissue <- tissues[i]
    tissue_targets <- targets_data[targets_data$Tissue.Cell.type == tissue, , drop = FALSE]
    # Rank CpGs independently within each experiment/tissue cohort.
    result <- create_simple_pca(
      betas_data, tissue_targets, ntop, paste0(output_prefix, "_", folder_names[i]),
      figure_dir = figure_dir, table_dir = table_dir,
      scope_label = paste0("Experiment ", experiments, " | ", tissue),
      plot_tissue = FALSE, color_palettes = color_palettes
    )
    results[i] <- list(result)
  }
  results
}
