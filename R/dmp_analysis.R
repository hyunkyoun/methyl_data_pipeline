# Fit each experiment/tissue independently, with genotype-defined treatment arms.
run_dmp_analysis <- function(betas_data, targets_data, fdr_cutoff = 0.05, delta_beta_cutoff = 0.1,
                              same_parent_experiments = c("1", "2")) {
  required <- c("Samples", "Experiment", "Tissue.Cell.type", "True.Genotype")
  if (!all(required %in% names(targets_data)) ||
      any(vapply(targets_data[required], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1)))) {
    stop("DMP analysis requires non-missing Samples, Experiment, Tissue.Cell.type, and True.Genotype.")
  }
  indices <- match(colnames(betas_data), targets_data$Samples)
  if (anyNA(indices) || anyDuplicated(targets_data$Samples)) stop("DMP analysis: invalid sample metadata.")
  targets_data <- targets_data[indices, , drop = FALSE]
  cohorts <- unique(targets_data[c("Experiment", "Tissue.Cell.type")])
  results <- weights <- list()
  wt_genotypes <- c(DP10 = "WT_DP10", DP16 = "WT_DP16", DP10_DP16 = "WT_DP10_WT_DP16")
  wt_crosses <- c(DP10 = "DP10", DP16 = "DP16", DP10_DP16 = "Double")
  for (i in seq_len(nrow(cohorts))) {
    experiment <- as.character(cohorts$Experiment[i])
    tissue <- as.character(cohorts$Tissue.Cell.type[i])
    samples <- targets_data[targets_data$Experiment == experiment & targets_data$Tissue.Cell.type == tissue, , drop = FALSE]
    if (sum(samples$True.Genotype == "WT") < 2L) {
      warning("Skipping Experiment ", experiment, " / ", tissue, ": fewer than two same-cohort WT samples.", call. = FALSE)
      next
    }
    genotypes <- setdiff(unique(as.character(samples$True.Genotype)), "WT")
    if (!length(genotypes)) next
    # Legacy cross/group names may repeat across experiments; do not use them as arms.
    samples$group <- make.names(samples$True.Genotype)
    if (anyDuplicated(make.names(c("WT", genotypes)))) stop("Genotype labels collide after conversion to model names.")
    contrasts <- paste0("group", make.names(genotypes), " - WT_", tissue)
    if (experiment %in% same_parent_experiments) {
      wt <- samples$True.Genotype == "WT"
      if (!"Genotype" %in% names(samples) ||
          any(!samples$Genotype[wt] %in% unname(wt_genotypes))) {
        stop("Same-parent WT selection requires Genotype WT_DP10, WT_DP16, or WT_DP10_WT_DP16 for every WT in Experiment ", experiment, ".")
      }
      if (any(!genotypes %in% names(wt_genotypes))) {
        stop("No same-parent WT rule for a genotype in Experiment ", experiment, ".")
      }
      fitted <- list(results = list())
      for (j in seq_along(genotypes)) {
        # All animals of this genotype stay in the Dp arm, regardless of their cross.
        selected <- samples$True.Genotype == genotypes[j] |
          (wt & samples$Genotype %in% wt_genotypes[[genotypes[j]]])
        pair <- samples[selected, , drop = FALSE]
        pair_fit <- run_dmp_cohort(betas_data[, pair$Samples, drop = FALSE], pair, contrasts[j],
                                   fdr_cutoff, delta_beta_cutoff)
        fitted$results <- c(fitted$results, pair_fit$results)
      }
    } else {
      fitted <- run_dmp_cohort(betas_data[, samples$Samples, drop = FALSE], samples, contrasts,
                               fdr_cutoff, delta_beta_cutoff)
    }
    for (result in fitted$results) {
      id <- paste(result$genotype, "vs_WT", paste0("Exp", experiment), tissue, sep = "_")
      id <- gsub("[^A-Za-z0-9_-]", "_", id)
      if (id %in% names(results)) stop("Duplicate comparison filename: ", id)
      result$experiment <- experiment
      result$contrast_name <- id
      result$contrast_expression <- paste(result$genotype, "- WT")
      result$control_group <- "WT"
      result$wt_cross <- if (experiment %in% same_parent_experiments) wt_crosses[[result$genotype]] else "All"
      for (field in c("table", "dmps")) {
        result[[field]]$Contrast <- rep(id, nrow(result[[field]]))
        result[[field]]$Comparison <- rep(result$contrast_expression, nrow(result[[field]]))
        result[[field]]$Experiment <- rep(experiment, nrow(result[[field]]))
        result[[field]]$Tissue <- rep(tissue, nrow(result[[field]]))
        result[[field]]$WT.Cross <- rep(result$wt_cross, nrow(result[[field]]))
      }
      results[[id]] <- result
      weights[[id]] <- data.frame(Comparison = id, Experiment = experiment, Tissue = tissue,
                                  WT.Cross = result$wt_cross,
                                  Genotype = c(result$genotype, "WT"),
                                  Sample.Count = c(length(result$experimental_samples), length(result$control_samples)),
                                  Coefficient = c(1, -1))
    }
  }
  list(results = results, contrast_weights = if (length(weights)) do.call(rbind, weights) else data.frame())
}

# Existing limma fit, applied only to one experiment/tissue by run_dmp_analysis.
run_dmp_cohort <- function(betas_data, targets_data, contrasts_list, fdr_cutoff = 0.05, delta_beta_cutoff = 0.1) {
  # QC and ComBat have checked the input data; keep sample order aligned here.
  sample_indices <- match(colnames(betas_data), targets_data$Samples)
  if (anyNA(sample_indices)) stop("DMP analysis: missing sample metadata.")
  targets_data <- targets_data[sample_indices, , drop = FALSE]

  # Pool WT before fitting so every retained WT sample has equal influence.
  is_wt <- targets_data$True.Genotype == "WT"
  group <- factor(ifelse(is_wt, paste0("WT_", targets_data$Tissue.Cell.type),
                         paste0("group", targets_data$group)))
  analyses <- list()
  valid_contrasts <- character()

  # Record each comparison's samples for the result tables and plots.
  for (contrast in contrasts_list) {
    terms <- trimws(strsplit(contrast, "-", fixed = TRUE)[[1]])
    if (length(terms) != 2L || !startsWith(terms[1], "group") || !startsWith(terms[2], "WT_")) {
      stop("Expected a contrast of the form 'groupNAME - WT_tissue': ", contrast)
    }
    exp_group <- sub("^group", "", terms[1])
    experimental <- group == terms[1]
    control <- group == terms[2]
    exp_samples <- targets_data$Samples[experimental]
    ctrl_samples <- targets_data$Samples[control]
    if (length(exp_samples) < 2 || length(ctrl_samples) < 2) {
      warning("Skipping ", contrast, ": fewer than two samples in the experimental or WT group.", call. = FALSE)
      next
    }
    comparison_tissue <- unique(targets_data$Tissue.Cell.type[experimental])
    comparison_genotype <- unique(targets_data$True.Genotype[experimental])
    if (length(comparison_tissue) != 1L || length(comparison_genotype) != 1L) {
      stop("Samples within each experimental group must share the same True.Genotype and tissue labels: ", contrast)
    }
    if (terms[2] != paste0("WT_", comparison_tissue)) {
      stop("WT controls must come from the experimental tissue: ", contrast)
    }
    contrast_name <- paste(exp_group, "vs", terms[2], sep = "_")

    analysis_id <- paste0(contrast_name, "__equal_mouse")
    valid_contrasts <- c(valid_contrasts, contrast)
    analyses[[analysis_id]] <- list(
      contrast_name = contrast_name,
      contrast_expression = paste(terms, collapse = " - "),
      tissue = comparison_tissue,
      genotype = comparison_genotype,
      experimental_group = exp_group,
      control_group = paste("All", comparison_tissue, "WT"),
      weighting = "equal_mouse",
      weighting_label = "Equal weight per sample",
      experimental_samples = exp_samples,
      control_samples = ctrl_samples
    )
  }

  contrast_weights <- data.frame(
    Group = levels(group),
    targets_data[match(levels(group), group), c("Tissue.Cell.type", "True.Genotype"), drop = FALSE],
    Sample.Count = as.integer(table(group)), row.names = NULL
  )
  if (length(analyses) == 0) {
    return(list(results = analyses, contrast_weights = contrast_weights))
  }
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  contrast_matrix <- limma::makeContrasts(contrasts = valid_contrasts, levels = design)
  colnames(contrast_matrix) <- names(analyses)
  contrast_weights <- cbind(contrast_weights, as.data.frame(contrast_matrix))
  fit <- limma::lmFit(betas_data, design)
  fit_contrasts <- limma::contrasts.fit(fit, contrast_matrix)
  fit_moderated <- limma::eBayes(fit_contrasts)

  for (analysis_id in names(analyses)) {
    result <- analyses[[analysis_id]]
    tt <- limma::topTable(fit_moderated, coef = analysis_id, number = Inf,
                          adjust.method = "BH", sort.by = "none")
    tt$Probe <- rownames(tt)
    tt$diff <- tt$logFC  # Beta values were fitted, so this is delta beta.
    tt$Contrast <- result$contrast_name
    tt$Comparison <- result$contrast_expression
    tt$Weighting <- result$weighting
    tt$Significant <- is.finite(tt$adj.P.Val) & is.finite(tt$diff) &
      tt$adj.P.Val < fdr_cutoff & abs(tt$diff) > delta_beta_cutoff
    result$table <- tt
    result$dmps <- tt[tt$Significant, , drop = FALSE]
    analyses[[analysis_id]] <- result
  }

  list(results = analyses, contrast_weights = contrast_weights)
}
