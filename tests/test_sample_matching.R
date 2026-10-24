# Run from the repository root: Rscript --vanilla tests/test_sample_matching.R
# These tests use base R and never read IDAT intensities or overwrite input files.
workflow_lines <- readLines("workflow_new.R", warn = FALSE)
workflow_block <- function(first_marker, last_marker) {
  first <- match(first_marker, workflow_lines)
  last <- match(last_marker, workflow_lines)
  stopifnot(!is.na(first), !is.na(last), first < last)
  parse(text = workflow_lines[seq.int(first, last - 1)])
}
load_samples <- workflow_block("# Load analysis samples before choosing IDAT pairs; brain is excluded throughout.",
                               "# Saving files")
before_qc <- workflow_block("common_sample_ids <- intersect(colnames(betas), targets_raw$Basename)", "# By run QC")
after_qc <- workflow_block("# Samplesheet matching", "# Technology between old vs new")

expect_error <- function(expr, text) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}

read_fixture <- function(sheet, betas) {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  write.csv(sheet, path, row.names = FALSE, fileEncoding = "UTF-8")
  # Include an Excel-style UTF-8 BOM to exercise the actual CSV reader.
  csv_lines <- readLines(path, encoding = "UTF-8")
  csv_lines[1] <- paste0("\ufeff", csv_lines[1])
  writeLines(csv_lines, path, useBytes = TRUE)
  env <- new.env(parent = globalenv())
  env$samplesheet_path <- path
  env$betas <- betas
  env$detP <- betas / 10
  env$idir <- tempdir()
  env$searchIDATprefixes <- function(idir) setNames(file.path(idir, colnames(betas)), colnames(betas))
  eval(load_samples, env)
  eval(before_qc, env)
  env
}

basenames <- c("205243950018_R01C01", "205243950019_R01C01", "205243950020_R02C01")
sheet <- data.frame(
  Sample.ID = c(101L, 102L, 103L), Experiment = c(1L, 1L, 2L),
  Beadchip = rep("2.05E+11", 3), Position = c("R01C01", "R01C01", "R02C01"),
  Basename = basenames, EPIC.Run = c("1_2", "3_4", "5"),
  Tissue.Cell.type = "Tcells",
  group = c("WT", "DP10", "DP16")
)
extra <- "205243950021_R01C01"
betas <- matrix(seq_len(12) / 100, nrow = 3,
                dimnames = list(paste0("cg", 1:3), c(basenames[2], extra, basenames[3], basenames[1])))

# Distinct full IDs survive even when rounded chip IDs and positions collide.
result <- read_fixture(sheet, betas)
expected_betas <- betas[, c(1, 3, 4), drop = FALSE]
stopifnot(identical(result$targets_raw$Basename, basenames),
          identical(result$betas, expected_betas),
          identical(result$detP, expected_betas / 10),
          identical(result$targets_matched$Samples, c("102_1", "103_2", "101_1")),
          identical(result$targets_matched$group, c("DP10", "DP16", "WT")),
          !file.exists(result$samplesheet_path))

# Mimic sample removal/reordering by QC. The deleted CSV cannot be read again.
result$betas <- result$betas[, c(3, 1), drop = FALSE]
result$detP <- result$detP[, c(3, 1), drop = FALSE]
invisible(capture.output(eval(after_qc, result)))
stopifnot(identical(result$targets$Basename, basenames[1:2]),
          identical(result$targets$Samples, c("101_1", "102_1")),
          identical(colnames(result$betas), c("101_1", "102_1")),
          identical(result$targets$group, c("WT", "DP10")),
          identical(unname(result$betas), unname(betas[, c(4, 1)])))

# Beadchip and Position are no longer needed for matching.
without_chip <- read_fixture(sheet[, !names(sheet) %in% c("Beadchip", "Position")], betas)
stopifnot(identical(without_chip$betas, expected_betas))
single <- read_fixture(sheet[1, ], betas[, 4, drop = FALSE])
stopifnot(is.matrix(single$betas), ncol(single$betas) == 1L,
          is.data.frame(single$targets_matched), nrow(single$targets_matched) == 1L)

expect_error(read_fixture(sheet[, names(sheet) != "Basename"], betas), "requires Basename")
for (bad_id in c(NA_character_, "", "  ")) {
  invalid <- sheet
  invalid$Basename[1] <- bad_id
  expect_error(read_fixture(invalid, betas), "missing or blank Basename")
}
duplicate <- sheet
duplicate$Basename[2] <- duplicate$Basename[1]
expect_error(read_fixture(duplicate, betas), "duplicate Basename")
unknown <- sheet
unknown$Basename[1] <- "999999999999_R01C01"
expect_error(read_fixture(unknown, betas), "absent from complete IDAT pairs")
colnames(single$betas) <- "unexpected_QC_sample"
expect_error(eval(after_qc, single), "without a matching samplesheet Basename")

# Brain must be removed before IDAT processing, QC, batch correction, or plotting.
with_brain <- sheet
with_brain$Tissue.Cell.type[2] <- " Brain "
without_brain <- read_fixture(with_brain, betas)
stopifnot(!basenames[2] %in% without_brain$targets_raw$Basename,
          !basenames[2] %in% names(without_brain$idat_prefixes),
          !basenames[2] %in% colnames(without_brain$betas),
          !basenames[2] %in% colnames(without_brain$detP),
          ncol(without_brain$betas) == 2L)

# Check the configured local samplesheet against both channels of the real IDATs.
config <- new.env(parent = globalenv())
eval(parse(text = workflow_lines[grep("^(idir|samplesheet_path) <-", workflow_lines)]), config)
if (file.exists(config$samplesheet_path)) {
  green <- sub("_Grn\\.idat$", "", list.files(config$idir, pattern = "_Grn\\.idat$"))
  red <- sub("_Red\\.idat$", "", list.files(config$idir, pattern = "_Red\\.idat$"))
  paired_ids <- intersect(green, red)
  env <- new.env(parent = globalenv())
  env$samplesheet_path <- config$samplesheet_path
  env$betas <- matrix(0.5, nrow = 3, ncol = length(paired_ids),
                     dimnames = list(paste0("cg", 1:3), paired_ids))
  env$detP <- env$betas / 100
  env$idir <- config$idir
  env$searchIDATprefixes <- function(idir) setNames(file.path(idir, paired_ids), paired_ids)
  eval(load_samples, env)
  eval(before_qc, env)
  saved_ids <- colnames(env$betas)
  invisible(capture.output(eval(after_qc, env)))
  stopifnot(identical(env$targets$Basename, saved_ids),
            identical(colnames(env$betas), env$targets$Samples),
            nrow(env$targets_raw) == ncol(env$betas), nrow(env$targets_raw) == 161L,
            !any(env$targets$Tissue.Cell.type == "brain"), !any(env$targets$Experiment == 6),
            env$targets$True.Genotype[env$targets$Samples == "12926_8"] == "DP10_DP16")
  cat(sprintf("Local samplesheet: %d Basenames match complete IDAT pairs; %d unlisted pairs excluded.\n",
              ncol(env$betas), length(paired_ids) - ncol(env$betas)))
}
cat("Sample matching tests passed: exact Basenames, brain excluded before IDAT selection, corrected genotype, BOM, QC ordering, and invalid identifiers.\n")
