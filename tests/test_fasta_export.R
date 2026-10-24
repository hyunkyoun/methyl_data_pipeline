# Run from the repository root: Rscript --vanilla tests/test_fasta_export.R
test_library <- Sys.getenv("METHYL_DMP_TEST_LIBRARY")
if (nzchar(test_library)) .libPaths(c(test_library, .libPaths()))
source("R/output_paths.R")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
})

# Test the actual exporter with a tiny indexed reference and real output files.
workflow_lines <- readLines("workflow_new.R", warn = FALSE)
first <- grep("# FASTA sequences from DMPS", workflow_lines, fixed = TRUE)
last <- grep("# Extract FASTA from the stored DMP calls", workflow_lines, fixed = TRUE)
stopifnot(length(first) == 1, length(last) == 1, first < last)
exporter <- new.env(parent = globalenv())
eval(parse(text = workflow_lines[seq.int(first, last - 1)]), exporter)
exporter$getSeq <- function(genome, regions) {
  # An empty category must never reach sequence extraction.
  stopifnot(length(regions) > 0)
  exporter$extraction_calls <- exporter$extraction_calls + 1L
  Biostrings::getSeq(genome, regions)
}

run_tests <- function() {
  original_dir <- getwd()
  scratch <- tempfile(pattern = "methyl-fasta-test-")
  stopifnot(dir.create(scratch))
  on.exit(setwd(original_dir))
  setwd(scratch)
  output_dirs <- experiment_output_dirs(create_output_dirs("results", "analysis_test"), 1)

  reference <- paste(rep("ACGTTCGA", 500), collapse = "")
  writeXStringSet(DNAStringSet(c(chr1 = reference)), "reference.fasta")
  indexFa("reference.fasta")
  genome <- FaFile("reference.fasta")
  annotation <- data.frame(Name = c("cg0001_BC21", "cg0002_TC21", "cg0003_BC21"),
                           chr = "chr1", pos = c(300L, 350L, 1000L))
  gains <- data.frame(Probe = c("cg0001", "cg0002"), diff = c(0.2, 0.3))
  losses <- data.frame(Probe = "cg0003", diff = -0.25)
  both <- rbind(gains, losses)

  run_export <- function(dmps, contrast) {
    exporter$extraction_calls <- 0L
    messages <- capture.output(
      result <- exporter$extract_and_write_FASTA(
        dmps, contrast, annotation, genome,
        fasta_dir = output_dirs$fasta,
        details_dir = output_dirs$fasta_details
      )
    )
    list(summary = result, messages = messages, calls = exporter$extraction_calls)
  }
  output_path <- function(contrast, direction, suffix) {
    directory <- if (suffix == ".fasta") output_dirs$fasta else output_dirs$fasta_details
    file.path(directory, direction, paste0(contrast, suffix))
  }
  check_direction <- function(contrast, direction, probes, header = NULL, start = NULL, end = NULL) {
    path <- output_path(contrast, direction, ".fasta")
    stopifnot(file.exists(path))
    sequences <- readDNAStringSet(path)
    details <- read.csv(output_path(contrast, direction, ".csv"))
    stopifnot(identical(names(details), c("Probe", "diff", "chr", "pos")),
              identical(as.character(details$Probe), probes))
    if (length(probes) == 0) {
      stopifnot(length(sequences) == 0L, file.size(path) == 0, nrow(details) == 0L)
    } else {
      stopifnot(length(sequences) == 1L, identical(names(sequences), header),
                identical(unname(as.character(sequences)), substr(reference, start, end)))
    }
  }
  check_counts <- function(result, gains, losses, written, calls) {
    stopifnot(result$summary$GOM == gains, result$summary$LOM == losses,
              result$summary$written == written, result$calls == calls)
  }

  gain_only <- run_export(gains, "gain_only")
  # Overlapping gain windows still merge: two DMPs produce one 250-base sequence.
  check_counts(gain_only, 2, 0, 1, 1)
  check_direction("gain_only", "GOM", gains$Probe, "chr1:201-450", 201, 450)
  check_direction("gain_only", "LOM", character())
  stopifnot(any(grepl("No methylation losses (LOM)", gain_only$messages, fixed = TRUE)))

  loss_only <- run_export(losses, "loss_only")
  check_counts(loss_only, 0, 1, 1, 1)
  check_direction("loss_only", "GOM", character())
  check_direction("loss_only", "LOM", losses$Probe, "chr1:901-1100", 901, 1100)
  stopifnot(any(grepl("No methylation gains (GOM)", loss_only$messages, fixed = TRUE)))

  mixed <- run_export(both, "mixed")
  check_counts(mixed, 2, 1, 2, 2)
  check_direction("mixed", "GOM", gains$Probe, "chr1:201-450", 201, 450)
  check_direction("mixed", "LOM", losses$Probe, "chr1:901-1100", 901, 1100)

  # A direction can become empty only after annotation filtering.
  missing_annotation <- run_export(rbind(gains, data.frame(Probe = "cg9999", diff = -0.2)),
                                   "missing_annotation")
  check_counts(missing_annotation, 2, 0, 1, 1)
  check_direction("missing_annotation", "LOM", character())

  # Rerunning the same comparison must clear previously populated empty-side files.
  rerun <- run_export(both, "rerun")
  check_counts(rerun, 2, 1, 2, 2)
  rerun <- run_export(gains, "rerun")
  check_counts(rerun, 2, 0, 1, 1)
  check_direction("rerun", "GOM", gains$Probe, "chr1:201-450", 201, 450)
  check_direction("rerun", "LOM", character())
  rerun <- run_export(losses, "rerun")
  check_counts(rerun, 0, 1, 1, 1)
  check_direction("rerun", "GOM", character())
  check_direction("rerun", "LOM", losses$Probe, "chr1:901-1100", 901, 1100)

  cat("FASTA export tests passed: gains only, losses only, mixed directions, annotation filtering,",
      "merged windows, sequence content, counts, and clearing stale empty-category outputs.\n")
}

run_tests()
