# Run once from the methyl_data_pipeline working directory:
# source("install_packages.R")

if (!file.exists("workflow_new.R")) {
  stop("Set your working directory to methyl_data_pipeline before running this script.")
}

package_dir <- file.path(getwd(), "packages")
if (!dir.exists(package_dir)) {
  dir.create(package_dir, recursive = TRUE)
}
if (!dir.exists(package_dir)) {
  stop("Could not create the package library: ", package_dir)
}
package_dir <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)
.libPaths(c(package_dir, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
# Genome and annotation packages can take longer than R's default timeout.
options(timeout = max(600, getOption("timeout")))
message("Installing packages into: ", package_dir)

# Check this project's library, even if a package exists elsewhere on the computer.
missing_local_packages <- function(packages) {
  setdiff(packages, rownames(installed.packages(lib.loc = package_dir)))
}

cran_packages <- c(
  "BiocManager", "remotes", "ggsci", "snow", "RColorBrewer",
  "readxl", "openxlsx", "ggplot2", "gridExtra", "pheatmap"
)
to_install <- missing_local_packages(cran_packages)
if (length(to_install) > 0L) {
  install.packages(to_install, lib = package_dir)
}

bioc_packages <- c(
  "sesame", "minfi", "limma", "sva", "wateRmelon", "BiocParallel",
  "GenomicRanges", "GenomicFeatures", "BSgenome.Mmusculus.UCSC.mm10",
  "Biostrings"
)
to_install <- missing_local_packages(bioc_packages)
if (length(to_install) > 0L) {
  # force ensures packages found in other libraries are installed here too.
  BiocManager::install(to_install, lib = package_dir,
                       ask = FALSE, update = FALSE, force = TRUE)
}

github_packages <- c(
  IlluminaMouseMethylationmanifest =
    "chiaraherzog/IlluminaMouseMethylationmanifest",
  IlluminaMouseMethylationanno.12.v1.mm10 =
    "chiaraherzog/IlluminaMouseMethylationanno.12.v1.mm10"
)
to_install <- unname(github_packages[missing_local_packages(names(github_packages))])
if (length(to_install) > 0L) {
  BiocManager::install(to_install, lib = package_dir,
                       ask = FALSE, update = FALSE, force = TRUE)
}

# parallel is included with R and does not need installation.
still_missing <- missing_local_packages(c(cran_packages, bioc_packages,
                                          names(github_packages)))
if (length(still_missing) > 0L) {
  stop("Some packages were not installed. Check the messages above: ",
       paste(still_missing, collapse = ", "))
}
message("Package installation complete. The workflow will use: ", package_dir)
