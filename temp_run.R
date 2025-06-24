
    .libPaths("./packages")

    # List required packages
    packages <- c("ggplot2", "dplyr", "tidyr", "sesame")

    # Install any missing packages
    new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(new_packages) > 0) {
        install.packages(new_packages, lib = "./packages", repos = "https://cloud.r-project.org/")
    }

    # Now run the actual script
    source("./idat.r")
    