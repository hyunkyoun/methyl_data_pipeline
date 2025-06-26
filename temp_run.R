
    .libPaths("./packages")

    # List required packages
    packages <- c("readxl", "openxlsx", "ggplot2", "dplyr", "tidyr", "sesame", "RPMM")

    # Install any missing packages
    new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(new_packages) > 0) {
        install.packages(new_packages, lib = "./packages", repos = "https://cloud.r-project.org/")
    }

    # Load the required libraries
    library(readxl)
    library(openxlsx)

    # Source the script to load the function
    source("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/bmiq/DoBMIQ.r")
    
    # Call the function with the correct paths
    DoBMIQ("data/probesample.xlsx", "combat_normalized.xlsx")
    