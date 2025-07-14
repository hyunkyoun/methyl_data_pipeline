import health_check as hc
import subprocess
import os

def get_idat_files(source_directory):
    output_idat_directory = './idat'
    hc.collect_idat_files(source_directory, output_idat_directory)
    idat_count = hc.count_idat_files(output_idat_directory)
    print(f"Collected {idat_count} IDAT files to {output_idat_directory}")

def read_idat_files():
    # Path to your custom R library
    R_LIBS_PATH = "./packages"

    # The R script you want to run
    R_SCRIPT_PATH = "./idat.r"

    # List the required R packages
    REQUIRED_PACKAGES = ["ggplot2", "dplyr", "tidyr", "sesame"]

    # Build the R script
    R_COMMAND = f"""
    .libPaths("{R_LIBS_PATH}")

    # List required packages
    packages <- c({', '.join(f'"{pkg}"' for pkg in REQUIRED_PACKAGES)})

    # Install any missing packages
    new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(new_packages) > 0) {{
        install.packages(new_packages, lib = "{R_LIBS_PATH}", repos = "https://cloud.r-project.org/")
    }}

    # Now run the actual script
    source("{R_SCRIPT_PATH}")
    """

    # Save the combined R command
    temp_file = "temp_run.R"
    with open(temp_file, "w") as f:
        f.write(R_COMMAND)

    # Run R script
    subprocess.run(["Rscript", temp_file], check=True)

    # Optional: remove temporary script
    os.remove(temp_file)

    print(f"✅ Done running {R_SCRIPT_PATH} with R libraries at {R_LIBS_PATH}")


if __name__ == "__main__":
    source_directory = 'C:\\Users\\ekim2\\Desktop\\tycko research\\DS idat files\\Down_Syndrome_Mouse_Microarray_Samples'
    get_idat_files(source_directory)
    print("IDAT files collection completed.")

    read_idat_files()
    print("IDAT files reading completed.")