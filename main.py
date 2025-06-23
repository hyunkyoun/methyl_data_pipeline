from filter import remove_intensity_columns, combine_sample_files
from preprocessing.combat_norm import combat_normalize
from preprocessing.data_parsing import SampleReportParsing

# This script processes two sets of sample files, removes intensity columns, combines them, and applies combat normalization.
# It is the main entry point for the pipeline for data analysis. 
def main():
    file1 = './data/Mu EPIC Run 3 3-28-2022/Mu EPIC Run 3 FinalReport-ctrl-bkg Intensities_AVGbeta.xlsx'
    file2 = './data/Mu EPIC Run 4 10_2024/Mu EPIC Run 4 FinalReport ctrl-bkg AVGbeta _ Intensities.xlsx'

    removed_intensity_output_file1 = './data/removed_intensity_values/Mu EPIC Run 3 FinalReport_ctrl-bkg AVGbeta.xlsx'
    removed_intensity_output_file2 = './data/removed_intensity_values/Mu EPIC Run 4 FinalReport_ctrl_bkg AVGbeta.xlsx'

    df1 = remove_intensity_columns(file1, removed_intensity_output_file1)
    df2 = remove_intensity_columns(file2, removed_intensity_output_file2)

    final_output_path = './data/combined_run3_run4.xlsx'

    final_df = combine_sample_files(df1, df2, final_output_path)

    return final_df

def read_idat_files():
    import subprocess
    import os

    # Path to your custom R library
    R_LIBS_PATH = "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/packages"

    # The R script you want to run
    R_SCRIPT_PATH = "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/idat.r"

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

    def get_sample_table(input_file_paths, output_file_path):
        """
        This function reads multiple sample table files, combines them, and writes the combined data to a CSV file.
        :param input_file_paths: Dictionary with file paths as keys and run numbers as values.
        :param output_file_path: Path to save the combined CSV file.
        """
        parser = SampleReportParsing(input_file_paths)
        parser.parse(output_file_path)
        print(f"Sample table combined and saved to {output_file_path}")


    
# This function applies combat normalization to the combined data from two sets of runs.
# Note: the final prototype of the pipeline will take the DataFrame as an input instead of file paths.
def combat():
    file1_2 = './data/combined_run1_run2.xlsx'
    file3_4 = './data/combined_run3_run4.xlsx'
    output_file = './data/combat_normalized_output.xlsx'

    combat_normalize(file1_2, file3_4, output_file)
    

if __name__ == "__main__":
    # df = main()
    # combat()

    # input_file_paths = {
    #     "./data/Mu EPIC Run 1 5-24-2021/SamplesTableFinalReport.txt": 1,
    #     "./data/Mu EPIC Run 2 RQ-022275 FINAL_02042022/TableControl.txt": 2,
    #     "./data/Mu EPIC Run 3 3-28-2022/SamplesTable.txt": 3,
    #     "./data/Mu EPIC Run 4 10_2024/SamplesTable.txt": 4,
    # }
    # sample_table_output_path = './data/sample_table_combined.csv'

    # get_sample_table(input_file_paths, sample_table_output_path)
    # read_idat_files()
    
    read_idat_files()