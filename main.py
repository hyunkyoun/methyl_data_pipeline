from filter import remove_intensity_columns, combine_sample_files
from preprocessing.combat_norm import combat_normalize
from preprocessing.data_parsing import SampleReportParsing, CSV_parsing
import subprocess
import os
import pandas as pd

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

def get_sample_table(input_file_paths, output_file_path):
    """
    This function reads multiple sample table files, combines them, and writes the combined data to a CSV file.
    :param input_file_paths: Dictionary with file paths as keys and run numbers as values.
    :param output_file_path: Path to save the combined CSV file.
    """
    parser = SampleReportParsing(input_file_paths)
    parser.parse(output_file_path)
    print(f"Sample table combined and saved to {output_file_path}")

def bmiq():
    # Path to your custom R library
    R_LIBS_PATH = "./packages"

    # The R script you want to run
    R_SCRIPT_PATH = "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/bmiq/DoBMIQ.r"

    # List the required R packages (fixed to include actual dependencies)
    REQUIRED_PACKAGES = ["readxl", "openxlsx", "ggplot2", "dplyr", "tidyr", "sesame", "RPMM"]

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

    # Load the required libraries
    library(readxl)
    library(openxlsx)

    # Source the script to load the function
    source("{R_SCRIPT_PATH}")
    
    # Call the function with the correct paths
    DoBMIQ("data/probesample.xlsx", "combat_normalized.xlsx")
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

    
# This function applies combat normalization to the combined data from two sets of runs.
# Note: the final prototype of the pipeline will take the DataFrame as an input instead of file paths.
def combat(file1_2, file3_4):
    output_file = './data/combat_normalized_output.xlsx'

    combat_normalize(file1_2, file3_4, output_file)

def filter_and_split_idat_by_run(
    sample_file,
    idat_file,
    output_dir='./data/split_runs'
):
    # Load sample sheet
    sample_df = pd.read_csv(sample_file)
    
    # Extract run numbers from Index and create mapping
    sample_df['Run'] = sample_df['Index'].apply(lambda x: x.split('_')[0])
    sample_df['IDAT_Column'] = sample_df['Sentrix Barcode'].astype(str) + "_" + sample_df['Sample Section']
    
    # Create dictionary for renaming: {old_column_name: Sample ID}
    rename_dict = dict(zip(sample_df['IDAT_Column'], sample_df['Sample ID']))

    # Map from run number to list of Sample IDs
    run_to_sample_ids = sample_df.groupby('Run')['Sample ID'].apply(list).to_dict()
    
    # Load IDAT data
    idat_df = pd.read_csv(idat_file)

    # Rename columns using rename_dict
    idat_df.columns = [rename_dict.get(col, col) for col in idat_df.columns]

    # Rename first column to 'TargetID'
    idat_df.rename(columns={idat_df.columns[0]: 'TargetID'}, inplace=True)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Split by run and save
    for run, sample_ids in run_to_sample_ids.items():
        subset_cols = ['TargetID'] + sample_ids
        subset_df = idat_df[subset_cols]
        output_file = os.path.join(output_dir, f"run_{run}.csv")
        subset_df.to_csv(output_file, index=False)
        print(f"Saved: {output_file}")


def _extract_run_number(path):
    import re
    basename = os.path.basename(path)
    match = re.match(r"run_(\d+)\.csv", basename)
    if match:
        return match.group(1)
    else:
        raise ValueError(f"Could not extract run number from filename: {path}")
    
def _combine_filename(file1, file2):
    n1 = _extract_run_number(file1)
    n2 = _extract_run_number(file2)

    run_nums = sorted([n1, n2], key=int)
    combined = f"run_{run_nums[0]}_{run_nums[1]}.csv"
    return combined

def combine_by_run(file1, file2, output_dir='./data/split_runs'):
    output_file = os.path.join(output_dir, _combine_filename(file1, file2))

    print(_combine_filename(file1, file2))
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    key_column = df1.columns[0]

    combined_df = pd.merge(df1, df2, on=key_column, how='outer')
    combined_df.to_csv(output_file, index=False)

    print(f"Combined file saved as: {output_file}")

if __name__ == "__main__":
    # remove_intensity_columns("./data/parsed_output.csv", "./data/Mu Epic Run 5 FinalReport_ctrl-bkg AVGbeta.xlsx")
    # remove_intensity_columns("./data/Mu EPIC Run 4 10_2024/Mu EPIC Run 4 FinalReport ctrl-bkg AVGbeta _ Intensities.xlsx", "./data/Mu Epic Run 4 FinalReport_ctrl-bkg AVGbeta.xlsx")
    # remove_intensity_columns("./data/Mu EPIC Run 3 3-28-2022/Mu EPIC Run 3 FinalReport-ctrl-bkg Intensities_AVGbeta.xlsx", "./data/Mu Epic Run 3 FinalReport_ctrl-bkg AVGbeta.xlsx")
    # remove_intensity_columns("./data/Mu EPIC Run 2 RQ-022275 FINAL_02042022/Mu EPIC Run 2 FinalReport_ctrl_bkg Instensities_AVGbeta.xlsx", "./data/Mu Epic Run 2 FinalReport_ctrl-bkg AVGbeta.xlsx")
    # remove_intensity_columns("./data/Mu EPIC Run 1 5-24-2021/Mu EPIC Run 1 FinalReport_ctrl-bkg Intensities_AVGbeta.xlsx", "./data/Mu Epic Run 1 FinalReport_ctrl-bkg AVGbeta.xlsx")
    # df1 = pd.read_excel('./data/Mu Epic Run 1 FinalReport_ctrl-bkg AVGbeta.xlsx')
    # df2 = pd.read_excel('./data/Mu Epic Run 2 FinalReport_ctrl-bkg AVGbeta.xlsx')
    # df3 = pd.read_excel('./data/Mu Epic Run 3 FinalReport_ctrl-bkg AVGbeta.xlsx')
    # df4 = pd.read_excel('./data/Mu Epic Run 4 FinalReport_ctrl-bkg AVGbeta.xlsx')
    # df5 = pd.read_excel('./data/Mu Epic Run 5 FinalReport_ctrl-bkg AVGbeta.xlsx')


    # combined_df = combine_sample_files([df1, df2], "Runs_1_2.xlsx")
    # combined_df = combine_sample_files([df3, df4, df5], "Runs_3_4_5.xlsx")
    
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
        
    # # read_idat_files()
    # filter_and_split_idat_by_run('./data/sample_table_combined.csv', './data/filtered_beta_matrix.csv')

    # combine_by_run('./data/split_runs/run_1.csv', './data/split_runs/run_2.csv')
    # combine_by_run('./data/split_runs/run_3.csv', './data/split_runs/run_4.csv')


    # combat_normalize('./Runs_1_2.xlsx', './Runs_3_4_5.xlsx', 'combat_normalized.xlsx')
    bmiq()

    '''
    Combined raw shape: (285143, 135)
Beta data shape before NaN removal: (285143, 134)
NaN count: 1255340
Beta data shape after NaN removal: (21541, 134)
Removed 263602 TargetIDs due to NaN values'''