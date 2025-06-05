from filter import remove_intensity_columns, combine_sample_files
from preprocess import combat_normalize

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
    
# This function applies combat normalization to the combined data from two sets of runs.
# Note: the final prototype of the pipeline will take the DataFrame as an input instead of file paths.
def combat():
    file1_2 = './data/combined_run1_run2.xlsx'
    file3_4 = './data/combined_run3_run4.xlsx'
    output_file = './data/combat_normalized_output.xlsx'

    combat_normalize(file1_2, file3_4, output_file)
    

if __name__ == "__main__":
    df = main()
    combat()