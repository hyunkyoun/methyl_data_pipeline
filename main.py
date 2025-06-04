from filter import remove_intensity_columns, combine_sample_files

def main():
    file1 = './data/Mu EPIC Run 1 5-24-2021/Mu EPIC Run 1 FinalReport_ctrl-bkg Intensities_AVGbeta.xlsx'
    file2 = './data/Mu EPIC Run 2 RQ-022275 FINAL_02042022/Mu EPIC Run 2 FinalReport_ctrl_bkg Instensities_AVGbeta.xlsx'

    removed_intensity_output_file1 = './data/removed_intensity_values/Mu EPIC Run 1 FinalReport_ctrl-bkg AVGbeta.xlsx'
    removed_intensity_output_file2 = './data/removed_intensity_values/Mu EPIC Run 2 FinalReport_ctrl_bkg AVGbeta.xlsx'

    df1 = remove_intensity_columns(file1, removed_intensity_output_file1)
    df2 = remove_intensity_columns(file2, removed_intensity_output_file2)

    final_output_path = './data/combined_run1_run2.xlsx'

    final_df = combine_sample_files(df1, df2, final_output_path)

if __name__ == "__main__":
    main()