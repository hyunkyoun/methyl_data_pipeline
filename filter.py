import pandas as pd

def remove_intensity_columns(input_file, output_file):
    # Read the Excel file
    df = pd.read_excel(input_file)
    
    # Get columns that don't contain 'Intensity'
    columns_to_keep = [col for col in df.columns if 'Intensity' not in str(col)]
    
    # Filter the dataframe
    filtered_df = df[columns_to_keep]
    
    # Save to new Excel file
    filtered_df.to_excel(output_file, index=False)
    
    print(f"Original columns: {len(df.columns)}")
    print(f"Filtered columns: {len(filtered_df.columns)}")
    print(f"Removed {len(df.columns) - len(filtered_df.columns)} columns containing 'Intensity'")

def combine_sample_files(file1, file2, output_file, id_column='TargetID'):
    # Read both Excel files
    df1 = pd.read_excel(file1)
    df2 = pd.read_excel(file2)
    
    print(f"File 1 shape: {df1.shape} (TargetIDs: {len(df1)})")
    print(f"File 2 shape: {df2.shape} (TargetIDs: {len(df2)})")
    
    # Check if TargetID columns exist
    if id_column not in df1.columns:
        print(f"Error: '{id_column}' not found in file 1")
        return
    if id_column not in df2.columns:
        print(f"Error: '{id_column}' not found in file 2")
        return
    
    # Merge on TargetID, keeping ALL TargetIDs from both files
    combined_df = pd.merge(df1, df2, on=id_column, how='outer')
    
    # Count how many TargetIDs are in each category
    only_in_file1 = len(df1) - len(pd.merge(df1, df2, on=id_column, how='inner'))
    only_in_file2 = len(df2) - len(pd.merge(df1, df2, on=id_column, how='inner'))
    in_both = len(pd.merge(df1, df2, on=id_column, how='inner'))
    
    print(f"\nCombined shape: {combined_df.shape}")
    print(f"TargetIDs in both files: {in_both}")
    print(f"TargetIDs only in file 1: {only_in_file1}")
    print(f"TargetIDs only in file 2: {only_in_file2}")
    print(f"Total unique TargetIDs: {len(combined_df)}")
    
    # Save combined file
    combined_df.to_excel(output_file, index=False)
    
    print(f"\nCombined file saved as: {output_file}")
    print("Note: Empty cells (NaN) indicate that TargetID didn't exist in that file")
    
    return combined_df

# Usage:
# combined = combine_sample_files('file1.xlsx', 'file2.xlsx', 'combined_all_samples.xlsx')
    