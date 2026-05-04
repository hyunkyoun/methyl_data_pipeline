import pandas as pd
from functools import reduce
import os

def remove_intensity_columns(input_file, output_file):
    # Read the input file based on its extension
    _, ext = os.path.splitext(input_file)
    if ext.lower() == '.csv':
        df = pd.read_csv(input_file)
    elif ext.lower() in ['.xls', '.xlsx']:
        df = pd.read_excel(input_file)
    else:
        raise ValueError("Unsupported file format: must be .csv, .xls, or .xlsx")

    # Get columns that contain 'AVG_Beta' or 'TargetID'
    columns_to_keep = [col for col in df.columns if 'AVG_Beta' in str(col) or 'TargetID' in str(col)]

    # Filter the dataframe
    filtered_df = df[columns_to_keep]

    # Save to new Excel file
    filtered_df.to_excel(output_file, index=False)

    # Print stats
    print(f"Original columns: {len(df.columns)}")
    print(f"Filtered columns: {len(filtered_df.columns)}")
    print(f"Got {len(df.columns) - len(filtered_df.columns)} columns containing 'AVG_Beta'")

    return filtered_df

def combine_sample_files(dfs, output_file, id_column='TargetID'):
    """
    Combine multiple sample DataFrames based on the given ID column.

    Args:
        dfs (list): List of DataFrames to combine.
        output_file (str): Path to save the combined output.
        id_column (str): The column name for merging. Defaults to 'TargetID'.
    """
    for i, df in enumerate(dfs):
        if id_column not in df.columns:
            raise ValueError(f"Error: '{id_column}' not found in DataFrame {i + 1}")
        print(f"Input DF {i + 1}: {df.shape}, columns: {list(df.columns)}")

    combined_df = reduce(lambda left, right: pd.merge(left, right, on=id_column, how='outer'), dfs)
    print(f"\nFinal combined shape: {combined_df.shape}")

    combined_df.to_excel(output_file, index=False)
    print(f"Final combined file saved as: {output_file}")
    print("Note: Empty cells (NaN) indicate that TargetID didn't exist in that file")

    return combined_df

# Usage example:
# df1 = pd.read_excel("file1.xlsx")
# df2 = pd.read_excel("file2.xlsx")
# df3 = pd.read_excel("file3.xlsx")
# combined_df = combine_sample_files([df1, df2, df3], "merged_output.xlsx")


