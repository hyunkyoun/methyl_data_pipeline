import pandas as pd
import numpy as np
from combat.pycombat import pycombat

# Function to apply ComBat normalization to combined sample files
# input: file1_2 (str), file3_4 (str), output_file (str), id_column (str)
def combat_normalize(file1_2, file3_4, output_file, id_column='TargetID'):
    
    # Read both files
    df1_2 = pd.read_excel(file1_2)
    df3_4 = pd.read_excel(file3_4)
    
    print(f"Run 1_2 shape: {df1_2.shape}")
    print(f"Run 3_4 shape: {df3_4.shape}")
    
    # Combine the files (keeping all TargetIDs)
    combined_raw = pd.merge(df1_2, df3_4, on=id_column, how='outer')
    print(f"Combined raw shape: {combined_raw.shape}")
    
    # Prepare data
    target_ids = combined_raw[id_column]
    beta_data = combined_raw.drop(columns=[id_column])
    
    print(f"Beta data shape before NaN removal: {beta_data.shape}")
    print(f"NaN count: {beta_data.isna().sum().sum()}") # there were a few NaN values that I think we should relook into in the future!
    
    # Remove rows with ANY NaN values
    mask_no_nan = beta_data.notna().all(axis=1)
    beta_data_clean = beta_data[mask_no_nan]
    target_ids_clean = target_ids[mask_no_nan]
    
    print(f"Beta data shape after NaN removal: {beta_data_clean.shape}") # 1047 Cpgs removed due to NaN values (combat cannot handle NaNs)
    print(f"Removed {len(beta_data) - len(beta_data_clean)} TargetIDs due to NaN values") # However, less than 0.5% (287465 total Cpgs) of the data was removed
    
    # Final NaN check
    if beta_data_clean.isna().sum().sum() > 0:
        raise ValueError("Not all NaN values removed... please check")
    
    # Create batch vector
    # Get sample columns for each batch
    cols_1_2 = [col for col in df1_2.columns if col != id_column]


    # print(cols_1_2)
    # print("-=-" * 20)
    # print(beta_data_clean.columns)

    for col in beta_data_clean.columns:
        print(col)
    # Create batch vector in the same order as beta_data_clean.columns
    batch = []
    for col in beta_data_clean.columns:
        if col in cols_1_2:
            batch.append(1)
        else:
            batch.append(2)
    
    batch = np.array(batch)

    # print(batch)
    
    print(f"Batch distribution: {np.unique(batch, return_counts=True)}")
    print(f"Data shape: {beta_data_clean.shape}")
    print(f"Batch length: {len(batch)}")
    
    # Dimension check
    if len(batch) != beta_data_clean.shape[1]:
        raise ValueError(f"Batch length {len(batch)} does not match number of samples {beta_data_clean.shape[1]}")
    
    # print(f"Applying ComBat normalization...")
    
    # # Apply ComBat
    # corrected_data = pycombat(beta_data_clean, batch)
    
    # print(f"ComBat successful! Type: {type(corrected_data)}, Shape: {corrected_data.shape}")
    
    # # Add TargetID column back
    # corrected_data.insert(0, id_column, target_ids_clean.values)
    
    # # Save
    # corrected_data.to_excel(output_file, index=False)
    # print(f"Saved to: {output_file}")
    
    # return corrected_data