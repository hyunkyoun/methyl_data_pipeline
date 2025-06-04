import pandas as pd
import numpy as np
import scimap as sm

def combat_normalize_beta_values(file1_2, file3_4, output_file, id_column='TargetID'):
    
    # Read both files
    df1_2 = pd.read_excel(file1_2)
    df3_4 = pd.read_excel(file3_4)
    
    print(f"Run 1_2 shape: {df1_2.shape}")
    print(f"Run 3_4 shape: {df3_4.shape}")
    
    # Combine the files (keeping all TargetIDs)
    combined_raw = pd.merge(df1_2, df3_4, on=id_column, how='outer')
    print(f"Combined raw shape: {combined_raw.shape}")
    
    # Separate TargetID from beta values
    target_ids = combined_raw[id_column]
    beta_data = combined_raw.drop(columns=[id_column])
    
    # Get column names to identify which batch each sample belongs to
    cols_1_2 = [col for col in df1_2.columns if col != id_column]
    cols_3_4 = [col for col in df3_4.columns if col != id_column]
    
    # Create batch labels for ComBat
    batch_labels = []
    for col in beta_data.columns:
        if col in cols_1_2:
            batch_labels.append('Run_1_2')
        elif col in cols_3_4:
            batch_labels.append('Run_3_4')
        else:
            batch_labels.append('Unknown')
    
    print(f"Batch distribution: {pd.Series(batch_labels).value_counts().to_dict()}")
    
    # Prepare data for ComBat (transpose so samples are rows, genes/targets are columns)
    # ComBat expects: rows = samples, columns = features (TargetIDs)
    beta_transposed = beta_data.T
    beta_transposed.columns = target_ids
    
    # Remove rows/columns with too many NaN values
    # Remove samples (rows) with >50% missing data
    beta_transposed = beta_transposed.dropna(thresh=len(beta_transposed.columns)*0.5)
    
    # Remove targets (columns) with >50% missing data  
    beta_transposed = beta_transposed.dropna(axis=1, thresh=len(beta_transposed)*0.5)
    
    print(f"After removing sparse data: {beta_transposed.shape}")
    
    # Update batch labels to match remaining samples
    remaining_samples = beta_transposed.index
    batch_labels_filtered = [batch_labels[i] for i, col in enumerate(beta_data.columns) if col in remaining_samples]
    
    # Apply ComBat normalization
    print("Applying ComBat normalization...")
    try:
        # Create batch series
        batch_series = pd.Series(batch_labels_filtered, index=remaining_samples)
        
        # Apply ComBat
        normalized_data = sm.pp.combat(
            beta_transposed, 
            batch=batch_series,
            inplace=False
        )
        
        # Transpose back to original format (TargetIDs as rows, samples as columns)
        normalized_transposed = normalized_data.T
        
        # Add TargetID column back
        normalized_df = normalized_transposed.reset_index()
        normalized_df.rename(columns={'index': id_column}, inplace=True)
        
        # Save normalized data
        normalized_df.to_excel(output_file, index=False)
        
        print(f"ComBat normalization complete!")
        print(f"Normalized data saved to: {output_file}")
        print(f"Final shape: {normalized_df.shape}")
        
        return normalized_df
        
    except Exception as e:
        print(f"Error during ComBat normalization: {e}")
        print("This might be due to insufficient data or batch structure issues")
        return None

# Usage:
# normalized_data = combat_normalize_beta_values(
#     'Run_1_2_combined.xlsx', 
#     'Run_3_4_combined.xlsx', 
#     'All_Runs_ComBat_Normalized.xlsx'
# )