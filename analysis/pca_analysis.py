import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def perform_pca(data_file, fileType, output_file, n_components=2):
    if fileType == 'csv':
        df = pd.read_csv(data_file)
    elif fileType in ['xls', 'xlsx']:
        df = pd.read_excel(data_file)
    
    # Set CpG_ID as index and transpose so samples are rows
    df = df.set_index('CpG_ID').T  # This is the key fix!
    
    # Now each row is a sample, each column is a CpG site
    print(f"Data shape after transpose: {df.shape}")
    print(f"Samples (rows): {df.index.tolist()}")
    
    # Standardize the features (CpG sites)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df)
    
    # Run PCA
    pca = PCA(n_components=n_components)
    pca_components = pca.fit_transform(scaled_data)
    
    # Create a DataFrame for the PCA results with sample names
    pca_df = pd.DataFrame(
        data=pca_components, 
        columns=[f'PC{i+1}' for i in range(n_components)],
        index=df.index  # Keep sample names as index
    )
    
    # Add sample names as a column for the output
    pca_df.reset_index(inplace=True)
    pca_df.rename(columns={'index': 'Sample'}, inplace=True)
    pca_df.to_excel(output_file, index=False)
    
    # Plot the results with sample labels
    plt.figure(figsize=(10, 8))
    plt.scatter(pca_df['PC1'], pca_df['PC2'], s=100)
    
    # Add sample labels to points
    for i, sample in enumerate(pca_df['Sample']):
        plt.annotate(sample, (pca_df['PC1'][i], pca_df['PC2'][i]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
    plt.title('PCA: Sample Clustering')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    print("Explained variance ratio:", pca.explained_variance_ratio_)
    print("Total variance explained:", sum(pca.explained_variance_ratio_))

if __name__ == "__main__":
    data_file = '../bmiq.xlsx'
    output_file = './results/pca_results.xlsx'
    perform_pca(data_file, "xlsx", output_file)