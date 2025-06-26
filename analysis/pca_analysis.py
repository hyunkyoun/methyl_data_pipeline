import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def perform_pca(data_file, fileType, output_file, n_components=2):
    if fileType == 'csv':
        df = pd.read_csv(data_file)
    elif fileType in ['xls', 'xlsx']:
        df = pd.read_excel(data_file)

    numeric_df = df.select_dtypes(include=["number"])

    # Standardize the features
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(numeric_df)

    # Run PCA (e.g., reduce to 2 components)
    pca = PCA(n_components=2)
    pca_components = pca.fit_transform(scaled_data)

    # Create a DataFrame for the PCA results
    pca_df = pd.DataFrame(data=pca_components, columns=['PC1', 'PC2'])
    pca_df.to_excel(output_file, index=False)

    # Plot the results
    plt.figure(figsize=(8,6))
    plt.scatter(pca_df['PC1'], pca_df['PC2'])
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('2D PCA Plot')
    plt.grid()
    plt.show()

    # Print explained variance
    print("Explained variance ratio:", pca.explained_variance_ratio_)


if __name__ == "__main__":
    data_file = './data/bmiq'  # Replace with your data file path
    output_file = './data/bmiq/results/pca_results.csv'  # Output file for PCA results
    perform_pca(data_file, output_file)