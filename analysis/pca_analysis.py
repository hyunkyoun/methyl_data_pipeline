import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def read_file(file_path):
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.csv':
        return pd.read_csv(file_path)
    elif ext in ['.xls', '.xlsx']:
        return pd.read_excel(file_path)
    else:
        raise ValueError(f"Unsupported file format: {ext}")

def perform_pca(data_file, sample_info_file, output_file, plot_3d=False):
    # Load files
    df = read_file(data_file)
    sample_info = read_file(sample_info_file)

    # Clean and transpose beta matrix
    df = df.set_index('CpG_ID').T
    df.index.name = 'Sample'
    df.index = df.index.str.replace(r'^X', '', regex=True)
    df.index = df.index.str.replace(r'\.AVG_Beta$', '', regex=True)
    df.index = df.index.str.strip().astype(str)

    # Clean sample info
    sample_info.columns = sample_info.columns.str.strip().str.lower()
    if 'sample id' not in sample_info.columns or 'index' not in sample_info.columns:
        raise ValueError(f"Expected 'Sample ID' and 'index' columns. Got: {sample_info.columns.tolist()}")

    sample_info = sample_info.rename(columns={'sample id': 'Sample', 'index': 'Run'})
    sample_info['Sample'] = sample_info['Sample'].astype(str).str.strip()

    # Standardize beta values
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df)

    # Perform PCA
    n_components = 3 if plot_3d else 2
    pca = PCA(n_components=n_components)
    pca_components = pca.fit_transform(scaled_data)

    # Create PCA DataFrame
    pca_columns = [f'PC{i+1}' for i in range(n_components)]
    pca_df = pd.DataFrame(pca_components, columns=pca_columns, index=df.index)
    pca_df.reset_index(inplace=True)

    # Merge with sample info
    merged_df = pd.merge(pca_df, sample_info[['Sample', 'Run']], on='Sample', how='left')

    if merged_df['Run'].isnull().any():
        missing = merged_df[merged_df['Run'].isnull()]
        print("⚠️ Warning: Missing Run info for these samples:", missing['Sample'].tolist())

    # Simplify Run to RunGroup
    merged_df['RunGroup'] = merged_df['Run'].astype(str).str.extract(r'^(\d+)')

    # Save to Excel
    merged_df.to_excel(output_file, index=False)

    # Define highlight samples
    highlight_samples = {'11145', 'tb11145'}

    # Setup color map
    unique_runs = sorted(merged_df['RunGroup'].dropna().unique())
    cmap = plt.cm.get_cmap('tab10', len(unique_runs))
    color_map = {run: cmap(i) for i, run in enumerate(unique_runs)}

    # Plotting
    fig = plt.figure(figsize=(10, 8))
    if plot_3d:
        ax = fig.add_subplot(111, projection='3d')
        for run in unique_runs:
            run_df = merged_df[merged_df['RunGroup'] == run]
            ax.scatter(run_df['PC1'], run_df['PC2'], run_df['PC3'],
                       s=100, color=color_map[run], label=f'Run {run}')

            for j, sample in enumerate(run_df['Sample']):
                x, y, z = run_df['PC1'].values[j], run_df['PC2'].values[j], run_df['PC3'].values[j]
                if sample in highlight_samples:
                    ax.scatter(x, y, z, s=200, facecolors='none', edgecolors='black', linewidths=2)
                    ax.text(x, y, z, sample, fontsize=9, color='red', weight='bold')
    else:
        for run in unique_runs:
            run_df = merged_df[merged_df['RunGroup'] == run]
            plt.scatter(run_df['PC1'], run_df['PC2'], s=100, color=color_map[run], label=f'Run {run}')
            for j, sample in enumerate(run_df['Sample']):
                x, y = run_df['PC1'].values[j], run_df['PC2'].values[j]
                if sample in highlight_samples:
                    plt.scatter(x, y, s=200, facecolors='none', edgecolors='black', linewidths=2)
                    plt.annotate(sample, (x, y), xytext=(5, 5), textcoords='offset points',
                                 fontsize=9, color='red', weight='bold')
                else:
                    plt.annotate(sample, (x, y), xytext=(5, 5), textcoords='offset points', fontsize=8)

    # Labels
    if plot_3d:
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})')
        ax.set_zlabel(f'PC3 ({pca.explained_variance_ratio_[2]:.2%})')
        ax.set_title('3D PCA: Sample Clustering by Run')
        ax.legend(title='Run', loc='upper left', bbox_to_anchor=(1.1, 1))
    else:
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        plt.title('2D PCA: Sample Clustering by Run')
        plt.legend(title='Run', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    # Print variance info
    print("✅ Explained variance ratio:", pca.explained_variance_ratio_)
    print("✅ Total variance explained:", sum(pca.explained_variance_ratio_))

if __name__ == "__main__":
    data_file = '../bmiq.xlsx'                        # Adjust path
    sample_info_file = '../samplesTable.csv'  # Adjust path
    output_file = './results/pca_results.xlsx'
    perform_pca(data_file, sample_info_file, output_file, plot_3d=True)  # Change to False for 2D
