import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import os
from matplotlib.lines import Line2D

# === USER INPUT ===
print("Separate multiple file paths with commas")
group1_files = input("Enter path(s) for Group 1 CSV(s): ").strip().split(",")
group2_files = input("Enter path(s) for Group 2 CSV(s): ").strip().split(",")
save_dir = input("Enter directory to save volcano plot (leave blank to skip saving): ").strip()

# === LOAD & COMBINE GROUPS ===
def load_and_merge(files):
    dfs = []
    for f in files:
        f = f.strip()
        if f:
            df = pd.read_csv(f, index_col=0)
            dfs.append(df)
    if not dfs:
        raise ValueError("No valid files provided")
    return pd.concat(dfs, axis=1)

print("Loading data...")
group1 = load_and_merge(group1_files)
group2 = load_and_merge(group2_files)

# === ALIGN CpGs ===
common_cpgs = group1.index.intersection(group2.index)
group1 = group1.loc[common_cpgs]
group2 = group2.loc[common_cpgs]

print(f"Group 1: {group1.shape[1]} samples")
print(f"Group 2: {group2.shape[1]} samples") 
print(f"Common CpGs: {len(common_cpgs)}")

# === CONVERT TO M-VALUES FOR STATISTICS ===
# Clamp beta values to avoid infinite M-values
group1_clamped = group1.clip(lower=0.001, upper=0.999)
group2_clamped = group2.clip(lower=0.001, upper=0.999)

# Convert to M-values: M = log2(beta/(1-beta))
mvals1 = np.log2(group1_clamped / (1 - group1_clamped))
mvals2 = np.log2(group2_clamped / (1 - group2_clamped))

# === CALCULATE STATISTICS ON M-VALUES ===
# Calculate means on M-values
mean_mvals1 = mvals1.mean(axis=1)
mean_mvals2 = mvals2.mean(axis=1)

# Log2 fold change (difference in M-values)
log2_fold_change = mean_mvals1 - mean_mvals2

# Also calculate delta beta for reference
delta_beta = group1_clamped.mean(axis=1) - group2_clamped.mean(axis=1)
delta_beta_percent = delta_beta * 100

# === MODERATED T-TEST (limma-style) ===
print("Performing moderated t-test...")

# Simple t-test on M-values (proper for methylation analysis)
p_values = []
t_stats = []

for cpg in mvals1.index:
    try:
        # Welch's t-test on M-values
        stat, p_val = ttest_ind(
            mvals1.loc[cpg].dropna(), 
            mvals2.loc[cpg].dropna(), 
            equal_var=False
        )
        p_values.append(p_val)
        t_stats.append(stat)
    except:
        p_values.append(1.0)
        t_stats.append(0.0)

p_values = pd.Series(p_values, index=mvals1.index)
t_stats = pd.Series(t_stats, index=mvals1.index)

# === FDR CORRECTION ===
rejected, fdrs, _, _ = multipletests(p_values.values, method='fdr_bh', alpha=0.05)
fdr_series = pd.Series(fdrs, index=p_values.index)

# Calculate -log10(p-value)
neg_log10_p = -np.log10(p_values.replace(0, 1e-300))

# === BUILD RESULTS ===
results_df = pd.DataFrame({
    'log2_fold_change': log2_fold_change,
    'delta_beta': delta_beta,
    'delta_beta_percent': delta_beta_percent,
    'p_value': p_values,
    'FDR': fdr_series,
    'neg_log10_p': neg_log10_p,
    't_stat': t_stats
})

# === SIGNIFICANCE THRESHOLDS ===
log2fc_threshold = 1.0  # Log2 fold change threshold (equivalent to 2-fold change)
fdr_threshold = 0.05
p_threshold = 0.05

# Identify significant CpGs
sig_fdr = results_df['FDR'] < fdr_threshold
large_effect = np.abs(results_df['log2_fold_change']) >= log2fc_threshold
significant = sig_fdr & large_effect

print(f"\nSignificance Summary:")
print(f"Total CpGs: {len(results_df)}")
print(f"FDR < {fdr_threshold}: {sig_fdr.sum()}")
print(f"|Log2FC| >= {log2fc_threshold}: {large_effect.sum()}")
print(f"Significant: {significant.sum()}")

# === ASSIGN COLORS ===
def get_color(row):
    # Use raw p-value threshold to match the horizontal line
    if row['p_value'] < p_threshold and abs(row['log2_fold_change']) >= log2fc_threshold:
        return 'red' if row['log2_fold_change'] > 0 else 'blue'
    return 'lightgray'

colors = results_df.apply(get_color, axis=1)

# === CREATE VOLCANO PLOT ===
plt.figure(figsize=(10, 8))

plt.scatter(results_df['log2_fold_change'], results_df['neg_log10_p'], 
           c=colors, s=8, alpha=0.6, edgecolors='none')

# Threshold lines
plt.axvline(log2fc_threshold, color='black', linestyle='--', alpha=0.7)
plt.axvline(-log2fc_threshold, color='black', linestyle='--', alpha=0.7)
plt.axhline(-np.log10(p_threshold), color='black', linestyle='--', alpha=0.7)

# Labels
plt.xlabel('Log₂ Fold Change', fontsize=12)
plt.ylabel('-log₁₀(p-value)', fontsize=12)
plt.title('Volcano Plot: Differential Methylation Analysis', fontsize=14)

# Legend - now using log2 fold change thresholds
n_hyper = ((results_df['log2_fold_change'] > log2fc_threshold) & 
           (results_df['p_value'] < p_threshold)).sum()
n_hypo = ((results_df['log2_fold_change'] < -log2fc_threshold) & 
          (results_df['p_value'] < p_threshold)).sum()

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8,
           label=f'Upregulated (n={n_hyper})'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8,
           label=f'Downregulated (n={n_hypo})'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgray', markersize=8,
           label='Not significant')
]
plt.legend(handles=legend_elements, loc='upper right')

plt.grid(True, alpha=0.3)
plt.tight_layout()

# === SAVE ===
if save_dir:
    os.makedirs(save_dir, exist_ok=True)
    base1 = "_".join([os.path.splitext(os.path.basename(f.strip()))[0] for f in group1_files if f.strip()])
    base2 = "_".join([os.path.splitext(os.path.basename(f.strip()))[0] for f in group2_files if f.strip()])
    
    plot_path = os.path.join(save_dir, f"volcano_{base1}_vs_{base2}_log2fc.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {plot_path}")
else:
    plt.show()