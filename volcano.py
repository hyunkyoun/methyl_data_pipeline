import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os
from matplotlib.lines import Line2D

# === USER INPUT ===
print("Separate multiple file paths with commas")
group1_files = input("Enter path(s) for Group 1 CSV(s): ").strip().split(",")
group2_files = input("Enter path(s) for Group 2 CSV(s): ").strip().split(",")
save_dir = input("Enter directory to save volcano plot (leave blank to skip saving): ").strip()
use_mannwhitney = input("Use Mann–Whitney U test instead of t-test? (yes/no): ").strip().lower() == "yes"

# === LOAD & COMBINE GROUPS ===
def load_and_merge(files):
    dfs = []
    for f in files:
        f = f.strip()
        if f:
            df = pd.read_csv(f, index_col=0)
            dfs.append(df)
    return pd.concat(dfs, axis=1)

group1 = load_and_merge(group1_files)
group2 = load_and_merge(group2_files)

# === ALIGN CpGs ===
common_cpgs = group1.index.intersection(group2.index)
group1 = group1.loc[common_cpgs]
group2 = group2.loc[common_cpgs]

# === DIAGNOSTICS ===
print(f"\nGroup 1 samples: {group1.shape[1]}")
print(f"Group 2 samples: {group2.shape[1]}")

# === FILTER LOW VARIANCE CpGs ===
combined = pd.concat([group1, group2], axis=1)
variances = combined.var(axis=1)
high_var_cpgs = variances[variances > 1e-5].index
group1 = group1.loc[high_var_cpgs]
group2 = group2.loc[high_var_cpgs]

# === CALCULATE Δβ ===
mean1 = group1.mean(axis=1)
mean2 = group2.mean(axis=1)
delta_beta = mean1 - mean2

# === CALCULATE P-VALUES ===
if use_mannwhitney:
    p_values = pd.Series(index=group1.index, dtype='float64')
    for cpg in group1.index:
        try:
            _, p = mannwhitneyu(group1.loc[cpg], group2.loc[cpg], alternative='two-sided')
        except:
            p = 1.0
        p_values[cpg] = p
else:
    ttest = ttest_ind(group1.T, group2.T, axis=0, equal_var=False, nan_policy='omit')
    p_values = pd.Series(ttest.pvalue, index=group1.index)

# === FDR CORRECTION ===
_, fdrs, _, _ = multipletests(p_values, method='fdr_bh')
fdr_series = pd.Series(fdrs, index=p_values.index)
neg_log_p = -np.log10(p_values)

# === BUILD RESULTS DF ===
results_df = pd.DataFrame({
    'Delta_Beta': delta_beta,
    'Delta_Beta_x100': delta_beta * 100,
    '-log10(p-value)': neg_log_p,
    'p-value': p_values,
    'FDR': fdr_series
})

# === STATS SUMMARY ===
threshold_db = 0.5
threshold_fdr = 0.05
sig_mask = (abs(results_df["Delta_Beta_x100"]) > threshold_db) & (results_df["FDR"] < threshold_fdr)

print("\nΔβ range:", delta_beta.min(), "to", delta_beta.max())
print("Significant CpGs (|Δβ×100| > 5 & FDR < 0.05):", sig_mask.sum())

# === ASSIGN COLORS ===
colors = np.where(
    sig_mask & (results_df["Delta_Beta_x100"] > threshold_db), 'red',
    np.where(sig_mask & (results_df["Delta_Beta_x100"] < -threshold_db), 'blue', 'gray')
)

# === PLOT ===
plt.figure(figsize=(8, 8))
plt.scatter(results_df["Delta_Beta_x100"], results_df["-log10(p-value)"], s=10, alpha=0.7, c=colors)

# Threshold lines
plt.axvline(threshold_db, color='red', linestyle='--', linewidth=1)
plt.axvline(-threshold_db, color='red', linestyle='--', linewidth=1)
plt.axhline(-np.log10(threshold_fdr), color='green', linestyle='--', linewidth=1)

# Labels
plt.xlabel("Fract. Methyl. Difference × 100 (Δβ × 100)")
plt.ylabel("-log10(p-value)")
plt.title("Volcano Plot of Differential Methylation (FDR-adjusted)")

# Optional legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Group 1 > Group 2', markerfacecolor='blue', markersize=6),
    Line2D([0], [0], marker='o', color='w', label='Group 2 > Group 1', markerfacecolor='red', markersize=6),
    Line2D([0], [0], marker='o', color='w', label='Not significant', markerfacecolor='gray', markersize=6)
]
plt.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()

# === SAVE OR SHOW ===
if save_dir:
    os.makedirs(save_dir, exist_ok=True)
    base1 = "_".join([os.path.splitext(os.path.basename(f.strip()))[0] for f in group1_files])
    base2 = "_".join([os.path.splitext(os.path.basename(f.strip()))[0] for f in group2_files])
    out_path = os.path.join(save_dir, f"volcano_{base1}_vs_{base2}.png")
    plt.savefig(out_path, dpi=300)
    print(f"\n✅ Volcano plot saved to: {out_path}")
else:
    plt.show()
