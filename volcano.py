import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, mannwhitneyu
import os

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
p_values = pd.Series(index=group1.index, dtype='float64')

if use_mannwhitney:
    for cpg in group1.index:
        try:
            _, p = mannwhitneyu(group1.loc[cpg], group2.loc[cpg], alternative='two-sided')
        except:
            p = 1.0
        p_values[cpg] = p
else:
    ttest = ttest_ind(group1.T, group2.T, axis=0, equal_var=False, nan_policy='omit')
    p_values = pd.Series(ttest.pvalue, index=group1.index)

neg_log_p = -np.log10(p_values)

# === BUILD RESULTS DF ===
results_df = pd.DataFrame({
    'Delta_Beta': delta_beta,
    '-log10(p-value)': neg_log_p,
    'p-value': p_values
})

# === STATS SUMMARY ===
sig = (abs(results_df["Delta_Beta"]) > 0.125) & (results_df["p-value"] < 0.05)
print("\nΔβ range:", delta_beta.min(), "to", delta_beta.max())
print("Significant CpGs (Δβ > 0.2 & p < 0.05):", sig.sum())

# === PLOT ===
plt.figure(figsize=(10, 6))
plt.scatter(results_df["Delta_Beta"], results_df["-log10(p-value)"], s=10, alpha=0.7, c='gray')
plt.scatter(results_df[sig]["Delta_Beta"], results_df[sig]["-log10(p-value)"], s=10, c='red')

plt.axvline(0.2, color='blue', linestyle='--', linewidth=1)
plt.axvline(-0.2, color='blue', linestyle='--', linewidth=1)
plt.axhline(-np.log10(0.05), color='green', linestyle='--', linewidth=1)

plt.xlabel("Fract. Methyl. Difference (Δβ)")
plt.ylabel("-log10(p-value)")
plt.title("Volcano Plot of Differential Methylation")
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
