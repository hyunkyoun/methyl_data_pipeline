import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, normaltest, levene
from statsmodels.stats.multitest import multipletests
import seaborn as sns

# === USER INPUT ===
print("=== VOLCANO PLOT DEBUGGING SCRIPT ===")
print("Separate multiple file paths with commas")
group1_files = input("Enter path(s) for Group 1 CSV(s): ").strip().split(",")
group2_files = input("Enter path(s) for Group 2 CSV(s): ").strip().split(",")

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

print("\n=== LOADING DATA ===")
group1 = load_and_merge(group1_files)
group2 = load_and_merge(group2_files)

# === BASIC DATA INSPECTION ===
print(f"\n=== BASIC DATA INSPECTION ===")
print(f"Group 1 shape: {group1.shape}")
print(f"Group 2 shape: {group2.shape}")
print(f"Group 1 columns: {list(group1.columns)}")
print(f"Group 2 columns: {list(group2.columns)}")

# Align CpGs
common_cpgs = group1.index.intersection(group2.index)
group1 = group1.loc[common_cpgs]
group2 = group2.loc[common_cpgs]
print(f"Common CpGs: {len(common_cpgs)}")

# === DATA QUALITY CHECKS ===
print(f"\n=== DATA QUALITY CHECKS ===")

# Check for missing values
print(f"Group 1 missing values: {group1.isnull().sum().sum()}")
print(f"Group 2 missing values: {group2.isnull().sum().sum()}")

# Check beta value ranges
print(f"\nGroup 1 beta range: {group1.min().min():.6f} to {group1.max().max():.6f}")
print(f"Group 2 beta range: {group2.min().min():.6f} to {group2.max().max():.6f}")

# Check for values outside [0,1]
invalid_g1 = ((group1 < 0) | (group1 > 1)).sum().sum()
invalid_g2 = ((group2 < 0) | (group2 > 1)).sum().sum()
print(f"Group 1 invalid beta values (outside [0,1]): {invalid_g1}")
print(f"Group 2 invalid beta values (outside [0,1]): {invalid_g2}")

# Check for identical samples (potential duplicates)
print(f"\n=== CHECKING FOR IDENTICAL SAMPLES ===")
for i, col1 in enumerate(group1.columns):
    for j, col2 in enumerate(group1.columns):
        if i < j:
            correlation = group1[col1].corr(group1[col2])
            if correlation > 0.99:
                print(f"WARNING: {col1} and {col2} are highly correlated (r={correlation:.4f})")

for i, col1 in enumerate(group2.columns):
    for j, col2 in enumerate(group2.columns):
        if i < j:
            correlation = group2[col1].corr(group2[col2])
            if correlation > 0.99:
                print(f"WARNING: {col1} and {col2} are highly correlated (r={correlation:.4f})")

# Check cross-group correlations (samples in wrong groups?)
print(f"\n=== CHECKING CROSS-GROUP SAMPLE SIMILARITIES ===")
max_cross_corr = 0
for col1 in group1.columns:
    for col2 in group2.columns:
        correlation = group1[col1].corr(group2[col2])
        if correlation > max_cross_corr:
            max_cross_corr = correlation
            most_similar = (col1, col2)

print(f"Highest cross-group correlation: {max_cross_corr:.4f} between {most_similar[0]} and {most_similar[1]}")
if max_cross_corr > 0.95:
    print("WARNING: Very high cross-group correlation suggests samples might be mislabeled!")

# === STATISTICAL DISTRIBUTIONS ===
print(f"\n=== STATISTICAL DISTRIBUTION ANALYSIS ===")

# Sample a subset of CpGs for detailed analysis
np.random.seed(42)
sample_cpgs = np.random.choice(common_cpgs, size=min(1000, len(common_cpgs)), replace=False)

# Calculate means and variances
group1_means = group1.loc[sample_cpgs].mean(axis=1)
group2_means = group2.loc[sample_cpgs].mean(axis=1)
group1_vars = group1.loc[sample_cpgs].var(axis=1)
group2_vars = group2.loc[sample_cpgs].var(axis=1)

print(f"Group 1 mean methylation: {group1_means.mean():.4f} ± {group1_means.std():.4f}")
print(f"Group 2 mean methylation: {group2_means.mean():.4f} ± {group2_means.std():.4f}")
print(f"Group 1 variance range: {group1_vars.min():.6f} to {group1_vars.max():.6f}")
print(f"Group 2 variance range: {group2_vars.min():.6f} to {group2_vars.max():.6f}")

# Check for zero variance CpGs
zero_var_g1 = (group1_vars == 0).sum()
zero_var_g2 = (group2_vars == 0).sum()
print(f"Zero variance CpGs - Group 1: {zero_var_g1}, Group 2: {zero_var_g2}")

# === DETAILED P-VALUE ANALYSIS ===
print(f"\n=== P-VALUE DISTRIBUTION ANALYSIS ===")

# Clamp beta values
group1_clamped = group1.clip(lower=0.001, upper=0.999)
group2_clamped = group2.clip(lower=0.001, upper=0.999)

# Convert to M-values
mvals1 = np.log2(group1_clamped / (1 - group1_clamped))
mvals2 = np.log2(group2_clamped / (1 - group2_clamped))

# Calculate delta beta
delta_beta = group1_clamped.mean(axis=1) - group2_clamped.mean(axis=1)
delta_beta_percent = delta_beta * 100

# Perform t-tests on a sample
print("Performing t-tests on sample CpGs...")
p_values = []
t_stats = []
effect_sizes = []

for cpg in sample_cpgs:
    try:
        stat, p_val = ttest_ind(
            mvals1.loc[cpg].dropna(), 
            mvals2.loc[cpg].dropna(), 
            equal_var=False
        )
        p_values.append(p_val)
        t_stats.append(stat)
        effect_sizes.append(abs(delta_beta_percent.loc[cpg]))
    except Exception as e:
        print(f"Error with CpG {cpg}: {e}")
        p_values.append(1.0)
        t_stats.append(0.0)
        effect_sizes.append(0.0)

p_values = np.array(p_values)
t_stats = np.array(t_stats)
effect_sizes = np.array(effect_sizes)

# Analyze p-value distribution
print(f"\nP-value distribution summary:")
print(f"Min p-value: {p_values.min():.2e}")
print(f"Max p-value: {p_values.max():.6f}")
print(f"Mean p-value: {p_values.mean():.6f}")
print(f"Median p-value: {np.median(p_values):.6f}")

# Check for exact values
exact_005 = (p_values == 0.05).sum()
near_005 = ((p_values > 0.049) & (p_values < 0.051)).sum()
print(f"Exactly p=0.05: {exact_005}")
print(f"Near p=0.05 (0.049-0.051): {near_005}")

# Check for other common exact values
exact_001 = (p_values == 0.01).sum()
exact_010 = (p_values == 0.10).sum()
print(f"Exactly p=0.01: {exact_001}")
print(f"Exactly p=0.10: {exact_010}")

# === PROBLEMATIC P-VALUES INVESTIGATION ===
print(f"\n=== INVESTIGATING PROBLEMATIC P-VALUES ===")

# Find CpGs with suspicious p-values
boundary_mask = (p_values > 0.049) & (p_values < 0.051) & (effect_sizes > 10)
if boundary_mask.sum() > 0:
    print(f"Found {boundary_mask.sum()} CpGs with p≈0.05 and high effect size")
    
    problematic_cpgs = np.array(sample_cpgs)[boundary_mask]
    for cpg in problematic_cpgs[:5]:  # Show first 5
        g1_vals = mvals1.loc[cpg].dropna()
        g2_vals = mvals2.loc[cpg].dropna()
        
        print(f"\nCpG: {cpg}")
        print(f"  Group 1: mean={g1_vals.mean():.4f}, std={g1_vals.std():.4f}, n={len(g1_vals)}")
        print(f"  Group 2: mean={g2_vals.mean():.4f}, std={g2_vals.std():.4f}, n={len(g2_vals)}")
        print(f"  Delta beta %: {delta_beta_percent.loc[cpg]:.2f}%")
        
        # Check if groups have identical values
        if g1_vals.std() == 0 or g2_vals.std() == 0:
            print(f"  WARNING: Zero variance detected!")
        
        # Manual t-test
        stat, p_val = ttest_ind(g1_vals, g2_vals, equal_var=False)
        print(f"  T-statistic: {stat:.6f}, p-value: {p_val:.10f}")

# === VISUAL DIAGNOSTICS ===
print(f"\n=== GENERATING DIAGNOSTIC PLOTS ===")

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Plot 1: P-value histogram
axes[0,0].hist(p_values, bins=50, alpha=0.7, edgecolor='black')
axes[0,0].axvline(0.05, color='red', linestyle='--', label='p=0.05')
axes[0,0].set_xlabel('P-value')
axes[0,0].set_ylabel('Frequency')
axes[0,0].set_title('P-value Distribution')
axes[0,0].legend()

# Plot 2: P-value vs Effect Size
axes[0,1].scatter(effect_sizes, -np.log10(p_values), alpha=0.6, s=10)
axes[0,1].axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
axes[0,1].axvline(5, color='black', linestyle='--', label='5% effect')
axes[0,1].set_xlabel('Effect Size (|Δβ%|)')
axes[0,1].set_ylabel('-log10(p-value)')
axes[0,1].set_title('Mini Volcano Plot (Sample)')
axes[0,1].legend()

# Plot 3: Group means comparison
axes[0,2].scatter(group1_means, group2_means, alpha=0.6, s=10)
axes[0,2].plot([0, 1], [0, 1], 'r--', label='y=x')
axes[0,2].set_xlabel('Group 1 Mean Beta')
axes[0,2].set_ylabel('Group 2 Mean Beta')
axes[0,2].set_title('Group Means Comparison')
axes[0,2].legend()

# Plot 4: Variance comparison
axes[1,0].scatter(group1_vars, group2_vars, alpha=0.6, s=10)
axes[1,0].plot([group1_vars.min(), group1_vars.max()], 
               [group1_vars.min(), group1_vars.max()], 'r--', label='y=x')
axes[1,0].set_xlabel('Group 1 Variance')
axes[1,0].set_ylabel('Group 2 Variance')
axes[1,0].set_title('Variance Comparison')
axes[1,0].set_xscale('log')
axes[1,0].set_yscale('log')
axes[1,0].legend()

# Plot 5: Delta beta distribution
axes[1,1].hist(delta_beta_percent.loc[sample_cpgs], bins=50, alpha=0.7, edgecolor='black')
axes[1,1].axvline(0, color='red', linestyle='-', label='No difference')
axes[1,1].axvline(5, color='black', linestyle='--', label='±5%')
axes[1,1].axvline(-5, color='black', linestyle='--')
axes[1,1].set_xlabel('Delta Beta (%)')
axes[1,1].set_ylabel('Frequency')
axes[1,1].set_title('Effect Size Distribution')
axes[1,1].legend()

# Plot 6: P-value vs T-statistic
axes[1,2].scatter(t_stats, -np.log10(p_values), alpha=0.6, s=10)
axes[1,2].axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
axes[1,2].set_xlabel('T-statistic')
axes[1,2].set_ylabel('-log10(p-value)')
axes[1,2].set_title('T-stat vs P-value')
axes[1,2].legend()

plt.tight_layout()
plt.savefig('volcano_debug_plots.png', dpi=300, bbox_inches='tight')
plt.show()

# === RECOMMENDATIONS ===
print(f"\n=== DEBUGGING RECOMMENDATIONS ===")

if max_cross_corr > 0.95:
    print("🚨 HIGH PRIORITY: Check sample labeling - very similar samples across groups")

if exact_005 > 10:
    print("🚨 HIGH PRIORITY: Too many exact p=0.05 values - statistical method issue")

if zero_var_g1 > 0 or zero_var_g2 > 0:
    print("⚠️  MEDIUM PRIORITY: Zero variance CpGs detected - filter these out")

if invalid_g1 > 0 or invalid_g2 > 0:
    print("⚠️  MEDIUM PRIORITY: Beta values outside [0,1] detected - check data preprocessing")

print(f"\n=== SUMMARY ===")
print(f"Groups appear to be: {'DIFFERENT' if abs(group1_means.mean() - group2_means.mean()) > 0.05 else 'SIMILAR'}")
print(f"Cross-group similarity: {'HIGH' if max_cross_corr > 0.9 else 'NORMAL'}")
print(f"P-value distribution: {'SUSPICIOUS' if exact_005 > 10 else 'NORMAL'}")
print(f"Data quality: {'GOOD' if invalid_g1 + invalid_g2 == 0 else 'NEEDS ATTENTION'}")

print(f"\n=== DEBUG COMPLETE ===")
print("Check 'volcano_debug_plots.png' for visual diagnostics")