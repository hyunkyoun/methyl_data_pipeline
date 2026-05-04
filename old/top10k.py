import pandas as pd

# === Step 1: Load your BMIQ-normalized beta matrix ===
# Assumes rows = CpGs, columns = samples
beta_df = pd.read_excel("/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/processed_bmiq_data.xlsx", index_col=0)

# === Step 2: Calculate variance across each CpG ===
cpg_variances = beta_df.var(axis=1, skipna=True)

# === Step 3: Get the top 10,000 CpGs by variance ===
top_10k_cpgs = cpg_variances.sort_values(ascending=False).head(25000).index

# === Step 4: Subset the beta matrix ===
top10k_beta_df = beta_df.loc[top_10k_cpgs]

# === Step 5: Save to CSV ===
top10k_beta_df.to_csv("beta_25k.csv")

print("✅ Top 10,000 most variable CpGs saved to 'betas_25k.csv'")
