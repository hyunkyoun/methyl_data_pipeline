import pandas as pd
import os
import re

# === USER INPUT FILES ===
beta_file = "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/completed_bmiq_data.xlsx"     # Replace with actual filename
sample_sheet_file = "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/master_samplesheet.csv" # Replace with actual filename

# === Load data ===
beta_df = pd.read_excel(beta_file)
sample_sheet = pd.read_csv(sample_sheet_file)

# === Normalize sample IDs ===
sample_sheet["tb"] = sample_sheet["tb"].astype(str).str.strip().str.upper()
beta_df.columns = beta_df.columns.astype(str).str.strip().str.upper()

# === Filter samples ===
sample_sheet = sample_sheet.dropna(subset=["Experiment"])
sample_sheet = sample_sheet[sample_sheet["tb"].isin(beta_df.columns)]

# === Build group label from 4 columns ===
def make_group_label(row):
    parts = [
        str(row["strain"]).strip(),
        str(row["genotype_dp10"]).strip(),
        str(row["genotype_dp16"]).strip(),
        str(row["Runx1"]).strip()
    ]
    label = "_".join(parts)
    return re.sub(r'[^\w\-]', '_', label)

sample_sheet["group_label"] = sample_sheet.apply(make_group_label, axis=1)

# === Build filtered beta matrix ===
filtered_beta_df = beta_df.set_index(beta_df.columns[0])
filtered_beta_df = filtered_beta_df.loc[:, filtered_beta_df.columns.isin(sample_sheet["tb"])]

# === Create experiment → group → sample list structure ===
experiment_group_data = {}
for _, row in sample_sheet.iterrows():
    exp = row["Experiment"]
    group = row["group_label"]
    sample = row["tb"]

    if exp not in experiment_group_data:
        experiment_group_data[exp] = {}
    if group not in experiment_group_data[exp]:
        experiment_group_data[exp][group] = []
    experiment_group_data[exp][group].append(sample)

# === Save each group to plots/ folder ===
base_dir = "plots"
os.makedirs(base_dir, exist_ok=True)

for exp, groups in experiment_group_data.items():
    exp_dir = os.path.join(base_dir, re.sub(r'[^\w\-]', '_', exp))
    os.makedirs(exp_dir, exist_ok=True)

    for group, samples in groups.items():
        valid_samples = [s for s in samples if s in filtered_beta_df.columns]
        group_df = filtered_beta_df[valid_samples]

        if group_df.shape[1] == 0:
            continue

        out_path = os.path.join(exp_dir, f"{group}.csv")
        group_df.to_csv(out_path)

print("✅ All beta value group files saved under 'plots/'")
