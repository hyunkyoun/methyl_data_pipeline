import pandas as pd
import os
import re

# === USER INPUT FILES ===
beta_file = "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/betas_final_top_10000.csv"     # Replace with actual filename
sample_sheet_file = "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/data/master_samplesheet.csv" # Replace with actual filename

# === Load data ===
beta_df = pd.read_csv(beta_file)
sample_sheet = pd.read_csv(sample_sheet_file)

# === Normalize sample IDs ===
sample_sheet["tb"] = sample_sheet["tb"].astype(str).str.strip().str.lower()

# Extract base column names (before underscore) from beta_df columns
beta_columns_original = beta_df.columns.astype(str).str.strip()
beta_columns_normalized = beta_columns_original.str.lower()
beta_base_columns_normalized = [col.split('_')[0] for col in beta_columns_normalized]

# Create mapping from normalized base column names to original full column names
column_mapping = dict(zip(beta_base_columns_normalized, beta_columns_original))

# === Filter samples ===
sample_sheet = sample_sheet.dropna(subset=["Experiment"])
# Filter sample sheet to only include samples that have matching base names in beta_df
sample_sheet = sample_sheet[sample_sheet["tb"].isin(beta_base_columns_normalized)]

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

# Create a list of full column names that match our filtered samples
matching_columns = []
for sample_id in sample_sheet["tb"]:
    if sample_id in column_mapping:
        matching_columns.append(column_mapping[sample_id])

# Filter beta dataframe to only include matching columns
filtered_beta_df = filtered_beta_df.loc[:, matching_columns]

# === Create experiment → group → sample list structure ===
experiment_group_data = {}
for _, row in sample_sheet.iterrows():
    exp = row["Experiment"]
    group = row["group_label"]
    sample_id = row["tb"]
    
    # Get the full column name for this sample
    if sample_id in column_mapping:
        full_column_name = column_mapping[sample_id]
    else:
        continue  # Skip if no matching column found

    if exp not in experiment_group_data:
        experiment_group_data[exp] = {}
    if group not in experiment_group_data[exp]:
        experiment_group_data[exp][group] = []
    experiment_group_data[exp][group].append(full_column_name)

# === Save each group to plots/ folder ===
base_dir = "plots_10k_new"
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

print("✅ All beta value group files saved under 'plots_10k_new/'")
print(f"📊 Processed {len(sample_sheet)} samples across {len(experiment_group_data)} experiments")

# Print summary of matches
total_beta_columns = len(beta_df.columns) - 1  # Subtract 1 for index column
matched_samples = len(matching_columns)
print(f"🔍 Matched {matched_samples} samples out of {total_beta_columns} beta columns and {len(sample_sheet)} sample sheet entries")

# Load original sample sheet to show unmatched samples
original_sample_sheet = pd.read_csv(sample_sheet_file)
original_sample_sheet["tb"] = original_sample_sheet["tb"].astype(str).str.strip().str.lower()
original_sample_sheet = original_sample_sheet.dropna(subset=["Experiment"])

# Find unmatched samples from sample sheet
unmatched_from_sample_sheet = original_sample_sheet[~original_sample_sheet["tb"].isin(beta_base_columns_normalized)]
if len(unmatched_from_sample_sheet) > 0:
    print(f"\n❌ Sample sheet entries NOT found in beta table ({len(unmatched_from_sample_sheet)}):")
    for _, row in unmatched_from_sample_sheet.iterrows():
        print(f"   - {row['tb']} (Experiment: {row['Experiment']})")
else:
    print("\n✅ All sample sheet entries found in beta table")

# Find unmatched beta columns
# Get all sample sheet tb values (normalized)
sample_sheet_tb_set = set(original_sample_sheet["tb"])

# Find beta columns whose base names are NOT in sample sheet
unmatched_beta_columns = []
for i, base_name in enumerate(beta_base_columns_normalized):
    if base_name not in sample_sheet_tb_set:
        unmatched_beta_columns.append(beta_columns_original[i])

if len(unmatched_beta_columns) > 0:
    print(f"\n❌ Beta table columns NOT found in sample sheet ({len(unmatched_beta_columns)}):")
    for col in sorted(unmatched_beta_columns):
        print(f"   - {col}")
else:
    print("\n✅ All beta table columns found in sample sheet")