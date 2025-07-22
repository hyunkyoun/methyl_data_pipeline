import pandas as pd
import os
import re

# === USER INPUT FILES ===
beta_file = "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/betas_combat_adjusted_final.csv"     # Replace with actual filename
sample_sheet_file = "/Users/elliottseo/Documents/GitHub/methyl_data_pipeline/data/samplesheet.csv" # Replace with actual filename

# === Load data ===
beta_df = pd.read_csv(beta_file)
sample_sheet = pd.read_csv(sample_sheet_file)

# === Create sample mapping based on Sample ID ===
# The beta file columns are named like "TB{Sample ID}_{Experiment}"
# We need to match these with the "Sample ID" column in the sample sheet

# Extract sample information from beta column names
beta_columns_original = beta_df.columns.astype(str).str.strip()
# Skip the first column (CpG_ID column)
sample_columns = beta_columns_original[1:]

# Create mapping from Sample ID to full column names
sample_id_to_column = {}
for col in sample_columns:
    # Extract TB number and experiment from column names like "TB123_1"
    match = re.match(r'TB(\d+)_(\d+)', col)
    if match:
        sample_id = int(match.group(1))
        sample_id_to_column[sample_id] = col

# === Filter samples ===
sample_sheet = sample_sheet.dropna(subset=["Experiment", "Genotype"])
# Filter sample sheet to only include samples that have matching Sample ID in beta_df
sample_sheet = sample_sheet[sample_sheet["Sample ID"].isin(sample_id_to_column.keys())]

# === Build group label from Genotype column ===
def make_group_label(genotype_value):
    """Clean the genotype value to make it a valid filename"""
    genotype_str = str(genotype_value).strip()
    # Replace any non-alphanumeric characters with underscores
    clean_genotype = re.sub(r'[^\w\-]', '_', genotype_str)
    return clean_genotype

sample_sheet["group_label"] = sample_sheet["Genotype"].apply(make_group_label)

# === Build filtered beta matrix ===
filtered_beta_df = beta_df.set_index(beta_df.columns[0])

# Create a list of full column names that match our filtered samples
matching_columns = []
for _, row in sample_sheet.iterrows():
    sample_id = row["Sample ID"]
    if sample_id in sample_id_to_column:
        matching_columns.append(sample_id_to_column[sample_id])

# Filter beta dataframe to only include matching columns
filtered_beta_df = filtered_beta_df.loc[:, matching_columns]

# === Create experiment → group → sample list structure ===
experiment_group_data = {}
for _, row in sample_sheet.iterrows():
    exp = row["Experiment"]
    group = row["group_label"]
    sample_id = row["Sample ID"]
    
    # Get the full column name for this sample
    if sample_id in sample_id_to_column:
        full_column_name = sample_id_to_column[sample_id]
    else:
        continue  # Skip if no matching column found

    if exp not in experiment_group_data:
        experiment_group_data[exp] = {}
    if group not in experiment_group_data[exp]:
        experiment_group_data[exp][group] = []
    experiment_group_data[exp][group].append(full_column_name)

# === Save each group to plots/ folder ===
base_dir = "plots_genotype"
os.makedirs(base_dir, exist_ok=True)

for exp, groups in experiment_group_data.items():
    exp_dir = os.path.join(base_dir, f"Experiment_{exp}")
    os.makedirs(exp_dir, exist_ok=True)

    for group, samples in groups.items():
        valid_samples = [s for s in samples if s in filtered_beta_df.columns]
        group_df = filtered_beta_df[valid_samples]

        if group_df.shape[1] == 0:
            continue

        out_path = os.path.join(exp_dir, f"{group}.csv")
        group_df.to_csv(out_path)

print("✅ All beta value group files saved under 'plots_10k_genotype/'")
print(f"📊 Processed {len(sample_sheet)} samples across {len(experiment_group_data)} experiments")

# Print summary of matches
total_beta_columns = len(beta_df.columns) - 1  # Subtract 1 for index column
matched_samples = len(matching_columns)
print(f"🔍 Matched {matched_samples} samples out of {total_beta_columns} beta columns and {len(sample_sheet)} sample sheet entries")

# Print experiment and genotype breakdown
print("\n📋 Experiment and Genotype breakdown:")
for exp in sorted(experiment_group_data.keys()):
    print(f"  Experiment {exp}:")
    for group, samples in experiment_group_data[exp].items():
        print(f"    - {group}: {len(samples)} samples")

# Load original sample sheet to show unmatched samples
original_sample_sheet = pd.read_csv(sample_sheet_file)
original_sample_sheet = original_sample_sheet.dropna(subset=["Experiment"])

# Find unmatched samples from sample sheet
unmatched_from_sample_sheet = original_sample_sheet[~original_sample_sheet["Sample ID"].isin(sample_id_to_column.keys())]
if len(unmatched_from_sample_sheet) > 0:
    print(f"\n❌ Sample sheet entries NOT found in beta table ({len(unmatched_from_sample_sheet)}):")
    for _, row in unmatched_from_sample_sheet.iterrows():
        print(f"   - Sample ID {row['Sample ID']} (Experiment: {row['Experiment']})")
else:
    print("\n✅ All sample sheet entries found in beta table")

# Find unmatched beta columns
sample_sheet_ids = set(original_sample_sheet["Sample ID"])
unmatched_beta_columns = []
for sample_id, col_name in sample_id_to_column.items():
    if sample_id not in sample_sheet_ids:
        unmatched_beta_columns.append(col_name)

if len(unmatched_beta_columns) > 0:
    print(f"\n❌ Beta table columns NOT found in sample sheet ({len(unmatched_beta_columns)}):")
    for col in sorted(unmatched_beta_columns):
        print(f"   - {col}")
else:
    print("\n✅ All beta table columns found in sample sheet")

# Show unique genotypes found
unique_genotypes = sample_sheet["Genotype"].unique()
print(f"\n🧬 Found {len(unique_genotypes)} unique genotypes:")
for genotype in sorted(unique_genotypes):
    count = len(sample_sheet[sample_sheet["Genotype"] == genotype])
    print(f"   - {genotype}: {count} samples")