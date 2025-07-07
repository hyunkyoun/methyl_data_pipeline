import pandas as pd

def process_bmiq_data(bmiq_data_path, sample_sheet_path):
    print("Reading sample sheet")
    sample_sheet = pd.read_csv(sample_sheet_path)
    sample_sheet['beadchip_position'] = sample_sheet['BEADCHIP'].astype(str) + "_" + sample_sheet['POSITION'].astype(str)

    print("Reading BMIQ data")
    bmiq = pd.read_excel(bmiq_data_path)

    # Create mapping from beadchip_position to TB
    mapping = dict(zip(sample_sheet['beadchip_position'], sample_sheet['tb']))

    print("Mapping column headers in BMIQ data")

    # Function to clean and map R-style column names (e.g., X205730000017_R06C02)
    def rename_column(col):
        stripped_col = col.lstrip("X")  # Remove leading 'X' added by R
        return mapping.get(stripped_col, col)  # Replace if match found, else keep original

    # Apply renaming
    bmiq.columns = [rename_column(col) for col in bmiq.columns]

    print("Saving to excel file")
    bmiq.to_excel("processed_bmiq_data.xlsx", index=False)

    print("Column headers updated successfully.")

if __name__ == "__main__":
    bmiq_data_path = "./data/bmiq_processed_values.xlsx"
    sample_sheet_path = "./data/master_samplesheet.csv"

    process_bmiq_data(bmiq_data_path, sample_sheet_path)
    print("BMIQ data processing completed.")
