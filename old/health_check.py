import pandas as pd

def count_empty_and_total_cells(csv_file_path):
    """
    Counts the number of empty cells and total cells in a CSV file.

    Parameters:
    csv_file_path (str): The path to the CSV file.

    Returns:
    tuple: (number of empty cells, total number of cells)
    """
    df = pd.read_csv(csv_file_path)

    # Replace empty strings with NaN
    df.replace('', pd.NA, inplace=True)

    # Count total and empty cells
    total_cells = df.size
    empty_cells = df.isna().sum().sum()

    return empty_cells, total_cells


# file_path1 = './data/split_runs/run_1_2.csv'
# empty, total = count_empty_and_total_cells(file_path1)
# print(f"Empty cells: {empty}")
# print(f"Total cells: {total}")
# print(f"Percentage empty: {100 * empty / total:.2f}%")

# file_path2 = './data/split_runs/run_3_4.csv'
# empty, total = count_empty_and_total_cells(file_path2)
# print(f"Empty cells: {empty}")
# print(f"Total cells: {total}")
# print(f"Percentage empty: {100 * empty / total:.2f}%")

import os
import zipfile
import shutil

def unzip_all_in_directory(directory):

    if not os.path.isdir(directory):
        print(f"Directory {directory} does not exist.")
        return
    
    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.lower().endswith('.zip'):
            zip_path = os.path.join(directory, filename)
            extract_folder = os.path.join(directory, os.path.splitext(filename)[0])

            # Create a folder for extracted files (optional)
            os.makedirs(extract_folder, exist_ok=True)

            # Extract the ZIP file
            try:
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    zip_ref.extractall(extract_folder)
                    print(f"Extracted: {filename} → {extract_folder}")
            except zipfile.BadZipFile:
                print(f"Error: Bad ZIP file - {filename}")

# Example usage
# unzip_all_in_directory('./data/Mu EPIC Zip files incl image files')

def find_all_idat_files(directory):
    idat_files = []

    # Walk through the directory tree
    for root, _, files in os.walk(directory):
        for file in files:
            if file.lower().endswith('.idat'):
                idat_files.append(os.path.join(root, file))

    return idat_files

def count_idat_files(directory):
    count = 0

    for root, _, files in os.walk(directory):
        for file in files:
            if file.lower().endswith('.idat'):
                count += 1

    return count

def collect_idat_files(source_directory, destination_directory):
    # Create the destination directory if it doesn't exist
    os.makedirs(destination_directory, exist_ok=True)

    # Walk through the source directory to find .idat files
    for root, _, files in os.walk(source_directory):
        for file in files:
            if file.lower().endswith('.idat'):
                source_path = os.path.join(root, file)
                destination_path = os.path.join(destination_directory, file)

                # Handle name collisions by appending a number
                base, ext = os.path.splitext(file)
                counter = 1
                while os.path.exists(destination_path):
                    destination_path = os.path.join(destination_directory, f"{base}_{counter}{ext}")
                    counter += 1

                shutil.copy2(source_path, destination_path)
                print(f"Copied: {source_path} → {destination_path}")


# Example usage
idat_file_paths = find_all_idat_files('./data/Mu EPIC Zip files incl image files')
for path in idat_file_paths:
    print(path)

collect_idat_files('./data/Down_Syndrome_Mouse_Microarray_Samples', './idat')
new_idat_count = count_idat_files('./idat')
print(new_idat_count)


# file_path1 = './data/filtered_beta_matrix.csv'
# empty, total = count_empty_and_total_cells(file_path1)
# print(f"Empty cells: {empty}")
# print(f"Total cells: {total}")
# print(f"Percentage empty: {100 * empty / total:.2f}%")

# file_path2 = './data/split_runs/run_3_4.csv'
# empty, total = count_empty_and_total_cells(file_path2)
# print(f"Empty cells: {empty}")
# print(f"Total cells: {total}")
# print(f"Percentage empty: {100 * empty / total:.2f}%")