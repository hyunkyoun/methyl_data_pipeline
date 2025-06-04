import csv

class CSV_parsing:
    def __init__(self, file_path):
        self.file_path = file_path

    def parse(self, output_file_path="parsed_output.csv"):
        with open(self.file_path, 'r', encoding='utf-8') as infile:
            lines = infile.readlines()

        # Find the start of the data table
        start_index = None
        for i, line in enumerate(lines):
            if line.strip() == '[Sample Methylation Profile]':
                start_index = i + 1
                break

        if start_index is None:
            raise ValueError("No '[Sample Methylation Profile]' section found.")

        # Extract header and data
        header = lines[start_index].strip().split('\t')
        data_lines = lines[start_index + 1:]

        # Write to CSV
        with open(output_file_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
            for line in data_lines:
                row = line.strip().split('\t')
                writer.writerow(row)

        print(f"Conversion completed. Saved to {output_file_path}")

if __name__ == "__main__":
    input_file_path = "./data/Mu EPIC Run 1 5-24-2021/SampleMethFinalReport_nonorm.txt"

    parser = CSV_parsing(input_file_path)
    parser.parse("./data/parsed_output.csv")
