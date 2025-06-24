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

class SampleReportParsing:
    def __init__(self, sample_tables_dict):
        self.sample_tables_dict = sample_tables_dict
        self.header_keywords = {"Index",
                                "Sample ID",
                                "Sample Group",
                                "Sentrix Barcode",
                                "Sample Section",
                                "Detected CpG (0.01)",
                                "Detected CpG (0.05)",
                                "Signal Average GRN",
                                "Signal Average RED",
                                "Signal P05 GRN",
                                "Signal P05 RED",
                                "Signal P25 GRN",
                                "Signal P25 RED",
                                "Signal P50 GRN",
                                "Signal P50 RED",
                                "Signal P75 GRN",
                                "Signal P75 RED",
                                "Signal P95 GRN",
                                "Signal P95 RED",
                                "Sample_Well",
                                "Sample_Plate",
                                "Pool_ID"}
        
    def is_header_line(self, line):
        columns = line.split('\t')
        return any(keyword in columns for keyword in self.header_keywords)
        
    def parse(self, output_file_path="sample_report_combined.csv"):
        headers = []
        all_data_rows = []

        # for path in self.file_paths:
        for path, run_number in self.sample_tables_dict.items():
            with open(path, 'r', encoding='utf-8') as f:
                found_headers = False
                for line in f:
                    line = line.strip()
                    if not found_headers and self.is_header_line(line):
                        current_headers = line.split('\t')
                        if not headers:
                            headers = current_headers
                        found_headers = True
                        continue
                    if found_headers:
                        if line:
                            row = line.split('\t')
                            row[0] = f"{run_number}_{row[0]}"
                            all_data_rows.append(row)

        if not headers:
            raise ValueError("No header line detected in input files.")

        with open(output_file_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(all_data_rows)

        print(f"Combined data written to {output_file_path}")
    
            

if __name__ == "__main__":
    # input_file_path = "./data/Mu EPIC Run 1 5-24-2021/SampleMethFinalReport_nonorm.txt"

    # parser = CSV_parsing(input_file_path)
    # parser.parse("./data/parsed_output.csv")
    input_file_path = {
        "./data/Mu EPIC Run 1 5-24-2021/SamplesTableFinalReport.txt": 1,
        "./data/Mu EPIC Run 2 RQ-022275 FINAL_02042022/TableControl.txt": 2,
    }

    parser = SampleReportParsing(input_file_path)
    parser.parse("./data/samplesTable.csv")

