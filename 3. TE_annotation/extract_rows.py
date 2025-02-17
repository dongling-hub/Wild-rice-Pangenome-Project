; This script extracts rows from a table file based on a list of target values.
; The table file and list file are provided as input, and the matching rows are written to an output file.
; The script can optionally keep the header row from the table file in the output.

; Function: extract_rows
; Parameters:
; - table_file (str): Path to the input table file (tab-delimited).
; - list_file (str): Path to the input list file containing target values (one per line).
; - output_file (str): Path to the output file where matching rows will be written (tab-delimited).
; - keep_header (bool): Flag indicating whether to keep the header row in the output file (default is True).

; The function reads the list file and converts it to a set of target values.
; It then reads the table file, and if keep_header is True, it adds the header row to the list of matching rows.
; For each subsequent row in the table file, if the first column value is in the target set, the row is added to the list of matching rows.
; Finally, the matching rows are written to the output file.

; Usage:
; The script is executed from the command line with the following arguments:
; python extract_rows.py <table_file> <list_file> <output_file> <keep_header>
; - <table_file>: Path to the input table file.
; - <list_file>: Path to the input list file.
; - <output_file>: Path to the output file.
; - <keep_header>: Boolean flag ('true' or 'false') indicating whether to keep the header row in the output file.

import csv
import sys

def extract_rows(table_file, list_file, output_file, keep_header=True):
    # Read the list file and convert it to a set of target values
    with open(list_file, 'r') as f:
        target_set = set(line.strip() for line in f)

    matching_rows = []

    # Read and filter the table data
    with open(table_file, 'r') as f:
        table_reader = csv.reader(f, delimiter='\t')
        header = next(table_reader)  # Read the header row

        if keep_header:
            matching_rows.append(header)

        # Start filtering other rows
        for row in table_reader:
            if row[0] in target_set:
                matching_rows.append(row)

    # Output the matching rows
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(matching_rows)

    print(f"Matching rows written to {output_file}")


# Pass file names and the flag indicating whether to keep the header row through command line arguments
if len(sys.argv) < 5:
    print("Usage: python extract_rows.py <table_file> <list_file> <output_file> <keep_header>")
else:
    table_file = sys.argv[1]
    list_file = sys.argv[2]
    output_file = sys.argv[3]
    keep_header = sys.argv[4].lower() == 'true'

    extract_rows(table_file, list_file, output_file, keep_header)