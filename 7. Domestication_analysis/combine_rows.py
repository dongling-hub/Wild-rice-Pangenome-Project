import argparse

# Create a command-line argument parser
parser = argparse.ArgumentParser(description='Combine rows with the same first three columns')
parser.add_argument('input_files', nargs='+', help='input file(s)')
parser.add_argument('--output', '-o', default='output.txt', help='output file')
# python combine_rows.py SNP DEL INS -o output

args = parser.parse_args()

# Open the output file
with open(args.output, 'w') as out_file:
    # Create a dictionary to store rows with the same first column
    row_dict = {}

    # Loop through each input file
    for input_file in args.input_files:
        with open(input_file, 'r') as in_file:
            for line in in_file:
                # Extract the first X columns from each line as the key
                key = tuple(line.strip().split()[:3])

                if key in row_dict:
                    # If the key already exists, append the current line to the existing value
                    row_dict[key].append(line.strip())
                else:
                    # If the key does not exist, create a new key-value pair
                    row_dict[key] = [line.strip()]

    # Write each key-value pair to the output file
    for key, values in row_dict.items():
        out_file.write('\t'.join(key) + '\t' + '\t'.join(values) + '\n')
