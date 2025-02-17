import sys

# Dictionary to store positions and corresponding TE family names from files
position_dict = {}

# Read bed file paths from the list file
list_file = sys.argv[1]
output_file = sys.argv[2]

with open(list_file, 'r') as file:
    bed_files = [line.strip() for line in file]

# Initialize the dictionary, listing all positions from all files
for bed_file in bed_files:
    with open(bed_file, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            position = '\t'.join(line[:3])

            # Set the file list corresponding to the position as None
            if position not in position_dict:
                position_dict[position] = {bed_file: None}

# Populate the dictionary, marking the TE family names present in the files
for bed_file in bed_files:
    with open(bed_file, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            position = '\t'.join(line[:3])
            te_family = line[3]  # Can be modified to extract superfamily or TE type information

            if position in position_dict:
                position_dict[position][bed_file] = te_family

# Sort positions
sorted_positions = sorted(position_dict.keys(), key=lambda x: (x.split('\t')[0], int(x.split('\t')[1])))

# Write results to the specified output file
with open(output_file, 'w') as file:
    for position in sorted_positions:
        file_dict = position_dict[position]
        te_families = [file_dict.get(bed_file, '-') for bed_file in bed_files]
        te_families_str = '\t'.join(te_families)
        file.write(f'{position}\t{te_families_str}\n')
