import sys

input_file = sys.argv[1]
max_distance = float(sys.argv[2])
output_file = sys.argv[3]

# Read input file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Merge adjacent lines
merged_lines = []
previous_position = None
previous_line = None

for line in lines:
    line = line.strip()
    fields = line.split('\t')
    position = fields[:3]
    te_families = fields[3:]

    if previous_position is not None:
        if max_distance >= 1:
            if (
                position[0] == previous_position[0]
                and abs(int(position[1]) - int(previous_position[1])) <= max_distance
            ):
                merged_families = []
                for i in range(len(te_families)):
                    if te_families[i] != '-' or previous_line.split('\t')[3 + i] != '-':
                        merged_families.append(te_families[i] if te_families[i] != '-' else previous_line.split('\t')[3 + i])
                    else:
                        merged_families.append('-')

                merged_line = '\t'.join(previous_position + merged_families)
                merged_lines[-1] = merged_line
                previous_line = merged_line
            else:
                merged_lines.append(line)
                previous_line = line
        elif 0 < max_distance < 1:
            distance_condition = (
                abs(int(position[1]) - int(previous_position[1])) + abs(int(position[2]) - int(previous_position[2]))
                <= max_distance * abs(int(position[2]) - int(position[1]))
            )
            if position[0] == previous_position[0] and distance_condition:
                merged_families = []
                for i in range(len(te_families)):
                    if te_families[i] != '-' or previous_line.split('\t')[3 + i] != '-':
                        merged_families.append(te_families[i] if te_families[i] != '-' else previous_line.split('\t')[3 + i])
                    else:
                        merged_families.append('-')

                merged_line = '\t'.join(previous_position + merged_families)
                merged_lines[-1] = merged_line
                previous_line = merged_line
            else:
                merged_lines.append(line)
                previous_line = line
    else:
        merged_lines.append(line)
        previous_line = line

    previous_position = position

# Write merged data to output file
with open(output_file, 'w') as file:
    for merged_line in merged_lines:
        file.write(f'{merged_line}\n')
