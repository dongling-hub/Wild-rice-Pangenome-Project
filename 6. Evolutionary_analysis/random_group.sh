#!/bin/bash

# Check if enough arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 file1 file2 num_combinations"
    exit 1
fi

# Read the file names and the number of combinations
file1=$1
file2=$2
num_combinations=$3

# Check if the files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
    echo "One or both files do not exist."
    exit 1
fi

# Generate unique combinations and store them in an associative array
declare -A combinations
while [ "${#combinations[@]}" -lt "$num_combinations" ]; do
    # Randomly shuffle and select 2 lines from each file
    combo="$(shuf -n 2 "$file1" | tr '\n' '\t')$(shuf -n 2 "$file2" | tr '\n' '\t')"
    
    # Add the combination to the associative array to ensure uniqueness
    combinations["$combo"]=1
done

# Output the results with group labels
for combo in "${!combinations[@]}"; do
    # Print each combination with a group number label
    echo "$combo" | awk '{print "group"NR"\t"$0}'
done
