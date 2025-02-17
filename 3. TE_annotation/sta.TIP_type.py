import sys
from collections import Counter

# Read the input file name from command line arguments
input_file = sys.argv[1]

# Count the most frequent string in each line, excluding "-"
with open(input_file, 'r') as file:
    for line in file:
        # Ignore the first 3 columns
        columns = line.strip().split('\t')[3:]
        # Count strings that are not "-"
        words = [word for word in columns if word != '-']
        if words:
            counter = Counter(words)
            most_common = counter.most_common(1)
            if most_common:
                most_common_word = most_common[0][0]
                print(most_common_word)
            else:
                print('-')  # If no non-dash words, print dash
        else:
            print('-')  # If no words, print dash