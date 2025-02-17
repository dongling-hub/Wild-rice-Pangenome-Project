import csv
from scipy.stats import fisher_exact
import sys
from statsmodels.stats.multitest import fdrcorrection

# Get command line arguments, the first argument is the input file path
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the CSV file, using tab as the delimiter
with open(input_file, 'r') as input_file:
    csv_reader = csv.DictReader(input_file, delimiter='\t')
    tables = [{'ID': row['ID'], 'Table': list(map(int, [row['Or_Ref'], row['Or_Alt'], row['Os_Ref'], row['Os_Alt']]))} for row in csv_reader]

# Perform Fisher's exact test and store the results in a list
results = []
for table in tables:
    odds_ratio, p_value = fisher_exact([table['Table'][:2], table['Table'][2:]])
    results.append({'ID': table['ID'], 'Odds ratio': odds_ratio, 'P-value': p_value})

# Perform FDR correction on the P-values
p_values = [result['P-value'] for result in results]
rejected, fdr_corrected_p_values = fdrcorrection(p_values)
for i, result in enumerate(results):
    result['FDR-corrected P-value'] = fdr_corrected_p_values[i]
    result['Rejected'] = rejected[i]

# Write the results to a CSV file, using tab as the delimiter
with open(output_file, 'w', newline='') as output_file:
    fieldnames = ['ID', 'Odds ratio', 'P-value', 'FDR-corrected P-value', 'Rejected']
    csv_writer = csv.DictWriter(output_file, delimiter='\t', fieldnames=fieldnames)
    csv_writer.writeheader()
    for result in results:
        csv_writer.writerow(result)