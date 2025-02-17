#!/usr/bin/python
# -*- coding: utf-8 -*-
# Usage: python script.py lstFile seqFile outfile
import argparse
from collections import OrderedDict
import time

parser = argparse.ArgumentParser()
parser.add_argument("lstFile", help="input a list file", type=str)
parser.add_argument("seqFile", help="input a sequence file", type=str)
parser.add_argument("outfile", help="output a file", type=str)
args = parser.parse_args()

lst = args.lstFile
seq = args.seqFile
out = args.outfile

start_time = time.time()

# Store query IDs in a set instead of a list
query = set(line.strip() for line in open(lst))

# Use OrderedDict to store database sequences in order to preserve the order of sequences in the input file
database = OrderedDict()
header = None
with open(seq) as f:
    for line in f:
        if line.startswith('>'):
            header = line.strip('>\n')
            if header in query:
                database[header] = []
        else:
            if header in query:
                database[header].append(line.strip())

# Write output directly to file instead of buffering in memory
with open(out, 'w') as f:
    for header, sequence in database.items():
        f.write('>{}\n{}\n'.format(header, ''.join(sequence)))

end_time = time.time()
print('Elapsed time: {:.2f} seconds'.format(end_time - start_time))
