#! /usr/bin/env python

"""
Create table of implemented genes for latex for science paper
"""

from __future__ import division

import csv
import os

import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
OUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUT_DIR):
	os.mkdir(OUT_DIR)
TSV_FILE = os.path.join(DATA_DIR, 'functional_genes.tsv')
OUTPUT_FILE = os.path.join(OUT_DIR, 'genes.txt')

# latex constant strings
TABLE = '\\begin{{table}}[H]\n' \
	'{}' \
	'\\end{{table}}\n'
PARBOX = '\\parbox{{{}\\linewidth}}{{\n' \
	'\\centering\n' \
	'\\scriptsize\n' \
	'\\begin{{tabular}}{{ c c }}\n' \
	'\\hline\n' \
	'{}' \
	'\\hline\n' \
	'{}' \
	'\\hline\n' \
	'\\end{{tabular}}\n' \
	'}}\n'

# Adjustable for dataset
n_parbox = 3
parbox_width = 0.3
max_lines = 55
continued_str = '\n\\newpage Table S6 continued.\n'
headers = ['Gene', 'Function']
header_line = ' & '.join(headers) + '\\\\\n'

# Read data to format
with open(TSV_FILE) as f:
	reader = csv.reader(f, delimiter='\t')
	data = list(reader)[2:]
lines = ['{} & {} \\\\\n'.format(l[0].replace('_', '\\_'), l[-1]) for l in data]
n_lines = len(lines)

# Assemble tables
tables = []
current = 0
while current < n_lines:
	lines_per = int(min(max_lines, np.ceil((n_lines - current) / n_parbox)))
	parboxes = []
	for col in range(n_parbox):
		parboxes += [PARBOX.format(parbox_width, header_line, ''.join(lines[current:current+lines_per]))]
		current += lines_per
	tables += [TABLE.format('\\hfill\n'.join(parboxes))]
output = continued_str.join(tables)


# Save output
with open(OUTPUT_FILE, 'w') as f:
	f.write(output)
