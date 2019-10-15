#! /usr/bin/env python

"""
Create table of data sources for latex for science paper
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
TSV_FILE = os.path.join(DATA_DIR, 'movie_data_sources.tsv')
OUTPUT_FILE = os.path.join(OUT_DIR, 'movie.txt')

# latex constant strings
TABLE = '\\begin{{table}}[H]\n' \
	'{}' \
	'\\end{{table}}\n'
PARBOX = '\\centering\n' \
	'\\scriptsize\n' \
	'\\begin{{tabular}}{{ c c c c c c }}\n' \
	'\\hline\n' \
	'{}' \
	'\\hline\n' \
	'{}' \
	'\\hline\n' \
	'\\end{{tabular}}\n' \

# Adjustable for dataset
first_lines = 50
max_lines = 55
continued_str = '\n\\newpage Table S9 continued.\n'

# Read data to format
with open(TSV_FILE) as f:
	reader = csv.reader(f, delimiter='\t')
	headers = reader.next()
	data = list(reader)
header_line = ' & '.join(headers[:-1]) + '\\\\\n'
lines = [' & '.join(l[:-1]) + '\\\\\n' for l in data]
n_lines = len(lines)

# Assemble tables
tables = []
current = 0
while current < n_lines:
	if current == 0:
		lines_per = first_lines
	else:
		lines_per = int(min(max_lines, np.ceil((n_lines - current))))
	content = PARBOX.format(header_line, ''.join(lines[current:current+lines_per]))
	current += lines_per
	tables += [TABLE.format(content)]
output = continued_str.join(tables)

# Save output
with open(OUTPUT_FILE, 'w') as f:
	f.write(output)
