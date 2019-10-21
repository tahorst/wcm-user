#! /usr/bin/env python

"""
Create table of implemented metabolites for latex for science paper
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
METABOLITES_FILE = os.path.join(DATA_DIR, 'metabolite_pools.tsv')
COMMON_NAMES_FILE = os.path.join(DATA_DIR, 'ecocyc-compounds.tsv')
OUTPUT_FILE = os.path.join(OUT_DIR, 'metabolites.txt')

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
n_parbox = 2
parbox_width = 0.45
max_lines = 45
continued_str = '\n\\newpage Table S7 continued.\n'
headers = ['Metabolite', 'Source']
header_line = ' & '.join(headers) + '\\\\\n'

# Read data to format
with open(METABOLITES_FILE) as f:
	reader = csv.reader(f, delimiter='\t')
	metabolites = list(reader)[2:]

# Read common names
with open(COMMON_NAMES_FILE) as f:
	reader = csv.reader(f, delimiter='\t')
	common_names = {row[0].lower(): row[1] for row in list(reader)[2:]}

lines = []
sources = {}
n_sources = 1
for name, source in metabolites:
	if source not in sources:
		sources[source] = n_sources
		n_sources += 1

	lines += ['{} & {} \\\\\n'.format(name.replace('_', '\\_'), sources[source])]
	# lines += ['{} & {} & {} \\\\\n'.format(name, common_names.get(name[:-3].lower()), sources[source])]
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
