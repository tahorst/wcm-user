#! /usr/bin/env python

"""
Create table of implemented kinetic constraints for latex for science paper
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
TSV_FILE = os.path.join(DATA_DIR, 'kinetic_constraints.tsv')
OUTPUT_FILE = os.path.join(OUT_DIR, 'kinetics.txt')

# latex constant strings
LANDSCAPE = '\\begin{{landscape}}\n' \
	'{{\\scriptsize\n' \
	'\\topmargin 0.5in\n' \
	'{}' \
	'}}\n' \
	'\\end{{landscape}}\n'
TABLE = '\\begin{{longtable}}{{ p{{0.75in}} p{{1.5in}} p{{0.3in}} p{{0.3in}} p{{0.3in}} p{{1.5in}} p{{1.5in}} p{{1.5in}} p{{0.3in}} p{{0.3in}} }}\n' \
	'\\hline\n' \
	'{}' \
	'\\hline\n' \
	'{}' \
	'\\hline\n' \
	'\\end{{longtable}}\n'

# Read data to format
with open(TSV_FILE) as f:
	reader = csv.reader(f, delimiter='\t')
	reader.next()
	headers = reader.next()
	data = list(reader)
header_line = ' & '.join(headers) + '\\\\\n'
lines = [' & '.join(l) + '\\\\\n' for l in data]

# Assemble tables
content = TABLE.format(header_line, ''.join(lines))
output = LANDSCAPE.format(content)
output = output.replace('_', '\\_')

# Save output
with open(OUTPUT_FILE, 'w') as f:
	f.write(output)
