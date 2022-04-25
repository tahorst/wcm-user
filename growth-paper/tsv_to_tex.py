#! /usr/bin/env python

# Script to convert tsv files from reflect/model_inspection.py into latex tables

import sys


with open(sys.argv[1]) as f:
    data = f.read()

data = data.replace('\n', ' \\\\\n')
data = data.replace('\t', ' & ')
print(data)
