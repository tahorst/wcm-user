'''
Uses raw_data and sim_data files to identify the metabolites with concentrations
in the model.

Requires:
	raw_data: cPickle object, specify path with RAW_DATA_FILE
	sim_data: cPickle object, specify path with SIM_DATA_FILE

Outputs:
	metabolites.tsv: tsv file with a list of metabolites with concentrations
'''

import cPickle
import csv
import os


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
BASE_DIR = os.path.dirname(FILE_LOCATION)
RAW_DATA_FILE = os.path.join(BASE_DIR, 'scripts', 'raw_data.cp')
SIM_DATA_FILE = os.path.join(BASE_DIR, 'scripts', 'sim_data.cp')
OUTPUT_FILE = os.path.join(FILE_LOCATION, 'metabolites.tsv')


with open(RAW_DATA_FILE, 'rb') as f:
	raw_data = cPickle.load(f)
with open(SIM_DATA_FILE, 'rb') as f:
	sim_data = cPickle.load(f)

metabolites = sim_data.process.metabolism.concDict.keys()
bennett = [m['Metabolite'] for m in raw_data.metaboliteConcentrations]

with open(OUTPUT_FILE, 'w') as f:
	writer = csv.writer(f, delimiter='\t')
	writer.writerow(['Metabolite', 'Source'])

	for m in sorted(metabolites):
		if m[:-3] in bennett:
			source = 'Bennett et al. 2009'
		else:
			source = 'Biomass'
		writer.writerow([m, source])
