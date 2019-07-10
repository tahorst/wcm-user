'''
Uses raw_data and sim_data files to identify the metabolites with concentrations
in the model and kinetic constraints used.

Requires:
	raw_data: cPickle object, specify path with RAW_DATA_FILE
	sim_data: cPickle object, specify path with SIM_DATA_FILE

Outputs:
	metabolites.tsv: tsv file with a list of metabolites with concentrations
	kinetics.tsv: tsv file with information about kinetic constraints used in the model
'''

import cPickle
import csv
import os

from wholecell.utils import units


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
BASE_DIR = os.path.dirname(FILE_LOCATION)
RAW_DATA_FILE = os.path.join(BASE_DIR, 'scripts', 'raw_data.cp')
SIM_DATA_FILE = os.path.join(BASE_DIR, 'scripts', 'sim_data.cp')
METABOLITES_FILE = os.path.join(FILE_LOCATION, 'metabolites.tsv')
KINETICS_FILE = os.path.join(FILE_LOCATION, 'kinetics.tsv')


# Load data
with open(RAW_DATA_FILE, 'rb') as f:
	raw_data = cPickle.load(f)
with open(SIM_DATA_FILE, 'rb') as f:
	sim_data = cPickle.load(f)

# Save metabolites
metabolites = sim_data.process.metabolism.concDict.keys()
bennett = [m['Metabolite'] for m in raw_data.metaboliteConcentrations]
with open(METABOLITES_FILE, 'w') as f:
	writer = csv.writer(f, delimiter='\t')
	writer.writerow(['Metabolite', 'Source'])

	for m in sorted(metabolites):
		if m[:-3] in bennett:
			source = 'Bennett et al. 2009'
		else:
			source = 'Biomass (EcoCyc GEM, Bremer and Dennis. 1996., and Neidhardt. 2006.)'
		writer.writerow([m, source])

# Save kinetic constraints
kinetic_constraints = sim_data.process.metabolism.constraintDict
with open(KINETICS_FILE, 'w') as f:
	writer = csv.writer(f, delimiter='\t')
	writer.writerow(['Pubmed ID', 'Reaction ID', 'kcat (1/s)',
		'KM (uM)', 'Adjusted kcat (1/s)', 'Enzyme',
		'Substrates', 'Substrates for KM', 'Temperature (C)'])

	for rxn in sorted(kinetic_constraints):
		c = kinetic_constraints[rxn]
		writer.writerow([c['Pubmed ID'], c['reactionID'], c['kcat'].asNumber(1 / units.s),
			c['kM'], c['kcatAdjusted'].asNumber(1 / units.s), c['enzymeIDs'],
			c['substrateIDs'], c['Concentration Substrates'], c['Temp']])
