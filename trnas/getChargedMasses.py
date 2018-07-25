# print id, mw7.2, location for charged tRNAs to add to modifiedForms.tsv

import cPickle
import os
import csv
import numpy as np

from reconstruction.spreadsheets import JsonWriter
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

raw_data = KnowledgeBaseEcoli()

# suppress scientific notation output
np.set_printoptions(suppress = True)

fileLoc = os.path.dirname(__file__)
# raw_data = cPickle.load(open(os.path.join(fileLoc, "raw_data.cp"), "rb"))
sim_data = cPickle.load(open(os.path.join(fileLoc, "sim_data.cp"), "rb"))

out = open(os.path.join(fileLoc, 'chargedData.tsv'), 'w')
writer = JsonWriter(out, ["id", "mw7.2", "location"], dialect = "excel-tab")

trnas = sim_data.process.transcription.rnaData['id'][sim_data.process.transcription.rnaData['isTRna']]
charged = [x['modifiedForms'] for x in raw_data.rnas if x['id']+'[c]' in trnas]
filteredCharged = []
for c1 in charged:
	for c2 in c1:
		if 'FMET' in c2 or 'modified' in c2:
			continue
		filteredCharged += [c2 + '[c]']

molNames = sim_data.state.bulkMolecules.bulkData['id']
mws = sim_data.state.bulkMolecules.bulkData['mass']
for rxn in raw_data.modificationReactions:
	reactants = []
	products = []
	for mol in rxn['stoichiometry']:
		if mol['coeff'] == -1:
			reactants += ['%s[%s]' % (mol['molecule'], mol['location'])]
		else:
			products += ['%s[%s]' % (mol['molecule'], mol['location'])]

	for trna, ctrna in zip(trnas, filteredCharged):
		if trna in reactants and ctrna in products:
			mass = 0
			for reactant in reactants:
				if reactant in molNames:
					mass += mws[np.where(molNames == reactant)[0][0]].asNumber()
				else:
					print 'could not get mass for %s' % (reactant)
			for product in products:
				if product == ctrna:
					continue

				if product in molNames:
					mass -= mws[np.where(molNames == product)[0][0]].asNumber()
				else:
					print 'could not get mass for %s' % (product)

			writer.writerow({
				'id': ctrna[:-3],
				'mw7.2': mass,
				'location': [ctrna[-2]],
				})

import ipdb; ipdb.set_trace()