'''
Script to identify proteins that catalyze reactions (catalysts) and catalysts that have enzyme kinetics (enzymes) and are also identified as essential
'''

import cPickle
import os
import csv
import numpy as np
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

fileDir = os.path.dirname(os.path.abspath(__file__))
home = os.path.dirname(os.path.dirname(fileDir))

with open(os.path.join(home, 'cached', 'simData_Fit_1.cPickle'), 'r') as f:
	sim_data = cPickle.load(f)

raw_data = KnowledgeBaseEcoli()

with open(os.path.join(home, 'validation', 'ecoli', 'flat', 'essentialGenes_536.tsv'), 'r') as f:
	reader = csv.reader(f, delimiter='\t')
	essentialData = np.array(list(reader))

essentialGenes = np.array([x[2] + '[' + x[3] + ']' for x in essentialData])

catalysts = [x for x in essentialGenes if x in sim_data.process.metabolism.catalystsList]
enzymes = [x for x in essentialGenes if x in sim_data.process.metabolism.enzymeIdList]

# TODO: add check for complexes from monomers

monomerToGeneSymbol = {x['monomerId'] : x['symbol'] for x in raw_data.genes}
catalystSymbols = [monomerToGeneSymbol[x[:-3]] for x in catalysts]
enzymeSymbols = [monomerToGeneSymbol[x[:-3]] for x in enzymes]

print 'Essential catalysts:'
for p in catalystSymbols:
	print p

print '\nEssential enzymes:'
for p in enzymeSymbols:
	print p

lagTimes = {}
with open(os.path.join(fileDir, 'lagTimes.csv')) as f:
	reader = csv.reader(f)
	headers = reader.next()
	for line in reader:
		lagTimes[line[0]] = float(line[1])

eStr = ''
print '\nEssential enzymes with a lag time measured:'
for gene in (set(lagTimes.keys()) & set(enzymeSymbols)):
	enzyme = enzymes[enzymeSymbols.index(gene)]
	eStr += "'%s', " % (enzyme)
	print '%s (%s): %f' % (gene, enzyme, lagTimes[gene])
print eStr

import ipdb; ipdb.set_trace()
