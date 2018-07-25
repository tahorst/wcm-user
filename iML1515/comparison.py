'''
Run from wcEcoli home directory
'''

import json
import cPickle
import os
import re
import csv
import urllib
import numpy as np

directory = os.path.join('user', 'iML1515')

def getWCMetabolismMonomers(sim_data):
	metMonomers = [x for x in sim_data.process.metabolism.catalystsList if x not in sim_data.process.complexation.complexNames and x not in sim_data.process.equilibrium.complexNameToRxnIdx]
	metComplexes = [x for x in sim_data.process.metabolism.catalystsList if x in sim_data.process.complexation.complexNames or x in sim_data.process.equilibrium.complexNameToRxnIdx]
	validMonomers = sim_data.process.translation.monomerData["id"]

	assert len(metMonomers) + len(metComplexes) == len(sim_data.process.metabolism.catalystsList)

	# find monomers that are accounted for in reactions in the model
	for metComplex in metComplexes:
		if metComplex in sim_data.process.complexation.complexNames:
			metMonomers += sim_data.process.complexation.getMonomers(metComplex)["subunitIds"].tolist()
		elif metComplex in sim_data.process.equilibrium.complexNameToRxnIdx:
			for subunit in sim_data.process.equilibrium.getMonomers(metComplex)["subunitIds"].tolist():
				if subunit in sim_data.process.complexation.complexNames:
					metMonomers += sim_data.process.complexation.getMonomers(subunit)["subunitIds"].tolist()
				elif subunit in validMonomers:
					metMonomers += [subunit]
		else:
			raise Exception

	import ipdb; ipdb.set_trace()
	return metMonomers

# load files
with open(os.path.join(directory, 'iML1515.json'), 'r') as f:
	iml1515 = json.load(f)

with open(os.path.join(directory, 'genes.json'), 'r') as f:
	wcGenes = json.load(f)

with open(os.path.join('cached', 'simData_Fit_1.cPickle'), 'r') as f:
	sim_data = cPickle.load(f)

with open(os.path.join(directory, 'All_genes_of_E._coli_K-12_substr._MG1655.txt'), 'r') as f:
	reader = csv.reader(f, delimiter='\t')
	ecocycHeaders = reader.next()
	ecocycGenes = np.array(list(reader))

with open(os.path.join(directory, 'universal_model.json'), 'r') as f:
	universal = json.load(f)

# create dict of possible gene names from ecocyc data, includes gene name (abcZ), BiGG (b####) and ECK (ECK####)
header = ecocycHeaders.index('Names')
allGenes = {}
for names in ecocycGenes[:, header]:
	genes = [x for x in re.split(' // ', re.sub('[",|]', '', names))]
	allGenes.update({gene : genes for gene in genes})

	# handle special cases where data from ecocyc doesn't match it's site for alternative names
	if 'ftm' in genes:
		allGenes.update({'fmt' : genes})
	if 'fms' in genes:
		allGenes.update({'def' : genes})
	if 'yihE' in genes:
		allGenes.update({'rdoA' : genes})
	if 'yadF' in genes:
		allGenes.update({'can' : genes})

# get gene ids (BiGG id - b####)
mlGeneNames = [gene['id'] for gene in iml1515['genes'] if gene['id'] != 's0001']

# get whole cell metabolism monomers
wcMetMonomers = np.unique([x[:-3] for x in getWCMetabolismMonomers(sim_data)])

# create dict with whole cell metabolism monomers to gene name
monomerToGene = {gene['monomerId'] : gene['symbol'] for gene in wcGenes}
wcGeneNames = [monomerToGene[monomer] for monomer in wcMetMonomers]

genesInBoth = []
genesInWC = []
genesInML = []
notInEC = []

# check each gene in iML1515 to see if it is included in whole cell model
for gene in mlGeneNames:
	if gene in wcGeneNames:
		genesInBoth.append(gene)
	else:
		both = False
		if gene in allGenes:
			for alt in allGenes[gene]:
				if alt in wcGeneNames:
					genesInBoth.append(gene)
					both = True
					break
		else:
			notInEC.append(gene)
		if not both:
			genesInML.append(gene)

# check each gene in whole cell model to see if it's not in iML1515
for gene in wcGeneNames:
	if gene not in mlGeneNames:
		both = False
		if gene in allGenes:
			for alt in allGenes[gene]:
				if alt in mlGeneNames:
					both = True
					break
		else:
			notInEC.append(gene)
		if not both:
			genesInWC.append(gene)

print 'Genes in both: %i' % len(genesInBoth)
print 'Genes only in ml: %i' % len(genesInML)
print 'Genes only in wc: %i' % len(genesInWC)

# check which genes that are present only in iML1515 reaction network are expressed in the wcm
mlGenesToWCM = {}
notIncluded = []
for mlGene in genesInML:
	included = False
	for gene in allGenes[mlGene]:
		if gene in monomerToGene.values():
			mlGenesToWCM[mlGene] = gene
			included = True
	if not included:
		notIncluded.append(mlGene)


# find reactions in iML1515 that have an EC number annotated in the universal reactions file from UCSD
mlRxns = {rxn['id'] : rxn for rxn in iml1515['reactions']}
universalRxns = {rxn['id'] : rxn for rxn in universal['reactions']}
ec = 0
ecMnx = 0
genesInMLWithECRxn = []
genesInMLWithECRxnMnx = []
geneFunction = {}
for rxn in mlRxns:
	oldRxn = rxn
	if rxn not in universalRxns:
		rxn = re.sub('_copy[0-9]', '', rxn)

	if 'EC Number' in universalRxns[rxn].get('annotation', {}):
		ec += 1
		for gene in re.findall('b[0-9]*', mlRxns[oldRxn]['gene_reaction_rule']):
			if gene in genesInML and gene not in genesInMLWithECRxn:
				genesInMLWithECRxn.append(gene)

	# get EC number from MNX site - takes a long time and gives same as EC annotation
	# if 'MetaNetX (MNX) Equation' in universalRxns[rxn].get('annotation', {}):
	# 	url = urllib.urlopen(universalRxns[rxn]['annotation']['MetaNetX (MNX) Equation'])
	# 	ecNumber = re.findall('EC number</td><td>(.*?)<', url.read(), re.DOTALL)
	# 	if 'NA' not in ecNumber:
	# 		ecMnx += 1
	# 		for gene in re.findall('b[0-9]*', mlRxns[oldRxn]['gene_reaction_rule']):
	# 			if gene in genesInML and gene not in genesInMLWithECRxnMnx:
	# 				genesInMLWithECRxnMnx.append(gene)

	# find function of reactions catalyzed by each gene not included in wcm
	for gene in re.findall('b[0-9]*', mlRxns[oldRxn]['gene_reaction_rule']):
		if gene in genesInML:
			if gene not in geneFunction:
				geneFunction[gene] = []
			geneFunction[gene].append([mlRxns[oldRxn]['subsystem']])

print 'Genes only in ml with an EC reaction: %i' % len(genesInMLWithECRxn)
print 'EC reactions: %i/%i' %(ec, len(mlRxns))

# count and print the functional categories for each gene
functionCounts = {}
for gene in geneFunction:
	for function in np.unique(geneFunction[gene]):
		if function not in functionCounts:
			functionCounts[function] = 0
		functionCounts[function] += 1

print 'Counts for each function (%i total):' % (sum(functionCounts.values()))
for function in functionCounts:
	print '\t%i for %s' % (functionCounts[function], function)


import ipdb; ipdb.set_trace()