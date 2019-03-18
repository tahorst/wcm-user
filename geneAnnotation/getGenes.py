'''
Uses raw_data and sim_data files to identify the genes that are functionally implemented
in the model.

Requires:
	raw_data: cPickle object, specify path with RAW_DATA_FILE
	sim_data: cPickle object, specify path with SIM_DATA_FILE

Outputs:
	genes.tsv: tsv file with a list of genes that are functionally implemented
'''

import cPickle
import csv
import os


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
BASE_DIR = os.path.dirname(FILE_LOCATION)
RAW_DATA_FILE = os.path.join(BASE_DIR, 'scripts', 'raw_data.cp')
SIM_DATA_FILE = os.path.join(BASE_DIR, 'scripts', 'sim_data.cp')
OUTPUT_FILE = os.path.join(FILE_LOCATION, 'genes.tsv')


with open(RAW_DATA_FILE, 'rb') as f:
	raw_data = cPickle.load(f)
with open(SIM_DATA_FILE, 'rb') as f:
	sim_data = cPickle.load(f)

validMonomers = sim_data.process.translation.monomerData["id"]

monomers = []

##### Metabolism #####
metMonomers = [
	x for x in sim_data.process.metabolism.catalystsList
	if x not in sim_data.process.complexation.complexNames
	and x not in sim_data.process.equilibrium.complexNameToRxnIdx
	]
metComplexes = [
	x for x in sim_data.process.metabolism.catalystsList
	if x in sim_data.process.complexation.complexNames
	or x in sim_data.process.equilibrium.complexNameToRxnIdx
	]

assert len(metMonomers) + len(metComplexes) == len(sim_data.process.metabolism.catalystsList)

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

monomers += metMonomers

##### Translation #####
translationMonomers = []
translationMonomers += sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s30_fullComplex)["subunitIds"].tolist()
translationMonomers += sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s50_fullComplex)["subunitIds"].tolist()

monomers += translationMonomers

##### Transcription #####
transcriptionMonomers = []
transcriptionMonomers += sim_data.process.complexation.getMonomers("APORNAP-CPLX[c]")["subunitIds"].tolist()

monomers += transcriptionMonomers

##### RNA Decay #####
rnaDecayMonomers = []
rnaDecayMonomers += sim_data.process.rna_decay.endoRnaseIds
rnaDecayMonomers += sim_data.moleculeGroups.exoRnaseIds

monomers += rnaDecayMonomers

##### Transcriptional Regulation #####
tfMonomers = []
tfIds = [x + "[c]" for x in sim_data.process.transcription_regulation.tfToTfType]
tfComplexes = [
	x for x in tfIds
	if x in sim_data.process.complexation.complexNames
	or x in sim_data.process.equilibrium.complexNameToRxnIdx
	or x in sim_data.process.two_component_system.complexToMonomer
	]
tfMonomers += [
	x for x in tfIds
	if x not in sim_data.process.complexation.complexNames
	and x not in sim_data.process.equilibrium.complexNameToRxnIdx
	and x not in sim_data.process.two_component_system.complexToMonomer
	]

assert len(tfMonomers) + len(tfComplexes) == len(tfIds)

for tfComplex in tfComplexes:
	if tfComplex in sim_data.process.complexation.complexNames:
		tfMonomers += sim_data.process.complexation.getMonomers(tfComplex)["subunitIds"].tolist()
	elif tfComplex in sim_data.process.equilibrium.complexNameToRxnIdx:
		for subunit in sim_data.process.equilibrium.getMonomers(tfComplex)["subunitIds"].tolist():
			if subunit in sim_data.process.complexation.complexNames:
				tfMonomers += sim_data.process.complexation.getMonomers(subunit)["subunitIds"].tolist()
			elif subunit in validMonomers:
				tfMonomers += [subunit]
	elif tfComplex in sim_data.process.two_component_system.complexToMonomer:
		for subunit in sim_data.process.two_component_system.complexToMonomer[tfComplex]:
			if subunit in sim_data.process.complexation.complexNames:
				tfMonomers += sim_data.process.complexation.getMonomers(subunit)["subunitIds"].tolist()
			elif subunit in validMonomers:
				tfMonomers += [subunit]
	else:
		raise Exception

tfMonomers = [x for x in tfMonomers if x in validMonomers]

monomers += tfMonomers

monomers = [x for x in monomers if x in validMonomers]

# Get gene names for each monomer implemented
rnaIdToSymbol = {x['rnaId']: x['symbol'] for x in raw_data.genes}
monomerToRna = {x['id'][:-3]: x['rnaId'][:-3] for x in sim_data.process.translation.monomerData}
geneNames = [rnaIdToSymbol[monomerToRna[monomer[:-3]]] for monomer in monomers]

# Save data to output tsv file
functional_monomers = {
	'Metabolism': metMonomers,
	'Translation': translationMonomers,
	'Transcription': transcriptionMonomers,
	'RNA Decay': rnaDecayMonomers,
	'Transcription Regulation': tfMonomers,
	}
with open(OUTPUT_FILE, 'w') as f:
	writer = csv.writer(f, delimiter='\t')
	writer.writerow(['Gene', 'Monomer', 'Process'])

	monomers_added = set()
	for gene, monomer in zip(geneNames, monomers):
		process = ''
		for function, subset in functional_monomers.items():
			if monomer in subset:
				process = function
				break

		# Prevent duplicates
		if monomer not in monomers_added:
			writer.writerow([gene, monomer, process])
			monomers_added.add(monomer)

nGenes = len(set(monomers))
print "Number of genes: %d" % nGenes
