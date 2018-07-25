import numpy as np
import cPickle
import os

fileLoc = os.path.dirname(os.path.dirname(__file__))
raw_data = cPickle.load(open(os.path.join(fileLoc, "scripts", "raw_data.cp"), "rb"))
sim_data = cPickle.load(open(os.path.join(fileLoc, "scripts", "sim_data.cp"), "rb"))

validMonomers = sim_data.process.translation.monomerData["id"]

monomers = []

##### Metabolism #####
metMonomers = [x for x in sim_data.process.metabolism.catalystsList if x not in sim_data.process.complexation.complexNames and x not in sim_data.process.equilibrium.complexNameToRxnIdx]
metComplexes = [x for x in sim_data.process.metabolism.catalystsList if x in sim_data.process.complexation.complexNames or x in sim_data.process.equilibrium.complexNameToRxnIdx]

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
translationMonomers += sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s30_fullComplex[0])["subunitIds"].tolist()
translationMonomers += sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s50_fullComplex[0])["subunitIds"].tolist()

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
tfComplexes = [x for x in tfIds if x in sim_data.process.complexation.complexNames or x in sim_data.process.equilibrium.complexNameToRxnIdx or x in sim_data.process.two_component_system.complexToMonomer]
tfMonomers += [x for x in tfIds if x not in sim_data.process.complexation.complexNames and x not in sim_data.process.equilibrium.complexNameToRxnIdx and x not in sim_data.process.two_component_system.complexToMonomer]

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
			print subunit, subunit in sim_data.process.complexation.complexNames, subunit in validMonomers
			if subunit in sim_data.process.complexation.complexNames:
				tfMonomers += sim_data.process.complexation.getMonomers(subunit)["subunitIds"].tolist()
			elif subunit in validMonomers:
				tfMonomers += [subunit]
	else:
		raise Exception

tfMonomers = [x for x in tfMonomers if x in validMonomers]

monomers += tfMonomers

monomers = [x for x in monomers if x in validMonomers]

# get gene names for each monomer implemented
rnaIdToSymbol = {x['rnaId']: x['symbol'] for x in raw_data.genes}
monomerToRna = {x['id'][:-3]: x['rnaId'][:-3] for x in sim_data.process.translation.monomerData}
geneNames = [rnaIdToSymbol[monomerToRna[monomer[:-3]]] for monomer in monomers]

# for monomer in sorted(set(monomers)):
# 	print monomer

nGenes = len(set(monomers))
print "Number of genes: %d" % nGenes
import ipdb; ipdb.set_trace()