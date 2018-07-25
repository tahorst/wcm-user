import csv
import cPickle
import os
import numpy as np
import re

# function to take a reaction and return the reaction name for ecocyc
# strips away molecule identifying information
def cleanRxn(reaction):
	rxn = re.findall(".+RXN", reaction)
	if len(rxn) == 0:
		rxn = re.findall("RXN[^-]*-[0-9A-Z]+", reaction)
	elif rxn[0] == "TRANS-RXN":
		rxn = re.findall("TRANS-RXN[^-]*-[0-9A-Z]+", reaction)
	if len(rxn) == 0:
		rxn = [reaction]

	return rxn[0]

fileLoc = os.path.dirname(os.path.dirname(__file__))
raw_data = cPickle.load(open(os.path.join(fileLoc, "scripts", "raw_data.cp"), "rb"))
sim_data = cPickle.load(open(os.path.join(fileLoc, "scripts", "sim_data.cp"), "rb"))

# load reactions downloaded from ecocyc
reader = csv.reader(open(os.path.join(fileLoc, "metabolism", "ecocycRxns051917.tsv"), "rb"), delimiter='\t')
reader.next()
ecocyc = np.array([[x[1], x[3]] for x in reader])
enzymes = ecocyc[:,0].tolist()
erxns = ecocyc[:,1].tolist()

# load reactions included in the model from raw_data
mrxns = np.unique([cleanRxn(x['reaction id']) for x in raw_data.reactions]).tolist()

# determine reactions from ecocyc that are not included in the model
notInModel = [x for x in erxns if x not in mrxns]
print "new reactions not in model: %i" % len(notInModel)

# find enzymes that could be included in the model with the new reactions that are not in the model
total = 0
newEnz = []
for rxn in notInModel:
	idx = erxns.index(rxn)
	if len(enzymes[idx]):
		for e in re.split(" // ", enzymes[idx]):
			newEnz.append(e)
		total += 1

newEnz = np.unique(newEnz)
print "new reactions not in model with an enzyme (# unique enzymes): %i (%i, )" % (total, newEnz.shape[0], )

# find all enzymes that can be added from the new reactions, might include catalysts for reactions already existing in the model
ecocycEnz = []
for enz in enzymes:
	if len(enz):
		for e in re.split(" // ", enz):
			ecocycEnz.append(e)

ecocycEnz = np.unique(ecocycEnz)


# get monomers already included in metabolism
validMonomers = sim_data.process.translation.monomerData["id"]
metMonomers = [x[:-3] for x in sim_data.process.metabolism.catalystsList if x not in sim_data.process.complexation.complexNames and x not in sim_data.process.equilibrium.complexNameToRxnIdx]
metComplexes = [x for x in sim_data.process.metabolism.catalystsList if x in sim_data.process.complexation.complexNames or x in sim_data.process.equilibrium.complexNameToRxnIdx]

for metComplex in metComplexes:
	if metComplex in sim_data.process.complexation.complexNames:
		metMonomers += [x[:-3] for x in sim_data.process.complexation.getMonomers(metComplex)["subunitIds"]]
	elif metComplex in sim_data.process.equilibrium.complexNameToRxnIdx:
		for subunit in sim_data.process.equilibrium.getMonomers(metComplex)["subunitIds"].tolist():
			if subunit in sim_data.process.complexation.complexNames:
				metMonomers += [x[:-3] for x in sim_data.process.complexation.getMonomers(subunit)["subunitIds"]]
			elif subunit in validMonomers:
				metMonomers += [subunit[:-3]]

# get monomers of enzymes that could be added to the model with new reactions
complexNames = list(sim_data.process.complexation.complexNames)
complexNamesStripped = [x[:-3] for x in complexNames]
equilComplex = sim_data.process.equilibrium.complexNameToRxnIdx.keys()
equilComplexStripped = [x[:-3] for x in equilComplex]
newMonomers = []
for metComplex in ecocycEnz:
	if metComplex in complexNamesStripped:
		idx = complexNamesStripped.index(metComplex) # need to add compartment id for getMonomers so pull from complexNames
		newMonomers += [x[:-3] for x in sim_data.process.complexation.getMonomers(complexNames[idx])["subunitIds"]]
	elif metComplex in equilComplex:
		idx = equilComplexStripped.index(metComplex)
		for subunit in sim_data.process.equilibrium.getMonomers(equilComplex[idx])["subunitIds"].tolist():
			if subunit in complexNames:
				newMonomers += [x[:-3] for x in sim_data.process.complexation.getMonomers(subunit)["subunitIds"]]
			elif subunit in validMonomers:
				newMonomers += [subunit[:-3]]
	else:
		newMonomers += [metComplex]

metMonomers = np.unique(metMonomers)
newMonomers = np.unique(newMonomers)
monomersToBeAdded = [x for x in newMonomers if x not in metMonomers]
print "new monomers to be added with additional reactions: %i" % len(monomersToBeAdded)

notInNew = [x for x in mrxns if x not in erxns]
print "model reactions not in new set: %i" % len(notInNew)


import ipdb; ipdb.set_trace()
