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

# get genes from Orth 2011 paper
def getOrthModelGenes():
	orthGenes = []
	reactionFile = csv.reader(open(os.path.join(fileLoc, 'reactions', 'orthReactions062317.tsv'), 'r'), delimiter='\t')
	reactionFile.next()

	# strip genes from the protein-reaction association column of the orth reactions
	for line in reactionFile:
		genes = line[5]
		genes = re.sub('[\(\)]', '', genes)
		genes = re.sub(' and ', ' ', genes)
		genes = re.sub(' or ', ' ', genes)
		for gene in re.findall('[A-z]+', genes):
			if len(gene) > 4:
				if gene == 'SPONTANEOUS':
					continue

				# some have ec tacked to the end
				if gene[-2:] == 'ec':
					gene = gene[:-2]
				# multiple genes in same family are joined together eg TrpABC
				if len(gene) >= 4:
					for c in gene[3:]:
						orthGenes += [gene[:3] + c]
				else:
					orthGenes += [gene]
			elif len(gene) > 1:
				orthGenes += [gene]

	orthGenes = np.unique(orthGenes)
	return [x.lower() for x in orthGenes]

# load sim data
fileLoc = os.path.dirname(os.path.dirname(__file__))
raw_data = cPickle.load(open(os.path.join(fileLoc, "scripts", "raw_data.cp"), "rb"))
sim_data = cPickle.load(open(os.path.join(fileLoc, "scripts", "sim_data.cp"), "rb"))

validMonomers = sim_data.process.translation.monomerData["id"]

metMonomers = [x for x in sim_data.process.metabolism.catalystsList if x not in sim_data.process.complexation.complexNames and x not in sim_data.process.equilibrium.complexNameToRxnIdx]
metComplexes = [x for x in sim_data.process.metabolism.catalystsList if x in sim_data.process.complexation.complexNames or x in sim_data.process.equilibrium.complexNameToRxnIdx]

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

monomers = [x[:-3] for x in metMonomers]

# load reactions downloaded from ecocyc
reader = csv.reader(open(os.path.join(fileLoc, "reactions", "ecocycRxns051917.tsv"), "rb"), delimiter='\t')
reader.next()
ecocyc = np.array([[x[1], x[3], x[5]] for x in reader])
enzymes = ecocyc[:,0].tolist()
erxns = ecocyc[:,1].tolist()
substrates = ecocyc[:,2].tolist()

# load reactions included in the model from raw_data
mrxns = np.unique([cleanRxn(x['reaction id']) for x in raw_data.reactions]).tolist()

# determine reactions from ecocyc that are not included in the model
notInModel = [x for x in erxns if x not in mrxns]
inModel = [x for x in erxns if x in mrxns]

# find enzymes that could be included in the model with the new reactions that are not in the model
newRxnEnz = []
for rxn in notInModel:
	idx = erxns.index(rxn)
	if len(enzymes[idx]):
		for e in re.split(" // ", enzymes[idx]):
			newRxnEnz.append(e)

# find enzymes that could be included in the model by adding the enzyme to an existing reaction in the model
newAddedEnz = []
for rxn in inModel:
	idx = erxns.index(rxn)
	if len(enzymes[idx]):
		for e in re.split(" // ", enzymes[idx]):
			newAddedEnz.append(e)

newRxnEnz = np.unique(newRxnEnz)
newAddedEnz = np.unique(newAddedEnz)

# find all enzymes that can be added from the new reactions, might include catalysts for reactions already existing in the model
ecocycEnz = []
ecocycSub = []
for i, enz in enumerate(enzymes):
	if len(enz):
		for e in re.split(" // ", enz):
			ecocycEnz.append(e)
			ecocycSub.append(substrates[i])

ecocycEnz = np.array(ecocycEnz)
ecocycSubList = np.array(ecocycSub)

# find monomers of enzymes that could be added to the model with new reactions
complexNames = list(sim_data.process.complexation.complexNames)
complexNamesStripped = [x[:-3] for x in complexNames]
equilComplex = sim_data.process.equilibrium.complexNameToRxnIdx.keys()
equilComplexStripped = [x[:-3] for x in equilComplex]
newRxnMonomers = []
for metComplex in newRxnEnz:
	if metComplex in complexNamesStripped:
		idx = complexNamesStripped.index(metComplex) # need to add compartment id for getMonomers so pull from complexNames
		newRxnMonomers += [x[:-3] for x in sim_data.process.complexation.getMonomers(complexNames[idx])["subunitIds"]]
	elif metComplex in equilComplex:
		idx = equilComplexStripped.index(metComplex)
		for subunit in sim_data.process.equilibrium.getMonomers(equilComplex[idx])["subunitIds"].tolist():
			if subunit in complexNames:
				newRxnMonomers += [x[:-3] for x in sim_data.process.complexation.getMonomers(subunit)["subunitIds"]]
			elif subunit in validMonomers:
				newRxnMonomers += [subunit[:-3]]
	else:
		newRxnMonomers += [metComplex]

newEnzMonomers = []
for metComplex in newAddedEnz:
	if metComplex in complexNamesStripped:
		idx = complexNamesStripped.index(metComplex) # need to add compartment id for getMonomers so pull from complexNames
		newEnzMonomers += [x[:-3] for x in sim_data.process.complexation.getMonomers(complexNames[idx])["subunitIds"]]
	elif metComplex in equilComplex:
		idx = equilComplexStripped.index(metComplex)
		for subunit in sim_data.process.equilibrium.getMonomers(equilComplex[idx])["subunitIds"].tolist():
			if subunit in complexNames:
				newEnzMonomers += [x[:-3] for x in sim_data.process.complexation.getMonomers(subunit)["subunitIds"]]
			elif subunit in validMonomers:
				newEnzMonomers += [subunit[:-3]]
	else:
		newEnzMonomers += [metComplex]

ecocycMonomers = []
ecocycSubDict = {}
for metComplex in np.unique(ecocycEnz):
	if metComplex in complexNamesStripped:
		idx = complexNamesStripped.index(metComplex) # need to add compartment id for getMonomers so pull from complexNames
		proteins = [x[:-3] for x in sim_data.process.complexation.getMonomers(complexNames[idx])["subunitIds"]]
		for p in proteins:
			ecocycMonomers += [p]
			if p not in ecocycSubDict:
				ecocycSubDict[p] = []
			ecocycSubDict[p] += ecocycSubList[ecocycEnz == metComplex]
	elif metComplex in equilComplex:
		idx = equilComplexStripped.index(metComplex)
		for subunit in sim_data.process.equilibrium.getMonomers(equilComplex[idx])["subunitIds"].tolist():
			if subunit in complexNames:
				proteins = [x[:-3] for x in sim_data.process.complexation.getMonomers(subunit)["subunitIds"]]
				for p in proteins:
					ecocycMonomers += [p]
					if p not in ecocycSubDict:
						ecocycSubDict[p] = []
					ecocycSubDict[p] += ecocycSubList[ecocycEnz == metComplex]
			elif subunit in validMonomers:
				p = [subunit[:-3]]
				ecocycMonomers += [p]
				if p not in ecocycSubDict:
					ecocycSubDict[p] = []
				ecocycSubDict[p] += ecocycSubList[ecocycEnz == metComplex]
	else:
		ecocycMonomers += [metComplex]
		if metComplex not in ecocycSubDict:
			ecocycSubDict[metComplex] = []
		ecocycSubDict[metComplex] += ecocycSubList[ecocycEnz == metComplex]

monomersToBeAdded = np.unique([x for x in ecocycMonomers if x not in monomers])

reader = csv.reader(open('user/geneAnnotation/annotations.csv', 'r'), delimiter = '\t')
writer = csv.writer(open('user/geneAnnotation/annotationsEdited.csv', 'w'), delimiter = '\t')

# write headers
writer.writerow([''] + reader.next())
writer.writerow(['Included'] + ['In Orth 2011'] + ['Substrates'] + reader.next())

genesMet = []
genesOrth = []
genesOther = []
predicted = []
ygene = []
orthModelGenes = getOrthModelGenes()

# only include genes associated with Met or Met Orth or a metabolism monomer
for line in reader:
	included = "No"
	inOrth = "No"
	gene = line[0]
	protein = [x['monomerId'] for x in raw_data.genes if x['id'] == gene][0]
	if line[-2] == 'Met ':
		genesMet.append(gene)
	elif line[-2] == "Met Orth ":
		genesOrth.append(gene)
	elif protein in monomers or protein in ecocycMonomers:
		genesOther.append(gene)
	else:
		continue

	if re.findall("predicted", line[2]):
		predicted.append(gene)
	if line[1][0] == 'y':
		ygene.append(gene)

	if protein in monomers:
		included = "Yes"
	elif protein in newRxnMonomers and protein in newEnzMonomers:
		included = "Can add reaction or catalyst"
	elif protein in newRxnMonomers:
		included = "Can add reaction"
	elif protein in newEnzMonomers:
		included = "Can add catalyst"
	elif protein in ecocycMonomers:
		included = "??"

	if line[1].lower() in orthModelGenes:
		inOrth = "Yes"

	subs = ecocycSubDict.get(protein, "")

	writer.writerow([included] + [inOrth] + [subs] + line)


countMet = 0
countOrth = 0
countOther = 0
rxnMet = 0
rxnOrth = 0
rxnOther = 0
catMet = 0
catOrth = 0
catOther = 0
predictedMet = 0
predictedOrth = 0
predictedOther = 0
yMet = 0
yOrth = 0
yOther = 0

proteins = []

# print "Met:"
for gene in genesMet:
	protein = [x['monomerId'] for x in raw_data.genes if x['id'] == gene][0]
	proteins.append(protein)
	if protein not in monomers:
		countMet += 1
		if protein in newRxnMonomers:
			rxnMet += 1
		elif protein in ecocycMonomers:
			catMet += 1
		else:
			if gene in predicted:
				predictedMet += 1
			if gene in ygene:
				yMet += 1


# print "\n\nOrth:"
for gene in genesOrth:
	protein = [x['monomerId'] for x in raw_data.genes if x['id'] == gene][0]
	proteins.append(protein)
	if protein not in monomers:
		countOrth += 1
		if protein in newRxnMonomers:
			rxnOrth += 1
		elif protein in ecocycMonomers:
			catOrth += 1
		else:
			if gene in predicted:
				predictedOrth += 1
			if gene in ygene:
				yOrth += 1

for gene in genesOther:
	protein = [x['monomerId'] for x in raw_data.genes if x['id'] == gene][0]
	proteins.append(protein)
	if protein not in monomers:
		if protein in newRxnMonomers:
			rxnOther += 1
		elif protein in ecocycMonomers:
			catOther += 1
		else:
			if gene in predicted:
				predictedOther += 1
			if gene in ygene:
				yOther += 1

print "Met:\t%i missing out of %i with %i that can be added with new reaction and %i with existing reaction \n\tOf remaining, %i are predicted function and %i are y genes" % (countMet, len(genesMet), rxnMet, catMet, predictedMet, yMet)
print "Orth:\t%i missing out of %i with %i that can be added with new reaction and %i with existing reaction \n\tOf remaining, %i are predicted function and %i are y genes" % (countOrth, len(genesOrth), rxnOrth, catOrth, predictedOrth, yOrth)
print "Other:\t%i genes included in reactions but not annotated Met or Orth with %i that can be added with new reaction and %i with existing reaction" % (len(genesOther) - rxnOther, rxnOther, catOther)
print "Unannotated:\t%i genes in new reactions are unannotated in file" % (monomersToBeAdded.shape[0] - rxnMet - rxnOrth - rxnOther - catMet - catOrth - catOther)

unannotated = [x for x in monomersToBeAdded if x not in proteins]

import ipdb; ipdb.set_trace()