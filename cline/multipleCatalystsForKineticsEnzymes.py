import cPickle
import os
import numpy as np

fileLoc = os.path.dirname(__file__)
raw_data = cPickle.load(open(os.path.join(fileLoc, "raw_data.cp"), "rb"))
sim_data = cPickle.load(open(os.path.join(fileLoc, "sim_data.cp"), "rb"))

for rxn in sim_data.process.metabolism.constrainedReactionList:
	noKinetics = False
	withKinetics = False
	for catalyst in sim_data.process.metabolism.reactionCatalysts[rxn]:
		if catalyst in sim_data.process.metabolism.enzymeIdList:
			withKinetics = True
		else:
			noKinetics = True

	if withKinetics and noKinetics:
		print rxn
