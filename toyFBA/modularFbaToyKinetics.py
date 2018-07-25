import numpy as np
import csv
from wholecell.utils.modular_fba import FluxBalanceAnalysis

MAX = np.inf

reactions = {
	"R1": {"A": -1, "ATP": -1, "B": 1},
	"R2a": {"B": -1, "C": 1, "ATP": 2, "NADH": 2},
	"R2b": {"B": 1, "C": -1, "ATP": -2, "NADH": -2},
	"R3": {"B": -1, "F": 1},
	"R4": {"C": -1, "G": 1},
	"R5a": {"C": 0.8, "G": -1, "NADH": 2},
	"R5b": {"C": 0.8, "G": -1, "NADH": 2},
	"R6": {"C": -1, "D": 3, "ATP": 2},
	"R7": {"C": -1, "E": 3, "NADH": -4},
	"R8a": {"G": -1, "H": 1, "ATP": -1, "NADH": -2},
	"R8b": {"G": 1, "H": -1, "ATP": 1, "NADH": 2},
	"Rres": {"O2": -1, "ATP": 1, "NADH": -1},
	"Tc1": {"A": 1, "Carbon1": -1},
	"Tc2": {"A": 1, "Carbon2": -1}
}

externalMolecules = {
	"Carbon1": 10.5,
	"Carbon2": 10.5,
	"D": -12.,
	"E": -12.,
	"F": 5.,
	"H": 5.,
	"O2": 15.
}

objective = {"A": 10, "B": 20, "C": 20, "D": 20, "E": 20, "F": 20, "G": 20, "H": 20, "O2": 10, "NADH": 25, "ATP": 50}
metaboliteConcentrations = [objective[x] for x in sorted(objective.keys())]

constrained = {'R2b': 2, 'R4': 1, 'R7': 0}

# Set up FBA solver
fbaObjectOptions = {
	"reactionStoich" : reactions,
	"externalExchangedMolecules" : externalMolecules.keys(),
	"objective" : objective,
	"objectiveType" : "homeostatic_kinetics_mixed",
	"objectiveParameters" : {
			"kineticObjectiveWeight":0.5,#self.metabolismKineticObjectiveWeight,
			"reactionRateTargets":{reaction: 1 for reaction in constrained}, #This target is arbitrary, it gets reset each timestep during evolveState
			"oneSidedReactionTargets":[], # TODO: deal with this dynamically (each time step)
			},
	"solver" : "glpk",
}

fba = FluxBalanceAnalysis(**fbaObjectOptions)

step = 1
import time
st = time.time()
while step <= 12:
	print "\nStep %i" % (step)
	fba.setKineticTarget(constrained.keys(), constrained.values())
	fba.setReactionFluxBounds('R7', lowerBounds=0, upperBounds=0)

	fba.internalMoleculeLevelsIs(metaboliteConcentrations)
	fba.externalMoleculeLevelsIs([externalMolecules[x] for x in externalMolecules])
	metaboliteConcentrations += fba.outputMoleculeLevelsChange()
	print fba.outputMoleculeLevelsChange()
	print fba.kineticTargetFluxes()
	print fba.reactionFluxes()
	print fba.objectiveValue()

	# fba.setInternalMoleculeLevels(metaboliteConcentrations)
	# fba.setExternalMoleculeLevels([externalMolecules[x] for x in externalMolecules])
	# metaboliteConcentrations += fba.getOutputMoleculeLevelsChange()
	# print fba.getOutputMoleculeLevelsChange()
	# print fba.kineticTargetFluxes()
	# print fba.getReactionFluxes()
	# print fba.getObjectiveValue()

	# consume some metabolites
	metaboliteConcentrations *= 0.5
	# consume all of ATP
	# metaboliteConcentrations[1] = 0

	step += 1
	# import ipdb; ipdb.set_trace()
print time.time() - st

model = fba.getArrayBasedModel()
mets = model['Metabolites']
rxns = model['Reactions']
stoich = model['S_matrix']
mat = np.vstack((np.array(['']+rxns).reshape(1,-1), np.hstack((np.array(mets).reshape(-1,1), stoich))))

with open('master.mat', 'w') as f:
	writer = csv.writer(f, delimiter='\t')
	for row in mat:
		row[row == '0.0'] = ''
		writer.writerow(row)


import ipdb; ipdb.set_trace()
