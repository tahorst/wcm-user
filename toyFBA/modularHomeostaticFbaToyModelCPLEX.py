import numpy as np
from wholecell.utils.modular_fba import FluxBalanceAnalysis

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


import cplex

prob = cplex.Cplex()

## set up rows
mets = {}
# normal metabolites
for rxn in reactions.values():
	mets.update(rxn)
mets = mets.keys()
# pseudo metabolites
for pseudo in objective:
	mets += [pseudo+'p']
# external metabolites
# for ext in externalMolecules:
# 	mets += [ext+'e']
nMets = len(mets)

## set up cols
rxns = reactions.keys()
# normal rxns
obj = list(np.zeros(len(rxns)))
# pseudo rxns
for pseudo in objective:
	rxns += [pseudo+'p', pseudo+'0', pseudo+'-1', pseudo+'1']
	obj += [0, 0, 1, 1]
# transport rxns
for ext in externalMolecules:
	# rxns += [ext+'e', ext+'t']
	# obj += [0, 0]
	rxns += [ext+'e']
	obj += [0]
nRxns = len(rxns)

## set S matrix
rows = []
cols = []
vals = []
# reaction stoich
for rxn in reactions:
	for met in reactions[rxn]:
		rows += [mets.index(met)]
		cols += [rxns.index(rxn)]
		vals += [reactions[rxn][met]]
for pseudo in objective:
	# pseudo rxn
	rows += [mets.index(pseudo), mets.index(pseudo+'p')]
	cols += [rxns.index(pseudo+'p'), rxns.index(pseudo+'p')]
	vals += [-objective[pseudo], 1]

	# pseudo target rxn
	rows += [mets.index(pseudo+'p')]
	cols += [rxns.index(pseudo+'0')]
	vals += [-1]

	# pseudo relax down
	rows += [mets.index(pseudo+'p')]
	cols += [rxns.index(pseudo+'-1')]
	vals += [-1]

	# pseudo relax up
	rows += [mets.index(pseudo+'p')]
	cols += [rxns.index(pseudo+'1')]
	vals += [1]
for ext in externalMolecules:
	# external source/sink
	rows += [mets.index(ext)]
	cols += [rxns.index(ext+'e')]
	vals += [np.sign(externalMolecules[ext])]

# set bounds
ub = list(10000*np.ones(nRxns))
lb = list(np.zeros(nRxns))
# transport constraints
for ext in externalMolecules:
	ub[rxns.index(ext+'e')] = np.abs(externalMolecules[ext])

# pseudo reactions
for pseudo in objective:
	idx = rxns.index(pseudo+'0')
	ub[idx] = 1
	lb[idx] = 1

prob.variables.add(obj=obj, ub=ub, lb=lb, names=rxns)
prob.linear_constraints.add(rhs=np.zeros(nMets), senses='E'*nMets, names=mets)
prob.linear_constraints.set_coefficients(zip(rows, cols, vals))
prob.objective.set_sense(prob.objective.sense.minimize)

import ipdb; ipdb.set_trace()

step = 1
while step <= 12:
	print "\nStep %i" % (step)
	fba = FluxBalanceAnalysis(reactions, externalMolecules, objective, objectiveType = "homeostatic", solver = "glpk")
	
	fba.setInternalMoleculeLevels(metaboliteConcentrations)
	fba.setExternalMoleculeLevels([externalMolecules[x] for x in externalMolecules])
	metaboliteConcentrations += fba.getOutputMoleculeLevelsChange()
	
	print "Metabolite changes:\n%s" % (fba.getOutputMoleculeLevelsChange())
	print "Metabolite conc:\n%s" % (metaboliteConcentrations)
	# print "Objective value: %f" % (fba.objectiveValue())

	# consume some metabolites
	metaboliteConcentrations *= 0.5
	# consume all of ATP
	# metaboliteConcentrations[1] = 0
	step += 1
	import ipdb; ipdb.set_trace()

import ipdb; ipdb.set_trace()
