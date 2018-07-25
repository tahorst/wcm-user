import numpy as np
from wholecell.utils.modular_fba import FluxBalanceAnalysis

reactions = {
	"R1": {"A": -1, "B": -1, "C": 1},
	"R2a": {"C": -2, "D": 1},
}

externalMolecules = {
	"A": 20,
	"B": 10,
}

objective = {"C": 10, "D": 10}
metaboliteConcentrations = [0 for x in sorted(objective.keys())]

fba = FluxBalanceAnalysis(reactions, externalMolecules, objective, objectiveType = "shift_pools", solver = "glpk")

fba.setInternalMoleculeLevels(metaboliteConcentrations)
fba.setExternalMoleculeLevels([externalMolecules[x] for x in externalMolecules])
metaboliteConcentrations += fba.getOutputMoleculeLevelsChange()

print "Metabolite changes:\n%s" % (fba.getOutputMoleculeLevelsChange())
print "Metabolite conc:\n%s" % (metaboliteConcentrations)
print "Objective value: %f" % (fba.getObjectiveValue())


