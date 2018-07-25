import numpy as np
import cvxpy
from wholecell.utils.modular_fba import FluxBalanceAnalysis

MAX = np.inf

# reactions = {
# 	"R1": {"A": -1, "ATP": -1, "B": 1},
# 	"R2a": {"B": -1, "C": 1, "ATP": 2, "NADH": 2},
# 	"R2b": {"B": 1, "C": -1, "ATP": -2, "NADH": -2},
# 	"R3": {"B": -1, "F": 1},
# 	"R4": {"C": -1, "G": 1},
# 	"R5a": {"C": 0.8, "G": -1, "NADH": 2},
# 	"R5b": {"C": 0.8, "G": -1, "NADH": 2},
# 	"R6": {"C": -1, "D": 3, "ATP": 2},
# 	"R7": {"C": -1, "E": 3, "NADH": -4},
# 	"R8a": {"G": -1, "H": 1, "ATP": -1, "NADH": -2},
# 	"R8b": {"G": 1, "H": -1, "ATP": 1, "NADH": 2},
# 	"Rres": {"O2": -1, "ATP": 1, "NADH": -1},
# 	"Tc1": {"A": 1, "Carbon1": -1},
# 	"Tc2": {"A": 1, "Carbon2": -1},
# 	"Biomass": {"C": -1, "F": -1, "H": -1, "ATP": -10, "Biomass": 1},
# 	"C1": {"A": -1, "B": 1},
# 	# "C2": {"A": 1, "B": -1},
# 	"C3": {"A": 1, "B": -1}
# }

reactions = {
	"C1": {"A": -1, "B": 1},
	"C2": {"B": -1, "C": 1},
	"C3": {"C": -1, "D": 1}
}

externalMolecules = {
	"A": 1
}

# externalMolecules = {
# 	"Carbon1": 10.5,
# 	"Carbon2": 10.5,
# 	# "D": -12.,
# 	# "E": -12.,
# 	"F": 5.,
# 	"H": 5.,
# 	"O2": 15.
# }

objective = {"C": 0.5, "A": 0.5}
constrained = ['C2']

## different outcome if C3 is constrained reaction instead
# constrained = ['C3']

# objective = {"C": 1, "F": 1, "H": 1, "ATP": 10}
# constrained = ['C1', 'C2', 'C3']

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
fba.setExternalMoleculeLevels([externalMolecules[x] for x in externalMolecules])
fba.setInternalMoleculeLevels([0, 0])
fba.setKineticTarget(constrained, [1e18])
fba.disableKineticTargets()
print fba.getOutputMoleculeLevelsChange()
print fba.kineticTargetFluxes()
print fba.getReactionFluxes()
print fba.getObjectiveValue()
import ipdb; ipdb.set_trace()
