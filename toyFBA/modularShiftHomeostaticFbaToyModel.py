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
	"Tc2": {"A": 1, "Carbon2": -1},
	"XY": {"X": -1, "Y": 1}
}

externalMolecules = {
	"Carbon1": 5,
	"Carbon2": 10.5,
	"D": -12.,
	"E": -12.,
	"F": 5.,
	"H": 5.,
	"O2": 15.,
	"X": 1.,
}

# objective = {"B": 2, "C": 2, "G": 2, "NADH": 5, "ATP": 30}
objective = {"Y": 12, "A": 10, "B": 20, "C": 20, "D": 20, "E": 20, "F": 20, "G": 20, "H": 20, "O2": 10, "NADH": 25, "ATP": 30}

fba = FluxBalanceAnalysis(reactions, externalMolecules, objective, objectiveType = "shift_pools", solver = "glpk")
fba.setExternalMoleculeLevels([7 * externalMolecules[x] for x in externalMolecules])
print fba.getOutputMoleculeLevelsChange()

