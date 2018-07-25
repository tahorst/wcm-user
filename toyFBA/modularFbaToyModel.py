import numpy as np
import cvxpy
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
	"Tc2": {"A": 1, "Carbon2": -1},
	# "Tf": {"F": 1, "Fext": -1},
	# "Td": {"D": -1, "Dext": 1},
	# "Te": {"E": -1, "Eext": 1},
	# "Th": {"H": 1, "Hext": -1},
	# "To2": {"O2": 1, "Oxygen": -1},
	# "C1ex": {"Carbon1": 1},
	# "C2ex": {"Carbon2": 1},
	# "Dex": {"Dext": -1},
	# "Eex": {"Eext": -1},
	# "Fex": {"Fext": 1},
	# "Hex": {"Hext": 1},
	# "Oex": {"Oxygen": 1},
	"Biomass": {"C": -1, "F": -1, "H": -1, "ATP": -10, "Biomass": 1}
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

objective = {"Biomass": 1}
# objective = {"C": 1, "F": 1, "H": 1, "ATP": 10}

fba = FluxBalanceAnalysis(reactions, externalMolecules, objective, solver = "glpk")
fba.setExternalMoleculeLevels([externalMolecules[x] for x in externalMolecules])
print fba.getBiomassReactionFlux()


# S_matrix = np.hstack((metabolic_reactions,  transport_rxns, source_reactions,  biomass_reaction))

# num_metabolites, num_reactions = S_matrix.shape

# # Construct the problem
# fluxes = cvxpy.Variable(num_reactions)

# # One-hot vector indicating location of biomass reaction
# c = np.zeros(num_reactions)
# c[-1] = 1

# # Maximize biomass reaction
# objective = cvxpy.Maximize(c*fluxes)

# min_fluxes = np.zeros(num_reactions)
# max_fluxes = np.hstack((np.ones(12) * MAX, 10.5, 10.5, 5, 12, 12, 5, 15, np.ones(9) * MAX))
# constraints = [
#     S_matrix*fluxes == 0,
#     min_fluxes <= fluxes,
#     fluxes <= max_fluxes,
#     ]

# prob = cvxpy.Problem(objective, constraints)

# # The optimal objective is returned by prob.solve().
# result = prob.solve(solver=cvxpy.GUROBI)

# print fluxes.value
