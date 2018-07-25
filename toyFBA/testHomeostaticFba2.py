import numpy as np
import cvxpy

MAX = np.inf

# A + B -> C
# 2C -> D
# A + 2D -> E
# B + F -> G
# Taret C, D, E, G conc
metabolic_reactions = [
[-1, 0, -1, 0, 0, 0, 0, 0],
[-1, 0, 0, -1, 0, 0, 0, 0],
[1, -2, 0, 0, -10, 0, 0, 0],
[0, 1, -2, 0, 0, -10, 0, 0],
[0, 0, 1, 0, 0, 0, -10, 0],
[0, 0, 0, -1, 0, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0, -10],
[0, 0, 0, 0, 1, 0, 0, 0],
[0, 0, 0, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 0, 1, 0],
[0, 0, 0, 0, 0, 0, 0, 1]
]

transport_reactions = [
[1, 0, 0],
[0, 1, 0],
[0, 0, 0],
[0, 0, 0],
[0, 0, 0],
[0, 0, 1],
[0, 0, 0],
[0, 0, 0],
[0, 0, 0],
[0, 0, 0],
[0, 0, 0]
]

pseudo_rxns = [
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[-1, 1, -1, 0, 0, 0, 1, 0, 0, 0],
[-1, 1, 0, -1, 0, 0, 0, 1, 0, 0],
[-1, 1, 0, 0, -1, 0, 0, 0, 1, 0],
[-1, 1, 0, 0, 0, -1, 0, 0, 0, 1],
]


S_matrix = np.hstack((metabolic_reactions, transport_reactions, pseudo_rxns))

num_metabolites, num_reactions = S_matrix.shape
nMetReactions = len(metabolic_reactions[0])
nPseudoReactions = len(pseudo_rxns[0])

# Construct the problem
fluxes = cvxpy.Variable(num_reactions)

# One-hot vector indicating location of biomass reaction
c = np.zeros(num_reactions)
c[-(nPseudoReactions-2):] = 1
c[-(nPseudoReactions-1)] = 1
# c[-4:] = 1

# Maximize biomass reaction
objective = cvxpy.Minimize(c*fluxes)

min_fluxes = np.zeros(num_reactions)
min_fluxes[-nPseudoReactions] = 1
max_fluxes = np.hstack((np.ones(nMetReactions) * MAX, 20, 10, 10, 1, np.ones(nPseudoReactions - 1) * MAX))
constraints = [
    S_matrix*fluxes == 0,
    min_fluxes <= fluxes,
    fluxes <= max_fluxes,
    ]

# Maximize biomass reaction
objective = cvxpy.Minimize(c*fluxes)

min_fluxes = np.zeros(num_reactions)
min_fluxes[6] = 1
max_fluxes = np.hstack((np.ones(4) * MAX, 20, 10, 1, 1, np.ones(4) * MAX))
constraints = [
    S_matrix*fluxes == 0,
    min_fluxes <= fluxes,
    fluxes <= max_fluxes,
    ]

prob = cvxpy.Problem(objective, constraints)

# The optimal objective is returned by prob.solve().
result = prob.solve(solver=cvxpy.GUROBI)

print fluxes.value
