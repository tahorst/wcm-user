import numpy as np
import cvxpy

MAX = np.inf

metabolic_reactions = [
[-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[1, -1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 1, -1, 0, -1, 0.8, 0.8, -1, -1, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 1, -1, -1, 0, 0, -1, 1, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0],
[-1, 2, -2, 0, 0, 0, 0, 2, 0, -1, 1, 1],
[0, 2, -2, 0, 0, 2, 2, 0, -4, -2, 2, -1],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

transport_rxns = [
[1, 1, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, -1, 0, 0, 0],
[0, 0, 0, 0, -1, 0, 0],
[0, 0, 1, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 1, 0],
[0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 1],
[-1, 0, 0, 0, 0, 0, 0],
[0, -1, 0, 0, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0],
[0, 0, 0, 0, 1, 0, 0],
[0, 0, -1, 0, 0, 0, 0],
[0, 0, 0, 0, 0, -1, 0],
[0, 0, 0, 0, 0, 0, -1],
[0, 0, 0, 0, 0, 0, 0]
]

source_reactions = [
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0],
[1, 0, 0, 0, 0, 0, 0, 0],
[0, 1, 0, 0, 0, 0, 0, 0],
[0, 0, -1, 0, 0, 0, 0, 0],
[0, 0, 0, -1, 0, 0, 0, 0],
[0, 0, 0, 0, 1, 0, 0, 0],
[0, 0, 0, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 0, 1, 0],
[0, 0, 0, 0, 0, 0, 0, -1]
]

biomass_reaction = [
[0],
[0],
[-1],
[0],
[0],
[-1],
[0],
[-1],
[-10],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[1]
]

S_matrix = np.hstack((metabolic_reactions,  transport_rxns, source_reactions,  biomass_reaction))

num_metabolites, num_reactions = S_matrix.shape

# Construct the problem
fluxes = cvxpy.Variable(num_reactions)

# One-hot vector indicating location of biomass reaction
c = np.zeros(num_reactions)
c[-1] = 1

# Maximize biomass reaction
objective = cvxpy.Maximize(c*fluxes)

min_fluxes = np.zeros(num_reactions)
max_fluxes = np.hstack((np.ones(12) * MAX, 10.5, 10.5, 5, 12, 12, 5, 15, np.ones(9) * MAX))
constraints = [
    S_matrix*fluxes == 0,
    min_fluxes <= fluxes,
    fluxes <= max_fluxes,
    ]

prob = cvxpy.Problem(objective, constraints)

# The optimal objective is returned by prob.solve().
result = prob.solve(solver=cvxpy.GUROBI)

print fluxes.value
