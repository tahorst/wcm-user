#! /usr/bin/env python

"""
ODE network to generate example data for fit_fold_changes.py.
"""

from __future__ import absolute_import, division, print_function

import argparse
import os
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUTPUT_DIR):
	os.makedirs(OUTPUT_DIR)

# Values to test with before real data
## Structure based on sim_data.process.metabolism.reactionStoich
REACTIONS = {
	'r1f': {'a': -1, 'b': 1},
	'r1r': {'b': -1, 'a': 1},
	'r2f': {'a': -1, 'c': 1},
	'r2r': {'c': -1, 'a': 1},
	'r3f': {'c': -1, 'd': -1, 'b': 1, 'e': 1},
	'r3r': {'b': -1, 'e': -1, 'c': 1, 'd': 1},
	'r4f': {'b': -1, 'd': 1},
	'r4r': {'d': -1, 'b': 1},
	'r5f': {'e': -1, 'd': 1},
	'r5r': {'d': -1, 'e': 1},
	}
## Structure based on sim_data.process.metabolism.reactionCatalysts
ENZYMES = {
	'r1f': ['e1'],
	'r1r': ['e2'],
	'r2f': ['e3'],
	'r2r': ['e3'],
	'r3f': ['e4'],
	'r3r': ['e4'],
	'r4f': ['e5'],
	'r4r': ['e5'],
	'r5f': ['e6'],
	'r5r': ['e7'],
	}
## Example K_M values for reactants in the reaction network
K_M = {
	'r1f': {'a': 10.},
	'r1r': {'b': 30.},
	'r2f': {'a': 5.},
	'r2r': {'c': 0.5},
	'r3f': {'c': 1., 'd': 5.},
	'r3r': {'b': 2., 'e': 10.},
	'r4f': {'b': 5.},
	'r4r': {'d': 3.},
	'r5f': {'e': 1.},
	'r5r': {'d': 20.},
	}
## Example k_cat values for reactions in the reaction network
K_CAT = {
	'r1f': 2.,
	'r1r': 1.,
	'r2f': 3.,
	'r2r': 0.3,
	'r3f': 10.,
	'r3r': 1.,
	'r4f': 12.,
	'r4r': 5.,
	'r5f': 8.,
	'r5r': 0.5,
	}
# Starting metabolite concentrations
METABOLITE_CONC = {
	'a': 1.,
	'b': 1.,
	'c': 1.,
	'd': 1.,
	'e': 1.,
	}
# Starting enzyme concentrations
ENZYME_CONC = {
	'e1': 5.,
	'e2': 2.,
	'e3': 1.,
	'e4': 10.,
	'e5': 0.1,
	'e6': 50.,
	'e7': 5.,
	}
METABOLITES = sorted(METABOLITE_CONC)
METABOLITE_INDEX = {m: i for i, m in enumerate(METABOLITES)}
N_METABOLITES = len(METABOLITES)


def dcdt(c, t, ko=None):
	"""
	Find change in metabolite concentrations.
	"""

	dc = np.zeros(N_METABOLITES)

	for rxn, stoich in REACTIONS.items():
		enzyme = ENZYMES[rxn][0]

		# Skip reaction if enzyme is knocked out
		if ko is not None and enzyme == ko:
			continue

		rate = K_CAT[rxn] * ENZYME_CONC[enzyme]
		for met, km in K_M[rxn].items():
			conc = c[METABOLITE_INDEX[met]]
			rate *= conc / (km + conc)

		for met, direction in stoich.items():
			dc[METABOLITE_INDEX[met]] += direction * rate

	return dc

def plot_solution(t, sol, filename):
	"""
	Plot concentration traces over time.
	"""

	plt.figure()
	plt.plot(t, sol)
	plt.savefig(os.path.join(OUTPUT_DIR, filename + '.png'))

def parse_args():
	# type: () -> argparse.Namespace
	"""
	Parses arguments from the command line.

	Returns:
		values of variables parsed from the command line
	"""

	parser = argparse.ArgumentParser()

	parser.add_argument('-o', '--output',
		default='solution',
		help='Output plot filename.')

	return parser.parse_args()

if __name__ == '__main__':
	start = time.time()
	args = parse_args()

	# Solve reaction network
	co = np.array([METABOLITE_CONC[c] for c in METABOLITES])
	t = np.linspace(0, 100, 10001)

	for ko in ENZYME_CONC:
		print('Knockout: {}'.format(ko))
		sol = odeint(dcdt, co, t, args=(ko,))

		# Plot results
		filename = '{}_{}'.format(args.output, ko)
		plot_solution(t, sol, filename)

		# Report results
		c_final = sol[-1, :]
		print(c_final)

	print('Completed in {:.2f} min'.format((time.time() - start) / 60))
