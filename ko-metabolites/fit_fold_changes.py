#! /usr/bin/env python

"""
Use metabolite fold change data from gene KOs to fit kinetic parameters.

TODO:
	- handle multiple enzymes for same reaction from model
	- handle metabolites in reaction network but without fold change data
	- handle KOs not in reaction network
"""

from __future__ import absolute_import, division, print_function

import argparse
import os
import time

import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

# Values to test with before real data
## Structure based on sim_data.process.metabolism.reactionStoich
TEST_REACTIONS = {
	'r1f': {'a': -1, 'b': 1},
	'r1r': {'b': -1, 'a': 1},
	'r2': {'a': -1, 'c': 1},
	'r3': {'c': -1, 'd': -1, 'b': 1, 'e': 1},
	'r4f': {'b': -1, 'd': 1},
	'r4r': {'d': -1, 'b': 1},
	'r5': {'e': -1, 'd': 1},
	}
## Structure based on sim_data.process.metabolism.reactionCatalysts
TEST_ENZYMES = {
	'r1f': ['e1'],
	'r1r': ['e2'],
	'r2': ['e3'],
	'r3': ['e4'],
	'r4f': ['e5'],
	'r4r': ['e5'],
	'r5': ['e6'],
	}
N_TEST_METABOLITES = len({m for stoich in TEST_REACTIONS.values() for m in stoich})
N_TEST_ENZYMES = len({e for enzymes in TEST_ENZYMES.values() for e in enzymes})
N_TEST_REACTIONS = len(TEST_REACTIONS)
N_TEST_SAMPLES = N_TEST_ENZYMES + 1


def load_data():
	"""
	Get fold change and KO data.

	TODO:
		- sample example data from actual ODE solution with noise
		- load full dataset
		- align ids with order in init_network()
	"""

	# N_TEST_SAMPLES x N_TEST_METABOLITES
	fcs = np.array([
		[1.01, 0.99, 1.03, 0.95, 1.02],
		[1.2, 0.5, 1.1, 0.8, 0.9],
		[0.8, 1.3, 0.95, 1.1, 1.05],
		[1.2, 1.1, 0.01, 1.3, 0.01],
		[1.1, 0.8, 1.6, 0.7, 0.01],
		[1.3, 1.2, 1.4, 0.3, 0.5],
		[1.05, 0.95, 1.1, 2.1, 0.4],
		])

	# N_TEST_SAMPLES x N_TEST_ENZYMES
	kos = np.vstack((
		np.ones(N_TEST_ENZYMES),
		-1 * (np.eye(N_TEST_ENZYMES) - 1)
		))

	return fcs, kos

def init_network(reactions, enzymes):
	"""
	Create W1, W2, b2, W3, K.

	TODO:
		- initialize values for W1, W2, b2
		- add metabolite inhibition for W2 (or equivalent matrix)
	"""

	# IDs for matrix setup
	reaction_ids = reactions.keys()
	metabolite_ids = sorted({m for stoich in reactions.values() for m in stoich})
	enzyme_ids = sorted({e for es in enzymes.values() for e in es})

	# Quick index lookups
	metabolite_idx = {m: i for i, m in enumerate(metabolite_ids)}
	enzyme_idx = {e: i for i, e in enumerate(enzyme_ids)}

	# Dimensions for matrix setup
	n_reactions = len(reaction_ids)
	n_metabolites = len(metabolite_ids)
	n_enzymes = len(enzyme_ids)

	# W1
	W1 = np.eye(n_metabolites)  # TODO: initialize to better value

	# W2 and W3
	W2 = np.zeros((n_reactions, n_metabolites))
	W3 = np.zeros((n_metabolites, n_reactions))
	for i, rxn in enumerate(reaction_ids):
		for met, stoich in reactions[rxn].items():
			j = metabolite_idx[met]

			# Stoich for mass balance layer
			W3[j, i] = stoich

			# Only include reactants in rate equations
			if stoich < 0:
				W2[i, j] = 1  # TODO: initialize to better value

	# b2
	b2 = np.ones((n_reactions))  # TODO: initialize to better value

	# K
	K = np.zeros((n_reactions, n_enzymes))
	for i, rxn in enumerate(reaction_ids):
		for enzyme in enzymes[rxn]:
			j = enzyme_idx[enzyme]
			K[i, j] = 1

	return W1, W2, b2, W3, K

def apply_ko(k, K, W):
	# type: (np.ndarray, np.ndarray, np.ndarray) -> np.ndarray
	"""
	Applies a knockout vector from a sample to the reaction mapping matrix.

	Args:
		k: sample KO, 1 if enzyme is present, 0 if KO'd (n enzymes)
		K: KO to reaction mapping, 1 if reaction is catalyzed by enzyme
			(r reactions, n enzymes)
		W: reaction to metabolite mapping, 1 if reaction produces metabolite,
			-1 if reaction consumes metabolite, 0 if metabolite not present
			(m metabolites, r reactions)

	Returns:
		W_ko: updated mapping that removes reactions catalyzed by KO'd enzymes
	"""

	W_ko = np.zeros_like(W)
	active_rxns = K.dot(k) != 0
	W_ko[:, active_rxns] = W[:, active_rxns]

	return W_ko

def sample_saturation(a, W):
	"""
	Calculate saturation terms for each reaction.
	"""

	return np.prod(a / (W + a), axis=1)


def forward(fc, W1, W2, b2, W3):
	"""
	Forward propagation through network to find loss function.
	"""

	a1 = W1.dot(fc)
	s2 = sample_saturation(a1, W2)
	a2 = b2 * s2
	a3 = W3.dot(a2)
	loss = 0.5 * a3**2

	return loss, a1, s2, a2, a3

def backward(fc, W1, W2, b2, W3, loss, a1, s2, a2, a3):
	"""
	Back propagation to find gradients.

	TODO:
		- copy down forward code except for loss calc
		- verify derivatives with test input to double check correct dimension math
		- optimize tile functions?
	"""

	# dL/da2 and other values used in multiple calculations
	n_metabolites = len(a1)
	n_reactions = len(s2)
	da2 = a3.dot(W3)
	b2s22 = b2 * s2**2

	# dL/db2
	db2 = da2 * s2

	# dL/dW2
	dW2 = np.tile((da2 * b2s22)[:, None], (1, n_metabolites)) * np.tile(a1, (n_reactions, 1)) * (W2 != 0)

	# dL/dW1
	denom = (np.tile(a1**2, (n_reactions, 1)) * W2)
	num = np.zeros_like(denom)
	active_mask = denom != 0
	num[active_mask] = 1 / denom[active_mask]
	da1 = da2.dot(np.tile(b2s22[:, None], (1, n_metabolites)) * num)
	dW1 = da1 * np.diag(fc)

	return dW1, dW2, db2

def parse_args():
	# type: () -> argparse.Namespace
	"""
	Parses arguments from the command line.

	Returns:
		values of variables parsed from the command line
	"""

	parser = argparse.ArgumentParser()

	parser.add_argument('-a', '--abc',
		type=int,
		default=1,
		help='help')
	parser.add_argument('--true',
		action='store_true',
		help='help')

	return parser.parse_args()

if __name__ == '__main__':
	start = time.time()

	args = parse_args()

	fcs, kos = load_data()
	W1, W2, b2, W3, K = init_network(TEST_REACTIONS, TEST_ENZYMES)

	lr = 0.01
	n_epochs = 500
	for epoch in range(n_epochs):
		for fc, ko in zip(fcs, kos):
			W3_ko = apply_ko(ko, K, W3)
			loss, a1, s2, a2, a3 = forward(fc, W1, W2, b2, W3_ko)
			dW1, dW2, db2 = backward(fc, W1, W2, b2, W3_ko, loss, a1, s2, a2, a3)

			# TODO: check values non-negative after update?
			W1 -= lr * dW1
			W2 -= lr * dW2
			b2 -= lr * db2

		print('Epoch {}: total loss = {:.4f}'.format(epoch, loss.sum()))

	print('Concentrations:')
	print(np.diag(W1))

	print('vmax:')
	print(b2)

	print('KM:')
	print(W2)

	print('Completed in {:.2f} min'.format((time.time() - start) / 60))
