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

from ode_network import REACTIONS, REACTION_ENZYMES


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))


def load_test_data():
	"""
	Load made up test data.
	"""

	# Test reaction network
	reactions = REACTIONS

	# Test reaction catalysts
	enzymes = REACTION_ENZYMES

	# n samples x m metabolites
	fcs = np.array([
		[1.01, 0.99, 1.03, 0.95, 1.02],
		[1.2, 0.5, 1.1, 0.8, 0.9],
		[0.8, 1.3, 0.95, 1.1, 1.05],
		[1.2, 1.1, 0.01, 1.3, 0.01],
		[1.1, 0.8, 1.6, 0.7, 0.01],
		[1.3, 1.2, 1.4, 0.3, 0.5],
		[1.05, 0.95, 1.1, 2.1, 0.4],
		])

	# n samples x m enzymes
	n_enzymes = len({e for es in enzymes.values() for e in es})
	kos = np.vstack((
		np.ones(n_enzymes),
		-1 * (np.eye(n_enzymes) - 1)
		))

	return reactions, enzymes, fcs, kos

def load_data(test=True):
	"""
	Get reaction network, fold change data and KO data.

	TODO:
		- sample example data from actual ODE solution with noise
		- load full dataset
		- align ids with order in init_network()
	"""

	if test:
		return load_test_data()
	else:
		raise NotImplementedError('Need to implement real data. Try running with test data option: -t')

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
	b2s2 = b2 * s2

	# dL/db2
	db2 = da2 * s2

	# dL/dW2
	dW2 = -np.tile((da2 * b2s2)[:, None], (1, n_metabolites)) / (W2 + np.tile(a1, (n_reactions, 1))) * (W2 != 0)

	# dL/dW1
	a1_tile = np.tile(a1, (n_reactions, 1))
	denom = a1_tile * (W2 + a1_tile)
	num = np.zeros_like(denom)
	active_mask = denom != 0
	num[active_mask] = (W2 / denom)[active_mask]
	da1 = da2.dot(np.tile(b2s2[:, None], (1, n_metabolites)) * num)
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

	# Data parameters
	parser.add_argument('-t', '--test_data',
		action='store_true',
		help='If set, loads test data, otherwise loads actual data.')

	# Learning parameters
	parser.add_argument('-l', '--learning-rate',
		type=float,
		default=0.01,
		help='Specify learning rate for training.')
	parser.add_argument('-e', '--epochs',
		type=int,
		default=500,
		help='Number of epochs to run.')
	parser.add_argument('-u', '--update-every',
		type=int,
		default=10,
		help='Number of epochs to run before printing an update.')

	return parser.parse_args()

if __name__ == '__main__':
	start = time.time()

	args = parse_args()

	reactions, enzymes, fcs, kos = load_data(args.test_data)
	W1, W2, b2, W3, K = init_network(reactions, enzymes)

	# Learning parameters
	lr = args.learning_rate
	n_epochs = args.epochs
	update_every = args.update_every

	# Find starting loss
	loss, _, _, _, _ = forward(np.ones(fcs.shape[1]), W1, W2, b2, W3)
	print('Initial loss: {:.4f}'.format(loss.sum()))

	# Backpropagation to learn model parameters
	for epoch in range(1, n_epochs + 1):
		for fc, ko in zip(fcs, kos):
			W3_ko = apply_ko(ko, K, W3)
			loss, a1, s2, a2, a3 = forward(fc, W1, W2, b2, W3_ko)
			dW1, dW2, db2 = backward(fc, W1, W2, b2, W3_ko, loss, a1, s2, a2, a3)

			# TODO: check values non-negative after update?
			W1 -= lr * dW1
			W2 -= lr * dW2
			b2 -= lr * db2

		if epoch % update_every == 0:
			print('Epoch {}: total loss = {:.4f}'.format(epoch, loss.sum()))

	# Print final state of parameters
	print('Concentrations:')
	print(np.diag(W1))

	print('vmax:')
	print(b2)

	print('KM:')
	print(W2)

	print('Completed in {:.2f} min'.format((time.time() - start) / 60))
