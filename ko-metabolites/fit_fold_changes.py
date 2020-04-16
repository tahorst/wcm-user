#! /usr/bin/env python

"""
Use metabolite fold change data from gene KOs to fit kinetic parameters.
"""

from __future__ import absolute_import, division, print_function

import argparse
import os
import time

import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

# Values to test with before real data
N_TEST_SAMPLES = 10
N_TEST_METABOLITES = 5
N_TEST_ENZYMES = 6
N_TEST_REACTIONS = 7


def load_data():
	"""
	Get fold change and KO data.

	TODO:
		- update for example data
		- load full dataset
	"""

	fcs = np.ones((N_TEST_SAMPLES, N_TEST_METABOLITES))
	kos = np.ones((N_TEST_SAMPLES, N_TEST_ENZYMES))

	return fcs, kos

def init_network():
	"""
	Create W1, W2, b2, W3, K.

	TODO:
		- update for example data
		- load from reaction network
	"""

	W1 = np.eye(N_TEST_METABOLITES)
	W2 = np.zeros((N_TEST_REACTIONS, N_TEST_METABOLITES))
	b2 = np.ones((N_TEST_REACTIONS))
	W3 = np.zeros((N_TEST_METABOLITES, N_TEST_REACTIONS))
	K = np.zeros((N_TEST_REACTIONS, N_TEST_ENZYMES))

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
		- verify derivatives with test input to double check correct dimensions
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
	dW1 = da2.dot(np.tile(b2s22[:, None], (1, n_metabolites)) / (np.tile(a1**2, (n_reactions, 1)) * W2)) * np.diag(fc)

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
	W1, W2, b2, W3, K = init_network()

	lr = 0.01
	for fc, ko in zip(fcs, kos):
		W3_ko = apply_ko(ko, K, W3)
		loss, a1, s2, a2, a3 = forward(fc, W1, W2, b2, W3_ko)
		dW1, dW2, db2 = backward(fc, W1, W2, b2, W3_ko, loss, a1, s2, a2, a3)

		# TODO: check values non-negative after update?
		W1 -= lr * dW1
		W2 -= lr * dW2
		b2 -= lr * db2

	print('Completed in {:.2f} min'.format((time.time() - start) / 60))
