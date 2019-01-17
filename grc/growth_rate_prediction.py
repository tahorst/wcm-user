#! /usr/bin/env python
'''
Exploration of growth rate prediction from ML techniques.

Requires:
	proteins.tsv: data from Schmidt et al. 2015 for monomer expression in various conditions
'''

import argparse
import csv
import os
import time

import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import FactorAnalysis


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
if not os.path.exists(DATA_DIR):
	os.makedirs(DATA_DIR)


def load_data(path):
	'''
	Loads protein data from file

	Args:
		path (str): path to the tsv to load

	Returns:
		ndarray[float]: concentration of proteins, shape=(n conditions, m proteins)
			(in units of uM)
		ndarray[float]: growth rate in each condition, shape=(n conditions,)
		list[str]: gene names for all m proteins
		list[str]: condition name for all n conditions

	Notes:
		The file should have the following row constraints:
			2nd: header names
			3rd: volume in condition
			4th: growth rate in condition
			5th-end: proteins
		the file should have the following column constraints:
			'Gene': contains gene identifier
			'Glucose': starting column for conditions
	'''

	# Load rows from file
	with open(path) as f:
		reader = csv.reader(f, delimiter='\t')

		reader.next()
		header = reader.next()
		volume = reader.next()
		growth_rate = reader.next()
		proteins = np.array(list(reader))

	condition_start_index = header.index('Glucose')
	gene_index = header.index('Gene')

	conditions = header[condition_start_index:]
	genes = list(proteins[:, gene_index])
	volume = np.array(volume[condition_start_index:], float)
	growth_rate = np.array(growth_rate[condition_start_index:], float)
	counts = proteins[:, condition_start_index:]

	counts[np.where(counts == 'NA')] = '0'
	conc = np.array(counts, float) / volume / 602.214  # Avogadro's number to get uM

	return conc.T, growth_rate, genes, conditions

def parse_args():
	'''
	Parses arguments from the command line.

	Returns:
		ArgumentParser namespace: values of variables parsed from the command line
	'''

	default_input = 'proteins.tsv'

	parser = argparse.ArgumentParser()

	parser.add_argument('-i', '--input',
		default=default_input,
		help='File in data/ containing protein counts and growth rates (default: {})'
		.format(default_input))

	return parser.parse_args()

if __name__ == '__main__':
	start = time.time()

	args = parse_args()

	protein_conc, growth_rate, gene_names, conditions = load_data(os.path.join(DATA_DIR, args.input))
	n_conditions = len(conditions)

	prediction = np.zeros(n_conditions)
	actual = np.zeros(n_conditions)

	for idx in range(n_conditions):
		mask = np.ones(n_conditions, bool)
		mask[idx] = False
		x_train = protein_conc[mask, :]
		x_valid = protein_conc[~mask, :]
		y_train = growth_rate[mask]
		y_valid = growth_rate[~mask]

		mean = y_train.mean()

		fa = FactorAnalysis(n_components=n_conditions)
		train_transform = fa.fit_transform(x_train)
		x, _, _, _ = np.linalg.lstsq(train_transform, y_train, rcond=None)
		transform = fa.transform(x_valid)
		y_pred = transform.dot(x)

		prediction[idx] = y_pred + mean
		actual[idx] = y_valid

	plt.figure()
	plt.plot(prediction, actual, 'x')
	plt.plot([1, 5], [1, 5], '--k')
	plt.show()

	print('Completed in {:.2f} min'.format((time.time() - start) / 60))
