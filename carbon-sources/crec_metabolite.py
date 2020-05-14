#! /usr/bin/env python

"""
Compare metabolite changes in carbon sources where creC is active.

Cariss et al. Defining the Growth Conditions and Promoter-Proximal DNA
Sequences Required for Activation of Gene Expression by CreBC in Escherichia
coli. 2008. shows activity of cre promoter in different carbon sources.

Results:
- PEP, Arg and maybe cAMP appear to increase with cre expression
"""

from __future__ import absolute_import, division, print_function

import csv
import os

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
FOLD_CHANGE_FILE = os.path.join(DATA_DIR, 'met-fc.tsv')
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUTPUT_DIR):
	os.makedirs(OUTPUT_DIR)

# Relative expression for active conditions estimated from Cariss et al. Fig 5
UP_CONDITIONS = {
	' M9 acetate 3.6g per L': 330,
	' M9 lactate 5g per L': 200,
	' M9 pyruvate 5g per L': 720,
	}

# Relative expression for inactive conditions estimated from Cariss et al. Fig 5
DOWN_CONDITIONS = {
	' M9 gluconate 2g per L': 20,
	' M9 glucose 2g per L first experiment': 30,
	' M9 glucose 2g per L plus CAA 2g per L': 25,
	}


def load_data():
	"""
	Load fold change data for molecules of interest
	"""

	with open(FOLD_CHANGE_FILE) as f:
		data = csv.reader(f, delimiter='\t')

		conditions = data.next()[1:]
		data = np.array(list(data))

	metabolites = data[:, 0]
	fold_changes = np.array(data[:, 1:], float)

	up_idx = np.array([conditions.index(condition) for condition in sorted(UP_CONDITIONS)])
	down_idx = np.array([conditions.index(condition) for condition in sorted(DOWN_CONDITIONS)])

	up_fold_changes = fold_changes[:, up_idx]
	down_fold_changes = fold_changes[:, down_idx]

	return metabolites, up_fold_changes, down_fold_changes

def plot_data(metabolites, up, down, label):
	"""
	Plot adjusted fold changes for each metabolite.
	"""

	x = np.arange(len(metabolites))
	ave_up = up.mean(1)
	ave_down = down.mean(1)
	diff = ave_up - ave_down
	sorted_idx = np.argsort(diff)

	plt.figure(figsize=(10, 10))

	plt.plot(x, up[sorted_idx, :], 'rx')
	plt.plot(x, down[sorted_idx, :], 'bo')

	plt.xticks(x, metabolites[sorted_idx], rotation=45, ha='right', fontsize=6)

	plt.tight_layout()
	plt.savefig(os.path.join(OUTPUT_DIR, '{}.png'.format(label)))
	plt.close('all')

def plot_relations(metabolites, up, down):
	"""
	Plot concentrations vs amount of expression for each metabolite.
	"""

	n_mets = len(metabolites)
	rows = int(np.ceil(np.sqrt(n_mets)))
	cols = int(np.ceil(n_mets / rows))

	up_exp = np.array([f for c, f in sorted(UP_CONDITIONS.items())])
	down_exp = np.array([f for c, f in sorted(DOWN_CONDITIONS.items())])

	fig = plt.figure(figsize=(10, 10))
	gs = gridspec.GridSpec(rows, cols)
	for i, (met, u, d) in enumerate(zip(metabolites, up, down)):
		row = i // cols
		col = i % cols
		ax = plt.subplot(gs[row, col])

		ax.plot(np.log(up_exp), u, 'rx')
		ax.plot(np.log(down_exp), d, 'bo')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		if row < rows - 1:
			ax.set_xticklabels([])
		else:
			ax.tick_params(axis='x', labelsize=6)
		ax.tick_params(axis='y', labelsize=6)

		ax.set_title(met, fontsize=8)

	fig.add_subplot(gs[:, :], frameon=False)
	plt.tick_params(labelcolor="none", bottom=False, left=False)
	plt.xlabel('log(cre expression)')
	plt.ylabel('Metabolite Fold Change')
	plt.tight_layout()
	plt.savefig(os.path.join(OUTPUT_DIR, 'correlation.png'))
	plt.close('all')


if __name__ == '__main__':
	metabolites, up_fold_changes, down_fold_changes = load_data()
	plot_data(metabolites, up_fold_changes, down_fold_changes, 'original')
	plot_relations(metabolites, up_fold_changes, down_fold_changes)

	up_factors = np.array([f for c, f in sorted(UP_CONDITIONS.items())])
	down_factors = np.array([f for c, f in sorted(DOWN_CONDITIONS.items())])
	up_adjusted = up_fold_changes / up_factors
	down_adjusted = down_fold_changes / down_factors
	plot_data(metabolites, up_adjusted, down_adjusted, 'adjusted')
