#! /usr/bin/env python
'''

'''

from __future__ import division

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data', 'wcm_params')
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUTPUT_DIR):
	os.makedirs(OUTPUT_DIR)

CONDITIONS = ['basal', 'with_aa', 'no_oxygen']
PARAMS = ['p', 'ppgpp', 'charging']


if __name__ == '__main__':
	n_rows = len(PARAMS)
	n_cols = len(CONDITIONS) + 1  # +1 for legend column
	plt.figure()

	y_max = np.zeros(n_rows)
	y_min = np.zeros(n_rows)

	for i, param in enumerate(PARAMS):
		print('{}:'.format(param))
		for j, condition in enumerate(CONDITIONS):
			print('  {}:'.format(condition))
			# Read parameter data
			filename = os.path.join(DATA_DIR, '{}_{}.tsv'.format(param, condition))
			with open(filename) as f:
				reader = csv.reader(f, delimiter='\t')
				headers = reader.next()
				data = np.array(list(reader), float)

			for h, d in zip(headers, data.mean(axis=0)):
				print('    {}: {}'.format(h, d))

			# Plot params
			plt.subplot(n_rows, n_cols, i*n_cols + j + 1)
			plt.plot(np.log(data))

			# Find common y limits
			lim = plt.ylim()
			if y_min[i] > lim[0]:
				y_min[i] = lim[0]
			if y_max[i] < lim[1]:
				y_max[i] = lim[1]

			# Label axes
			if j == 0:
				plt.ylabel(param)
			if i == n_rows - 1:
				plt.xlabel(condition)

			# Legend for parameters
			if j == n_cols - 2:
				plt.subplot(n_rows, n_cols, i*n_cols + j + 2)
				plt.plot([np.zeros(len(headers)), np.zeros(len(headers))], [np.zeros(len(headers)), np.zeros(len(headers))])
				plt.legend(headers)
				plt.axis('off')

	# Set common y limits
	for i in range(n_rows):
		for j in range(n_cols):
			plt.subplot(n_rows, n_cols, i*n_cols + j + 1)
			plt.ylim([y_min[i], y_max[i]])

	plt.tight_layout()
	plt.savefig(os.path.join(OUTPUT_DIR, 'wcm_params.png'))
