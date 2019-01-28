#! /usr/bin/env python

'''
Plots metabolite traces from Link et al. Real-time metabolome profiling of the
metabolic switch between starvation and growth. 2015.

Output plots in out/
'''

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
OUT_DIR = os.path.join(FILE_LOCATION, 'out')

# Data from Fig. 1 from Link et al. (shift from aerobic to anaerobic at 10 min)
ANAEROBIC_INFO_FILE = os.path.join(DATA_DIR, 'ecoli-metabolite-info.tsv')
ANAEROBIC_DATA1_FILE = os.path.join(DATA_DIR, 'ecoli-metabolite-data1.tsv')
ANAEROBIC_DATA2_FILE = os.path.join(DATA_DIR, 'ecoli-metabolite-data2.tsv')
ANAEROBIC_DATA3_FILE = os.path.join(DATA_DIR, 'ecoli-metabolite-data3.tsv')

# Amino acid names in data file
# Does not have data for cysteine, glycine, methionine, serine, tyrosine, or selenocysteine
# L-Leucine and L-Isoleucine will have the same trace with no mass difference
AA_NAMES = [
	'L-Alanine', 'L-Arginine', 'L-Asparagine', 'L-Aspartate', 'L-Glutamate',
	'L-Glutamine', 'L-Histidine', 'L-Leucine', 'L-Lysine', 'L-Phenylalanine',
	'L-Proline', 'L-Threonine', 'L-Tryptophan', 'L-Valine',
	]


def get_aa_indices():
	'''
	Find amino acid indices for the data files.

	Returns:
		dict (aa (str): index (int)): mapping of amino acid name to index in data files
	'''

	indices = {aa: None for aa in AA_NAMES}

	with open(ANAEROBIC_INFO_FILE) as f:
		reader = csv.reader(f, delimiter='\t')
		for line in reader:
			metabolite = line[-1]
			if metabolite in indices:
				indices[metabolite] = int(line[0])

	return indices

def plot_trace(indices, filename, output):
	'''
	Plots traces from data files for amino acids.

	Args:
		indices (dict): mapping of amino acid (str) to data file index (int)
		filename (str): input data filename
		output (str): output plot path
	'''

	aa_names = indices.keys()
	n_aa = len(aa_names)
	index_mapping = {indices[aa]: i for i, aa in enumerate(indices)}

	# Extract data traces for amino acids
	with open(filename) as f:
		reader = csv.reader(f, delimiter='\t')

		time = np.array(reader.next()[1:], float)

		aa_traces = np.zeros((n_aa, len(time)))

		for line in reader:
			index = int(line[0])
			if index in index_mapping:
				aa_traces[index_mapping[index], :] = line[1:]

	# Plot traces
	plt.figure()
	n_rows = np.ceil(np.sqrt(n_aa))
	n_cols = np.ceil(n_aa / n_rows)
	ma_size = 3  # points from either side to include

	for idx, (aa, trace) in enumerate(zip(aa_names, aa_traces)):
		ax = plt.subplot(n_rows, n_cols, idx + 1)

		ma = np.convolve(trace, np.ones(2*ma_size + 1), mode='valid')
		ax.plot(time[ma_size:-ma_size], ma / ma[0])
		ax.axvline(10, color='r', linestyle='--')
		ax.set_title(aa, fontsize=8)
		ax.tick_params(labelsize=7)

	plt.tight_layout()
	plt.savefig(output)


if __name__ == '__main__':
	aa_indices = get_aa_indices()
	plot_trace(aa_indices, ANAEROBIC_DATA1_FILE, os.path.join(OUT_DIR, 'data1_ma.png'))
	plot_trace(aa_indices, ANAEROBIC_DATA2_FILE, os.path.join(OUT_DIR, 'data2_ma.png'))
	plot_trace(aa_indices, ANAEROBIC_DATA3_FILE, os.path.join(OUT_DIR, 'data3_ma.png'))
