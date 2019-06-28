'''
Compares datasets for tRNA changes under different conditions to see if they
are consistent with each other.

Output:
	plot in comparison.png

TODO:
- adjust for absolute expression for rRNA ratio data - should shift all data
points the same amount so it shouldn't affect the correlation.
'''

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


FILE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(FILE_DIR, 'data')

RRNA_FILE = os.path.join(DATA_DIR, 'rrna_ratio.tsv')
PPGPP_FILE = os.path.join(DATA_DIR, 'ppgpp_fold_changes.tsv')
PLOT_FILE = os.path.join(FILE_DIR, 'comparison.png')

# Convert certain WCM keys to gene ID
RNA_CONVERSION = {
	'RNA0-300[c]': 'valZ',
	'RNA0-301[c]': 'lysY',
	'RNA0-302[c]': 'lysZ',
	'RNA0-303[c]': 'lysQ',
	'RNA0-304[c]': 'asnW',
	'RNA0-305[c]': 'ileY',
	'RNA0-306[c]': 'metV',
	}


def load_rrna_data(filename=RRNA_FILE):
	'''
	Load data from tRNA:rRNA ratio file.

	Args:
		filename (str): path to data file to load

	Returns:
		dict {RNA ID (str): ratio (ndarray[float])}: ratios for different
			growth rates for each tRNA

	Notes:
		- Only includes data for tRNA
		- gene ID is the key (eg. aspT)
	'''

	data = {}

	with open(filename) as f:
		reader = csv.reader(f, delimiter='\t')

		# Skip header
		reader.next()

		# Read data into dictionary
		for line in reader:
			rna_id = RNA_CONVERSION.get(line[0], line[0].split('-tRNA')[0])
			ratios = np.array(line[1:], float)
			data[rna_id] = ratios

	return data

def load_ppgpp_data(filename=PPGPP_FILE):
	'''
	Load data from ppGpp fold change file.

	Args:
		filename (str): path to data file to load

	Returns:
		dict {RNA ID (str): fold change (float)}: log2 fold changes for each gene

	Notes:
		- Includes all genes
		- gene ID is the key (eg. aspT)
	'''

	data = {}

	with open(filename) as f:
		reader = csv.reader(f, delimiter='\t')

		# Skip comments and header
		for line in reader:
			if line[0][0] != '#':
				break

		# Read data into dictionary
		for line in reader:
			rna_id = line[0]
			fold_change = float(line[2])
			data[rna_id] = fold_change

	return data

def convert_ratio_to_fold_change(rrna_data, slow_idx=0, fast_idx=4):
	'''
	Calculate fold change from relative expression in fast and slow
	conditions from tRNA:rRNA ratios.

	Args:
		rrna_data (dict): data from load_rrna_data
		slow_idx (int): index for slow growth condition in rrna_data
		fast_idx (int): index for fast growth condition in rrna_data

	Returns:
		dict {RNA ID (str): fold change (float)}: log2 fold changes for each tRNA
	'''

	total_from_rrna = np.sum(rrna_data.values(), axis=0)
	normalization = total_from_rrna[fast_idx] / total_from_rrna[slow_idx]

	rrna_fold_change = {
		key: np.log2(value[slow_idx] / value[fast_idx] * normalization)
		for key, value in rrna_data.items()
		}

	return rrna_fold_change

def plot_comparison(rrna_fc, ppgpp_fc, output=PLOT_FILE):
	'''
	Plots scatter of fold changes from both datasets and saves to a file.

	Args:
		rrna_fc (dict): data from convert_ratio_to_fold_change
		ppgpp_fc (dict): data from load_ppgpp_data
		output (str): path to plot file
	'''

	plt.figure()

	for trna_id in rrna_fc:
		plt.plot(rrna_fc[trna_id], ppgpp_fc[trna_id], 'bo')

	plt.xlabel('rRNA ratio fold change')
	plt.ylabel('ppGpp fold change')

	plt.savefig(output)


if __name__ == '__main__':
	# Load and manipulate data
	rrna_data = load_rrna_data()
	ppgpp_fc = load_ppgpp_data()
	rrna_fc = convert_ratio_to_fold_change(rrna_data)

	# Plot data
	plot_comparison(rrna_fc, ppgpp_fc)
