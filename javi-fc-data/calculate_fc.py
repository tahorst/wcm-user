#! /usr/bin/env python

"""
Aggregate data from Javi's repo to compare to foldChanges.tsv included in wcm.
"""

import csv
import os

import numpy as np
from typing import Dict, List


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
SRC_FILE = os.path.join(FILE_LOCATION, 'kb.tsv')
WCM_FILE = os.path.join(FILE_LOCATION, 'fold_changes.tsv')


def load_file(filename):
	# type: (str) -> List[str]
	"""Load a tsv file."""

	with open(filename) as f:
		reader = csv.reader(f, delimiter='\t')
		return list(reader)

def load_src():
	# type: () -> Dict[str, Dict[str, Dict[str, float]]]
	"""
	Loads and extracts gene regulation data from the source data.

	Returns:
		data: mean and standard deviation for each regulatory pair
			{TF gene name: {regulated gene: {'mean': mean, 'std': std}}}
	"""

	raw_data = load_file(SRC_FILE)

	# Extract fold changes from data
	data = {}
	for line in raw_data:
		# Columns of interest
		tf = line[7]
		regulated_gene = line[9]
		magnitude = float(line[10])

		if tf not in data:
			data[tf] = {}

		data[tf][regulated_gene] = data[tf].get(regulated_gene, []) + [magnitude]

	# Calculate mean and std from data
	processed_data = {}
	for tf, regulated in data.items():
		tf_data = {}
		for gene, fcs in regulated.items():
			fcs = np.abs(fcs)  # Bad!! - ignores annotated condition comparison regulation direction
			tf_data[gene] = {
				'mean': np.mean(fcs),
				'std': np.std(fcs, ddof=1),
				}
		processed_data[tf] = tf_data

	return processed_data

def load_wcm():
	# type: () -> Dict[str, Dict[str, Dict[str, float]]]
	"""
	Loads and extracts gene regulation data from the whole-cell model.

	Returns:
		data: mean and standard deviation for each regulatory pair
			{TF gene name: {regulated gene: {'mean': mean, 'std': std}}}
	"""

	raw_data = load_file(WCM_FILE)[1:]

	# Extract mean and std from data
	data = {}
	for line in raw_data:
		if line[0].startswith('#'):
			continue

		# Columns of interest
		tf = line[0].strip()
		regulated_gene = line[1].strip()
		mean = float(line[2])
		std = float(line[3])
		sign = 1  # np.sign(float(line[5]))

		if tf not in data:
			data[tf] = {}

		data[tf][regulated_gene] = {
			'mean': sign * mean,
			'std': std,
			}

	return data

def compare_data(src_data, wcm_data):
	# type: (Dict[str, Dict[str, Dict[str, float]]], Dict[str, Dict[str, Dict[str, float]]]) -> None
	"""
	Compares regulation from source data to wcm data.

	Args:
		src_data: mean and standard deviation for each regulatory pair in
			the source data
		wcm_data: mean and standard deviation for each regulatory pair in
			the whole-cell model

	Notes:
		dictionary structure for both inputs:
			{TF gene name: {regulated gene: {'mean': mean, 'std': std}}}

	TODO:
		check regulation not included in both wcm and src data
	"""

	# Track statistics
	total_regulation = {}
	discrepancies = {}
	direction_discrepancies = {}

	# Print regulation that is different from source and wcm
	print('TF -> regulated gene: source vs wcm')
	for tf, regulation in src_data.items():
		total_regulation[tf] = len(regulation)
		discrepancies[tf] = 0
		direction_discrepancies[tf] = 0

		for gene, data in regulation.items():
			mean1 = data['mean']
			mean2 = wcm_data.get(tf, {}).get(gene, {}).get('mean', 0)
			if np.abs(mean1 - mean2) > 0.01 and mean2 != 0:
				if np.sign(mean1) != np.sign(mean2):
					direction_discrepancies[tf] += 1
				discrepancies[tf] += 1
				print('{} -> {}: {:.2f} vs {:.2f}'.format(tf, gene, mean1, mean2))

	# Print summary statistics for each TF
	total_direction_discrepancies = 0
	total_discrepancies = 0
	total_interactions = 0
	print('\nTF: opposite direction, different mean, total')
	for tf, total in total_regulation.items():
		tf_dir_disc = direction_discrepancies[tf]
		tf_disc = discrepancies[tf]

		total_direction_discrepancies += tf_dir_disc
		total_discrepancies += tf_disc
		total_interactions += total

		print('{:5s}: {:3} {:3} {:3}'.format(tf, tf_dir_disc, tf_disc, total))
	print('Total: {:3} {:3} {:3}'.format(total_direction_discrepancies, total_discrepancies, total_interactions))


if __name__ == '__main__':
	src_data = load_src()
	wcm_data = load_wcm()

	compare_data(src_data, wcm_data)
