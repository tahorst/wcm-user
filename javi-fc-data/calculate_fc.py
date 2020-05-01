#! /usr/bin/env python

"""
Aggregate data from Javi's repo to compare to foldChanges.tsv included in wcm.
"""

import argparse
import csv
import os

import numpy as np
from typing import Dict, List, Tuple


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
SRC_FILE = os.path.join(FILE_LOCATION, 'kb.tsv')
WCM_FILE = os.path.join(FILE_LOCATION, 'fold_changes.tsv')
NEW_FILE = os.path.join(FILE_LOCATION, 'fc_single_shift.tsv')
SHIFTS_FILE = os.path.join(FILE_LOCATION, 'shifts.tsv')
GENES_FILES = os.path.join(FILE_LOCATION, 'gene_names.tsv')


def load_file(filename):
	# type: (str) -> List[str]
	"""Load a tsv file."""

	with open(filename) as f:
		reader = csv.reader(f, delimiter='\t')
		return list(reader)

def load_src(attempt_match):
	# type: (bool) -> Dict[str, Dict[str, Dict[str, float]]]
	"""
	Loads and extracts gene regulation data from the source data.

	Args:
		attempt_match: if True, handles data in way that was most likely the
			original processing, otherwise uses expected directionality

	Returns:
		data: mean and standard deviation for each regulatory pair
			{TF gene name: {regulated gene: {'mean': mean, 'std': std}}}
	"""

	raw_data = load_file(SRC_FILE)
	shifts = load_shifts()

	# Extract fold changes from data
	data = {}
	for line in raw_data:
		# Columns of interest
		condition = (line[4], line[5])
		tf = line[7]
		regulated_gene = line[9]
		magnitude = float(line[10])

		if tf not in data:
			data[tf] = {}
		direction = shifts[condition][tf]

		data[tf][regulated_gene] = data[tf].get(regulated_gene, []) + [direction * magnitude]

	# Calculate mean and std from data
	processed_data = {}
	for tf, regulated in data.items():
		tf_data = {}
		for gene, fcs in regulated.items():
			if attempt_match:
				fcs = np.abs(fcs)  # Bad!! - ignores annotated condition comparison regulation direction
			tf_data[gene] = {
				'mean': np.mean(fcs),
				'std': np.std(fcs, ddof=1),
				}
		processed_data[tf] = tf_data

	return processed_data

def load_wcm(attempt_match):
	# type: (bool) -> Dict[str, Dict[str, Dict[str, float]]]
	"""
	Loads and extracts gene regulation data from the whole-cell model.

	Args:
		attempt_match: if True, handles data in way that was most likely the
			original processing, otherwise uses expected directionality

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
		if attempt_match:
			sign = 1
		else:
			sign = np.sign(float(line[5]))

		if tf not in data:
			data[tf] = {}

		data[tf][regulated_gene] = {
			'mean': sign * mean,
			'std': std,
			}

	return data

def load_new():
	# type: () -> Dict[str, Dict[str, Dict[str, float]]]
	"""
	Loads and extracts gene regulation data from new processing of the raw data.

	Returns:
		data: mean for each regulatory pair
			{TF gene name: {regulated gene: {'mean': mean, 'std': std}}}
	"""

	raw_data = load_file(NEW_FILE)[1:]

	data = {}
	for line in raw_data:
		# Columns of interest
		tf = line[0]
		regulated_gene = line[1]
		fc = float(line[2])

		if tf not in data:
			data[tf] = {}

		data[tf][regulated_gene] = {
			'mean': fc,
			'std': 0,
			}

	return data

def load_shifts():
	# type: () -> Dict[Tuple[str], Dict[str, int]]
	"""
	Load TF regulation information for each shift.

	Returns:
		shifts: gene regulation direction for each condition
			1 for active TF in condition1 vs condition2
			-1 for active TF in condition2 vs condition1
			{(condition1, condition2): {TF1: 1, TF2: -1, ...}}
	"""

	genes = {int(line[0]): line[2] for line in load_file(GENES_FILES)}
	shifts = {}
	for line in load_file(SHIFTS_FILE):
		condition = (line[1], line[2])

		gene_dir = {genes[np.abs(int(v))]: np.sign(int(v)) for v in line[12:26] if v != '0'}
		shifts[condition] = gene_dir

	return shifts

def compare_data(data1, data2, verbose=True):
	# type: (Dict[str, Dict[str, Dict[str, float]]], Dict[str, Dict[str, Dict[str, float]]], bool) -> None
	"""
	Compares regulation from source data to wcm data.

	Args:
		data1: mean and standard deviation for each regulatory pair
		data2: mean and standard deviation for each regulatory pair
		verbose: if True, prints additional regulation information

	Notes:
		dictionary structure for both inputs:
			{TF gene name: {regulated gene: {'mean': mean, 'std': std}}}

	TODO:
		check regulation not included in both datasets
	"""

	# Track statistics
	total_regulation = {}
	discrepancies = {}
	direction_discrepancies = {}

	# Print regulation that is different in the datasets
	if verbose:
		print('TF -> regulated gene: data1 vs data2')
	for tf, regulation in data1.items():
		total_regulation[tf] = len(regulation)
		discrepancies[tf] = 0
		direction_discrepancies[tf] = 0

		for gene, data in regulation.items():
			mean1 = data['mean']
			mean2 = data2.get(tf, {}).get(gene, {}).get('mean', 0)
			if np.abs(mean1 - mean2) > 0.01 and mean2 != 0:
				if np.sign(mean1) != np.sign(mean2):
					direction_discrepancies[tf] += 1
				discrepancies[tf] += 1
				if verbose:
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

		if verbose:
			print('{:5s}: {:3} {:3} {:3}'.format(tf, tf_dir_disc, tf_disc, total))
	print('Total: {:3} {:3} {:3}'.format(total_direction_discrepancies, total_discrepancies, total_interactions))

def parse_args():
	# type: () -> argparse.Namespace
	"""
	Parses arguments from the command line.

	Returns:
		values of variables parsed from the command line
	"""

	parser = argparse.ArgumentParser()

	parser.add_argument('-m', '--match',
		action='store_true',
		help='If set, processes source data to best match wcm.')
	parser.add_argument('-v', '--verbose',
		action='store_true',
		help='If set, prints more information.')

	return parser.parse_args()

if __name__ == '__main__':
	args = parse_args()

	src_data = load_src(args.match)
	wcm_data = load_wcm(args.match)
	new_data = load_new()

	compare_data(src_data, wcm_data, args.verbose)
	compare_data(wcm_data, new_data, args.verbose)
