#! /usr/bin/env python

"""
Run FBA model from saved time points.
"""

from __future__ import absolute_import, division, print_function

import argparse
import cPickle
import os

from models.ecoli.processes.metabolism import FluxBalanceAnalysisModel
from wholecell.utils import units  # required for proper cPickle load


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'timepoints')
SIM_DATA_FILE = os.path.join(FILE_LOCATION, 'sim_data.cp')


def load_sim_data(path=SIM_DATA_FILE):
	"""

	"""

	with open(path) as f:
		sim_data = cPickle.load(f)
	return sim_data

def load_timepoints():
	"""

	"""

	all_timepoints = []
	for timepoint in os.listdir(DATA_DIR):
		data = {}
		timepoint_dir = os.path.join(DATA_DIR, timepoint)
		for filename in os.listdir(timepoint_dir):
			func = os.path.basename(filename).split('.')[0]
			with open(os.path.join(timepoint_dir, filename), 'rb') as f:
				data[func] = cPickle.load(f)

		all_timepoints.append(data)

	return all_timepoints

def parse_args():
	# type: () -> argparse.Namespace
	"""
	Parses arguments from the command line.

	Returns:
		values of variables parsed from the command line
	"""

	parser = argparse.ArgumentParser()

	parser.add_argument('-t', '--timepoint',
		type=int,
		default=0,
		help='Timepoint to analyze. Max value depends on number of timepoints saved in {}.'.format(DATA_DIR))
	parser.add_argument('-s', '--sim-data',
		default=SIM_DATA_FILE,
		help='Path to sim_data cPickle (default: {}).'.format(SIM_DATA_FILE))

	return parser.parse_args()

if __name__ == '__main__':
	args = parse_args()

	# Load data from files
	sim_data = load_sim_data(path=args.sim_data)
	timepoints = load_timepoints()
	timepoint = timepoints[args.timepoint]  # TODO: set up loop for all timepoints if desired

	# Create model and extract relevant data
	model = FluxBalanceAnalysisModel(sim_data)
	water_idx = model.fba.getOutputMoleculeIDs().index('WATER[c]')

	# TODO: setup loop to modify model to reach a desired output

	# Setup model for timepoint
	model.set_molecule_levels(*timepoint['set_molecule_levels'])
	model.set_reaction_bounds(*timepoint['set_reaction_bounds'])
	model.set_reaction_targets(*timepoint['set_reaction_targets'])

	# Solve model and check outputs
	print('Water change: {:.3f}'.format(model.fba.getOutputMoleculeLevelsChange()[water_idx]))
