#! /usr/bin/env python

"""
Run FBA model from saved time points.
"""

from __future__ import absolute_import, division, print_function

import argparse
import cPickle
import os

import numpy as np
from typing import Any, Dict, List

from models.ecoli.processes.metabolism import FluxBalanceAnalysisModel
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.utils import units  # required for proper cPickle load


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'timepoints')
SIM_DATA_FILE = os.path.join(FILE_LOCATION, 'sim_data.cp')


def load_sim_data(path=SIM_DATA_FILE):
	# type: (str) -> SimulationDataEcoli
	"""
	Load simulation data.

	Args:
		path: path to sim_data cPickle file to load

	Returns:
		sim_data: simulation data
	"""

	with open(path) as f:
		sim_data = cPickle.load(f)
	return sim_data

def load_timepoints():
	# type: () -> List[Dict[str, Any]]
	"""
	Load timepoint args for each FBA function from files.

	Returns:
		all_timepoints: args for FBA functions for each timepoint,
			dict keys are the same as the function they belong to
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

def setup_model(sim_data, timepoint):
	# type: (SimulationDataEcoli, Dict[str, Any]) -> FluxBalanceAnalysisModel
	"""
	Setup FBA model with loaded data.

	Args:
		sim_data: simulation data
		timepoint: args for FBA functions for a timepoint,
			dict keys are the same as the function they belong to

	Returns:
		model: FBA model initialized for the given timepoint
	"""

	model = FluxBalanceAnalysisModel(sim_data)
	model.set_molecule_levels(*timepoint['set_molecule_levels'])
	model.set_reaction_bounds(*timepoint['set_reaction_bounds'])
	model.set_reaction_targets(*timepoint['set_reaction_targets'])

	return model

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

	# Extract data from sim_data
	kinetic_constraint_reactions = set(sim_data.process.metabolism.kinetic_constraint_reactions)
	all_reactions = sorted(sim_data.process.metabolism.reactionStoich.keys())
	sim_data.process.metabolism.flux_regularization = 7.008e-4

	# Unmodified model
	original_model = setup_model(sim_data, timepoint)

	# Extract data from original model
	mol_id = 'CPD-8260[c]'
	mol_idx = original_model.fba.getOutputMoleculeIDs().index(mol_id)
	original_mol_change = original_model.fba.getOutputMoleculeLevelsChange()[mol_idx]

	# Iterate model for desired outcomes
	for rxn in all_reactions:
		if rxn in kinetic_constraint_reactions:
			continue

		# Create model
		model = setup_model(sim_data, timepoint)

		# Iteration specific modifications
		model.fba._solver.setFlowObjectiveCoeff(rxn, 0)

		# TODO: setup loop to modify model to reach a desired output

		# Solve model and check outputs
		mol_change = model.fba.getOutputMoleculeLevelsChange()[mol_idx]
		if mol_change > original_mol_change:
			print('{}: {} change: {:.3f}'.format(rxn, mol_id, mol_change))
