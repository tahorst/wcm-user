#! /usr/bin/env python

"""
Run FBA model from saved time points.
"""

from __future__ import absolute_import, division, print_function

import argparse
import cPickle
import os

import numpy as np
from typing import Any, Dict, List, Tuple

from models.ecoli.processes.metabolism import FluxBalanceAnalysisModel
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.utils import units  # required for proper cPickle load


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'timepoints')
SIM_DATA_FILE = os.path.join(FILE_LOCATION, 'sim_data.cp')

ALL_ANALYSIS_OPTIONS = [
	'adjust-parameter',
	'increase-molecule',
	]


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

def extract_sim_data_args(analysis_type, path):
	# type: (str, str) -> (SimulationDataEcoli, Dict[str, Any])
	"""
	Load simulation data, modify any sim_data attributes and extract args for
	the desired analysis.

	Args:
		analysis_type: label from ALL_ANALYSIS_OPTIONS to indicate which
			analysis to perform
		path: path to sim_data cPickle file to load

	Returns:
		sim_data: modified simulation data
		kwargs: sim_data args to use in model analysis

	TODO:
		create new functions for each analysis type or make a class?
	"""

	# Load from file
	with open(path) as f:
		sim_data = cPickle.load(f)

	# Handle each analysis type
	kwargs = {}
	if analysis_type == 'adjust-parameter':
		kwargs = {
			'factors': np.logspace(0, 4, 100),
			'sim_data_attr': 'process.metabolism.flux_regularization',
			}
	elif analysis_type == 'increase-molecule':
		sim_data.process.metabolism.flux_regularization = 7.008e-4
		kwargs = {
			'kinetic_constraint_reactions': set(sim_data.process.metabolism.kinetic_constraint_reactions),
			'all_reactions': sorted(sim_data.process.metabolism.reactionStoich.keys()),
			}

	return sim_data, kwargs

def extract_model_args(analysis_type, model):
	# type: (str, FluxBalanceAnalysisModel) -> Dict[str, Any]
	"""
	Extract data from the FBA model for the desired analysis.

	Args:
		analysis_type: label from ALL_ANALYSIS_OPTIONS to indicate which
			analysis to perform
		model: FBA model

	Returns:
		kwargs: model args to use in model analysis
	"""

	kwargs = {}
	if analysis_type == 'adjust-parameter':
		mol_id = 'CPD-8260[c]'
		mol_idx = model.fba.getOutputMoleculeIDs().index(mol_id)
		kwargs = {
			'mol_id': mol_id,
			'mol_idx': mol_idx,
			}
	elif analysis_type == 'increase-molecule':
		mol_id = 'CPD-8260[c]'
		mol_idx = model.fba.getOutputMoleculeIDs().index(mol_id)
		original_mol_change = model.fba.getOutputMoleculeLevelsChange()[mol_idx]
		kwargs = {
			'mol_id': mol_id,
			'mol_idx': mol_idx,
			'original_mol_change': original_mol_change,
			}

	return kwargs

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

def solve_model(analysis_type, sim_data, timepoint, **kwargs):
	# type: (str, SimulationDataEcoli, Dict[str, Any], **Any) -> None
	"""
	Solves the FBA model for the desired analysis.

	Args:
		analysis_type: label from ALL_ANALYSIS_OPTIONS to indicate which
			analysis to perform
		sim_data: simulation data
		timepoint: args for FBA functions for a timepoint,
			dict keys are the same as the function they belong to

	TODO:
		setup loop to modify model to reach a desired output
	"""

	if analysis_type == 'adjust-parameter':
		solve_adjust_parameter(sim_data, timepoint, **kwargs)
	elif analysis_type == 'increase-molecule':
		solve_increase_molecule(sim_data, timepoint, **kwargs)

def solve_adjust_parameter(sim_data, timepoint,
		factors, sim_data_attr,
		mol_id, mol_idx):
	"""
	Solves the model for the 'adjust-parameter' option.  This adjusts the given
	sim_data attribute by different factors to see the effect on the FBA
	output molecule of interest.

	Args:
		sim_data (SimulationDataEcoli): simulation data
		timepoint (Dict[str, Any]): args for FBA functions for a timepoint,
			dict keys are the same as the function they belong to
		factors (Iterable[float]): factors to multiply the original attribute
			value by
		sim_data_attr (str): sim_data attribute to adjust, separated by '.' if
			nested attributes
		mol_id (str): molecule ID to check for an increase in the change
			for each iteration
		mol_idx (int): index of output molecules for mol_id
	"""

	# Get attribute to adjust
	attrs = sim_data_attr.split('.')
	parent_attr = sim_data
	attr = attrs[-1]
	if len(attrs) > 1:
		for a in attrs[:-1]:
			parent_attr = getattr(parent_attr, a)
	original_value = getattr(parent_attr, attr)

	# Adjust attribute by each factor
	for factor in factors:
		# Iteration specific modifications
		value = original_value * factor
		setattr(parent_attr, attr, value)

		# Create model
		model = setup_model(sim_data, timepoint)

		# Solve model and check outputs
		mol_change = model.fba.getOutputMoleculeLevelsChange()[mol_idx]
		print('{}={:.3e}: {} change: {:.3f}'.format(attr, value, mol_id, mol_change))

def solve_increase_molecule(sim_data, timepoint,
		kinetic_constraint_reactions, all_reactions,
		mol_id, mol_idx, original_mol_change):
	"""
	Solves the model for the 'increase-molecule' option.  This finds the
	reactions that cause an increase in the change in concentration for a
	limited metabolite when the flux regularization for the reaction is
	disabled.

	Args:
		sim_data (SimulationDataEcoli): simulation data
		timepoint (Dict[str, Any]): args for FBA functions for a timepoint,
			dict keys are the same as the function they belong to
		kinetic_constraint_reactions (Set[str]): reactions that have a
			kinetic constraint (already have no flux regularization)
		all_reactions (List[str]): all reaction IDs to iterate over
		mol_id (str): molecule ID to check for an increase in the change
			for each iteration
		mol_idx (int): index of output molecules for mol_id
		original_mol_change (float): change in mol_id for the original model
	"""

	for rxn in all_reactions:
		if rxn in kinetic_constraint_reactions:
			continue

		# Create model
		model = setup_model(sim_data, timepoint)

		# Iteration specific modifications
		model.fba._solver.setFlowObjectiveCoeff(rxn, 0)

		# Solve model and check outputs
		mol_change = model.fba.getOutputMoleculeLevelsChange()[mol_idx]
		if mol_change > original_mol_change:
			print('{}: {} change: {:.3f}'.format(rxn, mol_id, mol_change))

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
	parser.add_argument('-a', '--analysis',
		default=ALL_ANALYSIS_OPTIONS[0],
		help='Analysis type to perform. Possible values: {}'.format(ALL_ANALYSIS_OPTIONS))

	return parser.parse_args()

if __name__ == '__main__':
	args = parse_args()

	# Load data from files
	sim_data, sim_data_args = extract_sim_data_args(args.analysis, args.sim_data)
	timepoints = load_timepoints()
	timepoint = timepoints[args.timepoint]  # TODO: set up loop for all timepoints if desired

	# Unmodified model
	original_model = setup_model(sim_data, timepoint)
	model_args = extract_model_args(args.analysis, original_model)

	# Iterate model for desired outcomes
	solver_args = dict(sim_data_args, **model_args)
	solve_model(args.analysis, sim_data, timepoint, **solver_args)
