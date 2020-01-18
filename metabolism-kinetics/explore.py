#! /usr/bin/env python

"""
Explore incorporating more kinetics data into the model.
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
RAW_DATA = os.path.join(DATA_DIR, 'raw_data.cp')
SIM_DATA = os.path.join(DATA_DIR, 'sim_data.cp')


def load_data():
	"""
	Load raw and sim data.

	Returns:
		raw_data: KnowledgeBaseEcoli object
		sim_data: SimulationData object
	"""

	with open(RAW_DATA) as f:
		raw_data = cPickle.load(f)
	with open(SIM_DATA) as f:
		sim_data = cPickle.load(f)

	return raw_data, sim_data

def extract_new(metabolism):
	"""
	Extract data of interest from new kinetics from raw_data.

	Args:
		metabolism (List[Dict[str, Any]]): each row of the new kinetics table
			as a list entry and columns as keys in a dict (from raw_data)

	Returns:
		rxns (set[str]): all reaction IDs
		mets (set[str]): all metabolite IDs
		enzs (set[str]): all enzyme IDs
	"""

	rxns = set()
	mets = set()
	enzs = set()
	for row in metabolism:
		rxn = row['reactionID']
		substrates = row['substrateIDs']
		enzymes = row['enzymeIDs']

		# Check type formatting
		if not isinstance(rxn, basestring):
			print('Invalid reaction: {}'.format(rxn))
		if isinstance(substrates, basestring):
			print('Invalid substrate: {}'.format(substrates))
			substrates = [substrates]
		if isinstance(enzymes, basestring):
			print('Invalid enzyme: {}'.format(enzymes))
			enzymes = [enzymes]

		rxns.update([row['reactionID']])
		mets.update(substrates)
		enzs.update(enzymes)

	return rxns, mets, enzs

def extract_sim(reactions, metabolism):
	"""
	Extract data of interest from sim_data.

	Args:
		reactions (List[Dict[str, Any]]): each row of the reactions table
			as a list entry and columns as keys in a dict (from raw_data)
		metabolism: Metabolism sim_data object

	Returns:
		raw_rxns (set[str]): reaction IDs from flat file
		all_rxns (set[str]): reaction IDs used in model, can be reverse or
			enzyme specific
		kinetics_rxns (set[str]): reaction IDs from all_rxns that currently
			have at least one kinetic constraint in the model
	"""

	raw_rxns = {row['reaction id'] for row in reactions}
	all_rxns = set(metabolism.reactionStoich.keys())
	kinetics_rxns = set(metabolism.reactionsToConstraintsDict)

	return raw_rxns, all_rxns, kinetics_rxns

if __name__ == '__main__':
	# Load data
	rd, sd = load_data()
	metabolism = sd.process.metabolism

	# Extract data of interest
	new_rxns, new_mets, new_enzs = extract_new(rd.metabolism_kinetics)
	raw_rxns, all_rxns, kinetics_rxns = extract_sim(rd.reactions, metabolism)

	# Compare data
	## Remove duplicates that are reverse or multiple enzyme kinetics reactions
	unique_kinetics = {rxn.strip(' (reverse)').split('__')[0] for rxn in kinetics_rxns}
	## Find kinetic reactions in the new data that are not in the current reactions
	unknown_rxns = {r for r in new_rxns if r not in raw_rxns}
	## Find kinetic reactions in the new data that have a partial match to current reactions
	## Indicates the need to add more information like specific molecules
	partial_match_rxns = {rxn for rxn in unknown_rxns if np.any([rxn in r for r in raw_rxns])}

	# Summarize comparisons
	print('Current model:')
	print('\t{} total reactions'.format(len(all_rxns)))
	print('\t{} reactions with at least one kinetic constraint'.format(len(kinetics_rxns)))
	print('\t{} unique reactions (either direction) with kinetics'.format(len(unique_kinetics)))

	print('New data:')
	print('\t{} unique reactions with kinetics in new kinetics'.format(len(new_rxns)))
	print('\t{} unknown reactions in raw_data'.format(len(unknown_rxns)))
	print('\t{} of unknowns with a partial match'.format(len(partial_match_rxns)))

	# TODO
	# compare enzymes
	# compare metabolites
	# print discrepancies (save to file)
	# get count of reactions that will be affected by change (kcat only or multiple kcat/km values)
	import ipdb; ipdb.set_trace()
