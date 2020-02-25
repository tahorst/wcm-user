#! /usr/bin/env python

"""
Find valid reactions and metabolites in reaction network.
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from typing import Dict, List, Set

from reconstruction.ecoli.dataclasses.process.metabolism import Metabolism


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
SIM_DATA_FILE = os.path.join(FILE_LOCATION, 'sim_data.cp')


def get_boundaries(metabolism, media='minimal'):
	# type: (Metabolism, str) -> (Set[str], Set[str])
	"""
	Get source and sink metabolites. Imports are sources. Concentration targets
	and secretions are sinks.

	Args:
		metabolism: sim_data metabolism process class
		media: media label for potential import molecules

	Returns:
		sources: set of metabolite IDs with location tag that are sources
		sinks: set of metabolite IDs with location tag that are sinks
	"""

	sources = set(metabolism.boundary.exchange_data_dict['importExchangeMolecules'][media])
	sinks = set(metabolism.boundary.secretion_exchange_molecules)
	sinks.update(metabolism.concDict)

	return sources, sinks

def get_mappings(reactions, map_reactants):
	# type: (Dict[str, Dict[str, int]], bool) -> (Dict[str, List[str]], Dict[str, List[str]])
	"""
	Mappings of metabolites to reactions and reactions to metabolites for easy
	parsing of all reactions.

	Args:
		reactions: sim_data stoichiometry structure
		map_reactants: True if mapping for reactants, False if mapping for products

	Returns:
		met_to_rxn: map each metabolite (reactant or product) to reactions that
			contain it
		rxn_to_met: map each reaction to metabolites that are either reactants
			or products
	"""

	if map_reactants:
		direction = -1
	else:
		direction = 1

	met_to_rxn = {}
	rxn_to_met = {}
	for rxn, stoich in reactions.items():
		for met, factor in stoich.items():
			if factor * direction > 0:
				met_to_rxn[met] = met_to_rxn.get(met, []) + [rxn]
				rxn_to_met[rxn] = rxn_to_met.get(rxn, []) + [met]

	return met_to_rxn, rxn_to_met

def trace_possible_reactions(start, met_to_rxn, rxn_to_met, excluded_rxns):
	# type: (Set[str], Dict[str, List[str]], Dict[str, List[str]], Set[str]) -> Set[str]
	"""
	Iteratively trace metabolites through valid reactions.

	Args:
		start: starting set of metabolites (sources or sinks)
		met_to_rxn: mapping of metabolites to reactions
		rxn_to_met: mapping of reactions to metabolites
		excluded_rxns: invalid reactions that will not link to new metabolites

	Returns:
		mets: set of metabolites that can be traced from start through valid
			reactions (source -> product or reactant -> sink)
	"""

	mets = set()
	rxns = set(excluded_rxns)

	new_mets = set(start)
	while new_mets:
		# Get reactions that have a source metabolite as an input
		possible_rxns = set()
		for met in new_mets:
			possible_rxns.update({r for r in met_to_rxn.get(met, [])})
		new_rxns = possible_rxns.difference(rxns)

		# Get products from reactions with sources
		possible_mets = set()
		for rxn in new_rxns:
			possible_mets.update({m for m in rxn_to_met.get(rxn, [])})

		# Update metabolites for the next round
		mets.update(new_mets)
		new_mets = possible_mets.difference(mets)

	return mets

def prune_reactions(reactions, valid_metabolites):
	# type: (Dict[str, Dict[str, int]], Set[str]) -> Set[str]
	"""
	Identify reactions that are not possible because mass balance must exist.

	Args:
		reactions: sim_data stoichiometry structure
		valid_metabolites: metabolites that have a path to a source and sink
			so that flux can flow through

	Returns:
		excluded_rxns: reactions that are not valid because one or more
			metabolites can not have flux
	"""

	excluded_rxns = set()
	for rxn, stoich in reactions.items():
		if not all([m in valid_metabolites for m in stoich]):
			excluded_rxns.add(rxn)

	return excluded_rxns


if __name__ == '__main__':
	# Load data to analyze
	with open(SIM_DATA_FILE) as f:
		sim_data = cPickle.load(f)
	metabolism = sim_data.process.metabolism
	reactions = metabolism.reactionStoich

	# Extract necessary data
	sources, sinks = get_boundaries(metabolism)
	reactant_to_rxn, rxn_to_reactant = get_mappings(reactions, True)
	product_to_rxn, rxn_to_product = get_mappings(reactions, False)

	# Find reactions that are not possible due to metabolites not being balanced
	excluded_rxns = set()
	rxn_len = -1  # used for break condition checking
	while len(excluded_rxns) != rxn_len:
		# Store old length to find when no new reactions are found
		rxn_len = len(excluded_rxns)

		# Find metabolites that are source -> product or reactant -> sink
		potential_products = trace_possible_reactions(
			sources, reactant_to_rxn, rxn_to_product, excluded_rxns)
		potential_reactants = trace_possible_reactions(
			sinks, product_to_rxn, rxn_to_reactant, excluded_rxns)

		# Only metabolites with path from a source and to a sink can carry flux and be balanced
		valid_mets = potential_products.intersection(potential_reactants)

		# Exclude reactions that have invalid metabolites
		excluded_rxns = prune_reactions(reactions, valid_mets)

	# Analyze results
	all_mets = set(reactant_to_rxn.keys() + product_to_rxn.keys())
	all_rxns = set(rxn_to_reactant.keys() + rxn_to_product.keys())
	excluded_metabolites = all_mets.difference(valid_mets)
	valid_rxns = all_rxns.difference(excluded_rxns)

	# Print summary
	print('{}/{} metabolites are valid'.format(len(valid_mets), len(all_mets)))
	print('{}/{} reactions are valid'.format(len(valid_rxns), len(all_rxns)))
