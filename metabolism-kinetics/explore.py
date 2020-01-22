#! /usr/bin/env python

"""
Explore incorporating more kinetics data into the model.
"""

from __future__ import absolute_import, division, print_function

import cPickle
import csv
import os
import subprocess
import time

import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

# Input data
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
RAW_DATA = os.path.join(DATA_DIR, 'raw_data.cp')
SIM_DATA = os.path.join(DATA_DIR, 'sim_data.cp')

# Output files
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUTPUT_DIR):
	os.makedirs(OUTPUT_DIR)
REACTION_FILE = os.path.join(OUTPUT_DIR, 'reactions.tsv')
METABOLITE_FILE = os.path.join(OUTPUT_DIR, 'metabolites.tsv')
ENZYME_FILE = os.path.join(OUTPUT_DIR, 'enzymes.tsv')


def get_git_hash():
	"""
	Get short git hash from the command line.

	Returns:
		git_hash (str): git hash of user repo
	"""

	cwd = os.getcwd()
	os.chdir(FILE_LOCATION)  # Ensure running git command from proper repo if nested
	req = subprocess.Popen(['git', 'rev-parse', 'HEAD'], stdout=subprocess.PIPE)
	req.wait()
	git_hash = req.communicate()[0].strip()[:10]
	os.chdir(cwd)

	return git_hash

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
		rxn_to_met (dict[str, list[str]]): map of each reaction to metabolites
		rxn_to_enz (dict[str, str]): map of each reaction to enzyme catalyst
	"""

	rxns = set()
	mets = set()
	enzs = set()
	rxn_to_met = {}
	rxn_to_enz = {}
	for row in metabolism:
		rxn = row['reactionID']
		substrates = row['substrateIDs']
		enzyme = row['enzymeIDs']

		# Check type formatting
		if not isinstance(rxn, basestring):
			print('Invalid reaction: {}'.format(rxn))
		if isinstance(substrates, basestring):
			print('Invalid substrate: {}'.format(substrates))
			substrates = [substrates]
		if isinstance(enzyme, basestring):
			print('str not list for enzyme: {}'.format(enzyme))
		else:
			if len(enzyme) > 1:
				print('Multiple enzymes, only using first: {}'.format(enzyme))
			enzyme = enzyme[0]

		rxn_to_met[rxn] = substrates
		rxn_to_enz[rxn] = enzyme

		rxns.update([row['reactionID']])
		mets.update(substrates)
		enzs.update([enzyme])

	return rxns, mets, enzs, rxn_to_met, rxn_to_enz

def extract_sim(reactions, sim_data):
	"""
	Extract data of interest from sim_data.

	Args:
		reactions (List[Dict[str, Any]]): each row of the reactions table
			as a list entry and columns as keys in a dict (from raw_data)
		sim_data: SimulationDataEcoli object

	Returns:
		raw_rxns (set[str]): reaction IDs from flat file
		all_rxns (set[str]): reaction IDs used in model, can be reverse or
			enzyme specific
		kinetics_rxns (set[str]): reaction IDs from all_rxns that currently
			have at least one kinetic constraint in the model
		mols (set[str]): possible molecule IDs (enzymes, metabolites, etc)
			with no location tag
		mets (set[str]): metabolites linked to a reaction with no location tag
		enzs (set[str]): enzymes linked to a reaction with no location tag
		met_to_rxn (dict[str, list[str]]): map of each metabolites to reactions
			that contain it
		enz_to_rxn (dict[str, list[str]]): map of each enzyme to reactions it
			catalyzes
	"""

	metabolism = sim_data.process.metabolism

	# Raw reactions data
	raw_rxns = set()
	mets = set()
	enzs = set()
	met_to_rxn = {}
	enz_to_rxn = {}
	for row in reactions:
		rxn = row['reaction id']
		stoich = [m[:-3] for m in row['stoichiometry']]
		enz = row['catalyzed by']

		raw_rxns.add(rxn)
		mets.update([m for m in stoich])

		for m in stoich:
			met_to_rxn[m] = met_to_rxn.get(m, []) + [rxn]

		for e in enz:
			enzs.add(e)
			enz_to_rxn[e] = enz_to_rxn.get(e, []) + [rxn]

	# Simulation reactions
	all_rxns = set(metabolism.reactionStoich.keys())
	kinetics_rxns = set(metabolism.reactionsToConstraintsDict)

	# Simulation molecules (enzymes and metabolites)
	mols = set(sim_data.getter._all_mass.keys())

	return raw_rxns, all_rxns, kinetics_rxns, mols, mets, enzs, met_to_rxn, enz_to_rxn

def create_valid_constraints(raw_data, unknown_rxns, unknown_mets, unknown_enzs, stringent_matches):
	"""

	TODO:
		- handle enzyme location
		- handle NAD vs NADP etc double reactions
		- handle enzyme used for multiple reactions (need to adjust kcat?)
	"""

	constraints = []
	for row in raw_data.metabolism_kinetics:
		# Get relevant entries from row
		rxn = row['reactionID']
		mets = row['substrateIDs']
		enz = row['enzymeIDs']
		kcat = row['kcat']
		kms = row['kM']
		kis = row['kI']
		direction = row['direction']  # TODO: handle
		constraint_type = row['rateEquationType']

		# Handle constraint type
		if constraint_type == 'custom':
			# TODO: handle
			pass

		# Handle reaction IDs
		if rxn in unknown_rxns:
			if rxn not in stringent_matches:
				continue

			rxn = stringent_matches[rxn]
			if len(rxn) > 1:
				print('Invalid reaction: {}'.format(rxn))
			rxn = rxn[0]

		# Handle enzyme IDs
		if len(enz) != 1:
			print('Invalid enzyme: {}: {}'.format(rxn, enz))
			continue
		enz = enz[0]
		if enz in unknown_enzs:
			continue

		# Handle metabolite IDs for parameters
		ks = []
		## KM parameters
		for met, k in zip(mets, kms):
			if met in unknown_mets:
				continue
			ks.append(1. * k)

		## KI parameters
		for met, k in zip(mets[len(kms):], kis):
			if met in unknown_mets:
				continue
			ks.append(1. * k)

		# Should only be one kcat
		# TODO: handle empty if kcat in 'customParameters' with value in 'customParameterConstantValues'
		# TODO: split if multiple kcats
		# TODO: remove kcat if km not matched to metabolite above
		if len(kcat) != 1:
			print('Invalid kcat: {}: {}'.format(rxn, kcat))
			continue
		kcat = kcat[0]

		new_constraint = {
			'reaction': rxn,
			'enzyme': enz,
			'kcat': 1. * kcat,
			'k': ks,
			}
		constraints.append(new_constraint)

	return constraints


if __name__ == '__main__':
	# Load data
	raw_data, sim_data = load_data()

	# For easy troubleshooting
	rd = raw_data
	sd = sim_data
	metabolism = sd.process.metabolism
	translation = sd.process.translation

	# Extract data of interest
	(new_rxns, new_mets, new_enzs, rxn_to_met, rxn_to_enz
		) = extract_new(raw_data.metabolism_kinetics)
	(raw_rxns, all_rxns, kinetics_rxns, all_mols, raw_mets, raw_enzs, met_to_rxn, enz_to_rxn
		) = extract_sim(raw_data.reactions, sim_data)

	# Compare data
	## Remove duplicates that are reverse or multiple enzyme kinetics reactions
	unique_kinetics = {rxn.strip(' (reverse)').split('__')[0] for rxn in kinetics_rxns}
	## Find current metabolites that are not represented in the wcm
	unknown_current_mets = {met for met in raw_mets if met not in all_mols}
	## Find current enzymes that are not represented in the wcm
	unknown_current_enzs = {enz for enz in raw_enzs if enz not in all_mols}
	## Find kinetic reactions in the new data that are not in the current reactions
	unknown_rxns = {r for r in new_rxns if r not in raw_rxns}
	## Find kinetic reactions in the new data that have a partial match to current reactions
	## Indicates the need to add more information like specific molecules
	partial_match_rxns = {rxn for rxn in unknown_rxns if np.any([rxn in r for r in raw_rxns])}
	## Find partial reaction matches that also match metabolites for more stringent match
	stringent_matches = {}
	for rxn in sorted(unknown_rxns):
		mets = rxn_to_met[rxn]
		matches = {r for r in raw_rxns if rxn in r}
		for met in mets:
			matches = matches.intersection(met_to_rxn.get(met, []))

		if matches:
			stringent_matches[rxn] = sorted(matches)
	## Find metabolites that are not represented in the current wcm
	unknown_mets = {met for met in new_mets if met not in all_mols and met.upper() not in all_mols}
	## Find enzymes that are not represented in the current wcm
	unknown_enzs = {enz for enz in new_enzs if enz not in all_mols}


	# Summarize comparisons
	print('\nCurrent model:')
	print('\t{} total reactions'.format(len(all_rxns)))
	print('\t{} reactions with at least one kinetic constraint'.format(len(kinetics_rxns)))
	print('\t{} unique reactions (either direction) with kinetics'.format(len(unique_kinetics)))
	print('\t{}/{} current metabolites that are not in molecules'.format(len(unknown_current_mets), len(raw_mets)))
	print('\t{}/{} current enzymes that are not in molecules'.format(len(unknown_current_enzs), len(raw_enzs)))

	print('\nNew data:')
	print('\t{} unique reactions with kinetics in new kinetics'.format(len(new_rxns)))
	print('\t{} unknown reactions in new kinetics'.format(len(unknown_rxns)))
	print('\t{} of unknowns with a partial match'.format(len(partial_match_rxns)))
	print('\t{} of unknowns with a partial match with metabolites'.format(len(stringent_matches)))
	print('\t{}/{} unknown metabolites in new kinetics'.format(len(unknown_mets), len(new_mets)))
	print('\t{}/{} unknown enzymes in new kinetics'.format(len(unknown_enzs), len(new_enzs)))
	print('')

	# Print details of discrepancies
	## Get autogenerated header
	metadata = 'Generated by {} on git commit {} at {}'.format(
		__file__, get_git_hash(), time.ctime()
		)

	## Reactions
	with open(REACTION_FILE, 'w') as f:
		print('Writing data to {}'.format(f.name))
		writer = csv.writer(f, delimiter='\t')
		writer.writerow([metadata])

		writer.writerow(['Unknown Reaction ID', 'Possible Reaction Match',
			'Metabolites', 'Enzyme', 'Enzyme in model',
			'Reactions catalyzed by enzyme'])
		first = []
		second = []
		third = []
		for rxn in sorted(unknown_rxns):
			matches = stringent_matches.get(rxn, [])
			partial_match = np.any([rxn in r for r in raw_rxns])
			mets = rxn_to_met[rxn]
			enz = rxn_to_enz[rxn]
			enz_in_model = enz in all_mols
			rxn_with_enz = enz_to_rxn.get(enz, [])

			# Sort reactions that don't have a match first
			if matches:
				group = third
			elif partial_match:
				group = second
			else:
				group = first

			group.append([rxn, matches, partial_match, mets, enz, enz_in_model, rxn_with_enz])

		writer.writerows(first)
		writer.writerows(second)
		writer.writerows(third)

	## Metabolites
	with open(METABOLITE_FILE, 'w') as f:
		print('Writing data to {}'.format(f.name))
		writer = csv.writer(f, delimiter='\t')
		writer.writerow([metadata])

		writer.writerow(['Unknown Metabolite ID', 'Metabolite in reactions'])
		for met in sorted(unknown_mets):
			met_in_rxns = met in raw_mets
			writer.writerow([met, met_in_rxns])

	## Enzymes
	with open(ENZYME_FILE, 'w') as f:
		print('Writing data to {}'.format(f.name))
		writer = csv.writer(f, delimiter='\t')
		writer.writerow([metadata])

		writer.writerow(['Unknown Enzyme ID', 'Enzyme in reactions'])
		for enz in sorted(unknown_enzs):
			enz_in_rxns = enz in raw_enzs
			writer.writerow([enz, enz_in_rxns])

	new_constraints = create_valid_constraints(raw_data, unknown_rxns,
		unknown_mets, unknown_enzs, stringent_matches)
	old_constraints = None  # TODO: pull from sim_data

	# TODO
	# get count of reactions that will be affected by change (kcat only or multiple kcat/km values)
	# function to compare new constraints to current constraints
	import ipdb; ipdb.set_trace()
