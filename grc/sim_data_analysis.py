#! /usr/bin/env python
'''
Works with output from the fitter to determine parameters and feasibility
of using ppGpp for growth rate control.

Requires (cPickle files, defaults to .cp files in data directory):
	sim_data (SimulationData object): output from fitter and knowledgebase
		for a simulation, saved after running fitter
	cell_specs (dict): information about each condition that was fit, not saved
		after running fitter but can save during execution

Output:
	For the --sgd option, an output tsv will be saved by default in out/sgd.tsv
	For the --plot-synthetases option, an output plot will be saved by default in out/synthetases.png
	For the --plot-parameters option, an output plot will be saved by default in out/<input file>.png
'''

from __future__ import division

import argparse
import cPickle
import csv
import os
import time

import matplotlib.pyplot as plt
import numpy as np
import plotly
import plotly.graph_objs as go

from wholecell.utils import units


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
if not os.path.exists(DATA_DIR):
	os.makedirs(DATA_DIR)
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUTPUT_DIR):
	os.makedirs(OUTPUT_DIR)

SGD_OUT = 'sgd.tsv'
SYNTHETASE_PLOT_OUT = 'synthetases.png'

MICROMOLAR_UNITS = units.umol / units.L
PARAMS = [
	'synthetase_charging_rate',
	'Km_synthetase_amino_acid',
	'Km_synthetase_uncharged_trna',
	'Kdissociation_charged_trna_ribosome',
	'Kdissociation_uncharged_trna_ribosome',
	'k_RelA_ppGpp_synthesis',
	'KD_RelA_ribosome',
	'k_SpoT_ppGpp_synthesis',
	'k_SpoT_ppGpp_degradation',
	'rrn_vmax',
	'KI_ppgpp_rnap',
	'KM_rrn_rnap',
	]

# Concentrations for each AA from proteomics from Schmidt et al. 2015
SCHMIDT_CONC = {
	'basal': np.array([
		0.570, 0.550, 1.616, 0.664, 0.322, 1.497, 0.584, 0.485, 0.552, 1.267, 1.080,
		0.894, 0.609, 1.469, 1.264, 0.739, 0.275, 0.531, 0.666, 0.011, 1.005
		]),
	'with_aa': np.array([
		0.631, 0.564, 1.127, 0.844, 0.332, 1.680, 0.663, 0.713, 0.551, 1.118, 1.043,
		1.396, 0.498, 1.275, 1.248, 1.125, 0.375, 0.387, 0.651, 0.005, 1.204
		]),
	'no_oxygen': np.array([
		0.661, 0.422, 1.655, 0.890, 0.246, 1.201, 0.598, 0.461, 0.445, 1.258, 1.250,
		1.035, 0.562, 1.532, 1.371, 0.783, 0.522, 0.407, 0.585, 0.015, 0.948
		]),
	}


def output_location(arg, parent, default):
	'''
	Assembles the path to the output file from command line argument.

	Args:
		arg (str): filename parsed from command line, if None uses default
		parent (str): path to the parent directory
		default (str): filename to use if arg is None
	'''

	if arg is None:
		arg = default

	return os.path.join(parent, arg)

def add_parameter_noise(constants):
	'''
	Adds noise to parameters.
	Uses gamma distribution to prevent negative values with mean of 1 and variance of 0.1.

	Args:
		constants (class): constants from sim_data
	'''

	for param in PARAMS:
		setattr(constants, param, getattr(constants, param) * np.random.gamma(10, 0.1))

def get_volume(sim_data, bulk_container):
	'''
	Calculates the volume of a cell for a given composition assuming constant
	cell density.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule

	Returns:
		float: volume of the cell (in units of L)

	Note:
		Can be time intensive so get_counts_to_micromolar may be a better option.
	'''

	mol_names = list(bulk_container.objectNames())
	mol_masses = sim_data.getter.getMass(mol_names)
	mol_counts = bulk_container.counts(mol_names)

	mass = mol_counts.dot(mol_masses) / sim_data.constants.nAvogadro
	volume = mass / sim_data.constants.cellDensity

	return volume.asNumber(units.L)

def get_counts_to_micromolar(doubling_time):
	'''
	Estimate of cell size in different conditions to convert counts to concentration.

	Args:
		doubling_time (float with time units): doubling time of the condition
			to get the conversion factor for

	Returns:
		float: conversion from counts of a molecule to uM

	Note:
		More accurate to use get_volume and convert with Avogadro's number
		but it can be time intensive.
	'''

	if doubling_time.asNumber() > 50:
		conversion = 0.0040315516
	elif doubling_time.asNumber() < 40:
		conversion = 0.0006922573
	else:
		conversion = 0.0015064190

	return conversion

def get_rrna_counts(sim_data, doubling_time, rnap_activation_rate, synth_prob):
	'''
	Gets the counts for rRNA based on RNAP activation rate and synthesis probability
	of rRNA, which can be calculated from ppGpp concentration for a feedback loop.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		doubling_time (float with time units): expected cell doubling time
		rnap_activation_rate (float): number of RNAP activations/terminations per second
		synth_prob (float): synthesis probability for rRNA

	Returns:
		float: counts of rRNA at steady state
	'''

	rna_data = sim_data.process.transcription.rnaData
	rrna_deg_rate = np.mean(rna_data['degRate'][rna_data['isRRna']]).asNumber(1 / units.s)  # currently all equal
	dilution_rate = np.log(2) / doubling_time.asNumber(units.s)

	counts = rnap_activation_rate * synth_prob / (rrna_deg_rate + dilution_rate)

	return counts

def get_active_ribosome_counts(sim_data, bulk_container, doubling_time, rnap_activation_rate,
		rrna_synth_prob, ribosome_control):
	'''
	Gets the counts of active ribosomes based on fraction active and subunits.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule
		doubling_time (float with time units): expected cell doubling time
		rnap_activation_rate (float): number of RNAP activations/terminations per second
		rrna_synth_prob (float): synthesis probability for rRNA
		ribosome_control (bool): if True, updates ribosome concentration based on
			ppGpp regulation, otherwise uses bulk container counts

	Returns:
		float: number of active ribosomes
	'''

	active_fraction = sim_data.growthRateParameters.getFractionActiveRibosome(doubling_time)

	if ribosome_control:
		count_rrna = get_rrna_counts(sim_data, doubling_time, rnap_activation_rate, rrna_synth_prob)
		ribosome_counts = count_rrna
	else:
		molecule_ids = sim_data.moleculeIds
		count_30s = bulk_container.count(molecule_ids.s30_fullComplex)
		count_50s = bulk_container.count(molecule_ids.s50_fullComplex)
		ribosome_counts = min(count_30s, count_50s)

	return ribosome_counts * active_fraction

def get_free_rnap_counts(sim_data, bulk_container, doubling_time):
	'''
	Gets the counts of free RNAP based on fraction active.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule
		doubling_time (float with time units): expected cell doubling time

	Returns:
		float: number of free RNAP
	'''

	active_fraction = sim_data.growthRateParameters.getFractionActiveRnap(doubling_time)
	counts = bulk_container.count(sim_data.moleculeIds.rnapFull)

	free_counts = (1 - active_fraction) * counts

	return free_counts

def get_bound_rnap_counts(sim_data, bulk_container, doubling_time):
	'''
	Gets the counts of bound RNAP based on fraction active.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule
		doubling_time (float with time units): expected cell doubling time

	Returns:
		float: number of bound RNAP
	'''

	active_fraction = sim_data.growthRateParameters.getFractionActiveRnap(doubling_time)
	counts = bulk_container.count(sim_data.moleculeIds.rnapFull)

	bound_counts = active_fraction * counts

	return bound_counts

def get_rnap_activation(sim_data, bulk_container, doubling_time, synth_prob):
	'''
	Gets the number of RNAP activations within one second.  Assumes steady state
	so activations will be equal to terminations to maintain the appropriate
	active fraction.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule
		doubling_time (float with time units): expected cell doubling time
		synth_prob (ndarray[float]): synthesis probabilities for each RNA

	Returns:
		float: number of RNAP activations/terminations per second

	Math:
		B = bound RNAP
		U = unbound RNAP
		R = rate of binding (1/s)
		f = fraction of total RNAP bound (f = B / (U + B))

		Want to solve for U*R (activation_rate - the rate of RNAP activations)

		For steady state, the amount of bound RNAP will not change:
			dB/dt = synth_prob * U * R - elong_rate / rna_lengths * B = 0

		Using the fraction bound relationship:
			f = B / (U + B)
			R = f / (1 - f) * 1 / (synth_prob * rna_lengths / elong_rate)

		Solving for U*R (with U = (1 - f) (U + B)):
			U*R = U * f / (1 - f) * 1 / (synth_prob * rna_lengths / elong_rate)
			    = B / (synth_prob * rna_lengths / elong_rate)
	'''

	rnap = get_bound_rnap_counts(sim_data, bulk_container, doubling_time)
	rna_lengths = sim_data.process.transcription.rnaData['length'].asNumber()
	elong_rate = sim_data.growthRateParameters.rnaPolymeraseElongationRate.asNumber()

	activation_rate = rnap / (rna_lengths / elong_rate).dot(synth_prob)

	return activation_rate

def get_expected_ppgpp(sim_data, doubling_time):
	'''
	Gets the expected ppGpp concentration for a given doubling time.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		doubling_time (float with time units): expected cell doubling time

	Returns:
		float: concentration of ppGpp (in units of uM)
	'''

	ppgpp = sim_data.growthRateParameters.getppGppConc(doubling_time)
	ppgpp = ppgpp * sim_data.constants.cellDensity * sim_data.mass.cellDryMassFraction

	return ppgpp.asNumber(MICROMOLAR_UNITS)

def get_expected_v_rib(sim_data, nutrients):
	'''
	Gets the expected ribosome elongation rate for a given nutrient condition.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		nutrients (str): nutrient label

	Returns:
		float: rate of ribosome elongation (in units of uM/s)
	'''

	rate = sim_data.translationSupplyRate[nutrients]
	rate = rate * sim_data.constants.cellDensity * sim_data.mass.cellDryMassFraction

	return rate.asNumber(MICROMOLAR_UNITS / units.s).sum()

def get_concentrations(sim_data, bulk_container, doubling_time, rnap_activation_rate, rrna_synth_prob,
		schmidt, ribosome_control):
	'''
	Gets concentrations for relevant molecules

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule
		doubling_time (float with time units): expected cell doubling time
		rnap_activation_rate (float): number of RNAP activations/terminations per second
		rrna_synth_prob (float): synthesis probability for rRNA
		schmidt (bool): if True, uses synthetase concentrations from proteomics
			from Schmidt et al. 2015
		ribosome_control (bool): if True, updates ribosome concentration based on
			ppGpp regulation, otherwise uses bulk container counts

	Returns (in units of uM):
		float: concentration of RelA
		ndarray[float]: concentration of total tRNA for each amino acid
		float: concentration of active ribosomes
		ndarray[float]: concentration of synthetases for each amino acid
		ndarray[float]: concentration of each amino acid
		float: concentration of free RNAPs
	'''

	transcription = sim_data.process.transcription
	molecule_ids = sim_data.moleculeIds
	molecule_groups = sim_data.moleculeGroups

	# Data structures for charging
	aa_from_synthetase = transcription.aa_from_synthetase
	aa_from_trna = transcription.aa_from_trna

	# Names of molecules associated with tRNA charging
	rela_name = molecule_ids.RelA
	uncharged_trna_names = transcription.rnaData['id'][transcription.rnaData['isTRna']]
	synthetase_names = transcription.synthetase_names
	aa_names = molecule_groups.aaIDs

	counts_to_micromolar = get_counts_to_micromolar(doubling_time)

	# Concentrations for tRNA charging molecules
	rela_counts = bulk_container.count(rela_name)
	total_trna_counts = aa_from_trna.dot(
		bulk_container.counts(uncharged_trna_names))
	ribosome_counts = get_active_ribosome_counts(
		sim_data, bulk_container, doubling_time, rnap_activation_rate, rrna_synth_prob, ribosome_control)
	synthetase_counts = aa_from_synthetase.dot(
		bulk_container.counts(synthetase_names))
	aa_counts = bulk_container.counts(aa_names)
	rnap_counts = get_free_rnap_counts(sim_data, bulk_container, doubling_time)

	rela_conc = rela_counts * counts_to_micromolar
	total_trna_conc = total_trna_counts * counts_to_micromolar
	ribosome_conc = ribosome_counts * counts_to_micromolar
	synthetase_conc = synthetase_counts * counts_to_micromolar
	aa_conc = aa_counts * counts_to_micromolar
	rnap_conc = rnap_counts * counts_to_micromolar

	# Make adjustments to synthetase concentrations to match proteomics data
	if schmidt:
		if doubling_time.asNumber() < 40:
			synthetase_conc = SCHMIDT_CONC['with_aa']
		elif doubling_time.asNumber() > 50:
			synthetase_conc = SCHMIDT_CONC['no_oxygen']
		else:
			synthetase_conc = SCHMIDT_CONC['basal']

	return rela_conc, total_trna_conc, ribosome_conc, synthetase_conc, aa_conc, rnap_conc

def get_aa_fraction(sim_data, bulk_container):
	'''
	Gets the fraction of amino acids that would be expected to be translated
	based on mRNA expression and translation efficiencies.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule

	Returns:
		ndarray[float]: fraction of expected elongations for each amino acid
	'''

	transcription = sim_data.process.transcription
	translation = sim_data.process.translation
	mapping = sim_data.relation.rnaIndexToMonomerMapping

	translation_efficiencies = translation.translationEfficienciesByMonomer
	aa_content = translation.monomerData['aaCounts'].asNumber()

	rna_names = transcription.rnaData['id'][mapping]
	rna_counts = bulk_container.counts(rna_names) + 1  # +1 for smoothing int counts

	total_aa = (rna_counts * translation_efficiencies).dot(aa_content)

	return total_aa / total_aa.sum()

def get_objective_value(rrna_prob, expected_rrna_prob, ppgpp, expected_ppgpp, v_rib,
		expected_v_rib, params):
	'''
	Calculate an objective value for the difference between calculated values and expected
	values for rRNA synthesis probability, ppGpp concentration and ribosome elongation rate.

	Args:
		rrna_prob (float): synthesis probability of rRNA based on ppGpp regulation
		expected_rrna_prob (float): expected average synthesis probability rRNA
		ppgpp (float): concentration of ppGpp from state of the cell (in units of uM)
		expected_ppgpp (float): expected concentration of ppGpp (in units of uM)
		v_rib (float): ribosome elongation rate from state of the cell (in units of uM/s)
		expected_v_rib (float): expected ribosome elongation rate (in units of uM/s)
		params (tuple(int, ndarray[float])): objective function and objective
			weights for each component, if None, all weights will be 1

	Returns:
		float: objective value
	'''

	def squared_normalized(actual, expected):
		return ((actual - expected) / expected)**2

	def squared_absolute(actual, expected):
		return (actual - expected)**2

	def accuracy(actual, expected):
		return (actual - expected) / expected

	def get_weight(weights, idx):
		if weights is not None and idx < len(weights):
			return weights[idx]
		else:
			return 1.

	objective_index = params[0]
	weights = params[1]

	# Select objective function from params
	objective_functions = [squared_normalized, squared_absolute, accuracy]
	if objective_index >= len(objective_functions):
		objective_index = 0
	objective_function = objective_functions[objective_index]

	objective = 0
	objective += get_weight(weights, 0) * objective_function(rrna_prob, expected_rrna_prob)
	objective += get_weight(weights, 1) * objective_function(ppgpp, expected_ppgpp)
	objective += get_weight(weights, 2) * objective_function(v_rib, expected_v_rib)

	return objective

def get_growth_constants(constants, no_units=True):
	'''
	Gets a list of growth constants of interest from sim_data constants

	Args:
		constants (class): constants from sim_data
		no_units (bool): strips units from constants if True

	Returns:
		list[float]: constant values from PARAMS
	'''

	if no_units:
		constants_list = []
		for param in PARAMS:
			try:
				constants_list += [getattr(constants, param).asNumber()]
			except AttributeError:
				constants_list += [getattr(constants, param)]
	else:
		constants_list = [getattr(constants, param) for param in PARAMS]

	return constants_list

def charge_trna(total_trna, synthetase_conc, aa_conc, ribosome_conc, f, constants, charged_fraction, t_limit=10):
	'''
	Calculates the concentration of charged and uncharged tRNA from the composition of the cell.

	Args:
		total_trna (ndarray[float]): concentration of all tRNA for each amino acid
			(in units of uM)
		synthetase_conc (ndarray[float]): concentration of synthetases for each amino acid
			(in units of uM)
		aa_conc (ndarray[float]): concentration of amino acids (in units of uM)
		ribosome_conc (float): concentration of active ribosomes
		f (ndarray[float]): fraction of each amino acid in sequences to be translated
		constants (class): constants from sim_data
		charged_fraction (float): fraction of each tRNA species that are charged
		t_limit (float): time limit for charging to prevent long computation times

	Returns:
		ndarray[float]: concentration of charged tRNA for each amino acid (in units of uM)
		ndarray[float]: concentration of uncharged tRNA for each amino acid (in units of uM)
		float: ribosome elongation rate (in units of uM/s)
	'''

	# Parameters from Bosdriesz et al
	k_s = constants.synthetase_charging_rate
	KM_aa = constants.Km_synthetase_amino_acid
	KM_tf = constants.Km_synthetase_uncharged_trna
	k_rta = constants.Kdissociation_charged_trna_ribosome
	k_rtf = constants.Kdissociation_uncharged_trna_ribosome
	rib_elong_rate = 22

	# Initialize to approximate charged levels
	charged_trna_conc = total_trna * charged_fraction
	uncharged_trna_conc = total_trna * (1 - charged_fraction)

	# Solve to steady state with short time steps
	t = 0
	dt = 0.001
	diff = 1
	while diff > 1e-3:
		v_charging = (k_s * synthetase_conc * uncharged_trna_conc * aa_conc
			/ (KM_aa * KM_tf + KM_aa * uncharged_trna_conc + KM_tf * aa_conc
			+ uncharged_trna_conc * aa_conc))
		numerator_ribosome = 1 + np.sum(f * k_rta / charged_trna_conc * (1
			+ uncharged_trna_conc / k_rtf))
		v_rib = rib_elong_rate * ribosome_conc / numerator_ribosome

		# Handle case when f is 0 and charged_trna_conc is 0
		if not np.isfinite(v_rib):
			v_rib = 0

		delta_conc = (v_charging - v_rib * f) * dt
		uncharged_trna_conc -= delta_conc
		charged_trna_conc += delta_conc
		diff = np.sqrt(delta_conc.dot(delta_conc))  # quick norm vs np.linalg.norm

		t += dt

		if t > t_limit:
			print('** Time limit reached, diff: {} **'.format(diff))
			break

	return charged_trna_conc, uncharged_trna_conc, v_rib

def create_ppgpp(rela_conc, charged_trna_conc, uncharged_trna_conc, ribosome_conc, f, constants):
	'''
	Calculates the concentration of ppGpp in a cell.

	Args:
		rela_conc (float): concentration of RelA protein
		charged_trna_conc (ndarray[float]): concentrations of charged tRNA for each amino acid
		uncharged_trna_conc (ndarray[float]): concentrations of uncharged tRNA for each amino acid
		ribosome_conc (float): concentration of active ribosomes
		f (ndarray[float]): fraction of each amino acid in sequences to be translated
		constants (class): constants from sim_data

	Returns:
		float: concentration of ppGpp (in units of uM)
	'''

	# Parameters from Bosdriesz et al
	k_rela = constants.k_RelA_ppGpp_synthesis
	KD_rela = constants.KD_RelA_ribosome
	k_spot_syn = constants.k_SpoT_ppGpp_synthesis
	k_spot_deg = constants.k_SpoT_ppGpp_degradation
	k_rta = constants.Kdissociation_charged_trna_ribosome
	k_rtf = constants.Kdissociation_uncharged_trna_ribosome

	numerator_ribosome = 1 + np.sum(f * (k_rta / charged_trna_conc
		+ uncharged_trna_conc / charged_trna_conc * k_rta / k_rtf))
	ribosome_bound_uncharged = ribosome_conc * (f * uncharged_trna_conc / charged_trna_conc
		* k_rta / k_rtf) / numerator_ribosome
	frac_rela = 1 / (1 + KD_rela / ribosome_bound_uncharged.sum())

	v_rela = k_rela * rela_conc * frac_rela
	v_syn = v_rela + k_spot_syn

	# Steady state solution for ppGpp when synthesis (v_syn) = degradation (k_spot_deg * ppgpp_conc)
	ppgpp_conc = v_syn / k_spot_deg

	return ppgpp_conc

def regulate_rrna_expression(ppgpp, rnap_free, new_rna_rate, constants):
	'''
	Calculates the expression of rRNA based on regulation from ppGpp.

	Args:
		ppgpp (float): concentration of ppGpp (in units of uM)
		rnap_free (float): concentration of free RNAP (in units of uM)
		new_rna_rate (float): rate of RNAP initialization (in units of 1/s)
		constants (class): constants from sim_data

	Returns:
		float: synthesis probability for rRNA
	'''

	vmax = constants.rrn_vmax
	KI_ppgpp = constants.KI_ppgpp_rnap
	KM_rrn = constants.KM_rrn_rnap

	new_rrna_rate = vmax * rnap_free / (KM_rrn + rnap_free) / (1 + ppgpp / KI_ppgpp)
	synth_prob = new_rrna_rate / new_rna_rate

	return synth_prob

def summarize_state(condition, rela, charged_trna, uncharged_trna, ribosomes, synthetases, aas, ppgpp, f):
	'''
	Prints a summary of the state (concentrations) of the cell to stdout

	Args:
		condition (str): cell condition
		rela (float): concentration of RelA
		charged_trna (ndarray[float]): concentrations of charged tRNA for each amino acid
		uncharged_trna (ndarray[float]): concentrations of uncharged tRNA for each amino acid
		ribosomes (float): concentration of active ribosomes
		synthetases (ndarray[float]): concentration of synthetases for each amino acid
		aas (ndarray[float]): concentration for each amino acid
		ppgpp (float): concentration of ppGpp
		f (ndarray[float]): fraction of expected elongations for each amino acid
	'''

	def print_conc(label, conc):
		if isinstance(conc, float):
			print('\t{}: {:.2f}'.format(label, conc))
		else:
			print('\t{}: {}'.format(label, conc))

	np.set_printoptions(precision=1, suppress=True, linewidth=120)
	print('\nSummary for {} condition:'.format(condition))
	print_conc('RelA', rela)
	print_conc('Charged tRNA', charged_trna)
	print_conc('Uncharged tRNA', uncharged_trna)
	print_conc('Ribosomes', ribosomes)
	print_conc('Synthetases', synthetases)
	print_conc('Amino acids', aas)
	print_conc('ppGpp', ppgpp)
	np.set_printoptions(precision=2)
	print_conc('f', f)
	np.set_printoptions(precision=8, suppress=False, linewidth=75)  # reset to defaults

def error_summary(rrna_prob, expected_rrna_prob, ppgpp, expected_ppgpp, v_rib, expected_v_rib):
	'''
	Report the error between expected state and calculated state for rRNA synthesis
	probability, ppGpp concentration and ribosome elongation rate.

	Args:
		rrna_prob (float): synthesis probability of rRNA based on ppGpp regulation
		expected_rrna_prob (ndarray[float]): synthesis probabilities of each rRNA
			from the fitter
		ppgpp (float): concentration of ppGpp from state of the cell (in units of uM)
		expected_ppgpp (float): expected concentration of ppGpp (in units of uM)
		v_rib (float): ribosome elongation rate from state of the cell (in units of uM/s)
		expected_v_rib (float): expected ribosome elongation rate (in units of uM/s)
	'''

	def print_error(label, actual, expected, format):
		error = np.abs(actual - expected) / expected
		display_str = '\t{{}} error: {{:.1f}}% ({{{}}} vs {{{}}})'.format(format, format)
		print(display_str.format(label, error*100, actual, expected))

	expected_rrna_prob = np.mean(expected_rrna_prob[expected_rrna_prob > 0])

	print('')
	print_error('rRNA probability', rrna_prob, expected_rrna_prob, ':.3f')
	print_error('ppGpp concentration', ppgpp, expected_ppgpp, ':.1f')
	print_error('Ribosome elongation rate', v_rib, expected_v_rib, ':.1f')

def sensitivity(sim_data, cell_specs, conditions, schmidt, objective_params, ribosome_control):
	'''
	Performs sensitivity analysis for each of the parameters.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		cell_specs (dict): information about each condition that was fit
		conditions (list[str]): set of conditions to test
			(eg. ['basal', 'with_aa', 'no_oxygen'])
		schmidt (bool): if True, uses synthetase concentrations from proteomics
			from Schmidt et al. 2015
		objective_params (tuple(int, ndarray[float])): objective function and objective
			weights for each component, if None, all weights will be 1
		ribosome_control (bool): if True, updates ribosome concentration based on
			ppGpp regulation, otherwise uses bulk container counts
	'''

	for param in PARAMS:
		print('For {}:'.format(param))
		original_value = getattr(sim_data.constants, param)

		for magnitude in [0.1, 0.2, 0.5, 0.75, 0.9, 1.1, 1.5, 2, 5, 10]:
			setattr(sim_data.constants, param, original_value * magnitude)
			error = [
				main(sim_data, cell_specs, [condition], schmidt, objective_params, ribosome_control, verbose=False)
				for condition in conditions
				]

			error_str = '  '.join(format(e[0], '.2f') for e in error)
			print('\t x{} error: {}'.format(magnitude, error_str))

		setattr(sim_data.constants, param, original_value)

def coordinate_descent(sim_data, cell_specs, conditions, schmidt, objective_params,
		ribosome_control, update_factors=0, update_synthetases=False, update_aas=False):
	'''
	Stochastic coordinate descent to determine optimal parameters.  Updates one
	parameter at a time to minimize error for a given number of iterations or
	until the error does not change.  Displays the set of parameters and results
	in each condition after convergence.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		cell_specs (dict): information about each condition that was fit
		conditions (list[str]): set of conditions to test
			(eg. ['basal', 'with_aa', 'no_oxygen'])
		schmidt (bool): if True, uses synthetase concentrations from proteomics
			from Schmidt et al. 2015
		objective_params (tuple(int, ndarray[float])): objective function and objective
			weights for each component, if None, all weights will be 1
		ribosome_control (bool): if True, updates ribosome concentration based on
			ppGpp regulation, otherwise uses bulk container counts
		update_factors (int): if a positive value, synthetase or amino acid concentrations
			will be updated every update_factors time steps
		update_synthetases (bool): if True, find synthetase concentrations to minimize objective
		update_aas (bool): if True, find amino acid concentrations to minimize objective

	Returns:
		Constants object: class of all constants values from sim_data
		list[float]: factors to multiply synthetase or amino acid concentrations by in each condition
		float: objective value reached
		float: reference objective value
	'''

	max_it = 1000
	max_constants = 15
	it = 1
	n_params = len(PARAMS)
	propensity = np.ones(n_params)
	objective = 10000000  # high starting value
	n_constants = 0  # variable for tracking how many steps have been held constant
	objective_limit = 0.001  # below limit, change is assumed constant
	delta = 0.1  # rate of change at each step
	decay = 0.9  # rate of decay of delta
	factors = np.ones(len(conditions))
	charged_fraction = None

	try:
		while it < max_it:
			idx = np.where(np.random.multinomial(1, propensity / propensity.sum()))[0][0]
			param = PARAMS[idx]

			if it % 20 == 0:
				delta *= decay
				print('Iteration {} (delta = {:.5f})'.format(it, delta))

			original_value = getattr(sim_data.constants, param)

			# Change high
			setattr(sim_data.constants, param, original_value * (1 + delta))
			high_objective, high_charged, _ = main(sim_data, cell_specs, conditions, schmidt, objective_params,
				ribosome_control, charged=charged_fraction, factors=factors, update_synthetases=update_synthetases,
				update_aas=update_aas, verbose=False)

			# Change low
			setattr(sim_data.constants, param, original_value * (1 - delta))
			low_objective, low_charged, _ = main(sim_data, cell_specs, conditions, schmidt, objective_params,
				ribosome_control, charged=charged_fraction, factors=factors, update_synthetases=update_synthetases,
				update_aas=update_aas, verbose=False)

			# Update to new parameter value based on lowest error
			if low_objective < objective and low_objective < high_objective:
				delta_objective = objective - low_objective
				propensity[idx] = max(delta_objective, propensity.min())
				objective = low_objective
				charged_fraction = low_charged
				status = 'decreased'
			elif high_objective < objective:
				delta_objective = objective - high_objective
				propensity[idx] = max(delta_objective, propensity.min())
				objective = high_objective
				charged_fraction = high_charged
				setattr(sim_data.constants, param, original_value * (1 + delta))
				status = 'increased'
			else:
				delta_objective = 0
				propensity[idx] = max(propensity.min() * 0.5, propensity.max() * 0.01)
				setattr(sim_data.constants, param, original_value)
				status = 'held constant'

			if delta_objective < objective_limit:
				n_constants += 1
			else:
				n_constants = 0

			# Initialize to equal propensities
			if it == 1:
				propensity[:] = objective

			# Update synthetase concentrations to improve objective if set
			if update_factors > 0 and it % update_factors == 0:
				new_factors, new_objective = find_concentrations(
					sim_data, cell_specs, conditions, schmidt, objective_params, ribosome_control,
					charged=charged_fraction, factors=factors, update_synthetases=update_synthetases,
					update_aas=update_aas, max_it=1, verbose=False)

				if new_objective < objective:
					if objective - new_objective > objective_limit:
						n_constants = 0
					factors = new_factors
					objective = new_objective

			it += 1
			print('{} {}: {:.3f}'.format(param, status, objective))

			# Early break condition if nothing changing
			if n_constants > max_constants:
				break
	# Allow early interrupt
	except KeyboardInterrupt:
		pass

	# Summarize results
	for param in PARAMS:
		print('{}: {}'.format(param, getattr(sim_data.constants, param)))

	if update_synthetases:
		conc_changed = 'Synthetase'
	else:
		conc_changed = 'Amino acid'
	for factor, condition in zip(factors, conditions):
		print('{} change in {}: {:.3f}'.format(conc_changed, condition, factor))

	objective, _, ref_objective = main(sim_data, cell_specs, conditions, schmidt, objective_params,
		ribosome_control, factors=factors, update_synthetases=update_synthetases, update_aas=update_aas)

	return sim_data.constants, list(factors), objective, ref_objective

def find_concentrations(sim_data, cell_specs, conditions, schmidt, objective_params, ribosome_control,
		charged=None, factors=None, update_synthetases=False, update_aas=False, max_it=None, verbose=True):
	'''
	Use gradient descent to determine optimal synthetase or amino acid concentrations.
	Updates concentrations in one condition at a time to minimize error until
	convergence. Displays the concentrations for each condition after convergence.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		cell_specs (dict): information about each condition that was fit
		conditions (list[str]): set of conditions to test
			(eg. ['basal', 'with_aa', 'no_oxygen'])
		schmidt (bool): if True, uses synthetase concentrations from proteomics
			from Schmidt et al. 2015
		objective_params (tuple(int, ndarray[float])): objective function and objective
			weights for each component, if None, all weights will be 1
		ribosome_control (bool): if True, updates ribosome concentration based on
			ppGpp regulation, otherwise uses bulk container counts
		charged (list[ndarray[float]]): fraction of charged tRNA for each amino acid
			in each condition
		factors (ndarray[float]): factor to multiply synthetase or amino acid
			concentrations by for each condition
		update_synthetases (bool): if True, find synthetase concentrations to minimize objective
		update_aas (bool): if True, find amino acid concentrations to minimize objective
		max_it (int): maximum number of iterations to update before returning
		verbose (bool): if True, prints function specific information

	Returns:
		ndarray[float]: factors to adjust synthetases or amino acids by in each condition
		float: objective value reached
	'''

	total_objective = 0
	delta = 0.01  # rate of change at each step
	eps = 0.5
	if factors is None:
		factors = np.ones(len(conditions))
	it = 0

	for idx, condition in enumerate(conditions):
		factor = factors[idx]
		if charged is None:
			charged_fraction = None
		else:
			charged_fraction = charged[idx]

		# Change high
		high_factor = factor * (1 + delta)
		high_objective, high_charged, _ = main(sim_data, cell_specs, [condition], schmidt, objective_params,
			ribosome_control, charged=charged_fraction, factors=[high_factor],
			update_synthetases=update_synthetases, update_aas=update_aas, verbose=False)

		# Change low
		low_factor = factor * (1 - delta)
		low_objective, low_charged, _ = main(sim_data, cell_specs, [condition], schmidt, objective_params,
			ribosome_control, charged=charged_fraction, factors=[low_factor],
			update_synthetases=update_synthetases, update_aas=update_aas, verbose=False)

		# Move in direction of decreasing objective
		if high_objective < low_objective:
			direction = 1
			objective = high_objective
			factor = high_factor
			charged_fraction = high_charged
		else:
			direction = -1
			objective = low_objective
			factor = low_factor
			charged_fraction = low_charged

		if verbose:
			print(objective, factor)

		# Update in direction of decreasing objective until objective converges
		converged = 0
		best_objective = objective
		best_factor = factor
		while converged < 3:
			it += 1

			if max_it and it >= max_it:
				break

			factor *= 1 + direction * delta
			objective, charged_fraction, _ = main(sim_data, cell_specs, [condition], schmidt, objective_params,
				ribosome_control, charged=charged_fraction, factors=[factor], update_synthetases=update_synthetases,
				update_aas=update_aas, verbose=False)

			if objective > best_objective - eps:
				converged += 1
			if objective < best_objective:
				best_objective = objective
				best_factor = factor

			if verbose:
				print(objective, factor)

		factors[idx] = best_factor
		total_objective += best_objective

	# Summarize results
	if verbose:
		main(sim_data, cell_specs, conditions, schmidt, objective_params, ribosome_control,
			factors=factors, update_synthetases=update_synthetases, update_aas=update_aas)

	return factors, total_objective

def plot_synthetases(sim_data, cell_specs, conditions, objective_params, ribosome_control, path):
	'''
	Plot synthetase concentrations in the different conditions comparing model
	and proteomics data from Schmidt et al. 2015.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		cell_specs (dict): information about each condition that was fit
		conditions (list[str]): set of conditions to test
			(eg. ['basal', 'with_aa', 'no_oxygen'])
		objective_params (tuple(int, ndarray[float])): objective function and objective
			weights for each component, if None, all weights will be 1
		ribosome_control (bool): if True, updates ribosome concentration based on
			ppGpp regulation, otherwise uses bulk container counts
		path (str): path to output file

	Output:
		.png file specified by path containing the plot
	'''

	def format_axes(ax, title, ticks, x_labeled=False):
		ax.set_ylim(0, 5)
		ax.set_title(title, fontsize=10)
		ax.set_ylabel('Conc (uM)', fontsize=8)
		ax.set_xticks(ticks)
		ax.tick_params(labelbottom=x_labeled)

	is_rrna = sim_data.process.transcription.rnaData['isRRna']
	model_factors, _ = find_concentrations(sim_data, cell_specs, conditions, False,
		objective_params, ribosome_control, update_synthetases=True)
	schmidt_factors, _ = find_concentrations(sim_data, cell_specs, conditions, True,
		objective_params, ribosome_control, update_synthetases=True)

	# Setup plot
	plt.figure(figsize=(8.5, 11))
	n_subplots = 4
	bar_width = 0.3
	ax_model = plt.subplot(n_subplots, 1, 1)
	ax_schmidt = plt.subplot(n_subplots, 1, 2)
	ax_model_adj = plt.subplot(n_subplots, 1, 3)
	ax_schmidt_adj = plt.subplot(n_subplots, 1, 4)

	# Get and plot concentrations for each condition
	for i, condition in enumerate(conditions):
		bulk_container = cell_specs[condition]['bulkAverageContainer']
		doubling_time = sim_data.conditionToDoublingTime[condition]
		synth_prob = cell_specs[condition]['synthProb']

		rrna_synth_prob = np.mean(synth_prob[is_rrna][synth_prob[is_rrna] > 0])

		rnap_activation_rate = get_rnap_activation(
			sim_data, bulk_container, doubling_time, synth_prob)
		_, _, _, model_synthetases, _, _ = get_concentrations(
			sim_data, bulk_container, doubling_time, rnap_activation_rate, rrna_synth_prob,
			False, ribosome_control)
		schmidt_synthetases = SCHMIDT_CONC[condition]

		x = np.arange(len(model_synthetases))
		ax_model.bar(x + i*bar_width, model_synthetases, bar_width)
		ax_schmidt.bar(x + i*bar_width, schmidt_synthetases, bar_width)
		ax_model_adj.bar(x + i*bar_width, model_synthetases * model_factors[i], bar_width)
		ax_schmidt_adj.bar(x + i*bar_width, schmidt_synthetases * schmidt_factors[i], bar_width)

	# Format plot
	format_axes(ax_model, 'Model Concentrations', x + bar_width)
	format_axes(ax_schmidt, 'Proteomics Concentrations', x + bar_width)
	format_axes(ax_model_adj, 'Adjusted Model Concentrations', x + bar_width)
	format_axes(ax_schmidt_adj, 'Adjusted Proteomics Concentrations', x + bar_width, x_labeled=True)
	ax_schmidt_adj.legend(conditions, fontsize=8)
	ax_schmidt_adj.set_xticklabels(sim_data.moleculeGroups.aaIDs, rotation=40, ha='right', fontsize=8)
	plt.tight_layout()

	# Save output
	plt.savefig(path)

def plot_parameters(data, path):
	'''
	Plot histograms of parameter distributions for a set of coordinate descent runs.

	Args:
		data (ndarray[str]): 2D array of data loaded from .tsv file
		path (str): path to output file (without extension)

	Output:
		Two .png files (<path>_hist.png and <path>_splom.png) containing the plots
	'''

	def display_param(param, length=20):
		display = ''
		line = ''
		for term in param.split('_'):
			if len(line):
				if len(line) + len(term) > length:
					display += line + '\n'
					line = term
				else:
					line += ' ' + term
			else:
				line = term
		display += line

		return display

	# Parse data
	header = data[0, :]
	objective_col = np.where(header == 'Objective')[0]
	start_parameter_col = int(objective_col + 1)
	params = header[start_parameter_col:]
	params = [display_param(param) for param in params]
	params_plotly = [param.replace('\n', '<br>') for param in params]
	n_params = len(params)

	reference_parameters = np.array(data[1, start_parameter_col:], float)
	modified_parameters = np.array(data[2:, start_parameter_col:], float)
	objective = np.array(data[2:, objective_col], float).squeeze()

	# Plot distributions
	plt.figure(figsize=(8.5, 11))
	n_rows = np.ceil(np.sqrt(n_params))
	n_cols = np.ceil(n_params / n_rows)

	## Subplot for each parameter
	for idx, param in enumerate(params):
		ax = plt.subplot(n_rows, n_cols, idx + 1)

		ax.hist(modified_parameters[:, idx])
		ax.axvline(reference_parameters[idx], color='r', linestyle='--')
		ax.set_title(param, fontsize=8)
		ax.tick_params(labelsize=7)

	## Save histogram output
	plt.tight_layout()
	plt.savefig(path + '_hist.png')

	# Plot splom
	reference_trace = go.Splom(
		name='Reference',
		dimensions=[
			dict(label=param, values=[v])
			for param, v in zip(params_plotly, reference_parameters)
			],
		showupperhalf=False,
		marker=dict(
			color='rgb(100,180,100)',
			line=dict(
				width=0.5,
				color='rgb(150,150,150)',
				),
			),
		)
	modified_trace = go.Splom(
		dimensions=[
			dict(label=param, values=v)
			for param, v in zip(params_plotly, modified_parameters.T)
			],
		showupperhalf=False,
		showlegend=False,
		marker=dict(
			color=np.log(objective),
			colorscale='RdBu',
			line=dict(
				width=0.5,
				color='rgb(150,150,150)',
				),
			colorbar=dict(
				title='ln(Objective)',
				ypad=50,
				),
			),
		)

	## Formatting
	axis_format = dict(
		showline=True,
		zeroline=False,
		gridcolor='#fff',
		ticklen=4,
		title=dict(
			font=dict(
				size=8,
				),
			),
		tickfont=dict(
			size=6,
			),
		)
	layout = go.Layout(
		title='Parameters values and objective',
		width=1600,
		height=1600,
		plot_bgcolor='rgba(240,240,240,0.95)',
		)
	layout.update({
		'xaxis{}'.format(idx + 1): axis_format
		for idx in range(n_params)
		})
	layout.update({
		'yaxis{}'.format(idx + 1): axis_format
		for idx in range(n_params)
		})

	## Save splom output
	fig = dict(data=[reference_trace, modified_trace], layout=layout)
	plotly.io.write_image(fig, path + '_splom.png')

def main(sim_data, cell_specs, conditions, schmidt, objective_params, ribosome_control,
		charged=None, factors=None, update_synthetases=False, update_aas=False, verbose=True):
	'''
	Main function to perform analysis.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		cell_specs (dict): information about each condition that was fit
		conditions (list[str]): set of conditions to test
			(eg. ['basal', 'with_aa', 'no_oxygen'])
		schmidt (bool): if True, uses synthetase concentrations from proteomics
			from Schmidt et al. 2015
		objective_params (tuple(int, ndarray[float])): objective function and objective
			weights for each component, if None, all weights will be 1
		ribosome_control (bool): if True, updates ribosome concentration based on
			ppGpp regulation, otherwise uses bulk container counts
		charged (list[ndarray[float]]): fraction of charged tRNA for each amino acid
			in each condition, if None, charged fraction starts at 0
		factors (ndarray[float]): factor to multiply synthetase or amino acid
			concentrations by for each condition, if None, defaults to 1
		update_synthetases (bool): if True, find synthetase concentrations to minimize objective
		update_aas (bool): if True, find amino acid concentrations to minimize objective
		verbose (bool): if True, prints to command line

	Returns:
		float: objective value for all conditions
		list[ndarray[float]]: fraction of charged tRNA for each amino acid in each
			condition, for fast iteration
		float: reference objective value for all conditions
	'''

	objective = 0
	ref_objective = 0

	if ribosome_control:
		iter = 10
	else:
		iter = 1

	constants = sim_data.constants
	is_rrna = sim_data.process.transcription.rnaData['isRRna']
	if factors is None:
		factors = np.ones(len(conditions))
	if charged is None:
		charged = np.zeros(len(conditions))

	final_charged_fractions = []

	# Calculate ppGpp concentrations for each condition
	for condition, factor, charged_fraction in zip(conditions, factors, charged):
		# Get condition specific values
		bulk_container = cell_specs[condition]['bulkAverageContainer']
		synth_prob = cell_specs[condition]['synthProb']
		doubling_time = sim_data.conditionToDoublingTime[condition]
		nutrients = sim_data.conditions[condition]['nutrients']

		expected_rrna_synth_prob = np.mean(synth_prob[is_rrna][synth_prob[is_rrna] > 0])
		rrna_synth_prob = expected_rrna_synth_prob

		f = get_aa_fraction(sim_data, bulk_container)
		expected_ppgpp = get_expected_ppgpp(sim_data, doubling_time)
		expected_v_rib = get_expected_v_rib(sim_data, nutrients)

		for i in range(iter):
			updated_synth_prob = synth_prob.copy()
			updated_synth_prob[is_rrna & (synth_prob > 0)] = rrna_synth_prob
			updated_synth_prob /= updated_synth_prob.sum()
			rnap_activation_rate = get_rnap_activation(
				sim_data, bulk_container, doubling_time, updated_synth_prob)
			rela, total_trnas, ribosomes, synthetases, aas, rnaps = get_concentrations(
				sim_data, bulk_container, doubling_time, rnap_activation_rate, rrna_synth_prob,
				schmidt, ribosome_control)
			if update_synthetases:
				synthetases = synthetases * factor
			elif update_aas:
				aas = aas * factor

			# Calculate values from current state
			if i > 0:
				charged_fraction = charged_trna / total_trnas
			charged_trna, uncharged_trna, v_rib = charge_trna(
				total_trnas, synthetases, aas, ribosomes, f, constants, charged_fraction)
			ppgpp = create_ppgpp(rela, charged_trna, uncharged_trna, ribosomes, f, constants)
			new_rrna_synth_prob = regulate_rrna_expression(ppgpp, rnaps, rnap_activation_rate, constants)

			if ribosome_control:
				rrna_synth_prob += (new_rrna_synth_prob - rrna_synth_prob) * 0.1  # 0.1 for stability
			else:
				rrna_synth_prob = new_rrna_synth_prob

		# Print current state
		if verbose:
			summarize_state(condition, rela, charged_trna, uncharged_trna,
				ribosomes, synthetases, aas, ppgpp, f)
			error_summary(rrna_synth_prob, synth_prob[is_rrna], ppgpp,
				expected_ppgpp, v_rib, expected_v_rib)

		objective += get_objective_value(rrna_synth_prob, expected_rrna_synth_prob, ppgpp,
			expected_ppgpp, v_rib, expected_v_rib, objective_params)
		ref_objective += get_objective_value(rrna_synth_prob, expected_rrna_synth_prob, ppgpp,
			expected_ppgpp, v_rib, expected_v_rib, (2, None))

		final_charged_fractions += [charged_trna / total_trnas]

	return objective, final_charged_fractions, ref_objective

def parse_args():
	'''
	Parses arguments from the command line.

	Returns:
		ArgumentParser namespace: values of variables parsed from the command line
	'''

	default_sim_data = os.path.join(DATA_DIR, 'sim_data.cp')
	default_cell_specs = os.path.join(DATA_DIR, 'cell_specs.cp')
	default_seeds = 1
	default_update = 0

	parser = argparse.ArgumentParser()

	# General arguments
	parser.add_argument('-s', '--sim-data',
		default=default_sim_data,
		help='Path to sim_data object (default: {})'.format(default_sim_data))
	parser.add_argument('-c', '--cell-specs',
		default=default_cell_specs,
		help='Path to cell_specs object (default: {})'.format(default_cell_specs))
	parser.add_argument('--condition',
		type=int,
		default=None,
		help='Only runs for specified condition if set (values: 0-2)')
	parser.add_argument('--schmidt',
		action='store_true',
		help='Uses synthetase concentrations from proteomics from Schmidt et al 2015 if set')
	parser.add_argument('--no-ribosome-control',
		action='store_false', dest='ribosome_control',
		help='Uses bulk container ribosome counts if set, otherwise updates ribosome'
			+ ' concentrations based on ppGpp regulation until convergence.')
	parser.add_argument('-o', '--out',
		default=None,
		help='Output file name saved in out/, (default: varies depending on what is saved)')
	parser.add_argument('--objective',
		type=int,
		default=0,
		help='Objective to select (default: 0)')
	parser.add_argument('--objective-weights',
		nargs='+',
		help='Objective function weights')
	parser.add_argument('--random-init',
		action='store_true',
		help='Initializes parameters with random noise if set')

	# Sensitivity arguments
	parser.add_argument('--sensitivity',
		action='store_true',
		help='Perform parameter sensitivity analysis if set')

	# Coordinate descent arguments
	parser.add_argument('--sgd',
		action='store_true',
		help='Perform stochastic gradient descent to find parameters')
	parser.add_argument('--seeds',
		type=int,
		default=default_seeds,
		help='Number of seeds to perform sgd for (default: {})'.format(default_seeds))
	parser.add_argument('--seed',
		type=int,
		default=None,
		help='Random seed to set for sgd, choosen randomly if not set')
	parser.add_argument('--update-synthetases',
		type=int,
		default=default_update,
		help='Number of time steps to update synthetases during sgd (default: {})'
			.format(default_update))
	parser.add_argument('--update-aas',
		type=int,
		default=default_update,
		help='Number of time steps to update amino acids during sgd (default: {})'
			.format(default_update))

	# Concentration search
	parser.add_argument('--synthetases',
		action='store_true',
		help='Search for optimal synthetase concentrations')
	parser.add_argument('--aas',
		action='store_true',
		help='Search for optimal amino acid concentrations')

	# Plotting
	parser.add_argument('--plot-synthetases',
		action='store_true',
		help='Plot analysis of problem if set')
	parser.add_argument('--plot-parameters',
		default=None,
		help='Specify input .tsv file in out/ for plotting parameter distribution')

	return parser.parse_args()


if __name__ == '__main__':
	start = time.time()

	args = parse_args()

	# Check args for conditions that do not run together
	single_args = ['sensitivity', 'sgd', 'synthetases', 'plot_synthetases', 'plot_parameters']
	if sum([1 for a in single_args if getattr(args, a)]) > 1:
		raise Exception('Cannot have multiple options selected {} at the same time.'.format(
			tuple(('--' + a.replace('_', '-') for a in single_args))))

	# Specify conditions
	conditions = ['with_aa', 'basal', 'no_oxygen']
	if args.condition is not None:
		conditions = [conditions[args.condition]]

	# Convert objective to float if provided
	if args.objective_weights:
		objective_params = (args.objective, np.array(args.objective_weights, float))
	else:
		objective_params = (args.objective, None)

	# Load necessary files
	if not args.plot_parameters:
		with open(args.sim_data) as f:
			sim_data = cPickle.load(f)
		with open(args.cell_specs) as f:
			cell_specs = cPickle.load(f)

		# Update constants that are not yet in sim_data
		# TODO: add to sim_data.constants in a flat file
		constants = sim_data.constants
		constants.rrn_vmax = 2000  # 1/s, Bosdriesz
		constants.KI_ppgpp_rnap = 1  # uM, Bosdriesz
		constants.KM_rrn_rnap = 20  # uM, Bosdriesz

		# Strip units
		constants.synthetase_charging_rate = constants.synthetase_charging_rate.asNumber(1 / units.s)
		constants.Km_synthetase_amino_acid = constants.Km_synthetase_amino_acid.asNumber(MICROMOLAR_UNITS)
		constants.Km_synthetase_uncharged_trna = constants.Km_synthetase_uncharged_trna.asNumber(MICROMOLAR_UNITS)
		constants.Kdissociation_charged_trna_ribosome = constants.Kdissociation_charged_trna_ribosome.asNumber(MICROMOLAR_UNITS)
		constants.Kdissociation_uncharged_trna_ribosome = constants.Kdissociation_uncharged_trna_ribosome.asNumber(MICROMOLAR_UNITS)
		constants.k_RelA_ppGpp_synthesis = constants.k_RelA_ppGpp_synthesis.asNumber(1 / units.s)
		constants.KD_RelA_ribosome = constants.KD_RelA_ribosome.asNumber(MICROMOLAR_UNITS)
		constants.k_SpoT_ppGpp_synthesis = constants.k_SpoT_ppGpp_synthesis.asNumber(MICROMOLAR_UNITS / units.s)
		constants.k_SpoT_ppGpp_degradation = constants.k_SpoT_ppGpp_degradation.asNumber(1 / units.s)

	# Perform desired analysis
	if args.sensitivity:
		if args.random_init:
			add_parameter_noise(sim_data.constants)
		sensitivity(sim_data, cell_specs, conditions, args.schmidt, objective_params, args.ribosome_control)
	elif args.sgd:
		out = output_location(args.out, OUTPUT_DIR, SGD_OUT)

		update_factors = 0
		update_synthetases = False
		update_aas = False
		conc_updated = 'Synthetases'
		if args.update_synthetases:
			update_factors = args.update_synthetases
			update_synthetases = True
		elif args.update_aas:
			update_factors = args.update_aas
			update_aas = True
			conc_updated = 'Amino acids'

		with open(out, 'w') as f:
			csv_writer = csv.writer(f, delimiter='\t')
			csv_writer.writerow(['Seed', 'Reference Objective', 'Objective'] + PARAMS
				+ ['{} in {}'.format(conc_updated, c) for c in conditions])
			csv_writer.writerow(['Original', '', ''] + get_growth_constants(sim_data.constants)
				+ [1 for c in conditions])
			f.flush()

			original_constants = {param: getattr(sim_data.constants, param) for param in PARAMS}

			for seed in np.random.randint(0, 100000, args.seeds):
				if args.seed is not None:
					seed = args.seed
				if args.random_init:
					add_parameter_noise(sim_data.constants)
				print('Seed: {}'.format(seed))
				np.random.seed(seed)
				constants, factors, objective, ref_objective = coordinate_descent(
					sim_data, cell_specs, conditions, args.schmidt, objective_params,
					args.ribosome_control, update_factors=update_factors,
					update_synthetases=update_synthetases, update_aas=update_aas)
				csv_writer.writerow([seed, ref_objective, objective]
					+ get_growth_constants(constants) + factors)
				f.flush()

				# Reset constants to original values
				for param, value in original_constants.items():
					setattr(sim_data.constants, param, value)
	elif args.synthetases:
		find_concentrations(sim_data, cell_specs, conditions, args.schmidt, objective_params,
			args.ribosome_control, update_synthetases=True)
	elif args.aas:
		find_concentrations(sim_data, cell_specs, conditions, args.schmidt, objective_params,
			args.ribosome_control, update_aas=True)
	elif args.plot_synthetases:
		out = output_location(args.out, OUTPUT_DIR, SYNTHETASE_PLOT_OUT)

		plot_synthetases(sim_data, cell_specs, conditions, objective_params, args.ribosome_control, out)
	elif args.plot_parameters:
		out = output_location(args.out, OUTPUT_DIR, args.plot_parameters.split('.')[0])

		with open(os.path.join(OUTPUT_DIR, args.plot_parameters)) as f:
			reader = csv.reader(f, delimiter='\t')
			data = np.array(list(reader))

		plot_parameters(data, out)
	else:
		main(sim_data, cell_specs, conditions, args.schmidt, objective_params, args.ribosome_control)

	print('Completed in {:.2f} min'.format((time.time() - start) / 60))
