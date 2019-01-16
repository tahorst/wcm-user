'''
Works with output from the fitter to determine parameters and feasibility
of using ppGpp for growth rate control.

Requires (cPickle files, defaults to files in script directory):
	sim_data (SimulationData object): output from fitter and knowledgebase
		for a simulation, saved after running fitter
	cell_specs (dict): information about each condition that was fit, not saved
		after running fitter but can save during execution
'''

from __future__ import division

import argparse
import cPickle
import csv
import os
import time

import numpy as np

from wholecell.utils import units


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
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


def get_volume(sim_data, bulk_container):
	'''
	Calculates the volume of a cell for a given composition assuming constant
	cell density.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule

	Returns:
		float: volume of the cell (in units of L)
	'''

	mol_names = list(bulk_container.objectNames())
	mol_masses = sim_data.getter.getMass(mol_names)
	mol_counts = bulk_container.counts(mol_names)

	mass = mol_counts.dot(mol_masses) / sim_data.constants.nAvogadro
	volume = mass / sim_data.constants.cellDensity

	return volume.asNumber(units.L)

def get_active_ribosome_counts(sim_data, bulk_container, doubling_time):
	'''
	Gets the counts of active ribosomes based on fraction active and subunits.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule
		doubling_time (float with time units): expected cell doubling time

	Returns:
		float: number of active ribosomes
	'''

	molecule_ids = sim_data.moleculeIds
	active_fraction = sim_data.growthRateParameters.getFractionActiveRibosome(doubling_time)

	count_30s = bulk_container.count(molecule_ids.s30_fullComplex)
	count_50s = bulk_container.count(molecule_ids.s50_fullComplex)

	ribosome_counts = min(count_30s, count_50s) * active_fraction

	return ribosome_counts

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

def get_expected_v_rib(sim_data, doubling_time, nutrients):
	'''
	Gets the expected ribosome elongation rate for a given doubling time.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		doubling_time (float with time units): expected cell doubling time
		nutrients (str): nutrient label

	Returns:
		float: rate of ribosome elongation (in units of uM/s)
	'''

	rate = sim_data.translationSupplyRate[nutrients]
	rate = rate * sim_data.constants.cellDensity * sim_data.mass.cellDryMassFraction

	return rate.asNumber(MICROMOLAR_UNITS / units.s).sum()

def get_concentrations(sim_data, bulk_container, doubling_time):
	'''
	Gets concentrations for relevant molecules

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule
		doubling_time (float with time units): expected cell doubling time

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
	n_avogadro = sim_data.constants.nAvogadro

	# Data structures for charging
	aa_from_synthetase = transcription.aa_from_synthetase
	aa_from_trna = transcription.aa_from_trna

	# Names of molecules associated with tRNA charging
	rela_name = molecule_ids.RelA
	uncharged_trna_names = transcription.rnaData['id'][transcription.rnaData['isTRna']]
	synthetase_names = transcription.synthetase_names
	aa_names = molecule_groups.aaIDs

	volume = get_volume(sim_data, bulk_container)
	counts_to_micromolar = 1 / n_avogadro.asNumber(1 / units.umol) / volume

	# Concentrations for tRNA charging molecules
	rela_counts = bulk_container.count(rela_name)
	total_trna_counts = aa_from_trna.dot(
		bulk_container.counts(uncharged_trna_names))
	ribosome_counts = get_active_ribosome_counts(sim_data, bulk_container, doubling_time)
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

def get_objective_value(rrna_prob, expected_rrna_prob, ppgpp, expected_ppgpp, v_rib, expected_v_rib):
	'''
	Calculate an objective value for the difference between calculated values and expected
	values for rRNA synthesis probability, ppGpp concentration and ribosome elongation rate.

	Args:
		rrna_prob (float): synthesis probability of rRNA based on ppGpp regulation
		expected_rrna_prob (ndarray[float]): synthesis probabilities of each rRNA
			from the fitter
		ppgpp (float): concentration of ppGpp from state of the cell (in units of uM)
		expected_ppgpp (float): expected concentration of ppGpp (in units of uM)
		v_rib (float): ribosome elongation rate from state of the cell (in units of uM/s)
		expected_v_rib (float): expected ribosome elongation rate (in units of uM/s)

	Returns:
		float: objective value

	TODO: allow for different objectives selected with an argument
	'''

	def objective_function(actual, expected):
		return (actual - expected)**2

	expected_rrna_prob = np.mean(expected_rrna_prob[expected_rrna_prob > 0])

	objective = 0
	objective += objective_function(rrna_prob, expected_rrna_prob)
	objective += objective_function(ppgpp, expected_ppgpp)
	objective += objective_function(v_rib, expected_v_rib)

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

def charge_trna(total_trna, synthetase_conc, aa_conc, ribosome_conc, f, constants, t_limit=10):
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
		t_limit (float): time limit for charging to prevent long computation times

	Returns:
		ndarray[float]: concentration of charged tRNA for each amino acid (in units of uM)
		ndarray[float]: concentration of uncharged tRNA for each amino acid (in units of uM)
		float: ribosome elongation rate (in units of uM/s)
	'''

	# Parameters from Bosdriesz et al
	k_s = constants.synthetase_charging_rate.asNumber(1 / units.s)
	KM_aa = constants.Km_synthetase_amino_acid.asNumber(MICROMOLAR_UNITS)
	KM_tf = constants.Km_synthetase_uncharged_trna.asNumber(MICROMOLAR_UNITS)
	k_rta = constants.Kdissociation_charged_trna_ribosome.asNumber(MICROMOLAR_UNITS)
	k_rtf = constants.Kdissociation_uncharged_trna_ribosome.asNumber(MICROMOLAR_UNITS)
	rib_elong_rate = 22

	# Initialize to approximate charged levels to avoid dividing by 0
	charged_fraction = 0
	charged_trna_conc = total_trna * charged_fraction
	uncharged_trna_conc = total_trna * (1 - charged_fraction)

	# Solve to steady state with short time steps
	t = 0
	dt = 0.001
	diff = 1
	while diff > 1e-3:
		v_charging = (k_s * synthetase_conc * uncharged_trna_conc * aa_conc / (KM_aa * KM_tf
			* (1 + uncharged_trna_conc / KM_tf + aa_conc / KM_aa
			+ uncharged_trna_conc * aa_conc / KM_tf / KM_aa)))
		numerator_ribosome = 1 + np.sum(f * (k_rta / charged_trna_conc
			+ uncharged_trna_conc / charged_trna_conc * k_rta / k_rtf))
		v_rib = rib_elong_rate * ribosome_conc / numerator_ribosome

		# Handle case when f is 0 and charged_trna_conc is 0
		if not np.isfinite(v_rib):
			v_rib = 0

		delta_conc = (v_charging - v_rib * f) * dt
		uncharged_trna_conc -= delta_conc
		charged_trna_conc += delta_conc
		diff = np.linalg.norm(delta_conc)

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
	k_rela = constants.k_RelA_ppGpp_synthesis.asNumber(1 / units.s)
	KD_rela = constants.KD_RelA_ribosome.asNumber(MICROMOLAR_UNITS)
	k_spot_syn = constants.k_SpoT_ppGpp_synthesis.asNumber(MICROMOLAR_UNITS / units.s)
	k_spot_deg = constants.k_SpoT_ppGpp_degradation.asNumber(1 / units.s)
	k_rta = constants.Kdissociation_charged_trna_ribosome.asNumber(MICROMOLAR_UNITS)
	k_rtf = constants.Kdissociation_uncharged_trna_ribosome.asNumber(MICROMOLAR_UNITS)

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

def sensitivity(sim_data, cell_specs, conditions):
	'''
	Performs sensitivity analysis for each of the parameters.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		cell_specs (dict): information about each condition that was fit
		conditions (list[str]): set of conditions to test
			(eg. ['basal', 'with_aa', 'no_oxygen'])
	'''

	for param in PARAMS:
		print('For {}:'.format(param))
		original_value = getattr(sim_data.constants, param)

		for magnitude in [0.1, 0.2, 0.5, 0.75, 0.9, 1.1, 1.5, 2, 5, 10]:
			setattr(sim_data.constants, param, original_value * magnitude)
			error = [
				main(sim_data, cell_specs, [condition], verbose=False)
				for condition in conditions
				]

			error_str = '  '.join(format(e, '.2f') for e in error)
			print('\t x{} error: {}'.format(magnitude, error_str))

		setattr(sim_data.constants, param, original_value)

def coordinate_descent(sim_data, cell_specs, conditions):
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

	Returns:
		Constants object: class of all constants values from sim_data
		float: objective value reached
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
			high_objective = main(sim_data, cell_specs, conditions, verbose=False)

			# Change low
			setattr(sim_data.constants, param, original_value * (1 - delta))
			low_objective = main(sim_data, cell_specs, conditions, verbose=False)

			# Update to new parameter value based on lowest error
			if low_objective < objective and low_objective < high_objective:
				delta_objective = objective - low_objective
				propensity[idx] = max(delta_objective, propensity.min())
				objective = low_objective
				status = 'decreased'
			elif high_objective < objective:
				delta_objective = objective - high_objective
				propensity[idx] = max(delta_objective, propensity.min())
				objective = high_objective
				setattr(sim_data.constants, param, original_value * (1 + delta))
				status = 'increased'
			else:
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

			it += 1
			print('{} {}: {:.3f}'.format(param, status, objective))

			# Early break condition if nothing changing
			if n_constants > max_constants:
				break
	# Allow early interrupt
	except KeyboardInterrupt:
		pass

	for param in PARAMS:
		print('{}: {}'.format(param, getattr(sim_data.constants, param)))
	main(sim_data, cell_specs, conditions)

	return sim_data.constants, objective

def main(sim_data, cell_specs, conditions, verbose=True):
	'''
	Main function to perform analysis.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		cell_specs (dict): information about each condition that was fit
		conditions (list[str]): set of conditions to test
			(eg. ['basal', 'with_aa', 'no_oxygen'])
		verbose (bool): if True, prints to command line

	Returns:
		float: objective value for all conditions
	'''

	objective = 0

	constants = sim_data.constants
	is_rrna = sim_data.process.transcription.rnaData['isRRna']

	# Calculate ppGpp concentrations for each condition
	for condition in conditions:
		# Get condition specific values
		bulk_container = cell_specs[condition]['bulkAverageContainer']
		synth_prob = cell_specs[condition]['synthProb']
		doubling_time = sim_data.conditionToDoublingTime[condition]
		nutrients = sim_data.conditions[condition]['nutrients']

		# Get current state
		rela, total_trnas, ribosomes, synthetases, aas, rnaps = get_concentrations(
			sim_data, bulk_container, doubling_time)
		f = get_aa_fraction(sim_data, bulk_container)
		rnap_activation_rate = get_rnap_activation(
			sim_data, bulk_container, doubling_time, synth_prob)
		expected_ppgpp = get_expected_ppgpp(sim_data, doubling_time)
		expected_v_rib = get_expected_v_rib(sim_data, doubling_time, nutrients)

		# Calculate values from current state
		charged_trna, uncharged_trna, v_rib = charge_trna(
			total_trnas, synthetases, aas, ribosomes, f, constants)
		ppgpp = create_ppgpp(rela, charged_trna, uncharged_trna, ribosomes, f, constants)
		rrna_synth_prob = regulate_rrna_expression(ppgpp, rnaps, rnap_activation_rate, constants)

		# Print current state
		if verbose:
			summarize_state(condition, rela, charged_trna, uncharged_trna,
				ribosomes, synthetases, aas, ppgpp, f)
			error_summary(rrna_synth_prob, synth_prob[is_rrna], ppgpp,
				expected_ppgpp, v_rib, expected_v_rib)

		objective += get_objective_value(rrna_synth_prob, synth_prob[is_rrna], ppgpp,
			expected_ppgpp, v_rib, expected_v_rib)

	return objective

def parse_args():
	'''
	Parses arguments from the command line.

	Returns:
		ArgumentParser namespace: values of variables parsed from the command line
	'''

	default_sim_data = os.path.join(FILE_LOCATION, 'sim_data.cp')
	default_cell_specs = os.path.join(FILE_LOCATION, 'cell_specs.cp')
	default_seeds = 1
	default_output = os.path.join(OUTPUT_DIR, 'sgd.tsv')

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
	parser.add_argument('-o', '--output',
		default=default_output,
		help='Output file name for sgd, (default: {})'.format(default_output))

	return parser.parse_args()


if __name__ == '__main__':
	start = time.time()

	args = parse_args()

	# Check args
	if args.sensitivity and args.sgd:
		raise Exception('Cannot do both sensitivity (--sensitivity) and'
			' stochastic gradient descent (--sgd) at the same time.')

	# Load necessary files
	with open(args.sim_data) as f:
		sim_data = cPickle.load(f)
	with open(args.cell_specs) as f:
		cell_specs = cPickle.load(f)

	# Specify conditions
	conditions = ['basal', 'with_aa', 'no_oxygen']
	if args.condition is not None:
		conditions = [conditions[args.condition]]

	# Update constants that are not yet in sim_data
	# TODO: add to sim_data.constants in a flat file
	constants = sim_data.constants
	constants.rrn_vmax = 2000  # 1/s, Bosdriesz
	constants.KI_ppgpp_rnap = 1  # uM, Bosdriesz
	constants.KM_rrn_rnap = 20  # uM, Bosdriesz

	# Perform desired analysis
	if args.sensitivity:
		sensitivity(sim_data, cell_specs, conditions)
	elif args.sgd:
		with open(args.output, 'w') as f:
			csv_writer = csv.writer(f, delimiter='\t')
			csv_writer.writerow(['Seed', 'Objective'] + PARAMS)
			csv_writer.writerow(['Original', ''] + get_growth_constants(sim_data.constants))

			original_constants = {param: getattr(sim_data.constants, param) for param in PARAMS}

			for seed in np.random.randint(0, 10000, args.seeds):
				print('Seed: {}'.format(seed))
				np.random.seed(seed)
				constants, objective = coordinate_descent(sim_data, cell_specs, conditions)
				csv_writer.writerow([seed, objective] + get_growth_constants(constants))

				# Reset constants to original values
				for param, value in original_constants.items():
					setattr(sim_data.constants, param, value)
	else:
		main(sim_data, cell_specs, conditions)

	print('Completed in {:.2f} min'.format((time.time() - start) / 60))
