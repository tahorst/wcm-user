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
import os

import numpy as np

from wholecell.utils import units


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
MICROMOLAR_UNITS = units.umol / units.L


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

def get_ribosome_counts(sim_data, bulk_container, doubling_time):
	'''
	Gets the counts of active ribosomes based on subunits, stoichiometry and
	fraction active.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		bulk_container (BulkMoleculesContainer object): counts for each molecule
		doubling_time (float with time units): expected cell doubling time

	Returns:
		float: number of active ribosomes
	'''

	complexation = sim_data.process.complexation
	molecule_ids = sim_data.moleculeIds
	active_fraction = sim_data.growthRateParameters.getFractionActiveRibosome(doubling_time)

	info_30s = complexation.getMonomers(molecule_ids.s30_fullComplex)
	info_50s = complexation.getMonomers(molecule_ids.s50_fullComplex)

	subunits_30s = info_30s['subunitIds']
	subunits_50s = info_50s['subunitIds']
	stoich_30s = info_30s['subunitStoich']
	stoich_50s = info_50s['subunitStoich']

	# Counts of each subunit limited by stoichiometry
	count_30s = (bulk_container.counts(subunits_30s) / stoich_30s).min()
	count_50s = (bulk_container.counts(subunits_50s) / stoich_50s).min()

	ribosome_counts = min(count_30s, count_50s) * active_fraction

	return ribosome_counts

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
	rna_counts = bulk_container.counts(rna_names)

	total_aa = (rna_counts * translation_efficiencies).dot(aa_content)

	return total_aa / total_aa.sum()

def charge_trna(total_trna, synthetase_conc, aa_conc, ribosome_conc, f, constants):
	'''
	Charges tRNA from the composition of the cell.

	Args:
		total_trna (ndarray[float]): concentration of all tRNA for each amino acid (in units of uM)
		synthetase_conc (ndarray[float]): concentration of synthetases for each amino acid (in units of uM)
		aa_conc (ndarray[float]): concentration of amino acids (in units of uM)
		ribosome_conc (float): concentration of active ribosomes
		f (ndarray[float]): fraction of each amino acid in sequences to be translated
		constants (class): constants from sim_data

	Returns:
		ndarray[float]: concentration of charged tRNA for each amino acid (in units of uM)
		ndarray[float]: concentration of uncharged tRNA for each amino acid (in units of uM)
	'''

	# Parameters from Bosdriesz et al
	k_s = constants.synthetase_charging_rate.asNumber(1 / units.s)
	KM_aa = constants.Km_synthetase_amino_acid.asNumber(MICROMOLAR_UNITS)
	KM_tf = constants.Km_synthetase_uncharged_trna.asNumber(MICROMOLAR_UNITS)
	k_rta = constants.Kdissociation_charged_trna_ribosome.asNumber(MICROMOLAR_UNITS)
	k_rtf = constants.Kdissociation_uncharged_trna_ribosome.asNumber(MICROMOLAR_UNITS)
	rib_elong_rate = 22

	# Initialize to approximate charged levels to avoid dividing by 0
	charged_fraction = 0.95
	charged_trna_conc = total_trna * charged_fraction
	uncharged_trna_conc = total_trna * (1 - charged_fraction)

	# Solve to steady state with short time steps
	dt = 0.0001
	diff = 1
	while diff > 1e-3:
		v_charging = (k_s * synthetase_conc * uncharged_trna_conc * aa_conc / (KM_aa * KM_tf *
			(1 + uncharged_trna_conc / KM_tf + aa_conc / KM_aa + uncharged_trna_conc * aa_conc / KM_tf / KM_aa)))
		numerator_ribosome = 1 + np.sum(f * (k_rta / charged_trna_conc + uncharged_trna_conc / charged_trna_conc * k_rta / k_rtf))
		v_rib = rib_elong_rate * ribosome_conc / numerator_ribosome

		# Handle case when f is 0 and charged_trna_conc is 0
		if not np.isfinite(v_rib):
			v_rib = 0

		delta_conc = (v_charging - v_rib * f) * dt
		uncharged_trna_conc -= delta_conc
		charged_trna_conc += delta_conc
		diff = np.linalg.norm(delta_conc)

	return charged_trna_conc, uncharged_trna_conc

def calc_ppgpp(rela_conc, charged_trna_conc, uncharged_trna_conc, ribosome_conc, f, constants):
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

	numerator_ribosome = 1 + np.sum(f * (k_rta / charged_trna_conc + uncharged_trna_conc / charged_trna_conc * k_rta / k_rtf))
	ribosome_bound_uncharged = ribosome_conc * (f * uncharged_trna_conc / charged_trna_conc * k_rta / k_rtf) / numerator_ribosome
	frac_rela = 1 / (1 + KD_rela / ribosome_bound_uncharged.sum())

	v_rela = k_rela * rela_conc * frac_rela
	v_syn = v_rela + k_spot_syn

	# Steady state solution for ppGpp when synthesis (v_syn) = degradation (k_spot_deg * ppgpp_conc)
	ppgpp_conc = v_syn / k_spot_deg

	return ppgpp_conc

def main(sim_data, cell_specs):
	'''
	Main function to perform analysis.

	Args:
		sim_data (SimulationData object): knowledgebase for a simulation
		cell_specs (dict): information about each condition that was fit
	'''

	transcription = sim_data.process.transcription
	molecule_ids = sim_data.moleculeIds
	molecule_groups = sim_data.moleculeGroups
	constants = sim_data.constants
	n_avogadro = constants.nAvogadro

	# Data structures for charging
	aa_from_synthetase = transcription.aa_from_synthetase
	aa_from_trna = transcription.aa_from_trna

	# Names of molecules associated with tRNA charging
	uncharged_trna_names = transcription.rnaData['id'][transcription.rnaData['isTRna']]
	synthetase_names = transcription.synthetase_names
	aa_names = molecule_groups.aaIDs

	# charging_stoich_matrix = transcription.charging_stoich_matrix()
	# charging_molecule_names = transcription.charging_molecules
	# charging_molecule_counts = bulk_container.counts(charging_molecule_names)

	# Calculate ppGpp concentrations for each condition
	for condition in ['basal', 'with_aa', 'no_oxygen']:
		bulk_container = cell_specs[condition]['bulkContainer']
		doubling_time = sim_data.conditionToDoublingTime[condition]

		volume = get_volume(sim_data, bulk_container)
		counts_to_micromolar = 1 / n_avogadro.asNumber(1 / units.umol) / volume

		# Concentrations for tRNA charging molecules
		rela_counts = bulk_container.count(molecule_ids.RelA)
		total_trna_counts = aa_from_trna.dot(bulk_container.counts(uncharged_trna_names))
		ribosome_counts = get_ribosome_counts(sim_data, bulk_container, doubling_time)
		synthetase_counts = aa_from_synthetase.dot(bulk_container.counts(synthetase_names))
		aa_counts = bulk_container.counts(aa_names)

		rela_conc = rela_counts * counts_to_micromolar
		total_trna_conc = total_trna_counts * counts_to_micromolar
		ribosome_conc = ribosome_counts * counts_to_micromolar
		synthetase_conc = synthetase_counts * counts_to_micromolar
		aa_conc = aa_counts * counts_to_micromolar

		f = get_aa_fraction(sim_data, bulk_container)

		charged_trna_conc, uncharged_trna_conc = charge_trna(
			total_trna_conc, synthetase_conc, aa_conc, ribosome_conc, f, constants)
		ppgpp_conc = calc_ppgpp(rela_conc, charged_trna_conc, uncharged_trna_conc, ribosome_conc, f, constants)

		import ipdb; ipdb.set_trace()
		print('ppGpp conc in {}: {:.1f}'.format(condition, ppgpp_conc))

def parse_args():
	'''
	Parses arguments from the command line.

	Returns:
		ArgumentParser namespace: values of variables parsed from the command line
	'''

	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--sim-data',
		help='Path to sim_data object',
		default=os.path.join(FILE_LOCATION, 'sim_data.cp'))
	parser.add_argument('-c', '--cell-specs',
		help='Path to cell_specs object',
		default=os.path.join(FILE_LOCATION, 'cell_specs.cp'))

	return parser.parse_args()


if __name__ == '__main__':
	args = parse_args()

	with open(args.sim_data) as f:
		sim_data = cPickle.load(f)

	with open(args.cell_specs) as f:
		cell_specs = cPickle.load(f)

	main(sim_data, cell_specs)
