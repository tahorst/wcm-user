#! /usr/bin/env python
"""
Used to simulate ppGpp concentrations given a time series of wcm states.

Assumes no impact from ppGpp synthesis/degradation on the rest of the simulation
but it will require varying amounts of ATP and GDP if the parameters change.
"""

import json
import os

import numpy as np
from typing import Dict, List, Tuple


FILE_LOCATION = os.path.realpath(os.path.dirname(__file__))
DATA_FILE = os.path.join(FILE_LOCATION, 'cell_states.json')

PARAMS = {
	'KD_RelA': 0.26,
	'KI_SpoT': 50,
	'k_RelA': 75,
	'k_SpoT_deg': 0.231,
	'k_SpoT_syn': 1,
	}


def load_data(filename):
	# type: (str) -> Tuple[np.ndarray[float]]
	"""
	Load concentrations of interest from whole cell model simulations.

	Args:
		filename: path to data file to load

	Returns:
		3D arrays (n cells, m time points, o concentrations) of data for each
		condition (basal, anaerobic, with_aa)
	"""

	with open(filename) as f:
		data = json.load(f)

	return data

def ppgpp_metabolite_changes(params, ppgpp_conc, rela_conc, spot_conc, ribosomes_bound, uncharged_trna_conc):
	# type: (Dict[str, float], float, float, float, float, float) -> float
	"""
	Calculates the change ppGpp concentration for a given state of the cell.
	All concentrations are in uM.

	Args:
		params: parameters for ppGpp synthesis and degradation
		ppgpp_conc: concentration of ppGpp
		ribosomes_bound: concentration of ribosomes bound to uncharged tRNA
		uncharged_trna_conc: total concentration of all uncharged tRNA
		rela_conc: concentration of RelA
		spot_conc: concentration of SpoT

	Returns:
		change in ppGpp concentration

	TODO:
		add time step? - currently assumes 1 sec step
	"""

	frac_rela = 1 / (1 + params['KD_RelA'] / ribosomes_bound)
	v_rela = params['k_RelA'] * rela_conc * frac_rela
	v_syn = v_rela + params['k_SpoT_syn'] * spot_conc
	v_deg = params['k_SpoT_deg'] * spot_conc * ppgpp_conc / (1 + uncharged_trna_conc / params['KI_SpoT'])

	return v_syn - v_deg

def simulate(data, params):
	# type: (np.ndarray[float], Dict[str, float]) -> List[np.ndarray[float]]
	"""
	Simulate all cells in a single condition.

	Args:
		data: 3D arrays (n cells, m time points, o concentrations) of data for each
			condition (basal, anaerobic, with_aa)
		params: parameters for ppGpp synthesis and degradation

	Returns:
		2D array (n cells, m time points) of ppGpp concentration traces
	"""

	ppgpp = []
	for cell in data:
		ppgpp = 47.5  # TODO: vary for each condition
		series = [ppgpp]

		for timepoint in cell.T:
			ppgpp += ppgpp_metabolite_changes(params, ppgpp, *timepoint)
			series += [ppgpp]

		ppgpp.append(np.array(series))

	return ppgpp

def simulate_all(data, params):
	# type: (Tuple[np.ndarray[float]], Dict[str, float]) -> None
	"""
	Simulate all conditions
	"""

	for d in data:
		ppgpp = simulate(d, params)



if __name__ == '__main__':
	data = load_data(DATA_FILE)
	simulate_all(data, PARAMS)
