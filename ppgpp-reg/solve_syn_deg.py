#! /usr/bin/env python

"""
Use ppGpp concentrations and simulation outputs to solve for ppGpp synthesis
and degradation parameters.

sim_data and simulation output data comes from 191113-ppgpp-copy-num ee3a481bd.
"""

from __future__ import division

import os

import numpy as np
import sympy as sp


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUT_DIR):
	os.mkdir(OUT_DIR)


if __name__ == '__main__':
	# Data for different conditions (no_oxygen, basal, with_aa)
	## From sim_data.growthRateParameters.getppGppConc (100, 44, 25 min)
	ppgpp = np.array([104.3, 47.6, 20.8])  # uM
	## Rough averages from 3 seeds for each condition within ppgpp_metabolite_changes function in polypeptide_elongation
	rela_spot = np.array([1.48, 1.37, 0.61])  # unitless
	trna_uncharged = np.array([8.1, 11.5, 12.4])  # uM
	trna_charged = np.array([208, 250, 348])  # uM
	ribosome_bound = np.array([0.054, 0.030, 0.003])  # uM
	## Steady state approximation (ppGpp = ksyn / kdeg) from Murray and Bremer. JMB. 1996.
	kspot_kdeg = 2.6  # uM

	# Set initial conditions from other measurements
	k_rela = 75  # 1/s
	k_deg = 0.231  # 1/s
	initial_k = k_rela / k_deg
	initial_KI = 50.  # uM
	initial_KD = 0.26  # uM

	# Create objective (steady state should be equal to 0 in all 3 conditions)
	# Add additional targets for measured data for KD and k
	k, KI, KD = sp.symbols('k KI KD')
	ss = k * rela_spot * ribosome_bound / (KD + ribosome_bound) + kspot_kdeg - ppgpp * KI / (KI + trna_uncharged)
	J = sp.sqrt(ss.dot(ss)) + (KD - initial_KD)**2 + (k - initial_k)**2

	# Find derivatives and lambdify
	dJdk = sp.lambdify((k, KI, KD), J.diff(k))
	dJdKI = sp.lambdify((k, KI, KD), J.diff(KI))
	dJdKD = sp.lambdify((k, KI, KD), J.diff(KD))
	J = sp.lambdify((k, KI, KD), J)

	# Initial parameters
	k = initial_k
	KI = initial_KI
	KD = initial_KD
	step_size = 0.01

	obj = J(k, KI, KD)
	old_obj = 10000
	step = 0
	tol = 1e-6
	rel_tol = 1e-9
	while obj > tol and 1 - obj / old_obj > rel_tol:
		k -= dJdk(k, KI, KD) * step_size
		KI -= dJdKI(k, KI, KD) * step_size
		KD -= dJdKD(k, KI, KD) * step_size * 0.01

		old_obj = obj
		obj = J(k, KI, KD)

		if step % 10000 == 0:
			print(k * rela_spot * ribosome_bound / (KD + ribosome_bound) + kspot_kdeg - ppgpp * KI / (KI + trna_uncharged))
			print(obj)
			print(k, KI, KD)
		step += 1

	print('k: {:.1f}'.format(k))
	print('KI: {:.3f}'.format(KI))
	print('KD: {:.3f}'.format(KD))

	import ipdb; ipdb.set_trace()
