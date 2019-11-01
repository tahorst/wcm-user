#! /usr/bin/env python

"""
Use ppGpp concentrations and RNA mass fractions to solve for a KM value
for ppGpp binding to RNAP.
"""

from __future__ import division

import os

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUT_DIR):
	os.mkdir(OUT_DIR)


if __name__ == '__main__':
	# Data for different doubling times (100, 60, 40, 30, 24 min)
	## From growthRateDependentParameters.tsv (pmol / ug)
	ppgpp = np.array([0.316, 0.219, 0.127, 0.0866, 0.0578])
	## From dryMassComposition.tsv
	rna = np.array([0.135135135, 0.151162791, 0.177829099, 0.205928237, 0.243930636])
	## From average fold change for negative regulation in sim_data.processs.transcription.ppgpp_fold_changes (511dcff41)
	fc_target = 2**-1.4037

	a1s, a2s, kms = sp.symbols('a1 a2 km')

	# Use sp.exp to prevent negative parameter values, also improves stability for larger step size
	relation = sp.exp(a1s)*(1 - ppgpp/(sp.exp(kms) + ppgpp)) + sp.exp(a2s)*ppgpp/(sp.exp(kms) + ppgpp) - rna
	f_low = ppgpp[-1] / (sp.exp(kms) + ppgpp[-1])
	fc = sp.exp(a2s) / (sp.exp(a1s)*(1 - f_low) + sp.exp(a2s)*f_low) - fc_target
	J = relation.dot(relation) + fc**2

	# Significantly faster than leaving in symbolic form and using subs at each iteration
	# 0.006 sec vs 38 sec for 1000 iterations
	dJda1 = sp.lambdify((a1s, a2s, kms), J.diff(a1s))
	dJda2 = sp.lambdify((a1s, a2s, kms), J.diff(a2s))
	dJdkm = sp.lambdify((a1s, a2s, kms), J.diff(kms))
	J = sp.lambdify((a1s, a2s, kms), J)

	# Initial parameters
	a1 = np.log(0.5)
	a2 = np.log(0.1)
	km = np.log(0.04)
	step_size = 0.1

	obj = J(a1, a2, km)
	old_obj = 100
	step = 0
	tol = 1e-6
	rel_tol = 1e-9
	while obj > tol and 1 - obj / old_obj > rel_tol:
		a1 -= dJda1(a1, a2, km) * step_size
		a2 -= dJda2(a1, a2, km) * step_size
		km -= dJdkm(a1, a2, km) * step_size

		old_obj = obj
		obj = J(a1, a2, km)

		step += 1
		if step % 1000 == 0:
			print(obj)
			print(np.exp(a1), np.exp(a2), np.exp(km))

	a1 = np.exp(a1)
	a2 = np.exp(a2)
	km = np.exp(km)
	rna_fit = a1*(1 - ppgpp/(km + ppgpp)) + a2*ppgpp/(km + ppgpp)
	f_low = ppgpp[-1] / (km + ppgpp[-1])
	fc = np.log2(a2 / (a1*(1 - f_low) + a2*f_low))
	print('fc: {}'.format(fc))

	# Plot results of fit to data
	plt.figure()

	plt.plot(ppgpp, rna, 'o')
	plt.plot(ppgpp, rna_fit, 'x')
	plt.axvline(km, color='k', linestyle='--')

	plt.xlabel('ppGpp (pmol / ng)')
	plt.ylabel('RNA Mass Fraction')
	plt.legend(['Measured', 'Fit', 'KM'])
	plt.title('a1: {:.3f}, a2: {:.3f}, KM: {:.3f}'.format(a1, a2, km))

	plt.savefig(os.path.join(OUT_DIR, 'ppgpp-rna.png'))
