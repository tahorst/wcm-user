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
	## From growthRateDependentParameters.tsv
	ppgpp = np.array([0.316, 0.219, 0.127, 0.0866, 0.0578])
	## From dryMassComposition.tsv
	rna = np.array([0.135135135, 0.151162791, 0.177829099, 0.205928237, 0.243930636])

	a1s, a2s, kms = sp.symbols('a1 a2 km')

	relation = a1s*(1 - ppgpp/(kms + ppgpp)) + a2s*ppgpp/(kms + ppgpp) - rna
	J = relation.dot(relation)

	# Significantly faster than leaving in symbolic form and using subs at each iteration
	# 0.006 sec vs 38 sec for 1000 iterations
	dJda1 = sp.lambdify((a1s, a2s, kms), J.diff(a1s))
	dJda2 = sp.lambdify((a1s, a2s, kms), J.diff(a2s))
	dJdkm = sp.lambdify((a1s, a2s, kms), J.diff(kms))
	J = sp.lambdify((a1s, a2s, kms), J)

	# Initial parameters
	a1 = 0.5
	a2 = 0.1
	km = 0.025
	step_size = 0.001

	obj = J(a1, a2, km)
	step = 0
	while obj > 5e-6:
		a1 -= dJda1(a1, a2, km) * step_size
		a2 -= dJda2(a1, a2, km) * step_size
		km -= dJdkm(a1, a2, km) * step_size

		obj = J(a1, a2, km)

		step += 1
		if step % 100000 == 0:
			print(obj)
			print(a1, a2, km)

	rna_fit = a1*(1 - ppgpp/(km + ppgpp)) + a2*ppgpp/(km + ppgpp)

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
