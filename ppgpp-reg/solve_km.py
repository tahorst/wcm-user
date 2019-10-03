#! /usr/bin/env python

"""
Use ppGpp concentrations and RNA mass fractions to solve for a KM value
for ppGpp binding to RNAP.
"""

from __future__ import division

import numpy as np


if __name__ == '__main__':
	# Data for different doubling times (100, 60, 40, 30, 24 min)
	## From growthRateDependentParameters.tsv
	ppgpp = np.array([0.316, 0.219, 0.127, 0.0866, 0.0578])
	## From dryMassComposition.tsv
	rna = np.array([0.135135135, 0.151162791, 0.177829099, 0.205928237, 0.243930636])

	# Derived values
	n_times = len(rna)
	p_sum = ppgpp.sum()
	r_sum = rna.sum()
	comb = ppgpp * rna
	c_sum = comb.sum()

	# Initial parameters
	a1 = 0.3
	a2 = 0.03
	km = 0.2
	step_size = 0.01

	J = np.linalg.norm(a1*km*np.ones(n_times) + a2*ppgpp - km*rna - comb)
	while J > 1e-5:
		dJda1 = 2*a1*km**2 + 2*a2*km*p_sum - 2*km**2*r_sum - 2*km*c_sum
		dJda2 = 2*a2*ppgpp.dot(ppgpp) + 2*a1*km*p_sum - 2*km*ppgpp.dot(rna) - 2*ppgpp.dot(comb)
		dJdkm = 2*a1**2*km + 2*a1*a2*p_sum - 4*a1*km*r_sum - 2*a1*c_sum - 2*a2*ppgpp.dot(rna) + 2*km*rna.dot(rna) + 2*rna.dot(comb)

		a1 -= dJda1 * step_size
		a2 -= dJda2 * step_size
		km -= dJdkm * step_size

		J = np.linalg.norm(a1*km*np.ones(n_times) + a2*ppgpp - km*rna - comb)
		new_rna = a1*(1-ppgpp/(km+ppgpp)) + a2*ppgpp/(km+ppgpp)
		print(new_rna)
		print(J)
