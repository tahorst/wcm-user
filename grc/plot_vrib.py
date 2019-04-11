'''
Plots to compare impacts from different version of v_rib calculation.

Uses two amino acids - one is varied and the other can be considered aggregation of
all other amino acids that are not varied.

Varies two parameters:
	f_i: fraction of AA
	sigma_i: fraction saturated with charged tRNA

Compares two v_rib equations:
	From Bosdriesz Eq S7:
		v_rib = k_rib * r * 1 / sum_i(f_i / sigma_i)

	Updated equation (can give higher v_rib and better AA condition growth):
		v_rib = k_riv * r * sum_i(f_i * sigma_i)

Output:
	two plots to out/v_rib_*.png

Results:
- For f: Bosdriesz requires higher charged fraction - could be good approximation
for one ribosome where future elongation depends on the limiting AA but less
appropriate in the whole cell model where f is calculated from many ribosomes
that are independent of each other and only short sequences for a timestep are
dependent on each other.
- For sigma: Bosdriesz equation performs poorly near 0 - if one amino acid is not
charged, then no elongation can occur.
'''

from __future__ import division

import os

import matplotlib.pyplot as plt
import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_LOCATION = os.path.join(FILE_LOCATION, 'out')


if __name__ == '__main__':
	k_rib = 22

	# Vary f
	f1 = np.linspace(0, 1, 101)
	f = np.array([f1, 1-f1]).T
	sigma = np.array([1, 0.01])

	v1 = k_rib / np.sum(f / sigma, axis=1)
	v2 = k_rib * np.sum(f * sigma, axis=1)

	plt.figure()
	plt.plot(f1, v1)
	plt.plot(f1, v2)
	plt.xlabel('Fraction of High Saturated Charged tRNA Species')
	plt.ylabel('v_rib')
	plt.legend(['Bosdriesz', 'Updated'])
	plt.tight_layout()
	plt.savefig(os.path.join(OUTPUT_LOCATION, 'v_rib_f.png'))

	# Vary sigma
	sigma1 = np.linspace(0, 1, 101)
	sigma = np.array([sigma1, np.ones(101)]).T
	f = np.array([0.05, 0.95])

	v1 = k_rib / np.sum(f / sigma, axis=1)
	v2 = k_rib * np.sum(f * sigma, axis=1)

	plt.figure()
	plt.plot(sigma1, v1)
	plt.plot(sigma1, v2)
	plt.xlabel('Fraction Saturated with Charged tRNA for one AA')
	plt.ylabel('v_rib')
	plt.legend(['Bosdriesz', 'Updated'])
	plt.tight_layout()
	plt.savefig(os.path.join(OUTPUT_LOCATION, 'v_rib_sigma.png'))
