'''
Plots to compare impacts from different version of v_rib calculation.

Uses two amino acids - one is varied and the other can be considered aggregation of
all other amino acids that are not varied.

Varies two parameters:
	f_i: fraction of AA
	sigma_i: fraction saturated with charged tRNA

Compares four v_rib equations:
	From Bosdriesz Eq S7:
		v_rib = k_rib * r * 1 / sum_i(f_i / sigma_i)

	Updated equation (can give higher v_rib and better AA condition growth):
		v_rib = k_rib * r * sum_i(f_i * sigma_i)

	Combination of both (Bosdriesz over each ribosome):
		v_rib = k_rib * sum_r(1 / sum_i(f_i,r / sigma_i))

	Using elongation times based on independent v_rib:
		rate_i = k_rib * sigma_i
		t_i = 1 / rate_i
		v_rib = elongations in t time (cumsum of t_i) / n_rib

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
- For combined v_rib, most similar to Bosdriesz even with a high number of independent
ribosomes that would be similar to conditions in a cell.
- Elongation time is affected by frequency of each amino acid but not by fraction
charged - matches combined v_rib for most fraction charged.
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

	n_rib = 15000
	n_seq = 25
	n_total = n_rib * n_seq
	idx = np.arange(n_total)
	f_mat = np.zeros((101, n_rib, 2))
	sequences = np.zeros((101, n_rib, n_seq))
	for i, prob in enumerate(f):
		np.random.shuffle(idx)
		n_to_change = int(n_total * prob[1])
		new_mat = np.zeros(n_total)
		new_mat[idx[:n_to_change]] = 1
		new_mat = new_mat.reshape((n_seq, n_rib))
		sequences[i] = new_mat.T

		counts = new_mat.sum(axis=0) / n_seq
		f_mat[i] = np.array([1 - counts, counts]).T

	v1 = k_rib / np.sum(f / sigma, axis=1)
	v2 = k_rib * np.sum(f * sigma, axis=1)
	v3 = k_rib * np.sum(1 / f_mat.dot(1 / sigma), axis=1) / n_rib

	# TODO: find way to vectorize to speed up
	# TODO: update f based on length that can be elongated
	v4 = np.zeros(101)
	rate = k_rib * sigma
	for v_idx, seq in enumerate(sequences):
		times = np.zeros_like(seq)
		elong_time = 1 / rate
		for i, t in enumerate(elong_time):
			times[seq == i] = t

		total_time = np.cumsum(times, axis=1)
		elongations = np.sum(total_time <= 1)
		rib_idx = 0
		rib, end = np.where(total_time > 1)
		for i, j in zip(rib, end):
			if rib_idx == i:
				rib_idx += 1
				if j == 0:
					start = 0
				else:
					start = total_time[i, j-1]
				frac = (1 - start) / (total_time[i, j] - start)
				elongations += frac

		v4[v_idx] = elongations / n_rib

	plt.figure()
	plt.plot(f1, v1)
	plt.plot(f1, v2)
	plt.plot(f1, v3)
	plt.plot(f1, v4)
	plt.xlabel('Fraction of High Saturated Charged tRNA Species')
	plt.ylabel('v_rib')
	plt.legend(['Bosdriesz', 'Updated', 'Combined', 'Actual Time'])
	plt.tight_layout()
	plt.savefig(os.path.join(OUTPUT_LOCATION, 'v_rib_f.png'))

	# Vary sigma
	sigma1 = np.linspace(0, 1, 101)
	sigma = np.array([sigma1, np.ones(101)]).T
	f = np.array([0.05, 0.95])

	np.random.shuffle(idx)
	n_to_change = int(n_total * f[1])
	new_mat = np.zeros(n_total)
	new_mat[idx[:n_to_change]] = 1
	new_mat = new_mat.reshape((n_seq, n_rib))
	counts = new_mat.sum(axis=0) / n_seq
	f_mat = np.array([1 - counts, counts]).T
	seq = new_mat.T

	v1 = k_rib / np.sum(f / sigma, axis=1)
	v2 = k_rib * np.sum(f * sigma, axis=1)
	v3 = k_rib * np.sum(1 / f_mat.dot(1 / sigma.T), axis=0) / n_rib

	# TODO: find way to vectorize to speed up
	# TODO: update f based on length that can be elongated
	v4 = np.zeros(101)
	rates = k_rib * sigma
	for v_idx, rate in enumerate(rates):
		times = np.zeros_like(seq)
		elong_time = 1 / rate
		for i, t in enumerate(elong_time):
			times[seq == i] = t

		total_time = np.cumsum(times, axis=1)
		elongations = np.sum(total_time <= 1)
		rib_idx = 0
		rib, end = np.where(total_time > 1)
		for i, j in zip(rib, end):
			if rib_idx == i:
				rib_idx += 1
				if j == 0:
					start = 0
				else:
					start = total_time[i, j-1]
				frac = (1 - start) / (total_time[i, j] - start)
				elongations += frac

		v4[v_idx] = elongations / n_rib

	plt.figure()
	plt.plot(sigma1, v1)
	plt.plot(sigma1, v2)
	plt.plot(sigma1, v3)
	plt.plot(sigma1, v4)
	plt.xlabel('Fraction Saturated with Charged tRNA for one AA')
	plt.ylabel('v_rib')
	plt.legend(['Bosdriesz', 'Updated', 'Combined', 'Actual Time'])
	plt.tight_layout()
	plt.savefig(os.path.join(OUTPUT_LOCATION, 'v_rib_sigma.png'))
