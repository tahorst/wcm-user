'''
Script to test this line of parca:
	synthProb = normalize(rnaLossRate.asNumber(1 / units.min))

Assumption at steady state:
	dm/dt = exp - loss = 0
	synth_prob = exp.normalize()

Expected length and time to produce a transcript should account for relationship
between expression and synth_prob but results showed that average counts are
the same even with different lengths and same prob.
'''

import numpy as np

rna_lengths = np.array([1, 5])
synth_prob = np.array([0.5, 0.5])
deg_rate = np.array([0.1, 0.1])
n_rna = 2
total_rnap = 2
rnap_positions = np.zeros(total_rnap)
rnap_bound = -np.ones(total_rnap)
rnap_rate = 1

t_limit = 100000
dt = 1

rna_counts = np.zeros(n_rna)
total = np.zeros(n_rna)

for t in np.arange(0, t_limit, dt):
	for i, (pos, bound) in enumerate(zip(rnap_positions, rnap_bound)):
		bound = int(bound)
		if bound >= 0:
			pos += rnap_rate * dt
		else:
			bound = np.random.choice(n_rna, p=synth_prob)
			pos = 1

		if pos >= rna_lengths[bound]:
			rna_counts[bound] += 1
			bound = -1

		rnap_positions[i] = pos
		rnap_bound[i] = bound

	rna_counts -= rna_counts * deg_rate
	total += rna_counts
	print(total / t)

import ipdb; ipdb.set_trace()
