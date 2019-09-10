#! /usr/bin/env python
"""
Explore parameters for ppGpp synthesis in the model from tsv output from
ppgpp_parameter_optimization.
"""

import csv
import os

import numpy as np


FILE_LOCATION = os.path.realpath(os.path.dirname(__file__))
DATA_FILE = os.path.join(FILE_LOCATION, 'new_ppgpp_parameter_optimization.tsv')
N_PARAMS = 5
N_CONDITIONS = 3

with open(DATA_FILE) as f:
	reader = csv.reader(f, delimiter='\t')
	headers = reader.next()
	data = np.array(list(reader), float)

param_names = headers[:N_PARAMS]
params = {param: data[:, i] for i, param in enumerate(param_names)}
param_min = {param: v.min(axis=0) for param, v in params.items()}
param_max = {param: v.max(axis=0) for param, v in params.items()}
ppgpp = data[:, N_PARAMS:N_PARAMS+N_CONDITIONS]
std = data[:, N_PARAMS+N_CONDITIONS:]


# Summarize parameter effects
## TODO: plot average ppGpp and stdev for each param value
for param, param_vals in params.items():
	print('\nParam: {}'.format(param))

	mask = param_vals == param_max[param]
	print('\tHigh ({}): {:.1f}  {:.1f}  {:.1f}'.format(
		param_max[param],
		*ppgpp[mask, :].mean(axis=0)
		))
	normalized_std = std[mask, :].mean(axis=0) / ppgpp[mask, :].mean(axis=0)
	print('\tStd: {:.1f}  {:.1f}  {:.1f}'.format(*normalized_std))

	mask = param_vals == param_min[param]
	print('\tLow ({}): {:.1f}  {:.1f}  {:.1f}'.format(
		param_min[param],
		*ppgpp[mask, :].mean(axis=0)
		))
	normalized_std = std[mask, :].mean(axis=0) / ppgpp[mask, :].mean(axis=0)
	print('\tStd: {:.1f}  {:.1f}  {:.1f}'.format(*normalized_std))


# Filter parameters of interest
## Ranges for each condition
lows = np.array([40, 90, 15])
highs = np.array([55, 110, 25])
max_normalized_std = np.array([1., 1., 1.])

mask = np.ones(ppgpp.shape[0], bool)
for low, high, norm_std, p, s in zip(lows, highs, max_normalized_std, ppgpp.T, std.T):
	mask = mask & (p > low) & (p < high) & (s / p < norm_std)
filtered_params = {param: v[mask] for param, v in params.items()}

import ipdb; ipdb.set_trace()

