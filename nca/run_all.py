#! /usr/bin/env python

"""
Script to run all options and models for NCA.  Need to specify options and
models of interest with the OPTIONS and MODELS variables.  Saves log files
to log/ and uses the default output location for NCA results.
"""

import multiprocessing
import os
import subprocess
from typing import List


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
LOG_DIR = os.path.join(FILE_LOCATION, 'log')
os.makedirs(LOG_DIR, exist_ok=True)

OPTIONS = [
	'average-seq',
	'linear',
	'global-expression',
	'split',
	'sigma-factors',
	]
MODELS = [
	'robust_nca',
	'constrained_nca',
	]


def solve_nca(model: str, options: List[str]):
	label = '-'.join([model] + [option.split('-')[0] for option in options])
	flags = ' '.join([f'--{option}' for option in options])
	cmd = f'./fold_changes.py -l {label} -m {model} {flags}'

	print(f'Running: {cmd}')
	with open(os.path.join(LOG_DIR, f'{label}.log'), 'w') as f:
		subprocess.run(cmd.split(), stdout=f)


if __name__ == '__main__':
	all_args = [
		(model, [option for j, option in enumerate(OPTIONS) if i//2**j % 2 == 0])
		for model in MODELS
		for i in range(2**len(OPTIONS))
		]

	pool = multiprocessing.Pool(8)
	results = [pool.apply_async(solve_nca, args) for args in all_args]
	pool.close()
	pool.join()
