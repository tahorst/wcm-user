#! /usr/bin/env python

"""
Compares profiles for functions using LineProfiler.

Used to compare processes and execution time with a single sim and multiple sims
running in parallel.

Sims run with debug-process-profile a6d9b206e:
	sim out/test-profiler --length-sec 5

Single run without other sims running.
Parallel run with 7 other sims running:
	for i in {1..7}; do (sim out/test-profiler/ --length-sec 1000 -v wildtype $i $i &); done
"""

import os

import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')

LINE_KEY = 'lines'
TIME_KEY = 'times'

CUTOFF = 500  # time cutoff to analyze more time intensive calls


def load_data():
	"""
	Extracts data from profile dumps from LineProfiler.
	"""

	data = {}
	for root, dirs, files in os.walk(DATA_DIR):
		for name in files:
			# Read data
			with open(os.path.join(root, name)) as f:
				lines = f.readlines()

			# Extract total time for each line
			line_num = []
			times = []
			for line in lines[8:]:
				measures = line.strip().split('\t')[0].split()
				if len(measures) > 2:
					line_num.append(int(measures[0]))
					times.append(float(measures[2]))

			# Save data
			condition, process, func = name.split('.')[0].split('-')
			if process not in data:
				data[process] = {}
			if func not in data[process]:
				data[process][func] = {}

			new = {
				LINE_KEY: np.array(line_num),
				TIME_KEY: np.array(times),
				}
			data[process][func][condition] = new

	return data

def compare_data(data, key1='single', key2='parallel'):
	"""
	Compares average and maximum relative ratio for each function.
	Positive values indicate key2 is slower than key1.
	"""

	for process, p in data.items():
		for func, f in p.items():
			time1 = f[key1][TIME_KEY]
			time2 = f[key2][TIME_KEY]
			filter = (time1 > CUTOFF) & (time2 > CUTOFF)
			filter1 = time1[filter]
			filter2 = time2[filter]
			ratio = (filter2 - filter1) / filter1

			if len(ratio):
				print('{:.2f}\t{:.2f}\t{} {}'.format(ratio.max(), ratio.mean(), process, func))

if __name__ == '__main__':
	data = load_data()
	compare_data(data)
