'''
Test for memory issues with parameter sensitivity analysis when doing parallel
analysis with 50000 variants (model/ecoli/analysis/variant/param_sensitivity.py).

Running script can freeze computer with quick allocation if SAMPLES set too high
unlike analysis plot which gets killed with a memory error.
'''

from future_builtins import zip

import gc
import multiprocessing as mp
import operator
import time

import numpy as np


SAMPLES = 5000


def get_results((x,)):
	return np.vstack((
		np.random.randint(0, 10, 20000),
		np.random.randint(0, 10, 20000),
		np.random.randint(0, 10, 20000),
		np.random.randint(0, 10, 20000),
		np.random.randint(0, 10, 20000),
		np.random.randint(0, 10, 20000),
		))

def mapped():
	for i in range(SAMPLES):
		yield get_results((i,))

def reduce_results(total, new):
	return (t + n for t,n in zip(total, new))

def reduce_results_gc(total, new):
	gc.collect()
	return (t + n for t,n in zip(total, new))

# Functions to test
def serial():
	reduce(reduce_results, mapped())

def serial_gc():
	reduce(reduce_results_gc, mapped())

def parallel():
	pool = mp.Pool(8)
	results = pool.imap_unordered(get_results, zip([None] * SAMPLES,))
	reduce(reduce_results, results)
	pool.close()
	pool.join()

def parallel_gc():
	pool = mp.Pool(8)
	results = pool.imap_unordered(get_results, zip([None] * SAMPLES,))
	reduce(reduce_results_gc, results)
	pool.close()
	pool.join()

def parallel_no_tuple():
	pool = mp.Pool(8)
	results = pool.imap_unordered(get_results, zip([None] * SAMPLES,))
	reduce(operator.add, results)
	pool.close()
	pool.join()

# Timing function
def time_fun(fun, desc):
	print('\n{}'.format(desc))
	st = time.time()
	fun()
	print('Completed: {:.1f} sec'.format(time.time() - st))

# Times with 5000 samples, tested individually to prevent memory buildup
time_fun(serial, 'Serial no gc')  # 8.7 sec
time_fun(serial_gc, 'Serial with gc')  # 22.1 sec - doesn't free memory
time_fun(parallel, 'Parallel no gc')  # 9.4 sec
time_fun(parallel_gc, 'Parallel with gc')  # 33.9 sec - much slower, doesn't free mem
time_fun(parallel_no_tuple, 'Parallel no tuple')  # 10.5 sec - Fixes memory build up
