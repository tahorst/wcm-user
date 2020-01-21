#! /usr/bin/env python

"""
Test read and write of sim_data objects with and without compression.
No compression appears to be best although file size can be reduced ~6-8x

Must supply path to a sim_data cPickle file as the first command line arg.

Sherlock results:
Write times for no compression:
[1.784773826599121, 1.921402931213379, 2.0703601837158203]
Read times for no compression:
[2.98009991645813, 1.109165906906128, 1.0727088451385498]
Write times for gzip compression:
[1.7210698127746582, 2.037745952606201, 2.0388450622558594]
Read times for gzip compression:
[3.704223871231079, 3.5736050605773926, 3.5095999240875244]
Write times for bz2 compression:
[5.744052171707153, 5.709091901779175, 5.744783163070679]
Read times for bz2 compression:
[3.1027019023895264, 2.7685089111328125, 2.869654893875122]

Local results:
Write times for no compression:
[0.2988240718841553, 0.4597010612487793, 0.5387649536132812]
Read times for no compression:
[0.9482810497283936, 0.7746589183807373, 0.8097898960113525]
Write times for gzip compression:
[1.154386043548584, 1.2513411045074463, 1.183440923690796]
Read times for gzip compression:
[2.7091009616851807, 2.6934940814971924, 2.751133918762207]
Write times for bz2 compression:
[4.494904041290283, 4.601887941360474, 4.744673013687134]
Read times for bz2 compression:
[2.6113109588623047, 2.74900484085083, 2.8047730922698975]
"""

from __future__ import absolute_import, division, print_function

import bz2
import cPickle
import gzip
import sys
import time


if __name__ == '__main__':
	path = sys.argv[1]

	# Load original data
	with open(path) as f:
		sim_data = cPickle.load(f)

	n = 3

	# No compress write
	print('Write times for no compression:')
	file_no_compress = path + '.none'
	times = []
	for i in range(n):
		start = time.time()
		with open(file_no_compress, 'w') as f:
			data = cPickle.dump(sim_data, f, protocol=cPickle.HIGHEST_PROTOCOL)
		end = time.time()
		times.append(end - start)
	print(times)

	# No compress read
	print('Read times for no compression:')
	times = []
	for i in range(n):
		start = time.time()
		with open(file_no_compress) as f:
			data = cPickle.load(f)
		end = time.time()
		times.append(end - start)
	print(times)

	# gzip write
	print('Write times for gzip compression:')
	file_gzip = path + '.gz'
	times = []
	for i in range(n):
		start = time.time()
		with gzip.GzipFile(file_gzip, 'w', compresslevel=1) as f:
			data = cPickle.dump(sim_data, f, protocol=cPickle.HIGHEST_PROTOCOL)
		end = time.time()
		times.append(end - start)
	print(times)

	# gzip read
	print('Read times for gzip compression:')
	times = []
	for i in range(n):
		start = time.time()
		with gzip.GzipFile(file_gzip) as f:
			data = cPickle.load(f)
		end = time.time()
		times.append(end - start)
	print(times)

	# bz2 write
	print('Write times for bz2 compression:')
	file_bz2 = path + '.bz'
	times = []
	for i in range(n):
		start = time.time()
		with bz2.BZ2File(file_bz2, 'w', compresslevel=1) as f:
			data = cPickle.dump(sim_data, f, protocol=cPickle.HIGHEST_PROTOCOL)
		end = time.time()
		times.append(end - start)
	print(times)

	# bz2 read
	print('Read times for bz2 compression:')
	times = []
	for i in range(n):
		start = time.time()
		with bz2.BZ2File(file_bz2) as f:
			data = cPickle.load(f)
		end = time.time()
		times.append(end - start)
	print(times)
