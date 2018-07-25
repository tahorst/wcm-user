import threading
import numpy as np
import time

def sumDict(mat, totals, thread):
	totals[thread] = mat.sum(axis = 1).cumsum(axis = 1)

def sum(mat, tm):
	tm += mat.sum(axis = 1).cumsum(axis = 1)

sequences = np.array(np.round(np.random.rand(21, 4000, 20)), np.bool)
count = 0
while count < 10:
	print "\n\n %i" % (count)
	count += 1

	tstart = time.time()
	nThreads = 4
	nActive = sequences.shape[1]
	totals = {}
	threads = []
	totalMonomers = np.zeros((sequences.shape[0], sequences.shape[2]))
	for thread in xrange(nThreads):
		totals[thread] = {}
		startInd = thread * nActive // nThreads
		endInd = (thread + 1) * nActive // nThreads

		t = threading.Thread(target = sumDict, args = (sequences[:, startInd:endInd, :], totals, thread,))
		# t = threading.Thread(target = sum, args = (sequences[:, startInd:endInd, :], totalMonomers,))
		threads.append(t)
		t.start()

	for thread in xrange(nThreads):
		threads[thread].join()
		if thread == 0:
			totalMonomers = totals[thread]
		else:
			totalMonomers += totals[thread]

	print "time to complete threaded: %f" % (time.time() - tstart)

	
	tstart = time.time()
	totalMonomers1 = sequences.sum(axis = 1).cumsum(axis = 1)
	print "time to complete single thread: %f" % (time.time() - tstart)

	# print totalMonomers