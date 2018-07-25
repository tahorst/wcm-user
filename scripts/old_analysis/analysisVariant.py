# runs all variant analysis plots for a given sim with data from one location to another
# only need parent directory in wcEcoli_out as argument

import cPickle
import time
import os
import sys

import models.ecoli.analysis.variant
import importlib
import multiprocessing as mp

def run_function(f, args, name):
	try:
		print "%s: Running %s" % (time.ctime(), name)
		f(*args)
	except KeyboardInterrupt:
		import sys; sys.exit(1)

BASE_DIR = "/scratch/users/thorst/wcEcoli_out/"

if len(sys.argv) < 2:
	for folder in sorted(os.listdir(BASE_DIR), reverse = True):
		if folder[0].isdigit():
			path = os.path.join(BASE_DIR, folder)
			break
else:
	arg = sys.argv[1]
	if arg.startswith("out/"):
		arg = arg[4:]
	path = os.path.join(BASE_DIR, arg)

inputDirectory = path
outputDir = os.path.join(path, "plotOut")
validationData = os.path.join(path, "kb/validationData.cPickle")
metadata = None

startTime = time.time()
print "%s: Running variant analysis from %s" % (time.ctime(startTime), inputDirectory)

directory = os.path.dirname(models.ecoli.analysis.variant.__file__)

# Run analysis scripts in order of modification, most recently edited first
fileList = os.listdir(directory)
fileList.sort(key=lambda x: os.stat(os.path.join(directory, x)).st_mtime, reverse=True)

if not os.path.exists(outputDir):
	os.makedirs(outputDir)

if "WC_ANALYZE_FAST" in os.environ:
	pool = mp.Pool(processes = 8)

for f in fileList:
	if f.endswith(".pyc") or f == "__init__.py":
		continue

	mod = importlib.import_module("models.ecoli.analysis.variant." + f[:-3])
	args = (
		inputDirectory,
		outputDir,
		f[:-3],
		validationData,
		metadata,
		)

	if "WC_ANALYZE_FAST" in os.environ:
		pool.apply_async(run_function, args = (mod.main, args, f))
	else:
		print "%s: Running %s" % (time.ctime(), f)
		mod.main(*args)

if "WC_ANALYZE_FAST" in os.environ:
	pool.close()
	pool.join()
timeTotal = time.time() - startTime
print "Completed multiple generation analysis in %s" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))
