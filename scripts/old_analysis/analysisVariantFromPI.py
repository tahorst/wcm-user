# runs all variant analysis plots for a given sim with data from one location to another
# only need parent directory in wcEcoli_out as argument

import cPickle
import time
import os
import sys

import models.ecoli.analysis.variant
import importlib
import multiprocessing as mp

SET = "SET_A"
BASE_DIR = os.path.join("/scratch/PI/mcovert/wc_ecoli/paper/", SET)
OUTPUT = os.path.join("/scratch/users/thorst/wcEcoli_out/", SET)

for folder in sorted(os.listdir(BASE_DIR), reverse = True):
	if folder[0].isdigit():
		pathIn = os.path.join(BASE_DIR, folder)
		pathOut = os.path.join(OUTPUT, folder)
		break

for folder in os.listdir(pathIn):
	if "_" in folder:
		variant = folder

inputDirectory = os.path.join(pathIn)
outputDir = os.path.join(pathOut, "plotOut")
validationData = os.path.join(pathIn, "kb/validationData.cPickle")
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