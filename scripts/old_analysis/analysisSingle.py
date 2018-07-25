# runs all single analysis plots for a given sim
# if no folder is given it runs for most recent in wcEcoli_out otherwise picks first arg passed

import cPickle
import time
import os
import sys

import models.ecoli.analysis.single
import importlib
import multiprocessing as mp
import traceback as tb

BASE_DIR = "/scratch/users/thorst/wcEcoli_out/"
DIRS = "000000/generation_000000/000000"

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

for folder in os.listdir(path):
	if "_" in folder:
		variant = folder

# variant = "condition_000000"

resultsDir = os.path.join(path, variant, DIRS, "simOut")
outputDir = os.path.join(path, variant, DIRS, "plotOut")
simData = os.path.join(path, variant, "kb/simData_Modified.cPickle")
validationData = os.path.join(path, "kb/validationData.cPickle")
metadata = {}

startTime = time.time()
print "%s: Running single simulation analysis in\n%s" % (time.ctime(startTime), resultsDir)

directory = os.path.dirname(models.ecoli.analysis.single.__file__)

# Run analysis scripts in order of modification, most recently edited first
fileList = os.listdir(directory)
fileList.sort(key=lambda x: os.stat(os.path.join(directory, x)).st_mtime, reverse=True)

if "WC_ANALYZE_FAST" in os.environ:
	pool = mp.Pool(processes = 8)

for f in fileList:
	if f.endswith(".pyc") or f == "__init__.py":
		continue

	mod = importlib.import_module("models.ecoli.analysis.single." + f[:-3])
	args = (
		resultsDir,
		outputDir,
		f[:-3],
		simData,
		validationData,
		)

	if "WC_ANALYZE_FAST" in os.environ:
		pool.apply_async(run_function, args = (mod.main, args, f))
	else:
		print "%s: Running %s" % (time.ctime(), f)
		mod.main(*args)
		# try:
		# 	mod.main(*args)
		# except SystemExit:
		# 	raise SystemExit(1)
		# except Exception as exc:
		# 	tb.print_exc()
		# 	continue

if "WC_ANALYZE_FAST" in os.environ:
	pool.close()
	pool.join()
timeTotal = time.time() - startTime
print "Completed single simulation analysis in %s" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))