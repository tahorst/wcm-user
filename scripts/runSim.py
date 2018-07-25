# include path from SCRATCH up to wildtype to run given file
# only runs first gen at the moment

import time
import sys
import os

from models.ecoli.sim.simulation import EcoliSimulation
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS

options = {}
# options['lengthSec'] = 30

CACHED = True

print "%s: Running simulation" % time.ctime()

if len(sys.argv) == 2:
	location = os.path.join("/home/users/thorst/wcEcoli/out", sys.argv[1])

	options["simDataLocation"] = os.path.join(location, "kb", "simData_Modified.cPickle")
	options["outputDir"] = os.path.join(location, "sim")
	options["logToDisk"] = False
	options["overwriteExistingFiles"] = False
else:
	location = os.path.dirname(__file__)
	if CACHED:
		simDataFile = "/home/users/thorst/wcEcoli/cached/simData_Fit_1.cPickle"
	else:
		simDataFile = os.path.join(location, "sim_data.cp")
	if not os.path.exists(simDataFile):
		import runFitter

	options["simDataLocation"] = simDataFile
	options["outputDir"] = "/home/users/thorst/wcEcoli/out/test"

sim = EcoliSimulation(**options)
sim.run()