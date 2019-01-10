# include path from SCRATCH up to wildtype to run given file
# only runs first gen at the moment

import cPickle
import time
import sys
import os

from models.ecoli.sim.simulation import EcoliSimulation
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS

options = {}
# options['lengthSec'] = 30

CACHED = True
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print "%s: Running simulation" % time.ctime()

if len(sys.argv) == 2:
	location = os.path.join(BASE_DIR, "out", sys.argv[1])

	options["simDataLocation"] = os.path.join(location, "kb", "simData_Modified.cPickle")
	options["outputDir"] = os.path.join(location, "sim")
	options["logToDisk"] = False
	options["overwriteExistingFiles"] = False
else:
	location = os.path.dirname(__file__)
	if CACHED:
		simDataFile = os.path.join(BASE_DIR, "cached", "simData_Fit_1.cPickle")
	else:
		simDataFile = os.path.join(location, "sim_data.cp")
	if not os.path.exists(simDataFile):
		import runFitter

	options["simData"] = cPickle.load(open(simDataFile))
	options["outputDir"] = os.path.join(BASE_DIR, "out", "test")

# options["seed"] = 48
options['massDistribution'] = 1
options['growthRateNoise'] = 1
options['dPeriodDivision'] = 1


sim = EcoliSimulation(**options)
sim.run()
