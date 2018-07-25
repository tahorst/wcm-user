#!/usr/bin/env python
"""
Plot fluxes for metabolic map figure during a shift

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/13/17
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis

START = 8000
SHIFT = 11000
END = 14000
BURNIN = 15

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()

	sim_data = cPickle.load(open(simDataFile, "rb"))
	rxnStoich = sim_data.process.metabolism.reactionStoich

	reactants = [
		"GLT[c]",
		"ATP[c]",
		]

	products = [
		"TYR[c]",
		"DATP[c]",
		]

	plt.figure(figsize = (8, 8))

	firstGen = True
	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")

		mainListener = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = mainListener.readAttribute("initialTime")
		time = mainListener.readColumn("time")
		timeStepSec = mainListener.readColumn("timeStepSec")
		mainListener.close()

		if initialTime > END or time[-1] < START:
			continue

		massListener = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massListener.readColumn("cellMass")
		dryMass = massListener.readColumn("dryMass")
		massListener.close()

		coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(units.g / units.L) * timeStepSec # units - g.s/L

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionIDs = fbaResults.readAttribute("reactionIDs")
		flux = (units.dmol / units.g / units.s) * (fbaResults.readColumn("reactionFluxes").T / coefficient).T
		fbaResults.close()

		flux = flux.asNumber(units.mmol / units.g / units.h)

		for idx, (reactant, product) in enumerate(zip(reactants, products)):
			ax1 = plt.subplot(2, 2, 2*idx + 1)
			ax2 = plt.subplot(2, 2, 2*(idx + 1))
			totalFlux = np.zeros_like(flux[:, 0])

			for rxn in rxnStoich:
				if reactant in rxnStoich[rxn] and product in rxnStoich[rxn]:
					if rxnStoich[rxn][reactant] < 0 and rxnStoich[rxn][product] > 0:
						direction = 1
					elif rxnStoich[rxn][reactant] > 0 and rxnStoich[rxn][product] < 0:
						direction = -1
					else:
						continue

					if rxn in reactionIDs:
						totalFlux += flux[:, reactionIDs.index(rxn)] * direction
						if rxn + " (reverse)" in reactionIDs:
							totalFlux -= flux[:, reactionIDs.index(rxn + " (reverse)")] * direction

			timeIdx = np.logical_and(np.logical_or(np.logical_and(time >= START, time < SHIFT), np.logical_and(time > SHIFT + BURNIN, time <= END)), time > initialTime + BURNIN)

			if firstGen:
				ax1.axvline(SHIFT, color = "#aaaaaa")
				ax1.set_title("%s to %s" % (reactant, product), fontsize = 4)
				ax1.set_xlabel("Time (s)", fontsize = 8)
				ax1.set_ylabel("Flux (mmol / g DCW / hr)", fontsize = 8)
				ax1.tick_params(axis = "both", labelsize = 6)
				ax1.yaxis.get_offset_text().set_fontsize(6)
				ax1.plot([SHIFT, SHIFT], [0, 0], color = "k")

				ax2.axvline(SHIFT, color = "#aaaaaa")
				ax2.set_title("%s to %s" % (reactant, product), fontsize = 4)
				ax2.set_xlabel("Time (s)", fontsize = 8)
				ax2.set_ylabel("Flux (mmol / hr)", fontsize = 8)
				ax2.tick_params(axis = "both", labelsize = 6)
				ax2.yaxis.get_offset_text().set_fontsize(6)
				ax2.plot([SHIFT, SHIFT], [0, 0], color = "k")

			ax1.plot(time[timeIdx], totalFlux[timeIdx])
			ax2.plot(time[timeIdx], totalFlux[timeIdx] * dryMass[timeIdx] * 1e-15)
			whitePadSparklineAxis(ax1)
			whitePadSparklineAxis(ax2)

		firstGen = False

	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
