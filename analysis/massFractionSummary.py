#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants



from wholecell.utils.sparkline import whitePadSparklineAxis

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	# Get all cells
	allDir = ap.get_cells()

	massNames = [
				"dryMass",
				"proteinMass",
				#"tRnaMass",
				# "rRnaMass",
				'mRnaMass',
				"dnaMass"
				]

	cleanNames = [
				"Dry\nmass",
				"Protein\nmass",
				#"tRNA\nmass",
				# "rRNA\nmass",
				"mRNA\nmass",
				"DNA\nmass"
				]

	#plt.figure(figsize = (8.5, 11))
	fig, axesList = plt.subplots(len(massNames), sharex = True)

	currentMaxTime = 0

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")
		#initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		for idx, massType in enumerate(massNames):
			massToPlot = mass.readColumn(massNames[idx])
			# massToPlot = massToPlot / massToPlot[0]
			axesList[idx].plot(time / 60. / 60., massToPlot, linewidth = 2, color = "b")
			
			# set axes to size that shows all generations
			cellCycleTime = ((time[-1] - time[0]) / 60. / 60. )
			if time[-1] > currentMaxTime:
				currentMaxTime = time[-1] / 3600

			# axesList[idx].set_xlim(0, currentMaxTime*int(metadata["total_gens"])*1.1)
			axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")

	for axes in axesList:
		axes.get_ylim()
		axes.set_yticks(list(axes.get_ylim()))
		whitePadSparklineAxis(axes)
		axes.spines['top'].set_visible(False)
		axes.spines['bottom'].set_visible(False)
		axes.spines['right'].set_visible(False)
		axes.get_xaxis().set_visible(False)
		axes.xaxis.set_ticks_position('none')
		axes.yaxis.set_ticks_position('left')
		axes.tick_params(which = 'both', direction = 'out', labelsize = 6)
		axes.set_xlim((0, currentMaxTime))

	xmin, xmax = axesList[-1].get_xlim()
	axesList[-1].spines['bottom'].set_visible(True)
	axesList[-1].get_xaxis().set_visible(True)
	axesList[-1].set_xticks([xmin, currentMaxTime])
	axesList[-1].set_xticklabels(["%0.2f" % xmin, "%0.2f" % currentMaxTime])
	axesList[-1].xaxis.set_ticks_position('bottom')
	axesList[0].set_title("Cell mass fractions")
	axesList[len(massNames) - 1].set_xlabel("Time (hr)")

	plt.subplots_adjust(hspace = 0.2, wspace = 0.5)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
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
