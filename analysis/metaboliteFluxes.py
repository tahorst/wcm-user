#!/usr/bin/env python
"""
Compare fluxes in simulation to target fluxes

@date: Created 12/15/16
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os
import cPickle
import csv
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.plotting_tools import COLORS_LARGE

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

BURN_IN_STEPS = 20
METABOLITE = 'FRUCTOSE-16-DIPHOSPHATE[c]'

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	return
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile))
	reactionStoich = sim_data.process.metabolism.reactionStoich
	constrainedReactions = np.array(sim_data.process.metabolism.constrainedReactionList)
	enzymeNames = sim_data.process.metabolism.enzymeIdList

	mainListener = TableReader(os.path.join(simOutDir, "Main"))
	initialTime = mainListener.readAttribute("initialTime")
	time = mainListener.readColumn("time") - initialTime
	timeStepSec = mainListener.readColumn("timeStepSec")
	mainListener.close()

	massListener = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = massListener.readColumn("cellMass")
	dryMass = massListener.readColumn("dryMass")
	massListener.close()

	coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(units.g / units.L) # units - g/L

	enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	targetFluxes = (units.dmol / units.g / units.s) * (enzymeKineticsReader.readColumn("targetFluxes").T / coefficient).T
	metaboliteNames = enzymeKineticsReader.readAttribute("metaboliteNames")
	metaboliteCounts = enzymeKineticsReader.readColumn("metaboliteCountsFinal")
	enzymeCounts = enzymeKineticsReader.readColumn("enzymeCountsInit")
	enzymeKineticsReader.close()

	targetFluxes = targetFluxes.asNumber(units.mmol / units.g / units.h)
	targetAve = np.mean(targetFluxes[BURN_IN_STEPS:, :], axis = 0)

	# FBA reader
	fbaReader = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionFluxes = (units.dmol / units.g / units.s) * (fbaReader.readColumn("reactionFluxes").T / coefficient).T
	reactionNames = fbaReader.readAttribute("reactionIDs")
	fbaReader.close()

	reactionFluxes = reactionFluxes.asNumber(units.mmol / units.g / units.h)
	fluxAve = np.mean(reactionFluxes[BURN_IN_STEPS:, :], axis = 0)

	print "Conc target:\t%f"
	for rxn in reactionStoich:
		if METABOLITE in reactionStoich[rxn]:
			stoich = reactionStoich[rxn][METABOLITE]
			flux = fluxAve[reactionNames.index(rxn)]
			target = targetAve[np.where(constrainedReactions == rxn)[0]]

			if rxn in constrainedReactions:
				print "%s:\t%.2f\t%.2f" % (rxn, flux*stoich, target*stoich)
			else:
				print "%s:\t%.2f\t---" % (rxn, flux*stoich)
	import ipdb; ipdb.set_trace()
		
	# read constraint data
	enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	actualFluxes = (units.dmol / units.g / units.s) * (enzymeKineticsReader.readColumn("actualFluxes").T / coefficient).T
	reactionConstraint = enzymeKineticsReader.readColumn("reactionConstraint")
	enzymeKineticsReader.close()

	actualFluxes = actualFluxes.asNumber(units.mmol / units.g / units.h)

	actualAve = np.mean(actualFluxes[BURN_IN_STEPS:, :], axis = 0)

	relError = np.abs((actualFluxes[BURN_IN_STEPS:, :] - targetFluxes[BURN_IN_STEPS:, :]) / (targetFluxes[BURN_IN_STEPS:, :] + 1e-15))
	aveError = np.mean(relError, axis = 0)

	kcatOnlyReactions = np.all(constraintIsKcatOnly[reactionConstraint[BURN_IN_STEPS:,:]], axis = 0)
	kmAndKcatReactions = ~np.any(constraintIsKcatOnly[reactionConstraint[BURN_IN_STEPS:,:]], axis = 0)
	mixedReactions = ~(kcatOnlyReactions ^ kmAndKcatReactions)

	kmAndKcatThresholds = [2, 10]
	kmAndKcatCategorization = np.zeros(np.sum(kmAndKcatReactions))
	categorization = np.zeros(reactionConstraint.shape[1])
	categorization[actualAve == 0] = -1
	for i, threshold in enumerate(kmAndKcatThresholds):
		# kmAndKcatCategorization[targetAve[kmAndKcatReactions] / actualAve[kmAndKcatReactions] > threshold] = i + 1
		kmAndKcatCategorization[actualAve[kmAndKcatReactions] / targetAve[kmAndKcatReactions] > threshold] = i + 1
		categorization[actualAve / targetAve > threshold] = i + 1
	kmAndKcatCategorization[actualAve[kmAndKcatReactions] == 0] = -1

	kcatOnlyThresholds = [2, 10]
	kcatOnlyCategorization = np.zeros(np.sum(kcatOnlyReactions))
	for i, threshold in enumerate(kcatOnlyThresholds):
		kcatOnlyCategorization[actualAve[kcatOnlyReactions] / targetAve[kcatOnlyReactions] > threshold] = i + 1
	kcatOnlyCategorization[actualAve[kcatOnlyReactions] == 0] = -1

	# url for ecocyc to highlight fluxes that are 0 on metabolic network diagram
	siteStr = "https://ecocyc.org/overviewsWeb/celOv.shtml?zoomlevel=1&orgid=ECOLI"
	excluded = ['RXN0-2201', 'RXN-16000', 'RXN-12583', 'RXN-11496', 'DIMESULFREDUCT-RXN', '3.6.1.41-R[4/63051]5-NUCLEOTID-RXN'] # reactions not recognized by ecocyc
	rxns = []
	for i, reaction in enumerate(constrainedReactions):
		if categorization[i] != -1:
			continue
		if actualAve[i] == 0:
			rxn = re.findall(".+RXN", reaction)
			if len(rxn) == 0:
				rxn = re.findall("RXN[^-]*-[0-9]+", reaction)
			if rxn[0] not in excluded:
				siteStr += "&rnids=%s" % rxn[0]
			rxns.append(rxn[0])
	print siteStr
	import ipdb; ipdb.set_trace()

	csvFile = open(os.path.join(plotOutDir, plotOutFileName + ".tsv"), "wb")
	output = csv.writer(csvFile, delimiter = "\t")
	output.writerow(["ecocyc link:", siteStr])
	output.writerow(["kM and kcat", "Target", "Actual", "Category"])
	for reaction, target, flux, category in zip(constrainedReactions[kmAndKcatReactions], targetAve[kmAndKcatReactions], actualAve[kmAndKcatReactions], kmAndKcatCategorization):
		output.writerow([reaction, target, flux, category])

	output.writerow(["kcat only"])
	for reaction, target, flux, category in zip(constrainedReactions[kcatOnlyReactions], targetAve[kcatOnlyReactions], actualAve[kcatOnlyReactions], kcatOnlyCategorization):
		output.writerow([reaction, target, flux, category])

	csvFile.close()

	targetAve += 1e-6
	actualAve += 1e-6

	plt.figure(figsize = (8, 8))
	from scipy.stats import pearsonr
	targetPearson = targetAve[kmAndKcatReactions]
	actualPearson = actualAve[kmAndKcatReactions]
	# plt.title(pearsonr(np.log10(targetPearson[actualPearson > 0]), np.log10(actualPearson[actualPearson > 0])))
	# plt.loglog(targetAve[kmAndKcatReactions][kmAndKcatCategorization == 0], actualAve[kmAndKcatReactions][kmAndKcatCategorization == 0], "og")
	# plt.loglog(targetAve[kmAndKcatReactions][kmAndKcatCategorization == 1], actualAve[kmAndKcatReactions][kmAndKcatCategorization == 1], "o")
	# plt.loglog(targetAve[kmAndKcatReactions][kmAndKcatCategorization == 2], actualAve[kmAndKcatReactions][kmAndKcatCategorization == 2], "or")
	# plt.loglog(targetAve[kmAndKcatReactions][kmAndKcatCategorization == -1], actualAve[kmAndKcatReactions][kmAndKcatCategorization == -1], "og")
	# plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 0], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 0], "og")
	# plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 1], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 1], "o")
	# plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == 2], actualAve[kcatOnlyReactions][kcatOnlyCategorization == 2], "or")
	# plt.loglog(targetAve[kcatOnlyReactions][kcatOnlyCategorization == -1], actualAve[kcatOnlyReactions][kcatOnlyCategorization == -1], "og")
	plt.loglog(targetAve[categorization == 0], actualAve[categorization == 0], "og", markeredgewidth = 0.25)
	plt.loglog(targetAve[categorization == 1], actualAve[categorization == 1], "o", markeredgewidth = 0.25)
	plt.loglog(targetAve[categorization == 2], actualAve[categorization == 2], "or", markeredgewidth = 0.25)
	plt.loglog(targetAve[categorization == -1], actualAve[categorization == -1], "og", markeredgewidth = 0.25)
	# plt.loglog(targetAve[kmAndKcatReactions], actualAve[kmAndKcatReactions], "o")
	# plt.loglog(targetAve[kcatOnlyReactions], actualAve[kcatOnlyReactions], "ro")
	plt.loglog([1e-7, 1e4], [1e-7, 1e4], '--g')
	plt.loglog([1e-7, 1e3], [1e-6, 1e4], '--r')
	# plt.loglog([1e-13, 1], [1e-14, 0.1], '--r')
	plt.xlabel("Target Flux (mmol/g/hr)")
	plt.ylabel("Actual Flux (mmol/g/hr)")
	plt.minorticks_off()
	whitePadSparklineAxis(plt.axes())

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
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
