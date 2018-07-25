#!/usr/bin/env python
"""
Compare fluxes in simulation to target fluxes for supplemental figure 2

@date: Created 4/3/17
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

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.plotting_tools import COLORS_LARGE

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS

from multiprocessing import Pool

# ignore data from metabolism burnin period
BURN_IN_TIME = 1

def exportData(i, simDir, rnaNames, proteinNames, plotOutDir):
	seed = i // 32
	gen = i % 32
	simOutDir = os.path.join(simDir, "simOut")

	mainListener = TableReader(os.path.join(simOutDir, "Main"))
	time = mainListener.readColumn("time")
	mainListener.close()

	bulkListener = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeNames = bulkListener.readAttribute('objectNames')
	bulkMolecules = bulkListener.readColumn("counts")
	bulkListener.close()

	rnaIdx = [moleculeNames.index(x) for x in rnaNames]
	rnaCounts = bulkMolecules[:, rnaIdx]

	proteinIdx = [moleculeNames.index(x) for x in proteinNames]
	proteinCounts = bulkMolecules[:, proteinIdx]

	parent = os.path.join(plotOutDir, 'fathom', str(seed))
	if not os.path.exists(parent):
		os.makedirs(parent)

	with open(os.path.join(parent, 'mRNA%s.tsv' % gen), "wb") as csvFile:
		output = csv.writer(csvFile, delimiter = "\t")
		output.writerow(["Time"] + rnaNames.tolist())
		for t, c in zip(time, rnaCounts):
			output.writerow([t] + c.tolist())

	with open(os.path.join(parent, 'protein%s.tsv' % gen), "wb") as csvFile:
		output = csv.writer(csvFile, delimiter = "\t")
		output.writerow(["Time"] + proteinNames.tolist())
		for t, c in zip(time, proteinCounts):
			output.writerow([t] + c.tolist())

	print seed, gen

def exportMonomers(i, simDir, ids_complexation, ids_complexation_complexes, ids_equilibrium, ids_equilibrium_complexes, ids_translation,
				   ribosome_subunit_ids, rnap_subunit_ids, ribosome_subunit_stoich, rnap_subunit_stoich, complexationStoichMatrixMonomers,
				   equilibriumStoichMatrixMonomers, plotOutDir):
	seed = i // 32
	gen = i % 32
	simOutDir = os.path.join(simDir, "simOut")

	mainListener = TableReader(os.path.join(simOutDir, "Main"))
	time = mainListener.readColumn("time")
	mainListener.close()

	bulkListener = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkListener.readAttribute('objectNames')
	bulkCounts = bulkListener.readColumn("counts")
	bulkListener.close()

	complexationIdx = np.array([moleculeIds.index(x) for x in ids_complexation])  # Complexes of proteins, and protein monomers
	complexation_complexesIdx = np.array([moleculeIds.index(x) for x in ids_complexation_complexes])  # Only complexes
	equilibriumIdx = np.array([moleculeIds.index(x) for x in ids_equilibrium])  # Complexes of proteins + small molecules, small molecules, protein monomers
	equilibrium_complexesIdx = np.array([moleculeIds.index(x) for x in ids_equilibrium_complexes])  # Only complexes
	translationIdx = np.array([moleculeIds.index(x) for x in ids_translation])  # Only protein monomers

	ribosomeIdx = np.array([moleculeIds.index(x) for x in ribosome_subunit_ids])
	rnapIdx = np.array([moleculeIds.index(x) for x in rnap_subunit_ids])
	bulkCounts[:, complexationIdx] += np.dot(complexationStoichMatrixMonomers, bulkCounts[:, complexation_complexesIdx].transpose() * -1).transpose()

	# Dissociate protein-small molecule complexes
	bulkCounts[:, equilibriumIdx] += np.dot(equilibriumStoichMatrixMonomers, bulkCounts[:, equilibrium_complexesIdx].transpose() * -1).transpose()

	# Load unique molecule data for RNAP and ribosomes
	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
	ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
	rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
	nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
	nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
	uniqueMoleculeCounts.close()

	# Add subunits from RNAP and ribosomes
	ribosomeSubunitCounts = (nActiveRibosome.reshape((nActiveRibosome.size, 1)) * ribosome_subunit_stoich.reshape(
		(1, ribosome_subunit_stoich.size)))
	rnapSubunitCounts = (nActiveRnaPoly.reshape((nActiveRnaPoly.size, 1)) * rnap_subunit_stoich.reshape(
		(1, rnap_subunit_stoich.size)))

	bulkCounts[:, ribosomeIdx] += ribosomeSubunitCounts
	bulkCounts[:, rnapIdx] += rnapSubunitCounts

	# Get protein monomer counts for calculations now that all complexes are dissociated
	proteinMonomerCounts = bulkCounts[:, translationIdx]

	parent = os.path.join(plotOutDir, 'fathom', str(seed))
	if not os.path.exists(parent):
		os.makedirs(parent)

	with open(os.path.join(parent, 'monomers%s.tsv' % gen), "wb") as csvFile:
		output = csv.writer(csvFile, delimiter = "\t")
		output.writerow(["Time"] + ids_translation)
		for t, c in zip(time, proteinMonomerCounts):
			output.writerow([t] + c.tolist())

	print seed, gen

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile = None, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	allDir = ap.get_cells()

	sim_data = cPickle.load(open(simDataFile, "rb"))
	rnaNames = sim_data.process.transcription.rnaData['id'][sim_data.process.transcription.rnaData['isMRna']]
	monomerNames = sim_data.process.translation.monomerData['id']
	complexNames = np.array([x for x in sim_data.process.complexation.complexNames])
	proteinNames = np.hstack((monomerNames, complexNames))
	#
	# pool = Pool(processes=16)
	# results = [pool.apply_async(exportData, (i, simDir, rnaNames, proteinNames, plotOutDir)) for i, simDir in enumerate(sorted(allDir))]
	# pool.close()
	# pool.join()
	#
	# for result in results:
	# 	if not result.successful():
	# 		print result.get()

	# Get all ids required for getting monomers from complexes
	ids_complexation = sim_data.process.complexation.moleculeNames # Complexes of proteins, and protein monomers
	ids_complexation_complexes = [ids_complexation[i] for i in np.where((sim_data.process.complexation.stoichMatrix() == 1).sum(axis = 1))[0]] # Only complexes
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames # Complexes of proteins + small molecules, small molecules, protein monomers
	ids_equilibrium_complexes = [ids_equilibrium[i] for i in np.where((sim_data.process.equilibrium.stoichMatrix() == 1).sum(axis = 1))[0]] # Only complexes
	ids_translation = sim_data.process.translation.monomerData["id"].tolist() # Only protein monomers

	data_50s = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s50_fullComplex[0])
	data_30s = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s30_fullComplex[0])
	ribosome_subunit_ids = data_50s["subunitIds"].tolist() + data_30s["subunitIds"].tolist()
	ribosome_subunit_stoich = np.hstack((data_50s["subunitStoich"],data_30s["subunitStoich"]))
	data_rnap = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.rnapFull[0])
	rnap_subunit_ids = data_rnap["subunitIds"].tolist()
	rnap_subunit_stoich = data_rnap["subunitStoich"]
	complexationStoichMatrixMonomers = sim_data.process.complexation.stoichMatrixMonomers()
	equilibriumStoichMatrixMonomers = sim_data.process.equilibrium.stoichMatrixMonomers()

	pool = Pool(processes=16)
	results = [pool.apply_async(exportMonomers, (i, simDir, ids_complexation, ids_complexation_complexes, ids_equilibrium, ids_equilibrium_complexes,
						ids_translation, ribosome_subunit_ids, rnap_subunit_ids, ribosome_subunit_stoich, rnap_subunit_stoich,
						complexationStoichMatrixMonomers, equilibriumStoichMatrixMonomers, plotOutDir))
			   			for i, simDir in enumerate(sorted(allDir))]
	pool.close()
	pool.join()

	for result in results:
		if not result.successful():
			print result.get()

	# for i, simDir in enumerate(sorted(allDir)):
	# 	seed = i // 32
	# 	gen = i % 32
	# 	print seed, gen
	# 	simOutDir = os.path.join(simDir, "simOut")
	#
	# 	mainListener = TableReader(os.path.join(simOutDir, "Main"))
	# 	time = mainListener.readColumn("time")
	# 	mainListener.close()
	#
	# 	bulkListener = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	# 	moleculeNames = bulkListener.readAttribute('objectNames')
	# 	bulkMolecules = bulkListener.readColumn("counts")
	# 	bulkListener.close()
	#
	# 	rnaIdx = [moleculeNames.index(x) for x in rnaNames]
	# 	rnaCounts = bulkMolecules[:, rnaIdx]
	#
	# 	proteinIdx = [moleculeNames.index(x) for x in proteinNames]
	# 	proteinCounts = bulkMolecules[:, proteinIdx]
	#
	# 	parent = os.path.join(plotOutDir, 'fathom', str(seed))
	# 	if not os.path.exists(parent):
	# 		os.makedirs(parent)
	#
	# 	with open(os.path.join(parent, 'mRNA%s.tsv' % gen), "wb") as csvFile:
	# 		output = csv.writer(csvFile, delimiter = "\t")
	# 		output.writerow(["Time"] + rnaNames.tolist())
	# 		for t, c in zip(time, rnaCounts):
	# 			output.writerow([t] + c.tolist())
	#
	# 	with open(os.path.join(parent, 'protein%s.tsv' % gen), "wb") as csvFile:
	# 		output = csv.writer(csvFile, delimiter = "\t")
	# 		output.writerow(["Time"] + proteinNames.tolist())
	# 		for t, c in zip(time, proteinCounts):
	# 			output.writerow([t] + c.tolist())


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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
