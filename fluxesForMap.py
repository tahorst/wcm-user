#!/usr/bin/env python
"""
Plot amino acid counts

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from copy import copy

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

BURNIN = 25

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	# return

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	validation_data = cPickle.load(open(validationDataFile, "rb"))
	toya_reactions = validation_data.getReactionFlux.toya2010fluxes["reactionID"]

	sim_data = cPickle.load(open(simDataFile))
	kineticsReactions = sim_data.process.metabolism.constrainedReactionList
	enzymeIDs = [sim_data.process.metabolism.reactionCatalysts[rxn][0] for rxn in kineticsReactions]
	cellDensity = sim_data.constants.cellDensity
	nAvogadro = sim_data.constants.nAvogadro
	rxnStoich = sim_data.process.metabolism.reactionStoich

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	fbaReader = TableReader(os.path.join(simOutDir, "FBAResults"))
	flux = fbaReader.readColumn("reactionFluxes")
	reactionIDs = fbaReader.readAttribute("reactionIDs")
	fbaReader.close()

	reactants = [
		'GLC-6-P[c]',
		'FRUCTOSE-6P[c]',
		'FRUCTOSE-16-DIPHOSPHATE[c]',
		'DIHYDROXY-ACETONE-PHOSPHATE[c]',
		'GAP[c]',
		'DPG[c]',
		'G3P[c]',
		'2-PG[c]',
		'PHOSPHO-ENOL-PYRUVATE[c]',
		'PYRUVATE[c]',
		'ACETYL-COA[c]',
		'CIT[c]',
		'CIS-ACONITATE[c]',
		'THREO-DS-ISO-CITRATE[c]',
		'2-KETOGLUTARATE[c]',
		'SUC-COA[c]',
		'SUC[c]',
		'FUM[c]',
		'MAL[c]',
		'GLC-6-P[c]',
		'D-6-P-GLUCONO-DELTA-LACTONE[c]',
		'CPD-2961[c]',
		'RIBULOSE-5P[c]',
		'RIBULOSE-5P[c]',
		'XYLULOSE-5-PHOSPHATE[c]',
		'D-SEDOHEPTULOSE-7-P[c]',
		'FRUCTOSE-6P[c]',
		'3-P-SERINE[c]',
		'SER[c]',
		'ACETYLSERINE[c]',
		'HOMO-CYS[c]',
		'GLT[c]',
		'GLT[c]',
		'SER[c]',
		'HISTIDINAL[c]',
		'PYRUVATE[c]',
		'GLT[c]',
		'GLT[c]',
		'GLT[c]',
		'GLT[c]',
		'L-ASPARTATE[c]',
		'O-PHOSPHO-L-HOMOSERINE[c]',
		'MESO-DIAMINOPIMELATE[c]',
		'L-ARGININO-SUCCINATE[c]',
		'2-KETOGLUTARATE[c]',
		'GLT[c]',
		'L-DELTA1-PYRROLINE_5-CARBOXYLATE[c]',
		'DGMP[c]',
		'DGDP[c]',
		'DGMP[c]',
		'DAMP[c]',
		'DADP[c]',
		'DAMP[c]',
		'TMP[c]',
		'TDP[c]',
		'TMP[c]',
		'DUMP[c]',
		'DUDP[c]',
		'DUMP[c]',
		'DCMP[c]',
		'DCDP[c]',
		'DCMP[c]',
		'GMP[c]',
		'GDP[c]',
		'GMP[c]',
		'AMP[c]',
		'ADP[c]',
		'AMP[c]',
		'UMP[c]',
		'UDP[c]',
		'UMP[c]',
		'CMP[c]',
		'CDP[c]',
		'CMP[c]',
		'GDP[c]',
		'GTP[c]',
		'ADP[c]',
		'ATP[c]',
		'TMP[c]',
		'UDP[c]',
		'UTP[c]',
		'UTP[c]',
		'CDP[c]',
		'CTP[c]',
		'GUANINE[c]',
		'ADENINE[c]',
		'URACIL[c]',
		]

	products = [
		'FRUCTOSE-6P[c]',
		'FRUCTOSE-16-DIPHOSPHATE[c]',
		'DIHYDROXY-ACETONE-PHOSPHATE[c]',
		'GAP[c]',
		'DPG[c]',
		'G3P[c]',
		'2-PG[c]',
		'PHOSPHO-ENOL-PYRUVATE[c]',
		'PYRUVATE[c]',
		'ACETYL-COA[c]',
		'CIT[c]',
		'CIS-ACONITATE[c]',
		'THREO-DS-ISO-CITRATE[c]',
		'2-KETOGLUTARATE[c]',
		'SUC-COA[c]',
		'SUC[c]',
		'FUM[c]',
		'MAL[c]',
		'OXALACETIC_ACID[c]',
		'D-6-P-GLUCONO-DELTA-LACTONE[c]',
		'CPD-2961[c]',
		'RIBULOSE-5P[c]',
		'XYLULOSE-5-PHOSPHATE[c]',
		'RIBOSE-5P[c]',
		'D-SEDOHEPTULOSE-7-P[c]',
		'ERYTHROSE-4P[c]',
		'ERYTHROSE-4P[c]',
		'SER[c]',
		'GLY[c]',
		'CYS[c]',
		'MET[c]',
		'TYR[c]',
		'PHE[c]',
		'TRP[c]',
		'HIS[c]',
		'L-ALPHA-ALANINE[c]',
		'VAL[c]',
		'LEU[c]',
		'ILE[c]',
		'L-ASPARTATE[c]',
		'ASN[c]',
		'THR[c]',
		'LYS[c]',
		'ARG[c]',
		'GLT[c]',
		'GLN[c]',
		'PRO[c]',
		'DGDP[c]',
		'DGTP[c]',
		'DGTP[c]',
		'DADP[c]',
		'DATP[c]',
		'DATP[c]',
		'TDP[c]',
		'TTP[c]',
		'TTP[c]',
		'DUDP[c]',
		'DUTP[c]',
		'DUTP[c]',
		'DCDP[c]',
		'DCTP[c]',
		'DCTP[c]',
		'GDP[c]',
		'GTP[c]',
		'GTP[c]',
		'ADP[c]',
		'ATP[c]',
		'ATP[c]',
		'UDP[c]',
		'UTP[c]',
		'UTP[c]',
		'CDP[c]',
		'CTP[c]',
		'CTP[c]',
		'DGDP[c]',
		'DGTP[c]',
		'DADP[c]',
		'DATP[c]',
		'DUMP[c]',
		'DUDP[c]',
		'DUTP[c]',
		'CTP[c]',
		'DCDP[c]',
		'DCTP[c]',
		'GMP[c]',
		'AMP[c]',
		'UMP[c]',
		]


	plt.figure(figsize = (17, 22))

	for idx, (reactant, product) in enumerate(zip(reactants, products)):
		plt.subplot(10, 9, idx + 1)
		totalFlux = np.zeros_like(flux[:, 0])

		# rxns = [rxn for rxn in rxnStoich if product in rxnStoich[rxn] and reactant in rxnStoich[rxn] and rxnStoich[rxn][product] / rxnStoich[rxn][reactant] < 0]
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

		plt.plot(time, totalFlux)
		plt.title("%s to %s" % (reactant, product), fontsize = 4)
		plt.tick_params(axis = "both", labelsize = 4)

	# reactions = [
	# 	'PGLUCISOM-RXN',
	# 	'6PFRUCTPHOS-RXN',
	# 	'F16ALDOLASE-RXN',
	# 	'TRIOSEPISOMERIZATION-RXN',
	# 	'GAPOXNPHOSPHN-RXN',
	# 	'PHOSGLYPHOS-RXN',
	# 	'RXN-15513',
	# 	'2PGADEHYDRAT-RXN',
	# 	'PEPDEPHOS-RXN',
	# 	'PYRUVDEH-RXN',
	# 	'CITSYN-RXN',
	# 	'ACONITATEDEHYDR-RXN',
	# 	'ACONITATEHYDR-RXN',
	# 	'ISOCITDEH-RXN',
	# 	'2OXOGLUTARATEDEH-RXN',
	# 	'SUCCCOASYN-RXN',
	# 	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	# 	'FUMHYDR-RXN',
	# 	'MALATE-DEH-RXN',
	# 	'ISOCIT-CLEAV-RXN',
	# 	'MALSYN-RXN',
	# 	'GLU6PDEHYDROG-RXN',
	# 	'6PGLUCONOLACT-RXN',
	# 	'RXN-9952',
	# 	'RIBULP3EPIM-RXN',
	# 	'RIB5PISOM-RXN',
	# 	'1TRANSKETO-RXN',
	# 	'TRANSALDOL-RXN',
	# 	'2TRANSKETO-RXN',
	# 	'RXN0-5114',
	# 	'GLYOHMETRANS-RXN',
	# 	'ACSERLY-RXN',
	# 	'O-SUCCHOMOSERLYASE-RXN',
	# 	'TYROSINE-AMINOTRANSFERASE-RXN',
	# 	'PHEAMINOTRANS-RXN',
	# 	'RXN0-2382',
	# 	'HISTALDEHYD-RXN',
	# 	'ALANINE-AMINOTRANSFERASE-RXN',
	# 	'VALINE-PYRUVATE-AMINOTRANSFER-RXN',
	# 	'BRANCHED-CHAINAMINOTRANSFERVAL-RXN',
	# 	'BRANCHED-CHAINAMINOTRANSFERLEU-RXN',
	# 	'BRANCHED-CHAINAMINOTRANSFERILEU-RXN',
	# 	'ASPAMINOTRANS-RXN',
	# 	'ASNSYNA-RXN',
	# 	'THRESYN-RXN',
	# 	'DIAMINOPIMDECARB-RXN',
	# 	'ARGSUCCINLYA-RXN',
	# 	'GLUTDEHYD-RXN',
	# 	'GLUTAMINESYN-RXN',
	# 	'PYRROLINECARBREDUCT-RXN',
	# 	]

	# reverse = {
	# 	"6PFRUCTPHOS-RXN": "F16BDEPHOS-RXN",
	# 	"PEPDEPHOS-RXN": "PEPSYNTH-RXN",
	# 	"ASNSYNA-RXN": "ASPARAGHYD-RXN",
	# 	"GLUTAMINESYN-RXN": "GLUTAMATESYN-RXN",
	# 	}

	# branch = {
	# 	"RXN-15513": "3PGAREARR-RXN",
	# 	"MALATE-DEH-RXN": "MALATE-DEHYDROGENASE-ACCEPTOR-RXN",
	# 	"ASNSYNA-RXN": "ASNSYNB-RXN",
	# }

	# plt.figure(figsize = (17, 22))

	# for idx, rxn in enumerate(reactions):
	# 	plt.subplot(9, 6, idx + 1)

	# 	if rxn in reactionIDs:
	# 		totalFlux = copy(flux[:, reactionIDs.index(rxn)])

	# 		if rxn + " (reverse)" in reactionIDs:
	# 			totalFlux -= flux[:, reactionIDs.index(rxn + " (reverse)")]

	# 		if rxn in reverse and reverse[rxn] in reactionIDs:
	# 			totalFlux -= flux[:, reactionIDs.index(reverse[rxn])]
	# 	else:
	# 		continue

	# 	plt.plot(time, totalFlux)
	# 	plt.title(rxn, fontsize = 4)
	# 	plt.tick_params(axis = "both", labelsize = 4)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")
	import ipdb; ipdb.set_trace()

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
