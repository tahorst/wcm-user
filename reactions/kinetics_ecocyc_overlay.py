'''
Prints a url for ecocyc reaction network overlay with kinetic constrained reactions.
Specifically focused on AA synthesis with some reactions excluded.

Requires:
	user/scripts/sim_data.cp
'''

import cPickle
import os
import re


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
SIM_DATA_FILE = os.path.join(ROOT_DIR, 'scripts', 'sim_data.cp')


if __name__ == '__main__':
	with open(SIM_DATA_FILE) as f:
		sim_data = cPickle.load(f)
	constrained_rxns = sim_data.process.metabolism.constrainedReactionList

	# url for ecocyc to highlight fluxes that are 0 on metabolic network diagram
	siteStr = "https://ecocyc.org/overviewsWeb/celOv.shtml?zoomlevel=1&orgid=ECOLI"
	excluded = ['RXN0-2201', 'RXN-16000', 'RXN-12583', 'RXN-11496', 'DIMESULFREDUCT-RXN', '3.6.1.41-R[4/63051]5-NUCLEOTID-RXN',
		'1.1.1.262-RXN', 'FORMATEDEHYDROG-RXN', 'PHEAMINOTRANS-RXN', 'NADH-DEHYDROG-A-RXN', 'RXN0-3281', 'RXN-14325'] # reactions not recognized by ecocyc
	not_aa = ['1.1.1.39-RXN', 'DTDPGLUCDEHYDRAT-RXN', 'GSPSYN-RXN', 'GLYOHMETRANS-RXN', 'GLYOHMETRANS-RXN', 'DXS-RXN', 'ERYTH4PDEHYDROG-RXN', '2.5.1.64-RXN', 'CHORPYRLY-RXN',
		'GARTRANSFORMYL2-RXN', 'CTPSYN-RXN', 'DIHYDROOROT-RXN', 'ASPCARBTRANS-RXN']  # subset of reactions to display AA rxns (ecocyc limits number of reactions that can overlay)
	replace = {'PHEAMINOTRANS-RXN': 'RXN-10814'}  # AA synthesis rxn not recognized by ecocyc

	rxns = []
	for i, reaction in enumerate(constrained_rxns):
		# if actualAve[i] == 0:
		rxn = re.findall(".+RXN", reaction)
		if len(rxn) == 0:
			rxn = re.findall("RXN[^-]*-[0-9]+", reaction)
		rxn[0] = replace.get(rxn[0], rxn[0])
		if rxn[0] not in excluded and rxn[0] not in not_aa and rxn[0] not in rxns and '--TRNA' not in rxn[0]:
			siteStr += "&rnids=%s" % rxn[0]
		rxns.append(rxn[0])
	print siteStr
