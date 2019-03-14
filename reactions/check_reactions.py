'''
Checks the metabolites in reactions to see how many contain molecules not
in any other reactions.  With FBA, this reaction will not be able to occur.

Requires:
	user/scripts/sim_data.cp
'''

import cPickle
import os


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
SIM_DATA_FILE = os.path.join(ROOT_DIR, 'scripts', 'sim_data.cp')


if __name__ == '__main__':
	with open(SIM_DATA_FILE) as f:
		sim_data = cPickle.load(f)

	metabolism = sim_data.process.metabolism
	stoich = metabolism.reactionStoich
	reversible = metabolism.reversibleReactions

	# Identify number of times metabolites appear
	all_mets = set()
	reactants = set()
	products = set()
	for name, rxn in stoich.items():
		rev = name in reversible
		for met, direction in rxn.items():
			all_mets.add(met)

			if rev:
				reactants.add(met)
				products.add(met)
			elif direction > 0:
				products.add(met)
			elif direction < 0:
				reactants.add(met)

	# Add all metabolites that can be exchanged
	for met in metabolism.all_external_exchange_molecules:
		all_mets.add(met)
		reactants.add(met)
		products.add(met)

	single_mets = set([m for m in all_mets if not (m in reactants and m in products)])

	# Get reactions that will not occur
	no_flux_rxns = []
	for name, rxn in stoich.items():
		for met in rxn:
			if met in single_mets:
				no_flux_rxns.append(name)
				break
	no_flux_kinetics = [rxn for rxn in metabolism.constrainedReactionList if rxn in no_flux_rxns]

	print(no_flux_rxns)
	print(no_flux_kinetics)
	print('# of single metabolites: {}'.format(len(single_mets)))
	print('# of rxns that will not occur: {}'.format(len(no_flux_rxns)))
	print('# of rxns with kinetics that will not occur: {}'.format(len(no_flux_kinetics)))

	import ipdb; ipdb.set_trace()
