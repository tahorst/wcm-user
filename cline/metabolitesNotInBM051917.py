import csv
import cPickle
import os
import numpy as np
import re

fileLoc = os.path.dirname(os.path.dirname(__file__))
raw_data = cPickle.load(open(os.path.join(fileLoc, "scripts", "raw_data.cp"), "rb"))
sim_data = cPickle.load(open(os.path.join(fileLoc, "scripts", "sim_data.cp"), "rb"))

bm = [x[0] for x in sim_data.state.bulkMolecules.bulkData.fullArray()]

stoich = [x for x in sim_data.process.metabolism.reactionStoich.values()]
met = []
for m in stoich:
	[met.append(x) for x in m]

notin = np.unique([x for x in met if x not in bm])
print '%i metabolites not in reaction S matrix' % notin.shape[0]
