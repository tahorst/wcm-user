import os
import time
import cPickle

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1

print "%s: Loading raw" % time.ctime()
raw_data = KnowledgeBaseEcoli()
print "%s: Fitting" % time.ctime()
sim_data = fitSimData_1(raw_data)

fileLoc = os.path.dirname(__file__)
cPickle.dump(raw_data, open(os.path.join(fileLoc, "raw_data.cp"), "wb"))
cPickle.dump(sim_data, open(os.path.join(fileLoc, "sim_data.cp"), "wb"))