import cPickle
import os
import numpy as np

fileLoc = os.path.dirname(__file__)
raw_data = cPickle.load(open(os.path.join(fileLoc, "raw_data.cp"), "rb"))
sim_data = cPickle.load(open(os.path.join(fileLoc, "sim_data.cp"), "rb"))

rd = raw_data
sd = sim_data

import ipdb; ipdb.set_trace()