run user/scripts/loadSimData.py
import numpy as np

enk = [x['enzymeIDs'] for x in rd.enzymeKinetics]
enk = np.unique(enk)
e = [x[:-3] for x in enk]
enk = e
enk = np.unique(enk)

re = [x['catalyzed by'] for x in rd.reactions]
re = np.unique(re)

r = []
for x in re: [r.append(y) for y in x]
re = np.unique(r)

both = [x for x in enk if x in re]
nore = [x for x in enk if x not in re]
