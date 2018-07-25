run user/scripts/loadSimData.py
import numpy as np

su = [x[:-3] for x in sd.process.complexation.subunitNames]

re = [x['catalyzed by'] for x in rd.reactions]
re = np.unique(re)

r = []
for x in re: [r.append(y) for y in x]
re = np.unique(r)

