import numpy as np
import re

# compare metabolicEnzymes to reactionEnzymes
me = sd.process.metabolism.metabolicEnzymes	# len = 949, used in metabolism to +1
rxnE = np.unique([sd.process.metabolism.reactionEnzymes[x].keys() for x in sd.process.metabolism.reactionEnzymes]) # len = 900
rxne = []
for x in rxnE: rxne.append(x[0])
diff = [x for x in me if x not in rxne] # len = 165

subme = [re.sub('\[.\]', '', x) for x in me]
subrxne = [re.sub('\[.\]', '', x) for x in rxne]

subdiff = [x for x in subme if x not in subrxne] # len 152

