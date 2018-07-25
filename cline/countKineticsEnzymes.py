import re
import numpy as np

# count of reactions that are constrained (not counting multiple enzyme reactions multiple times)
fwd = [re.split('__', x)[0] for x in self.kineticsConstrainedReactions if not (re.findall('(reverse)', x) and re.findall('__', x))]
rev = [re.split('__', x)[0] + ' (reverse)' for x in self.kineticsConstrainedReactions if re.findall('(reverse)', x) and re.findall('__', x)]

print np.unique(fwd + rev).shape - len(self.constraintsToDisable)


# count of constraints that are used disregarding ones for disabled constraints
idx = [self.kineticsConstrainedReactions.index(x) for x in self.kineticsConstrainedReactions if x not in self.constraintsToDisable]
disabled = np.where(np.sum(self.constraintToReactionMatrix[idx], axis=0) == 0)[0].shape[0]
nConstraints = len(sim_data.process.metabolism.constraintIdList) - disabled
nKcatOnly = int(sum(self.constraintIsKcatOnly[np.where(np.sum(self.constraintToReactionMatrix[idx], axis = 0))[0]]))
nKmKcat = nConstraints - nKcatOnly

print "Number of constraints: %i" % nConstraints
print "Number of kcat only constraints: %i" % nKcatOnly
print "Number of Km and kcat constraints: %i" % nKmKcat