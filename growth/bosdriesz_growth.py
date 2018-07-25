'''
Recapitulating model in Bosdriesz 2015 FEBS paper
Uses ppGpp kinetics to model up/down shifts
'''

from __future__ import division

from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np

n_aa = 20
frac_ribosomes = 0.5
rho = 1/20

# Constants from paper
## amino acids, TODO - values for each AA
kn = np.random.lognormal(1, 0.15, 20)  # AA production rate, from text under S4 equations
mi = 1/20  # AA fraction
KI = 100  # uM - inhibition constant for synthesis

## tRNA synthetases, TODO - values for each synthetase
kS = 100  # 1/s - kcat for synthetase
Stot = 1  # uM - synthetase conc
KMt = 100  # uM - Michaelis constant for synthetases for AA
KMa = 1  # uM - Michaelis constant for synthetases for uncharged tRNA

## ribosome, TODO - adjust f for each AA
krib = 20  # 1/s - protein elongation rate
f = 1/20  # protein composition per AA
kta = 1  # uM - charged tRNA-ribosome dissociation constant
kt = 500  # uM - uncharged tRNA-ribosome dissociation constant

## ppGpp
kRelA = 75  # 1/s - RelA ppGpp synthesis
RelAtot = 1/15  # uM - RelA conc
KDRelA = 0.26  # uM - uncharged tRNA conc for half maximal ppGpp production
vSpoTsyn = 0.001  # uM/s - Spot ppGpp synthesis
kSpoTdeg = np.log(2) / 30  # 1/s - ppGpp half life

## synthesis
Vmaxrrn = 2000  # 1/s - max rrn rate
RNAPf = 1  # uM - free RNAP conc
KMrrn = 20  # uM - Michaelis constant for rrn-promoters
KippGpp = 1  # uM - inhibition constant for rrn promoters

## others
p = 1e6  # uM - protein concentration
Vcell = 2.5e-9  # uL - cell volume
Nav = 6.022e23  # Avogadro's number
Nar = 12307  # AA per ribosome
Nam = 300  # AA per protein

# indices for concentrations
aa_indices = range(n_aa)
ta_indices = range(n_aa, 2 * n_aa)
ppgpp_index = ta_indices[-1] + 1
r_index = ppgpp_index + 1


# Equations from the paper
# Note: ti has been changed to tu for uncharged tRNA
def dcdt(c, t):
	# amino acid synthesis shift
	kaa = 0.1
	if t > 4000:
		kaa = 1

	dc = np.zeros_like(c)

	# concentrations
	a = c[aa_indices]
	ta = c[ta_indices]
	ppGpp = c[ppgpp_index]
	r = c[r_index]

	# derived concentrations
	ttot = 0.5 * r  # table S1
	tu = ttot - ta  # eq S3c

	rti = r * rho * tu/kt / (1 + ta/kta + tu/kt)  # eq S6b
	rttot = np.sum(rti)

	# rates
	va = kaa * kn * mi / (1 + a/KI)  # eq S4a, modified from text in paragraph
	vta = kS * Stot * (tu/KMt * a/KMa) / (1 + tu/KMt + a/KMa + tu/KMt*a/KMa)  # eq S4b, adjusted
	vrib = krib * r / (1 + np.sum(f * (kta/ta + tu/kt*kta/ta)))  # eq S4c
	vRelA = kRelA * RelAtot * (rttot / KDRelA) / (1 + rttot/KDRelA)  # eq S9a
	vSpoTdeg = kSpoTdeg * ppGpp  # eq S9c
	vrnn = Vmaxrrn * RNAPf / (KMrrn + RNAPf) / (1 + ppGpp / KippGpp)  # eq S13a
	vr = min(vrnn/(Vcell*Nav), Nar*vrib)  # eq S12

	mu = vrib / p  # eq S2
	# print mu

	# derivatives
	dc[aa_indices] = va - vta  # eq S3a
	dc[ta_indices] = vta - f*vrib  # eq S3b
	dc[ppgpp_index] = vRelA + vSpoTsyn - vSpoTdeg  # eq S8
	dc[r_index] = vr - mu*r  # eq S10

	# import ipdb; ipdb.set_trace()

	return dc

co = np.zeros(2*n_aa + 3)
co[aa_indices] = kn  # aa
co[ta_indices] = kn  # charged tRNA
co[ppgpp_index] = 5  # ppGpp
co[r_index] = 15  # ribosome
t = range(10000)

sol = odeint(dcdt, co, t)

# import ipdb; ipdb.set_trace()
# print text
# print 'Final gene number: %s' % (np.dot(tu, sol[-1].reshape(-1,1)))
# print 'Final gene deg rate: %s' % (np.dot(tu, (sol[-1]*kd).reshape(-1,1)))
# print 'Final TU number: %s' % (sol[-1])

# print sol[-1]

plt.subplot(3,1,1)
plt.plot(t, sol[:,ppgpp_index])
plt.plot(t, sol[:,r_index])

plt.subplot(3,1,2)
plt.plot(t, sol[:,aa_indices])
plt.ylabel('[AA]')

plt.subplot(3,1,3)
plt.plot(t, sol[:,ta_indices])
plt.ylabel('[charged tRNA]')
plt.show()

# import ipdb; ipdb.set_trace()