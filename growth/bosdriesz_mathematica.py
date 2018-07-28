'''
Direct port of Bosdriesz from mathematica file
'''

from __future__ import division

from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
import numpy as np

nAA = 20
tmax = 10000
cellVolume = 2.5e-15
nAvogadro = 6e23
proteinContent = 15e8

# indices for concentrations
aa_indices = range(nAA)
ta_indices = range(nAA, 2 * nAA)
ppgpp_index = ta_indices[-1] + 1
r_index = ppgpp_index + 1

# pars
e = 0.05
kn = 0.2 * np.random.lognormal(np.log(1), 0.2, 20)
kIa = 100
sTot = 1
ks = 100
kMaa = 100
kMtf = 1
krib = 20
krta = 1
krt = 500
kRelA = 75
RelAtot = 100 / (nAvogadro * cellVolume / 1e6)
kDRelA = 0.26
vSpoTSynt = 0.001
kSpoTdeg = np.log(2) / 30
vInitMax = 2000
rnapF = 1
kMrrn = 20
kIppGpp = 1
nppGpp = 1
gammamax = 1
tau = 0.5
f = 0.05
bm = proteinContent / (cellVolume*nAvogadro/1e6)
nARib = 7459 * 1.65
nAmet = 300
rmax = proteinContent / (cellVolume*nAvogadro/1e6) / (7459*1.65)

def dcdt_ode(t, c):
	return dcdt(c, t)

def dcdt(c, t):
	dc = np.zeros_like(c)

	# shift - not in mathematica file
	shift = 1
	# if t > 2000:
	# 	shift = 0.1
	# shift = 0.5
	# if t > 2000:
	# 	shift = 2

	aa = c[aa_indices]
	taa = c[ta_indices]
	ppGpp = c[ppgpp_index]
	r = c[r_index]

	tf = tau * r - taa

	vAAsynt = shift * bm * e * kn * (1 - r/rmax) / (nAmet * (1 + aa / kIa))
	vtRNAchar = ks * sTot * tf * aa / (kMaa * kMtf * (tf / kMtf + aa / kMaa + tf * aa / kMaa / kMtf))
	numeratorRibosome = 1 + np.sum(f * (krta/taa + tf/taa*krta/krt))
	vR = krib*r / numeratorRibosome
	mu = vR / bm
	vrrnInit = vInitMax*rnapF/(kMrrn + rnapF) / (1 + ppGpp / kIppGpp * nppGpp) / (cellVolume*nAvogadro/1e6)
	vribosome = vR
	vRsynt = min(vrrnInit, gammamax*vR/nARib)
	vRdilution = r * mu
	rtfSolutions = r*(f*tf/taa*krta/krt)/numeratorRibosome
	rtfTot = np.sum(rtfSolutions)
	vRelA = kRelA * RelAtot / (1 + kDRelA/rtfTot)
	vSpoTdeg = kSpoTdeg * ppGpp

	odesAA = vAAsynt - vtRNAchar
	odesTAA = vtRNAchar - f*vribosome
	odesppGpp = vRelA + vSpoTSynt - vSpoTdeg
	odesMacromolComp = vRsynt - vRdilution

	# derivatives
	dc[aa_indices] = odesAA
	dc[ta_indices] = odesTAA
	dc[ppgpp_index] = odesppGpp
	dc[r_index] = odesMacromolComp

	return dc


co = np.zeros(2*nAA + 2)
co[aa_indices] = kIa  # aa (100)
co[ta_indices] = 0.1*tau*rmax  # charged tRNA (4.063)
co[ppgpp_index] = kIppGpp  # ppGpp (1)
co[r_index] = 0.2*rmax  # ribosome (16.25)
t = np.linspace(0,1000,10000)

# ode15s = ode(dcdt_ode)
# ode15s.set_integrator('vode', method='bdf', order=15)
# ode15s.set_initial_value(co, t[0])
#
# sol = np.zeros((len(t), len(co)))
# sol[0,:] = co
# for idx, next_time in enumerate(t[1:]):
# 	ode15s.integrate(next_time)
# 	sol[idx + 1, :] = ode15s.y

sol = odeint(dcdt, co, t)

plt.subplot(3,1,1)
plt.plot(t, sol[:,ppgpp_index])
plt.plot(t, sol[:,r_index])
plt.legend(['ppGpp', 'ribosomes'], fontsize=6)

plt.subplot(3,1,2)
plt.plot(t, sol[:,aa_indices])
plt.ylabel('[AA]')

plt.subplot(3,1,3)
plt.plot(t, sol[:,ta_indices])
plt.ylabel('[charged tRNA]')
plt.show()