'''
Trying to simulate control by ppGpp due to different time scales of small molecule,
RNA and protein production.  Need to implement a more complicated model and account
for dilution
'''

from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np

kg1 = 8.
kg2 = 10.
kgm = 5.

kp = 0.1
kpm1 = 1.
kpm2 = 1.

km = 0.25
kmm = 1.

kn = 1.

ka = 1.

def dcdt(c, t):
	# concentrations
	g = c[0]
	p = c[1]
	m = c[2]
	n = c[3]
	a = c[4]
	u = c[5]

	# derivatives
	dg = kg1 - kg2 * p / (kgm + p)
	dp = kp * a / (kpm1 + a) * m / (kpm2 + m)
	dm = km * n * g / (kmm + g)
	dn = kn - dm
	da = ka - dp

	return np.array([dg, dp, dm, dn, da])

co = np.ones(6)
t = range(1000)

sol = odeint(dcdt, co, t)

import ipdb; ipdb.set_trace()
# print text
# print 'Final gene number: %s' % (np.dot(tu, sol[-1].reshape(-1,1)))
# print 'Final gene deg rate: %s' % (np.dot(tu, (sol[-1]*kd).reshape(-1,1)))
# print 'Final TU number: %s' % (sol[-1])

# print sol[-1]

plt.plot(t, sol)
plt.show()
