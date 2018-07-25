from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np

## A, B single
text = 'single'
ks = np.array([5, 3])
kd = np.array([0.1, 0.6])
tu = np.array([[1, 0], [0, 1]])

## A, B, AB TUs (underdetermined)
text = 'A, B, AB'
# arbitrary
ks = np.array([4, 2, 1])
kd = np.array([0.087, 2.14, 0.246])

# set deg/exp for joint TU as geometric mean - no way to determine what exp should be
ks = np.array([5-0.246, 3-0.246, 0.246])
kd = np.array([0.097, 0.6885, 0.246])


# # ks = np.array([4.5, 1.5, 5])
# # kd = np.array([0.1, 0.3, 1])
#
# ks = np.array([4.835, 1, 1])
# kd = np.array([0.1075, 0.1, 0.1])
#
tu = np.array([[1, 0, 1], [0, 1, 1]])

# A, AB TUs (fully determined)
# text = 'A, AB'
# ks = np.array([2, 3])
# kd = np.array([0.0445, 0.6])
# tu = np.array([[1, 1], [0, 1]])

# AB TU (overdetermined)
# text = 'AB'
# ks = np.array([3.88])
# kd = np.array([0.246])
# tu = np.array([[1, 1]])

def dRdt(R, t):
	return ks - kd*R

Ro = np.zeros(ks.shape)
t = range(1000)

sol = odeint(dRdt, Ro, t)

# import ipdb; ipdb.set_trace()
print text
print 'Final gene number: %s' % (np.dot(tu, sol[-1].reshape(-1,1)))
print 'Final gene deg rate: %s' % (np.dot(tu, (sol[-1]*kd).reshape(-1,1)))
print 'Final TU number: %s' % (sol[-1])

# print sol[-1]

# plt.plot(t, sol)
# plt.show()
