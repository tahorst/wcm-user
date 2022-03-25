#! /usr/bin/env python

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

wt = np.load(os.path.join(FILE_LOCATION, 'leu-wt.npy'))
leua = np.load(os.path.join(FILE_LOCATION, 'leu-leua.npy'))


def plot_all(axes, data, splrep=False):
    plot(axes[0], data, splrep=splrep)
    plot(axes[1], data, input_ma=500, splrep=splrep)
    plot(axes[2], data, deriv_ma=500, splrep=splrep)
    plot(axes[3], data, input_ma=1000, deriv_ma=1000, splrep=splrep)

def plot(ax, data, input_ma=False, deriv_ma=False, splrep=True):
    if input_ma:
        data = np.convolve(data, np.ones(input_ma) / input_ma, 'valid')

    t = np.arange(len(data)) / 3600

    if splrep:
        spline = interpolate.splrep(t, data, k=5)
        y = interpolate.splev(t, spline)
        deriv = interpolate.splev(t, spline, der=2)
    else:
        spline = interpolate.CubicSpline(t, data)
        y = spline(t)
        deriv = spline.derivative(2)(t)

    if deriv_ma:
        deriv = np.convolve(deriv, np.ones(deriv_ma) / deriv_ma, 'same')

    sign = np.sign(deriv)

    print(np.sum(sign[:-1] != sign[1:]))

    ax.plot(t, data)
    ax.plot(t, y, alpha=0.5)
    ax.plot(t, deriv, alpha=0.2)

    ax.set_ylim([data.min(), data.max()])


_, axes = plt.subplots(nrows=4, ncols=4)

plot_all(axes[0, :], wt)
plot_all(axes[1, :], leua)
plot_all(axes[2, :], wt, splrep=True)
plot_all(axes[3, :], leua, splrep=True)

plt.savefig(os.path.join(FILE_LOCATION, 'interp.pdf'))
