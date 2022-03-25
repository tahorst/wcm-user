#! /usr/bin/env python

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

wt = np.load(os.path.join(FILE_LOCATION, 'leu-wt.npy'))
leua = np.load(os.path.join(FILE_LOCATION, 'leu-leua.npy'))


def plot_all(axes, data, input_ma=3600, deriv_ma=3600, **kwargs):
    plot(axes[0], data, **kwargs)
    plot(axes[1], data, input_ma=input_ma, **kwargs)
    plot(axes[2], data, deriv_ma=deriv_ma, **kwargs)
    plot(axes[3], data, input_ma=input_ma, deriv_ma=deriv_ma, **kwargs)

def plot(ax, data, input_drop=1, input_ma=False, deriv_ma=False, splrep=False):
    if input_drop > 1:
        n_points = len(data) // input_drop
        n_drop = len(data) % input_drop
        if n_drop == 0:
            index = slice(None)
        else:
            index = slice(-n_drop)
        data = data[index].reshape(n_points, -1).mean(1)

    if input_ma:
        input_ma = min(len(data) - 1, input_ma)
        data = np.convolve(data, np.ones(input_ma) / input_ma, 'valid')

    t = np.arange(len(data)) / 3600 * input_drop

    if splrep:
        spline = interpolate.splrep(t, data, k=5)
        y = interpolate.splev(t, spline)
        deriv = interpolate.splev(t, spline, der=2)
    else:
        spline = interpolate.CubicSpline(t, data)
        y = spline(t)
        deriv = spline.derivative(2)(t)

    if deriv_ma:
        deriv_ma = min(len(deriv) - 1, deriv_ma)
        deriv = np.convolve(deriv, np.ones(deriv_ma) / deriv_ma, 'same')

    sign = np.sign(deriv)

    print(np.sum(sign[:-1] != sign[1:]))
    data_range = data.max() - data.min()
    sign[sign > 0] = data.max() - 0.2 * data_range
    sign[sign < 0] = data.min() + 0.2 * data_range

    ax.plot(t, data)
    ax.plot(t, y, alpha=0.5)
    ax.plot(t, sign, alpha=0.2)

    ax.set_ylim([data.min(), data.max()])


_, axes = plt.subplots(nrows=6, ncols=4)

plot_all(axes[0, :], wt)
plot_all(axes[1, :], leua)
plot_all(axes[2, :], wt, splrep=True)
plot_all(axes[3, :], leua, splrep=True)
plot_all(axes[4, :], wt, input_ma=30, deriv_ma=30, input_drop=60)
plot_all(axes[5, :], leua, input_ma=30, deriv_ma=30, input_drop=60)

plt.savefig(os.path.join(FILE_LOCATION, 'interp.pdf'))
