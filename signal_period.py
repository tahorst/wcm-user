#! /usr/bin/env python

"""
Exploring how to extract period from a signal using FFT and periodogram analysis.
"""

import os

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_FILE = os.path.join(FILE_LOCATION, 'period.pdf')


def plot(axes, xmax):
    xsampled = 10000
    freq = xmax / xsampled

    x = np.linspace(0, xmax, xsampled)
    y = np.sin(0.4 * np.pi * x) + np.sin(2 * np.pi * x) + np.sin(10 * np.pi * x) + np.sin(20 * np.pi * x)  # freq of 0.2, 1, 5, 10

    f, Pxx = signal.periodogram(y, fs=1 / freq)
    fft_f = np.fft.fftfreq(xsampled)
    fft_mag = np.real(np.fft.fft(y))
    fft_mask = fft_f > 0

    ax_trace, ax_periodgram, ax_fft = axes

    ax_trace.plot(x, y)

    ax_periodgram.plot(f, Pxx)
    ax_periodgram.set_xscale('log')

    ax_fft.plot(fft_f[fft_mask] / freq, fft_mag[fft_mask], label='power')
    ax_fft.plot(fft_f[fft_mask] / freq, fft_mag[fft_mask] / fft_f[fft_mask] * freq, label='normalized')
    ax_fft.set_xscale('log')
    ax_fft.legend()

_, axes = plt.subplots(3, 2, figsize=(6, 9))

plot(axes[:, 0], 10)
plot(axes[:, 1], 8)

plt.tight_layout()
plt.savefig(OUTPUT_FILE)
