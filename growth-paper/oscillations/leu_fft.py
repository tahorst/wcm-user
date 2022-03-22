#! /usr/bin/env python

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

OPTIONS = [
    dict(real=True, normalize_input=True, taper=True, ma_window=None, normalize_output=True, ma_output=None, period=True, n_points=None),
    dict(real=False, normalize_input=True, taper=True, ma_window=None, normalize_output=True, ma_output=None, period=True, n_points=None),
    dict(real=True, normalize_input=False, taper=True, ma_window=None, normalize_output=True, ma_output=None, period=True, n_points=None),
    dict(real=True, normalize_input=True, taper=False, ma_window=None, normalize_output=True, ma_output=None, period=True, n_points=None),
    dict(real=True, normalize_input=True, taper=True, ma_window=1001, normalize_output=True, ma_output=None, period=True, n_points=None),
    dict(real=True, normalize_input=True, taper=True, ma_window=None, normalize_output=True, ma_output=201, period=True, n_points=None),
    dict(real=True, normalize_input=True, taper=True, ma_window=None, normalize_output=False, ma_output=None, period=True, n_points=None),
    dict(real=True, normalize_input=True, taper=True, ma_window=None, normalize_output=True, ma_output=None, period=False, n_points=None),
    dict(real=True, normalize_input=True, taper=True, ma_window=None, normalize_output=True, ma_output=None, period=True, n_points=7200),
    dict(real=True, normalize_input=True, taper=False, ma_window=None, normalize_output=True, ma_output=201, period=True, n_points=None),
    dict(real=True, normalize_input=True, taper=True, ma_window=1001, normalize_output=True, ma_output=201, period=True, n_points=None),
    dict(real=True, normalize_input=True, taper=True, ma_window=None, normalize_output=True, ma_output=201, period=True, n_points=7200),
    dict(real=True, normalize_input=False, taper=True, ma_window=None, normalize_output=True, ma_output=201, period=True, n_points=7200),
    ]


wt = np.load(os.path.join(FILE_LOCATION, 'leu-wt.npy'))
leua = np.load(os.path.join(FILE_LOCATION, 'leu-leua.npy'))

data = leua
data_backup = data.copy()

def fft(data, real=True, normalize_input=True, taper=True, ma_window=None, normalize_output=True, ma_output=None, period=True, n_points=None):
    if real:
        fun = np.fft.rfft
        freq_fun = np.fft.rfftfreq
    else:
        fun = np.fft.fft
        freq_fun = np.fft.fftfreq

    if normalize_input:
        data = (data - data.mean()) / data.std()

    if taper:
        data = data * signal.cosine(len(data))

    if ma_window:
        data = np.convolve(data, np.ones(ma_window) / ma_window, mode='valid')

    if n_points is None:
        n_points = len(data)

    f = freq_fun(n_points)
    Pxx = np.abs(fun(data, n_points))

    if normalize_output:
        Pxx *= f

    if ma_output:
        Pxx = np.convolve(Pxx, np.ones(ma_output) / ma_output, mode='same')

    if period:
        f = 1 / f / 3600

    return f, Pxx


n_plots = len(OPTIONS)
rows = int(np.ceil(np.sqrt(n_plots)))
cols = int(np.ceil(n_plots / rows))

_, axes = plt.subplots(rows, cols, figsize=(4*cols, 2*rows))

for ax, options in zip(axes.flatten(), OPTIONS):
    f, Pxx = fft(data, **options)
    ax.plot(f, Pxx, alpha=0.5)

    # twin = ax.twinx()
    f, Pxx = fft(wt, **options)
    ax.plot(f, Pxx, 'k', alpha=0.5)

    ax.set_xscale('log')
    ax.set_yscale('log')
    # twin.set_yscale('log')
    title = ', '.join([f'{k}={v}' for k, v in options.items() if OPTIONS[0][k] != v])
    ax.set_title(title, fontsize=6)

    if not np.all(data == data_backup):
        print(options)

plt.tight_layout()
plt.savefig(os.path.join(FILE_LOCATION, 'fft.pdf'))

## autocorrelation
wt_corr = signal.correlate(wt, wt, method='fft')
adjusted_wt = wt_corr[(len(wt_corr)-1)//2:] / wt_corr.mean()
leua_corr = signal.correlate(leua, leua, method='fft')
adjusted_leua = leua_corr[(len(leua_corr)-1)//2:] / leua_corr.mean()

plt.figure()
plt.plot(np.arange(len(adjusted_wt)) / 3600, adjusted_wt)
plt.plot(np.arange(len(adjusted_leua)) / 3600, adjusted_leua)
plt.savefig(os.path.join(FILE_LOCATION, 'autocorr.pdf'))
