#! /usr/bin/env python

# using general pyenv with pyts installed

import os

import matplotlib.pyplot as plt
import numpy as np
from pyts import decomposition


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))


leua = np.load('leu-leua.npy')
ssa = decomposition.ssa.SingularSpectrumAnalysis(window_size=5)

# Might be able to fit to multiple samples (seeds)
fit = ssa.fit_transform(leua.reshape(1, -1))

x = np.arange(len(leua))
plt.figure()

plt.plot(x, leua)
plt.plot(x, fit[0, :])
plt.plot(x, fit[1, :])

plt.savefig(os.path.join(FILE_LOCATION, 'ssa.pdf'))


plt.figure()

plt.plot(x, fit[1, :])
plt.plot(x, fit[2, :])

plt.savefig(os.path.join(FILE_LOCATION, 'ssa-1.pdf'))
