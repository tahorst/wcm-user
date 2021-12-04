#! /usr/bin/env python

"""
Plots to explore growth rate vs R/P ratio in 4 conditions (all AA, +12 AA, +6 AA, no AA).
Data is saved from generation 5+ from sims run for the growth paper (11da79b362df16f5730c5df199421962d7a6e5fc):

DESC="Minimal to rich media conditions" PARALLEL_PARCA=1 RUN_AGGREGATE_ANALYSIS=0 \
  N_GENS=24 N_INIT_SIMS=24 \
  VARIANT=remove_aas_shift FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=3 \
  TRNA_ATTENUATION=1 PPGPP_REGULATION=1 MECHANISTIC_TRANSLATION_SUPPLY=1 AA_SUPPLY_IN_CHARGING=1 \
  D_PERIOD_DIVISION=1 MECHANISTIC_REPLISOME=0 TIMESTEP_MAX=1 \
  python runscripts/fireworks/fw_queue.py

all_growth and all_ratio are saved from growth_trajectory.py variant plot
"""

import os
from typing import Iterable, List

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'plots')
os.makedirs(OUTPUT_DIR, exist_ok=True)


## Helper functions ##

def ma(xs: Iterable[np.ndarray], window=501) -> List[np.ndarray]:
    convolution_array = np.ones(window) / window
    return [
        np.convolve(x, convolution_array, mode='valid')
        if len(x) > window else np.array([0])
        for x in xs
        ]

def mean_std(xs):
    x = np.hstack(xs)
    return x.mean(), x.std()


## Disk functions ##

def load_files(attr: str):
    filenames = [os.path.join(FILE_LOCATION, f'{attr}-{i}.npy') for i in range(4)]
    return [np.load(f, allow_pickle=True) for f in filenames]

def save_fig(filename: str):
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, filename)
    print(f'Saving figure to {path}')
    plt.savefig(path)


## Plot helper functions ##

def plot_binned(xs, ys, bins=10):
    if isinstance(xs, list):
        x = np.hstack(np.hstack(xs))
        y = np.hstack(np.hstack(ys))
    else:
        x = np.hstack(xs)
        y = np.hstack(ys)

    sort_idx = np.argsort(x)
    size = len(sort_idx) // bins

    x_bin = np.array([np.mean(x[sort_idx[i * size:(i + 1) * size]]) for i in range(bins)])
    y_bin = np.array([np.mean(y[sort_idx[i * size:(i + 1) * size]]) for i in range(bins)])

    trace, = plt.plot(x_bin, y_bin, 'o')

    line = stats.linregress(x_bin[2:-2], y_bin[2:-2])
    x_regress = np.array([x_bin.min(), x_bin.max()])
    y_regress = line.slope * x_regress + line.intercept
    plt.plot(x_regress, y_regress, color=trace.get_color())


## Saved plot functions ##

def plot_all(xs: List[np.ndarray], ys: List[np.ndarray]):
    for x, y in zip(xs, ys):
        x_stack = np.hstack(ma(x))
        y_stack = np.hstack(ma(y))
        plt.plot(x_stack, y_stack, linewidth=1, alpha=0.2)
        plot_binned(x, y)

    plot_binned(xs, ys)

    save_fig('all.pdf')

def plot_traces(xs, ys, variant=0, length=5000):
    x = xs[variant]
    y = ys[variant]

    rows = x.shape[0]
    cols = np.max([int(len(x_) / length) for x_ in x])

    x_mean, x_std = mean_std(x)
    y_mean, y_std = mean_std(y)

    _, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(2*cols, 2*rows))

    for row in range(rows):
        for col in range(cols):
            ax = axes[row, col]
            start = col * length
            end = start + length
            if start > len(x[row]):
                continue
            x_split = x[row][start:end]
            y_split = y[row][start:end]

            ax.plot(x_split, y_split, linewidth=0.5)
            try:
                ax.plot(x_split[0], y_split[0], 'go')
            except:
                import ipdb; ipdb.set_trace()
            ax.plot(x_split[-1], y_split[-1], 'ro')
            ax.plot(x_mean, y_mean, 'ko')

            ax.tick_params(labelsize=6)
            ax.set_xlim([x_mean - 2*x_std, x_mean + 2*x_std])
            ax.set_ylim([y_mean - 2*y_std, y_mean + 2*y_std])

    save_fig(f'traces-{variant}.pdf')

def plot_quiver(xs, ys, variant=0, future=500, n_bins=10):
    x = ma(xs[variant])
    y = ma(ys[variant])
    x_mean, x_std = mean_std(x)
    y_mean, y_std = mean_std(y)

    x_diff = np.hstack([x_[future:] - x_[:-future] for x_ in x])
    y_diff = np.hstack([y_[future:] - y_[:-future] for y_ in y])
    x_stack = np.hstack([x_[:-future] for x_ in x])
    y_stack = np.hstack([y_[:-future] for y_ in y])

    mask = (x_stack > x_mean - 2*x_std) & (x_stack < x_mean + 2*x_std) & (y_stack > y_mean - 2*y_std) & (y_stack < y_mean + 2*y_std)
    x_diff = x_diff[mask]
    y_diff = y_diff[mask]
    x_stack = x_stack[mask]
    y_stack = y_stack[mask]

    x_bin = np.floor((x_stack - x_stack.min()) / (x_stack.max() - x_stack.min()) * n_bins).astype(int)
    x_bin[x_bin == n_bins] = n_bins - 1
    y_bin = np.floor((y_stack - y_stack.min()) / (y_stack.max() - y_stack.min()) * n_bins).astype(int)
    y_bin[y_bin == n_bins] = n_bins - 1
    bins = n_bins * x_bin + y_bin

    quiver_x = []
    quiver_y = []
    quiver_u = []
    quiver_v = []
    for bin in np.unique(bins):
        mask = bins == bin
        x_bin_mean = x_stack[mask].mean()
        y_bin_mean = y_stack[mask].mean()
        quiver_x.append(x_bin_mean)
        quiver_y.append(y_bin_mean)
        quiver_u.append(x_diff[mask].mean() / x_bin_mean)
        quiver_v.append(y_diff[mask].mean() / y_bin_mean)

    plt.figure()
    plt.quiver(quiver_x, quiver_y, quiver_u, quiver_v)

    save_fig(f'quiver-{variant}.pdf')


if __name__ == '__main__':
    ratio = load_files('ratio')
    growth = load_files('growth')

    plot_all(ratio, growth)
    plot_traces(ratio, growth)
    plot_traces(ratio, growth, variant=3)
    plot_quiver(ratio, growth)
    plot_quiver(ratio, growth, variant=3)
