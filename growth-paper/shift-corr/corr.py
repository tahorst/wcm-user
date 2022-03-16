#! /usr/bin/env python

"""
Data comes from growth_trajectory changes on user-shift-corr.

On paper-2022-long:
DESC="Long lead in down and up shifts with regulation" PARALLEL_PARCA=1 RUN_AGGREGATE_ANALYSIS=0 \
  N_GENS=28 N_INIT_SIMS=32 \
  VARIANT=timelines FIRST_VARIANT_INDEX=27 LAST_VARIANT_INDEX=27 \
  TRNA_ATTENUATION=1 PPGPP_REGULATION=1 MECHANISTIC_TRANSLATION_SUPPLY=1 MECHANISTIC_AA_TRANSPORT=1 AA_SUPPLY_IN_CHARGING=1 \
  D_PERIOD_DIVISION=1 MECHANISTIC_REPLISOME=0 TIMESTEP_MAX=1 \
  python runscripts/fireworks/fw_queue.py

no-ppgpp file: PPGPP_REGULATION=0
no-mech file: MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 AA_SUPPLY_IN_CHARGING=0
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
# DATA = os.path.join(FILE_LOCATION, 'no-ppgpp-growth_trajectory-shifts.tsv')
# DATA = os.path.join(FILE_LOCATION, 'no-mech-growth_trajectory-shifts.tsv')
DATA = os.path.join(FILE_LOCATION, 'growth_trajectory-shifts.tsv')
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)

REMOVE_DATA = ['variant', 'seed', 'pre time']


def load_data():
    with open(DATA) as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        data = np.array(list(reader), float).T

    return dict(zip(headers, data))

def modify_data(data):
    for removed in REMOVE_DATA:
        data.pop(removed)

    intercept = -0.5  # TODO: get this actual value for growth vs rp trend

    data['growth decrease'] = data['pre growth'] - data['min growth growth']
    data['ratio increase'] = data['max ratio ratio'] - data['pre ratio']
    data['ppgpp increase'] = data['max ppgpp ppgpp'] / data['pre ppgpp']
    data['pre ppgpp to ratio'] = data['pre ppgpp'] / data['pre ratio']
    data['pre slope'] = data['pre growth'] / data['pre ratio']
    data['pre slope, adjusted'] = (data['pre growth'] - intercept) / data['pre ratio']
    data['min growth slope'] = data['min growth growth'] / data['min growth ratio']
    data['min growth slope, adjusted'] = (data['min growth growth'] - intercept) / data['min growth ratio']
    data['max ratio slope'] = data['max ratio growth'] / data['max ratio ratio']
    data['max ratio slope, adjusted'] = (data['max ratio growth'] - intercept) / data['max ratio ratio']

def plot_single(data, var):
    n_data = len(data)
    n_rows = int(np.ceil(np.sqrt(n_data)))
    n_cols = int(np.ceil(n_data / n_rows))
    _, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(3*n_cols, 3*n_rows))

    ref_data = data[var]
    for i, (label, d) in enumerate(data.items()):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]

        ax.plot(d, ref_data, 'o')

        r, p = stats.pearsonr(d, ref_data)

        ax.tick_params(labelsize=8)
        ax.set_xlabel(f'{label} {r:.3f}')

def plot_corr(data):
    n_data = len(data)
    _, axes = plt.subplots(nrows=n_data, ncols=n_data, figsize=(3*n_data, 3*n_data))

    for i, (row_label, row_data) in enumerate(data.items()):
        for j, (col_label, col_data) in enumerate(data.items()):
            ax = axes[i, j]

            ax.plot(col_data, row_data, 'o')

            ax.tick_params(labelsize=8)

            if j == 0:
                ax.set_ylabel(row_label)
            if i == n_data - 1:
                ax.set_xlabel(col_label)

def save_fig(filename):
    plt.tight_layout()
    file = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(file)
    plt.close('all')
    print(f'Saved to {file}')


if __name__ == '__main__':
    data = load_data()
    modify_data(data)

    plot_single(data, 'growth recovery time')
    save_fig('growth-recovery-corr.pdf')

    plot_single(data, 'both recovery time')
    save_fig('both-recovery-corr.pdf')

    plot_single(data, 'pre ppgpp')
    save_fig('pre-ppgpp-corr.pdf')

    plot_single(data, 'max ppgpp ppgpp')
    save_fig('max-ppgpp-corr.pdf')

    plot_single(data, 'pre ppgpp to ratio')
    save_fig('pre-ppgpp-ratio-corr.pdf')

    plot_corr(data)
    save_fig('corr.pdf')
