#!/usr/bin/env python
"""
Plot data from ppgpp_limitations variants to show effects of adjustments to
ppGpp, ribosome expression and enzyme expression.

Data comes from paper runs (compiled from separate runs for low, normal and high ppGpp)
run on hash 38c67ccb5e0068b1550e891cace8d4d583d03f39:
DESC="ppGpp limitations" PARALLEL_PARCA=1 RUN_AGGREGATE_ANALYSIS=0 \
  N_GENS=8 N_INIT_SIMS=4 \
  VARIANT=ppgpp_limitations FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=110 \
  TRNA_ATTENUATION=1 PPGPP_REGULATION=1 MECHANISTIC_TRANSLATION_SUPPLY=1 MECHANISTIC_AA_TRANSPORT=1 AA_SUPPLY_IN_CHARGING=1 \
  D_PERIOD_DIVISION=1 MECHANISTIC_REPLISOME=0 TIMESTEP_MAX=1 \
  python runscripts/fireworks/fw_queue.py

Analysis from 4781718600d0e8af33a51609fac77caac205bbb1
(run for low, normal and high and output compiled into one tsv):
python /home/users/thorst/wcEcoli/runscripts/manual/analysisVariant.py out/20220120.060837__ppGpp_limitations_-
_high_ppGpp/ -p ppgpp_conc --generation-path-range 2 8 -o 2+gen-

TODO:
    add grouping and subgroup labels
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


# Paths
FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_FILE = os.path.join(FILE_LOCATION, '2+gen-ppgpp_conc.tsv')
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Key is ppGpp, ribosome, enzyme adjustment (-1: low, 0: no adjust, 1: high)
# Value is (condition, factor) key in data
ADJUSTMENT_KEYS = {
    # Controls
    (-1, 0, 0): (3, 0),
    (0, 0, 0): (15, 0),
    (1, 0, 0): (31, 0),
    # Low ppGpp adjustmnets
    (-1, -1, 0): (7, -3),
    (-1, 0, 1): (6, 3),
    # Normal ppGpp adjustmnets
    (0, -1, 0): (19, -3),
    (0, 1, 0): (19, 3),
    (0, 0, -1): (18, -3),
    (0, 0, 1): (18, 3),
    # High ppGpp adjustmnets
    (1, 1, 0): (35, 3),
    (1, 0, -1): (34, -3),
    }
PAIRWISE_COMPARISONS = [
    # Low ppGpp
    [(0, 1, 0), (-1, 0, 1)],  # and high ribosome
    [(0, 0, -1), (-1, -1, 0)],  # and low enzymes
    # High ppGpp
    [(0, -1, 0), (1, 0, -1)],  # and low ribosomes
    [(0, 0, 1), (1, 1, 0)],  # and high enzymes
    ]
# Control, rib, enz, ppGpp and rib, ppGpp and enz, all adjusted
GROUP_COMPARISONS = [
    [(0, 0, 0), (0, 1, 0), (0, 0, -1), (-1, 0, 1), (-1, -1, 0), (-1, 0, 0)],  # low ppGpp -> high rib, low enz
    [(0, 0, 0), (0, -1, 0), (0, 0, 1), (1, 0, -1), (1, 1, 0), (1, 0, 0)],  # high ppGpp -> low rib, high enz
    ]

def load_data():
    data = {}
    with open(DATA_FILE) as f:
        reader = csv.reader(f, delimiter='\t')

        conditions = np.array(next(reader)[1:], int)
        factors = np.array(next(reader)[1:], int)
        keys = list(zip(conditions, factors))

        for line in reader:
            label = line[0]
            values = np.array(line[1:], float)
            data[label] = dict(zip(keys, values))

    return data

def plot_data(data, labels, comparisons):
    n_plots = len(labels)
    n_rows = int(np.ceil(np.sqrt(n_plots)))
    n_cols = int(np.ceil(n_plots / n_rows))

    _, axes = plt.subplots(n_rows, n_cols, figsize=(2*n_cols, 2*n_rows))
    hide_axes = np.ones_like(axes, dtype=bool)

    for i, label in enumerate(labels):
        values = data[label]
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        hide_axes[row, col] = False

        extracted = np.array([[values[ADJUSTMENT_KEYS[k]] for k in comp] for comp in comparisons])
        x = np.arange(extracted.shape[0])
        width = 0.8 / extracted.shape[1]
        offsets = np.arange(extracted.shape[1]) * width - 0.4 + width/2

        for i, offset in enumerate(offsets):
            ax.bar(x + offset, extracted[:, i], width)

        ax.set_ylabel(label, fontsize=8)
        ax.tick_params(labelsize=6)

    # Hide unused axes
    for ax in axes[hide_axes]:
        ax.set_visible(False)

def save_fig(output_file):
    plt.tight_layout()

    filename = os.path.join(OUTPUT_DIR, output_file)
    plt.savefig(filename)
    print(f'Saved to {filename}')
    plt.close('all')


if __name__ == '__main__':
    data = load_data()
    labels = [k for k in data if 'mean' in k]

    plot_data(data, labels, PAIRWISE_COMPARISONS)
    save_fig('pairwise.pdf')

    plot_data(data, labels, GROUP_COMPARISONS)
    save_fig('group.pdf')
