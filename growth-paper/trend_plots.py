#!/usr/bin/env python
"""
Explore trends between growth related data to make a similar comparison like
growth vs R/P ratio.

Requires the directory of all the sims run for the paper as SIM_DIR with
certain descriptions assumed in SIM_DESC to select sets of runs.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


# TODO: highlight sim sets that are expected to be off
SIM_DIR = '/home/travis/scratch/wcEcoli_out/'
CONDITION_SIMS = [
    'Amino_acid_combinations_in_media',  # TODO: generate file based on the right number of gens (after shift adjustment)
    'Add_one_amino_acid_shift',
    'Remove_one_amino_acid_shift',
    'Conditions_with_regulation',
    ]
PERTURBATION_SIMS = [
    'ppGpp_sensitivity',
    'ppGpp_limitations_-_low_ppGpp',
    'Remove_amino_acid_inhibition',
    'Amino_acid_synthesis_network_sensitivity_-_glt',
    'ppGpp_limitations_-_high_ppGpp',
    'Amino_acid_synthesis_network_sensitivity_-_control',
    'ppGpp_limitations_-_normal_ppGpp',
    'ppGpp_sensitivity_-_no_mechanistic_transport',
    'Amino_acid_combinations_in_media_without_regulation_or_charging',
    'Conditions_without_regulation_or_charging',
    'Amino_acid_combinations_in_media_-_no_mechanistic_transport',
    'Remove_amino_acid_inhibition_-_no_ppgpp',
    ]
FILE_PATH = 'plotOut/{}.tsv'
OUTPUT_FILE = 'growth-trends.pdf'


def load_datasets(sims):
    return {sim_desc: load_data(sim_desc) for sim_desc in sims}

def load_data(desc, filename='rp_data'):
    dirs = os.listdir(SIM_DIR)
    for d in dirs:
        if d.endswith(desc):
            path = os.path.join(SIM_DIR, d, FILE_PATH.format(filename))
            break
    else:
        raise ValueError(f'{desc} not found in sim directory')

    data = {}
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')

        headers = next(reader)
        columns = np.array(list(reader), float).T

        for header, column in zip(headers, columns):
            data[header] = column

    return data

def plot(all_data):
    # TODO: take list of x and y and plot on subplots
    x_key = 'rp_ratio'
    y_keys = list(next(iter(all_data.values())).keys())
    n_keys = len(y_keys)
    rows = int(np.ceil(np.sqrt(n_keys)))
    cols = int(np.ceil(n_keys / rows))
    scale = 3

    _, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(cols*scale, rows*scale))

    for i, y_key in enumerate(y_keys):
        row = i % rows
        col = i // rows
        ax = axes[row, col]

        for sim, data in all_data.items():
            ax.plot(data[x_key], data[y_key], 'o', alpha=0.4, label=sim)

        format_ax(ax, x_key, y_key)

    save_fig(OUTPUT_FILE)

def format_ax(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def format_plot():
    plt.legend(fontsize=8, frameon=False)

def save_fig(output_file):
    plt.tight_layout()
    plt.savefig(output_file)
    print(f'Saved to {output_file}')
    plt.close('all')


if __name__ == '__main__':
    condition_data = load_datasets(CONDITION_SIMS)

    plot(condition_data)
