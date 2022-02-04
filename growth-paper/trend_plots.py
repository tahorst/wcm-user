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


# TODO: set back to sherlock dirs
# SIM_DIR = '/home/travis/scratch/wcEcoli_out/'
SIM_DIR = '/home/travis/wcEcoli/out/'
SIM_DESC = [
    'ecocyc',
    # 'Conditions_without_regulation_or_charging',
    # 'Conditions_with_regulation',
    # 'Add_one_amino_acid_shift',
    # 'Remove_one_amino_acid_shift',
    # 'ppGpp_sensitivity_-_no_mechanistic_transport',
    # 'Amino_acid_combinations_in_media',
    # 'Remove_amino_acid_inhibition_-_no_ppgpp',
    # 'Remove_amino_acid_inhibition',
    ]
FILE_PATH = 'plotOut/{}.tsv'
OUTPUT_FILE = 'growth-trends.pdf'


def load_all_data():
    return {sim_desc: load_data(sim_desc) for sim_desc in SIM_DESC}

def load_data(desc, filename='rp_data', ):
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

def plot_all(all_data):
    # TODO: take list of x and y and plot on subplots
    x_key = 'rp_ratio'
    y_key = 'growth_rate'

    plt.figure()
    ax = plt.gca()

    for sim, data in all_data.items():
        ax.plot(data[x_key], data[y_key], 'o', label=sim)

    format_plot(x_key, y_key)
    save_fig(OUTPUT_FILE)

def format_plot(xlabel, ylabel):
    plt.legend(fontsize=8, frameon=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()

def save_fig(output_file):
    plt.savefig(output_file)
    print(f'Saved to {output_file}')
    plt.close('all')


if __name__ == '__main__':
    data = load_all_data()
    plot_all(data)
