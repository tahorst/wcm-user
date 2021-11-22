#!/usr/bin/env python
"""
Plot the growth rate vs R/P ratio in multiple conditions to show the new
environmental capabilities of the model.

Requires the directory of all the sims run for the paper as SIM_DIR with
certain descriptions assumed to select sets of runs.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


SIM_DIR = '/home/travis/scratch/wcEcoli_out/growth-paper/'
BASE_SIM_DIR = 'Conditions_without_regulation_or_charging'
NEW_C_SOURCE_DIR = 'Conditions_with_regulation'
ADD_ONE_DIR = 'Add_one_amino_acid_shift'
REMOVE_ONE_DIR = 'Remove_one_amino_acid_shift'
FILE_PATH = 'plotOut/growth_trajectory.tsv'
OUTPUT_FILE = 'combined-growth-rp.pdf'

CONTROL_IDX = 19
GROWTH_HEADER = 'Growth'
RP_HEADER = 'R/P ratio'
MEAN_HEADER = ' mean'
STD_HEADER = ' std'



def load_data(desc):
    dirs = os.listdir(SIM_DIR)
    for d in dirs:
        if d.endswith(desc):
            path = os.path.join(SIM_DIR, d, FILE_PATH)
            break
    else:
        raise ValueError(f'{desc} not found in sim directory')

    data = {}
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')

        headers = next(reader)
        for row in reader:
            data[float(row[0])] = dict(zip(headers[1:], np.array(row[1:], float)))

    return data

def plot(data, variants=None, exclude=None, std=True, label=None, options=None):
    def extract_data(header, data_type):
        return np.array([data[v][header + data_type] for v in variants if v not in exclude])

    if variants is None:
        variants = list(data)
    if exclude is None:
        exclude = set()
    if options is None:
        options = {}

    rp_ratio = extract_data(RP_HEADER, MEAN_HEADER)
    growth = extract_data(GROWTH_HEADER, MEAN_HEADER)

    if std:
        rp_ratio_std = extract_data(RP_HEADER, STD_HEADER)
        growth_std = extract_data(GROWTH_HEADER, STD_HEADER)
    else:
        rp_ratio_std = None
        growth_std = None

    plt.errorbar(rp_ratio, growth, xerr=rp_ratio_std, yerr=growth_std, fmt='o', label=label, **options)


if __name__ == '__main__':
    no_regulation = load_data(BASE_SIM_DIR)
    regulation = load_data(NEW_C_SOURCE_DIR)
    add_one = load_data(ADD_ONE_DIR)
    remove_one = load_data(REMOVE_ONE_DIR)

    one_aa_options = dict(alpha=0.5, markersize=4)

    plt.figure()

    plot(no_regulation, label='Original conditions', variants=np.arange(3))
    plot(regulation, label='New carbon sources', variants=np.arange(3, 5))
    plot(add_one, variants=[CONTROL_IDX], label='Minimal + glc with growth regulation')
    plot(remove_one, variants=[CONTROL_IDX], label='Rich + glc with growth regulation')
    plot(add_one, std=False, exclude=[CONTROL_IDX], label='Add one AA to minimal with growth regulation', options=one_aa_options)
    plot(remove_one, std=False, exclude=[CONTROL_IDX], label='Remove one AA from rich with growth regulation', options=one_aa_options)

    plt.legend(fontsize=8, frameon=False)
    plt.xlabel('RNA/protein mass ratio')
    plt.ylabel('Growth rate (1/hr)')
    plt.xlim([0, 0.6])
    plt.ylim([0, 2])

    plt.tight_layout()
    plt.savefig(OUTPUT_FILE)
    print(f'Saved to {OUTPUT_FILE}')
