#!/usr/bin/env python
"""
Plot the growth rate from fixed ppGpp runs with expression adjustments to
amino acid synthesis enzymes or ribosomes to show limitations.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


# Input paths
SIM_DIR = '/home/travis/scratch/wcEcoli_out/'
SENSITIVITY_DIR = 'ppGpp_sensitivity_-_no_mechanistic_transport'
LIMITATIONS_LOW_DIR = 'ppGpp_limitations_-_low_ppGpp'
LIMITATIONS_HIGH_DIR = 'ppGpp_limitations_-_high_ppGpp'
# TODO: update these dirs
RIBOSOME_LIMITATIONS_INHIBITION_DIR = 'ppGpp_limitations_with_split_rRNA_rProtein_at_high_ppGpp_-_updated'
RIBOSOME_LIMITATIONS_NO_INHIBITION_DIR = 'ppGpp_limitations_with_split_rRNA_rProtein_at_high_ppGpp,_no_translation_inhibition_-_updated'
FILE_PATH = 'plotOut/{}.tsv'

# Output paths
FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)

GROWTH_HEADER = 'Growth'
RP_HEADER = 'R/P ratio'
MEAN_HEADER = ' mean'
STD_HEADER = ' std'


def load_data(desc, filename='2+gen-growth_trajectory'):
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
        for row in reader:
            data[float(row[0])] = dict(zip(headers[1:], np.array(row[1:], float)))

    return data

def plot_data(data):
    labels = list(data.keys())
    y = list(data.values())
    x = np.arange(len(y))

    plt.figure()
    plt.bar(x, y)

    plt.xticks(x, labels, rotation=45)
    plt.ylabel('Growth rate (1/hr)')

    # Remove axes borders
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def save_fig(output_file):
    plt.tight_layout()
    filename = os.path.join(OUTPUT_DIR, output_file + '.pdf')
    plt.savefig(filename)
    print(f'Saved to {filename}')
    plt.close('all')


if __name__ == '__main__':
    # Load data from disk
    sensitivity = load_data(SENSITIVITY_DIR)
    limitations_low = load_data(LIMITATIONS_LOW_DIR)
    limitations_high = load_data(LIMITATIONS_HIGH_DIR)
    ribosome_limit_inhibition = load_data(RIBOSOME_LIMITATIONS_INHIBITION_DIR)
    ribosome_limit_no_inhibition = load_data(RIBOSOME_LIMITATIONS_NO_INHIBITION_DIR)

    # Compile data to plot
    low_data = {
        'Control': sensitivity[4][GROWTH_HEADER + MEAN_HEADER],
        'Low ppGpp': limitations_low[0][GROWTH_HEADER + MEAN_HEADER],
        # 'Decrease enzymes': limitations_low[20][GROWTH_HEADER + MEAN_HEADER],
        'Increase enzymes': limitations_low[25][GROWTH_HEADER + MEAN_HEADER],
        # 'Decrease ribosomes': limitations_low[35][GROWTH_HEADER + MEAN_HEADER],
        'Increase ribosomes': limitations_low[30][GROWTH_HEADER + MEAN_HEADER],
        }
    high_data = {
        'Control': sensitivity[4][GROWTH_HEADER + MEAN_HEADER],
        'High ppGpp': limitations_high[74][GROWTH_HEADER + MEAN_HEADER],
        # 'Decrease enzymes': limitations_high[94][GROWTH_HEADER + MEAN_HEADER],
        'Increase enzymes': limitations_high[99][GROWTH_HEADER + MEAN_HEADER],
        # 'Decrease ribosomes': limitations_high[103][GROWTH_HEADER + MEAN_HEADER],  # TODO: run with low rRNA?
        'Increase ribosomes': ribosome_limit_inhibition[45][GROWTH_HEADER + MEAN_HEADER],
        'Increase ribosomes,\nno GTPase inhibition': ribosome_limit_no_inhibition[45][GROWTH_HEADER + MEAN_HEADER],
        }

    # Plot low ppGpp sim results
    plot_data(low_data)
    save_fig('low-ppgpp')

    # Plot high ppGpp sim results
    plot_data(high_data)
    save_fig('high-ppgpp')
