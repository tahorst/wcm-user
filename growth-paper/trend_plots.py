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
CONDITION_SIMS = {
    'Amino_acid_combinations_in_media': '6+',
    'Add_one_amino_acid_shift': '6+',
    'Remove_one_amino_acid_shift': '4+',
    'Conditions_with_regulation': '',
    }
PERTURBATION_SIMS = {
    'ppGpp_sensitivity': '4+',
    'ppGpp_limitations_-_low_ppGpp': '',
    'Remove_amino_acid_inhibition': '',
    'Amino_acid_synthesis_network_sensitivity_-_glt': '',
    'ppGpp_limitations_-_high_ppGpp': '',
    'Amino_acid_synthesis_network_sensitivity_-_control': '',
    'ppGpp_limitations_-_normal_ppGpp': '',
    'ppGpp_sensitivity_-_no_mechanistic_transport': '4+',
    'Amino_acid_combinations_in_media_without_regulation_or_charging': '',
    'Conditions_without_regulation_or_charging': '',
    'Amino_acid_combinations_in_media_-_no_mechanistic_transport': '',
    'Remove_amino_acid_inhibition_-_no_ppgpp': '',
    }
FILE_PATH = 'plotOut/{}{}.tsv'
OUTPUT_FILE = 'growth-trends.pdf'


def load_datasets(sims):
    all_data = {}
    all_headers = set()
    min_headers = None
    for sim_desc, prepend in sims.items():
        data, headers = load_data(sim_desc, prepend=prepend)
        all_data[sim_desc] = data

        header_set = set(headers)
        all_headers |= header_set
        if min_headers is None:
            min_headers = header_set
        else:
            min_headers &= header_set

    # Check to make sure headers are consistent across all sims
    if min_headers != all_headers:
        raise ValueError('Headers do not match for certain sims')

    return all_data, sorted(all_headers)

def load_data(desc, prepend='', filename='rp_data'):
    dirs = os.listdir(SIM_DIR)
    for d in dirs:
        if d.endswith(desc):
            path = os.path.join(SIM_DIR, d, FILE_PATH.format(prepend, filename))
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

    return data, headers

def get_comparisons(y_keys):
    # TODO: take list of x and y and plot on subplots
    n_keys = len(y_keys)
    rows = int(np.ceil(np.sqrt(n_keys)))
    cols = int(np.ceil(n_keys / rows))

    return rows, cols

def plot_setup(y_keys, scale=3):
    # TODO: take list of x and y to pass to get_comparisons
    # TODO: add subplot for legend
    rows, cols = get_comparisons(y_keys)
    _, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(cols*scale, rows*scale))

    return axes

def plot(axes, all_data, y_keys, fade=False):
    if fade:
        options = dict(color='k', alpha=0.1, markersize=2)
    else:
        options = dict(alpha=0.4)

    # TODO: multiple comparisons of x and y
    x_key = 'rp_ratio'
    rows, cols  = get_comparisons(y_keys)

    for i, y_key in enumerate(y_keys):
        row = i % rows
        col = i // rows
        ax = axes[row, col]

        for sim, data in all_data.items():
            ax.plot(data[x_key], data[y_key], 'o', label=sim, **options)

        format_ax(ax, x_key, y_key)

def format_ax(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def save_fig(output_file):
    plt.tight_layout()
    plt.savefig(output_file)
    print(f'Saved to {output_file}')
    plt.close('all')

def get_ax_lim(axes):
    return [(ax.get_xlim(), ax.get_ylim()) for ax in axes.flatten()]

def set_ax_lim(axes, lims):
    for ax, lim in zip(axes.flatten(), lims):
        ax.set_xlim(lim[0])
        ax.set_ylim(lim[1])


if __name__ == '__main__':
    condition_data, headers = load_datasets(CONDITION_SIMS)
    perturbation_data, _ = load_datasets(PERTURBATION_SIMS)

    # Condition sims
    axes = plot_setup(headers)
    plot(axes, condition_data, headers)
    save_fig('conditions-' + OUTPUT_FILE)

    # Condition sims with faded perturbations
    axes = plot_setup(headers)
    plot(axes, condition_data, headers)
    plot(axes, perturbation_data, headers, fade=True)
    save_fig('all-' + OUTPUT_FILE)

    # Condition sims with faded perturbations and fixed axes
    axes = plot_setup(headers)
    plot(axes, condition_data, headers)
    ax_lim = get_ax_lim(axes)
    plot(axes, perturbation_data, headers, fade=True)
    set_ax_lim(axes, ax_lim)
    save_fig('all-fixed-' + OUTPUT_FILE)
