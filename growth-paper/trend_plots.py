#!/usr/bin/env python
"""
Explore trends between growth related data to make a similar comparison like
growth vs R/P ratio.

Requires the directory of all the sims run for the paper as SIM_DIR with
certain descriptions assumed in SIM_DESC to select sets of runs.
"""

import csv
import os
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np


# Input paths
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

# Output paths
FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_FILE = 'growth-trends.pdf'

ADDED_DATA = {
    'aa / ribosome': lambda data: data['aa_mass'] / data['ribosome_mass'],
    'glt / aa': lambda data: data['glt'] / data['aa_mass'],
    'glt fraction / ppGpp': lambda data: data['glt'] / data['aa_mass'] / data['ppGpp'],
    'ribosome / cell mass fraction': lambda data: data['ribosome_mass'] / data['cell_mass'],
    'ribosome / protein mass fraction': lambda data: data['ribosome_mass'] / data['protein_mass'],
    'excess aa': lambda data: data['aa_mass'] / data['protein_mass'],
    'total excess': lambda data: data['aa_mass'] / data['protein_mass'] + data['excess rrna fraction'] + data['excess rprotein fraction'],
    }
ALL_X_COMPARISONS = ['growth_rate']
ALL_Y_COMPARISONS = ['rp_ratio', 'ppGpp']
OTHER_COMPARISONS = [
    ('glt / aa', 'aa / ribosome'),
    ('glt fraction / ppGpp', 'aa / ribosome'),
    ]


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

        for key, fun in ADDED_DATA.items():
            data[key] = fun(data)
            headers.append(key)

    return data, headers

def get_comparisons(
        keys: List[str],
        all_x: Optional[List[str]] = None,
        all_y: Optional[List[str]] = None,
        others: Optional[List[Tuple[str, str]]] = None,
        ) -> Tuple[Tuple[int, int], List[Tuple[str, str]]]:
    # Default to empty lists to skip iterations below
    if all_x is None:
        all_x = []
    if all_y is None:
        all_y = []
    if others is None:
        others = []

    x_keys = []
    y_keys = []

    # Add key pairs for x and y axes
    for y in all_x:
        x_keys += keys
        y_keys += [y] * len(keys)
    for x in all_y:
        x_keys += [x] * len(keys)
        y_keys += keys
    for x, y in others:
        x_keys.append(x)
        y_keys.append(y)
    key_pairs = list(zip(x_keys, y_keys))

    # Get near square dimensions for plotting
    n_keys = len(y_keys)
    rows = int(np.ceil(np.sqrt(n_keys)))
    cols = int(np.ceil(n_keys / rows))
    dims = (rows, cols)

    return dims, key_pairs

def plot_setup(rows, cols, scale=3):
    # TODO: add subplot for legend
    _, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(cols*scale, rows*scale))

    return axes

def plot(axes, all_data, dims, keys, fade=False, sims=None, mask_funs=None, **plot_options):
    if fade:
        options = dict(color='k', alpha=0.1, markersize=2)
    else:
        options = dict(alpha=0.4)
    options.update(plot_options)

    if sims is None:
        sims = all_data.keys()
    if mask_funs is None:
        mask_funs = [lambda data: slice(None)]

    rows, cols = dims
    for i, (x_key, y_key) in enumerate(keys):
        row = i % rows
        col = i // rows
        ax = axes[row, col]

        for sim in sims:
            data = all_data[sim]
            for mask_fun in mask_funs:
                mask = mask_fun(data)
                ax.plot(data[x_key][mask], data[y_key][mask], 'o', label=sim, **options)

        format_ax(ax, x_key, y_key)

def format_ax(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def save_fig(output_file):
    filename = os.path.join(OUTPUT_DIR, output_file)
    plt.tight_layout()
    plt.savefig(filename)
    print(f'Saved to {filename}')
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
    dims, keys = get_comparisons(headers, all_x=ALL_X_COMPARISONS,
        all_y=ALL_Y_COMPARISONS, others=OTHER_COMPARISONS)

    # Condition sims
    axes = plot_setup(*dims)
    plot(axes, condition_data, dims, keys)
    save_fig('conditions-' + OUTPUT_FILE)

    # Condition sims with faded perturbations
    axes = plot_setup(*dims)
    plot(axes, condition_data, dims, keys)
    plot(axes, perturbation_data, dims, keys, fade=True)
    save_fig('all-' + OUTPUT_FILE)

    # Condition sims with faded perturbations and fixed axes
    axes = plot_setup(*dims)
    plot(axes, condition_data, dims, keys)
    ax_lim = get_ax_lim(axes)
    plot(axes, perturbation_data, dims, keys, fade=True)
    set_ax_lim(axes, ax_lim)
    save_fig('all-fixed-' + OUTPUT_FILE)

    # Condition sims with ppGpp sensitivity high and low sims to see if any attributes are off for one but not the other
    sims = ['ppGpp_sensitivity']
    low_mask = [lambda data: (data['variant'] >= 0) & (data['variant'] < 4)]
    high_mask = [lambda data: (data['variant'] > 4) & (data['variant'] < 10)]
    axes = plot_setup(*dims)
    plot(axes, condition_data, dims, keys)
    ax_lim = get_ax_lim(axes)
    plot(axes, perturbation_data, dims, keys, sims=sims, mask_funs=low_mask, marker='X', color='k')
    plot(axes, perturbation_data, dims, keys, sims=sims, mask_funs=high_mask, marker='X')
    set_ax_lim(axes, ax_lim)
    save_fig('ppgpp-sensitivity-' + OUTPUT_FILE)

    # Low ppGpp limitation sims with enzymes and ribosomes plotted
    sims = ['ppGpp_limitations_-_low_ppGpp']
    masks = [
        lambda data: (data['variant'] > 18) & (data['variant'] < 28),  # enzymes
        lambda data: (data['variant'] > 27) & (data['variant'] < 37),  # ribosomes
        ]
    axes = plot_setup(*dims)
    plot(axes, condition_data, dims, keys)
    ax_lim = get_ax_lim(axes)
    plot(axes, perturbation_data, dims, keys, sims=sims, mask_funs=masks, marker='x')
    set_ax_lim(axes, ax_lim)
    save_fig('low-ppgpp-' + OUTPUT_FILE)

    # High ppGpp limitation sims with enzymes and ribosomes plotted
    sims = ['ppGpp_limitations_-_high_ppGpp']
    masks = [
        lambda data: (data['variant'] > 92) & (data['variant'] < 102),  # enzymes
        lambda data: (data['variant'] > 101) & (data['variant'] < 111),  # ribosomes
        ]
    axes = plot_setup(*dims)
    plot(axes, condition_data, dims, keys)
    ax_lim = get_ax_lim(axes)
    plot(axes, perturbation_data, dims, keys, sims=sims, mask_funs=masks, marker='x')
    set_ax_lim(axes, ax_lim)
    save_fig('high-ppgpp-' + OUTPUT_FILE)
