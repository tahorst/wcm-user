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


SIM_DIR = '/home/travis/scratch/wcEcoli_out/'
BASE_SIM_DIR = 'Conditions_without_regulation_or_charging'
NEW_C_SOURCE_DIR = 'Conditions_with_regulation'
ADD_ONE_DIR = 'Add_one_amino_acid_shift'
REMOVE_ONE_DIR = 'Remove_one_amino_acid_shift'
PPGPP_DIR = 'ppGpp_sensitivity_-_no_mechanistic_transport'
PPGPP_LIMITATION_LOW_DIR = 'ppGpp_limitations_-_low_ppGpp'
PPGPP_LIMITATION_HIGH_DIR = 'ppGpp_limitations_-_high_ppGpp'
NEW_AA_SOURCE_DIR = 'Amino_acid_combinations_in_media'
INHIBITION_NO_PPGPP_DIR = 'Remove_amino_acid_inhibition_-_no_ppgpp'
INHIBITION_DIR = 'Remove_amino_acid_inhibition'
FILE_PATH = 'plotOut/{}.tsv'
OUTPUT_FILE = 'combined-growth-rp.pdf'

CONTROL_IDX = 19
GROWTH_HEADER = 'Growth'
RP_HEADER = 'R/P ratio'
MEAN_HEADER = ' mean'
STD_HEADER = ' std'

DENNIS_BREMER_2021 = np.array([
[0.1691, 0.415888308335967],
[0.2056, 0.693147180559945],
[0.2576, 1.03972077083992],
[0.3307, 1.38629436111989],
[0.4176, 1.73286795139986],
[0.5023, 2.07944154167984],
])

ONE_AA_OPTIONS = dict(alpha=0.5, markersize=4)
PPGPP_OPTIONS = dict(alpha=0.5, markersize=6)
FADE_OPTIONS = dict(alpha=0.2, markersize=4, color='black')


def load_data(desc, filename='growth_trajectory'):
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

def plot_conditions(fade=False, options=None, std=True, label=True):
    if fade:
        options = FADE_OPTIONS
        std = False
        label = False

    plot(no_regulation, std=std, label='Original conditions' if label else '', variants=np.arange(3), options=options)
    plot(regulation, std=std, label='New carbon sources with growth regulation' if label else '', variants=np.arange(3, 5), options=options)
    plot(add_one, std=std, variants=[CONTROL_IDX], label='Minimal + glc with growth regulation' if label else '', options=options)
    plot(remove_one, std=std, variants=[CONTROL_IDX], label='Rich + glc with growth regulation' if label else '', options=options)
    plot(new_aa, std=std, variants=[1, 2], label='New amino acid media with growth regulation' if label else '', options=options)
    plot(add_one, std=False, exclude=[CONTROL_IDX], label='Add one AA to minimal with growth regulation' if label else '',
        options=options if options else ONE_AA_OPTIONS)
    plot(remove_one, std=False, exclude=[CONTROL_IDX], label='Remove one AA from rich with growth regulation' if label else '',
        options=options if options else ONE_AA_OPTIONS)

def plot_ppgpp():
    plot(ppgpp, std=False, label='Minimal lower ppGpp', variants=range(2, 4), options=PPGPP_OPTIONS)  # 0, 4 for all
    plot(ppgpp, std=False, label='Minimal higher ppGpp', variants=range(5, 8), options=PPGPP_OPTIONS)  # 5, 10 for all
    plot(add_one, variants=[CONTROL_IDX], label='Minimal + glc')
    plot(remove_one, variants=[CONTROL_IDX], label='Rich + glc')
    plot(ppgpp, std=False, label='Rich higher ppGpp', variants=range(12, 15), options=PPGPP_OPTIONS)  # 12, 20 for all
    plot(ppgpp_low, std=False, label='Low inhibition enzymes', variants=range(24, 28), options=PPGPP_OPTIONS)  # 19 starts low, 19 failed
    plot(ppgpp_low, std=False, label='Low inhibition ribosomes', variants=range(33, 37), options=PPGPP_OPTIONS)  # 28 starts low
    plot(ppgpp_high, std=False, label='High inhibition enzymes', variants=range(98, 102), options=PPGPP_OPTIONS)  # 93 starts low
    plot(ppgpp_high, std=False, label='High inhibition ribosomes', variants=range(107, 111), options=PPGPP_OPTIONS)  # 102 starts low

def plot_inhibition():
    plot(inhib_no_ppgpp, variants=range(1, 8), std=False, label='Removed allosteric inhibition without ppGpp', options=PPGPP_OPTIONS)
    plot(inhib, variants=range(1, 8), std=False, label='Removed allosteric inhibition with ppGpp', options=PPGPP_OPTIONS)
    plot(add_one, variants=[CONTROL_IDX], label='Minimal + glc')
    plot(remove_one, variants=[CONTROL_IDX], label='Rich + glc')

    # Control variants
    # plot(inhib_no_ppgpp, variants=[0], std=False, label='Removed allosteric inhibition without ppGpp', options=PPGPP_OPTIONS)
    # plot(inhib, variants=[0], std=False, label='Removed allosteric inhibition with ppGpp', options=PPGPP_OPTIONS)

def plot_trends(all=False):
    plt.plot([0.11, 0.52], [0, 2], '--k')  # Dennis and Bremer (dry_mass_composition.tsv)

    if all:
        plt.plot([0.07, 0.49], [0, 2], '--k')  # Zhu et al. Growth suppression by altered (p)ppGpp levels... 2019.
        plt.plot(DENNIS_BREMER_2021[:, 0], DENNIS_BREMER_2021[:, 1], 'x')

def format_plot():
    plt.legend(fontsize=8, frameon=False)
    plt.xlabel('RNA/protein mass ratio')
    plt.ylabel('Growth rate (1/hr)')
    plt.xlim([0, 0.6])
    plt.ylim([0, 2])
    plt.tight_layout()

def save_fig(output_file):
    plt.savefig(output_file)
    print(f'Saved to {output_file}')
    plt.close('all')


if __name__ == '__main__':
    no_regulation = load_data(BASE_SIM_DIR)
    regulation = load_data(NEW_C_SOURCE_DIR)
    add_one = load_data(ADD_ONE_DIR)
    remove_one = load_data(REMOVE_ONE_DIR)
    ppgpp = load_data(PPGPP_DIR)
    ppgpp_low = load_data(PPGPP_LIMITATION_LOW_DIR)
    ppgpp_high = load_data(PPGPP_LIMITATION_HIGH_DIR)
    new_aa = load_data(NEW_AA_SOURCE_DIR)
    inhib_no_ppgpp = load_data(INHIBITION_NO_PPGPP_DIR)
    inhib = load_data(INHIBITION_DIR)

    plt.figure()
    plot_conditions()
    plot_trends()
    format_plot()
    save_fig(OUTPUT_FILE)

    plt.figure()
    plot_ppgpp()
    plot_trends()
    format_plot()
    save_fig('ppgpp-' + OUTPUT_FILE)

    plt.figure()
    plot_inhibition()
    plot_trends()
    format_plot()
    save_fig('inhibition-' + OUTPUT_FILE)

    # Optional plots
    # ppgpp_aa = load_data(PPGPP_DIR, 'protein_aa-growth_trajectory')
    # plot(ppgpp, std=False, label='Rich lower ppGpp', variants=range(10, 11), options=PPGPP_OPTIONS)  # off window
    # plot(ppgpp_aa, std=False, label='Minimal lower ppGpp', variants=range(2, 4), options=PPGPP_OPTIONS)
    # plot(ppgpp_aa, std=False, label='Minimal higher ppGpp', variants=range(5, 8), options=PPGPP_OPTIONS)
    # plot(ppgpp_aa, std=False, label='Rich lower ppGpp', variants=range(10, 11), options=PPGPP_OPTIONS)
    # plot(ppgpp_aa, std=False, label='Rich higher ppGpp', variants=range(12, 17), options=PPGPP_OPTIONS)

