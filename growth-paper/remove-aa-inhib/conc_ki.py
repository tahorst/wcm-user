#!/usr/bin/env python
"""
Compare Sander et al data to output from the model at different effective KI
values.

Requires data from remove_aa_inhibition variant and saving numpy arrays from the
variant analysis plot:
    control_conc: remove-inhib-conc-control.npy
    conc: remove-inhib-conc.npy
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
VALIDATION_DATA = os.path.join(FILE_LOCATION, 'sander-allosteric-aa-conc.tsv')
MODEL_DATA = os.path.join(FILE_LOCATION, 'remove-inhib-conc')
OUTPUT_FILE = 'aa-ki.pdf'

# From variant analysis plot
AMINO_ACIDS = [
	'GLN',
	'GLT',
	'ARG',
	'PRO',
	'L-ASPARTATE',
	'ASN',
	'LYS',
	'MET',
	'THR',
	'VAL',
	'L-ALPHA-ALANINE',
	'SER',
	'GLY',
	'HIS',
	'PHE',
	'TRP',
	'TYR',
	'ILE',
	'LEU',
    ]
ENZYMES = [
	'argA',
	'trpE',
	'hisG',
	'leuA',
	'thrA',
	'ilvA',
	'proB',
    ]
AA_IDX = {aa: i for i, aa in enumerate(sorted(AMINO_ACIDS))}
ENZ_IDX = {enz: i for i, enz in enumerate(ENZYMES)}
KI_FACTORS = sorted([np.inf, 2, 5, 10, 100])[::-1]  # from sim variant file

# Data to plot
COMPARISONS = [
    ('ARG', 'argA'),
    ('TRP', 'trpE'),
    ('HIS', 'hisG'),
    ('LEU', 'leuA'),
    ('THR', 'thrA'),
    ('PRO', 'proB'),
    ('ILE', 'ilvA'),
    ]


def load_validation():
    data = {}
    with open(VALIDATION_DATA) as f:
        reader = csv.reader(f, delimiter='\t')

        for _ in range(3):
            next(reader)
        conditions = next(reader)[1:]
        for row in reader:
            aa = row[0]
            conc = np.array(row[1:], float)
            data[aa] = dict(zip(conditions, conc))

    return data

def load_model():
    control = np.load(MODEL_DATA + '-control.npy')
    conc = np.load(MODEL_DATA + '.npy')
    return control, conc

def get_validation(validation, aa, enz):
    return validation[aa]['WT'], np.array([validation[aa][enz]])

def get_model(model, aa, enz):
    aa_idx = AA_IDX[aa]
    control = model[0][aa_idx]
    conc = model[1][AA_IDX[aa], ENZ_IDX[enz], :]
    return control, conc

def plot_ki_range(validation, model):
    _, axes = plt.subplots(len(COMPARISONS), 2, figsize=(8, 20))
    for i, (aa, enz) in enumerate(COMPARISONS):
        val = get_validation(validation, aa, enz)
        wcm = get_model(model, aa, enz)

        ax = axes[i, 0]
        data = [val[1][0] / val[0]] + list(wcm[1] / wcm[0])
        x = np.arange(len(data))
        ax.bar(x, data)
        ax.set_yscale('log')
        ax.set_ylabel(f'Normalized {aa} conc to WT')
        ax.set_xticks(x)
        ax.set_xticklabels(['Val'] + KI_FACTORS, rotation=45, fontsize=6)

        ax = axes[i, 1]
        data = [val[0]] + list(val[1]) + [wcm[0]] + list(wcm[1])
        x = np.arange(len(data))
        ax.bar(x, data)
        ax.set_yscale('log')
        ax.set_ylabel(f'{aa} conc')
        ax.set_xticks(x)
        ax.set_xticklabels(['Val WT', 'Val mutant', 'WCM WT'] + KI_FACTORS, rotation=45, fontsize=6)

def plot_bars(data, fun, log=False):
    # TODO: label x and y axes
    # TODO: highlight enzyme
    cols = 3
    rows = int(np.ceil(len(COMPARISONS) / cols))
    _, axes = plt.subplots(rows, cols)

    for i, (aa, _) in enumerate(COMPARISONS):
        row = i // cols
        col = i % cols
        ax = axes[row, col]

        conc = [fun(data, aa, ENZYMES[0])[0]] + [fun(data, aa, enz)[1][0] for enz in ENZYMES]
        ax.bar(range(len(conc)), conc)

        if log:
            ax.set_yscale('log')

def plot_scatter(validation, model):
    val_control = []
    model_control = []
    val_mutants = []
    model_mutants = []
    val_other = []
    model_other = []

    for aa, enz in COMPARISONS:
        val = get_validation(validation, aa, enz)
        val_control.append(val[0])
        val_mutants.append(val[1][0])

        mod = get_model(model, aa, enz)
        model_control.append(mod[0])
        model_mutants.append(mod[1][0])

    for aa in AMINO_ACIDS:
        for enz in ENZYMES:
            if (aa, enz) in COMPARISONS:
                continue
            val_other.append(get_validation(validation, aa, enz)[1][0])
            model_other.append(get_model(model, aa, enz)[1][0])

    plt.figure()
    plt.loglog(val_control, model_control, 'or', alpha=0.5, label='Allosteric AA in WT')
    plt.loglog(val_mutants, model_mutants, 'ob', alpha=0.5, label='Allosteric AA in mutant')
    plt.loglog(val_other, model_other, 'ok', alpha=0.2, markersize=2, label='Other AA')
    plt.legend(fontsize=6, frameon=False)

    xlim = plt.xlim()
    ylim = plt.ylim()
    min_ax = min(xlim[0], ylim[0])
    max_ax = max(xlim[1], ylim[1])
    xy_line = [min_ax, max_ax]
    plt.loglog(xy_line, xy_line, '--k', alpha=0.2, linewidth=1)

    plt.xlabel('Amino acid conc\nin validation (mM)')
    plt.ylabel('Amino acid conc\nin model (mM)')

def save_fig(filename):
    plt.tight_layout()
    plt.savefig(filename)
    plt.close('all')
    print(f'Saved to {filename}')


if __name__ == '__main__':
    validation = load_validation()
    model = load_model()

    # Compare validation data to range of KI values produced in the model
    plot_ki_range(validation, model)
    save_fig(OUTPUT_FILE)

    # Validation bar plot to reproduce Fig 1B from paper
    plot_bars(validation, get_validation)
    save_fig('validation-bar.pdf')

    # Comparable bar plot from model output
    plot_bars(model, get_model, log=True)
    save_fig('model-bar.pdf')

    # Scatter plot between validation and model
    plot_scatter(validation, model)
    save_fig('scatter.pdf')
