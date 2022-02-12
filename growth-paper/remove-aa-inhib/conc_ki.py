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
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
VALIDATION_DATA = os.path.join(FILE_LOCATION, 'sander-allosteric-aa-conc.tsv')
MODEL_DATA = os.path.join(FILE_LOCATION, 'remove-inhib-conc')
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)

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
KIS = np.array([0.15, 0.17, 0.07327, 0.28, 0.167, 0.15, 0.06])  # KIs (mM) from sim_data.process.metabolism.aa_kis corresponding to AA in COMPARISONS

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

def plot_ki_prediction(validation, model):
    # TODO: decide on metric (need to do 1 - metric?)
    # TODO: fit curve to points and drop prediction vline
    cols = 3
    rows = int(np.ceil(len(COMPARISONS) / cols))
    _, axes = plt.subplots(rows, cols)
    hide_axes = np.ones_like(axes, dtype=bool)

    ki_factors = np.array(KI_FACTORS)
    for i, ((aa, enz), ki) in enumerate(zip(COMPARISONS, KIS)):
        row = i // cols
        col = i % cols
        ax = axes[row, col]
        hide_axes[row, col] = False

        val = get_validation(validation, aa, enz)
        wcm = get_model(model, aa, enz)

        val_increase = val[1][0] / val[0]
        wcm_increase = wcm[1] / wcm[0]

        kis = ki * ki_factors
        fraction_wt = 1 - 1 / (1 + wcm[0] / ki)
        fraction_inhibited = 1 - 1 / (1 + wcm[0] / kis)
        reduced_inhibition = fraction_inhibited / fraction_wt

        # Original data
        ax.plot(reduced_inhibition, wcm_increase, 'o', alpha=0.5)
        ax.plot(1 / ki_factors, wcm_increase, 'o', alpha=0.5)

        # Curve fit data
        fun = lambda x, a, d: a*np.exp(d*x)
        sol = curve_fit(fun, 1/ki_factors, np.log(wcm_increase))  # Fit in log space for closer fit at higher x values
        ax.plot(1 / ki_factors, np.exp(fun(1 / ki_factors, *sol[0])), 'x')

        # Interp function
        x = np.log(wcm_increase)  # Fit in log space for smoother fit
        interp = interp1d(x, 1 / ki_factors)
        val_match = interp(np.log(val_increase))
        y = np.linspace(x.min(), x.max(), 1000)
        ax.plot(interp(y), np.exp(y))

        ax.axhline(val_increase, linestyle='--', color='k', linewidth=0.5, alpha=0.5)
        ax.axvline(val_match, linestyle='--', color='k', linewidth=0.5, alpha=0.5)

        ax.set_yscale('log')

        ylabel = f'{aa} conc (mM)'
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=6)

    # Hide unused axes
    for ax in axes[hide_axes]:
        ax.set_visible(False)

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

def plot_bars(datasets, functions, log=False, normalize=False):
    # TODO: label x and y axes
    # TODO: highlight enzyme
    # TODO: option to normalize based on control
    cols = 3
    rows = int(np.ceil(len(COMPARISONS) / cols))
    _, axes = plt.subplots(rows, cols)
    hide_axes = np.ones_like(axes, dtype=bool)

    for i, (aa, _) in enumerate(COMPARISONS):
        row = i // cols
        col = i % cols
        ax = axes[row, col]
        hide_axes[row, col] = False

        all_conc = []
        for data, fun in zip(datasets, functions):
            conc = np.array([fun(data, aa, ENZYMES[0])[0]] + [fun(data, aa, enz)[1][0] for enz in ENZYMES])
            if normalize:
                conc = np.log10(conc / conc[0])
            all_conc.append(conc)
        if len(all_conc) == 0:
            continue

        n_mutants = len(all_conc[0])
        n_datasets = len(all_conc)
        x = np.arange(n_mutants)
        width = 0.8 / n_datasets
        offsets = np.arange(n_datasets) * width - 0.4 + width/2

        for i, offset in enumerate(offsets):
            ax.bar(x + offset, all_conc[i], width)

        if log:
            ax.set_yscale('log')

        if normalize:
            ax.axhline(0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)
            ylabel = f'{aa} log10 increase\nover WT'
        else:
            ylabel = f'{aa} conc (mM)'

        ax.set_xticks(x)
        ax.set_xticklabels(['WT'] + ENZYMES, rotation=45, fontsize=6)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=6)

    # Hide unused axes
    for ax in axes[hide_axes]:
        ax.set_visible(False)

def plot_sub_scatter(validation, model):
    cols = 3
    rows = int(np.ceil(len(COMPARISONS) / cols))
    _, axes = plt.subplots(rows, cols, figsize=(10, 10))
    hide_axes = np.ones_like(axes, dtype=bool)

    for i, (aa, _) in enumerate(COMPARISONS):
        row = i // cols
        col = i % cols
        ax = axes[row, col]
        hide_axes[row, col] = False

        plot_scatter(ax, validation, model, label=aa, amino_acids=[aa], plot_all=True)

    # Hide unused axes
    for ax in axes[hide_axes]:
        ax.set_visible(False)

def plot_scatter(ax, validation, model, label='Amino acid', amino_acids=None,
        enzymes=None, legend=True, plot_all=False):
    if amino_acids is None:
        amino_acids = AMINO_ACIDS
    if enzymes is None:
        enzymes = ENZYMES

    val_control = []
    model_control = []
    val_mutants = []
    model_mutants = []
    val_other = []
    model_other = []

    for aa, enz in COMPARISONS:
        if aa not in amino_acids:
            continue

        val = get_validation(validation, aa, enz)
        val_control.append(val[0])
        val_mutants.append(val[1][0])

        mod = get_model(model, aa, enz)
        model_control.append(mod[0])
        model_mutants.append(mod[1][0])

        if plot_all:
            for m in mod[1][1:]:
                val_mutants.append(val[1][0])
                model_mutants.append(m)

    for aa in amino_acids:
        for enz in enzymes:
            if (aa, enz) in COMPARISONS:
                continue
            val_other.append(get_validation(validation, aa, enz)[1][0])
            model_other.append(get_model(model, aa, enz)[1][0])

    ax.loglog(val_control, model_control, 'or', alpha=0.5, label='Allosteric AA in WT')
    ax.loglog(val_mutants, model_mutants, 'ob', alpha=0.5, label='Allosteric AA in mutant')
    ax.loglog(val_other, model_other, 'ok', alpha=0.2, markersize=2, label='Other AA')

    if legend:
        ax.legend(fontsize=6, frameon=False)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    min_ax = min(xlim[0], ylim[0])
    max_ax = max(xlim[1], ylim[1])
    xy_line = [min_ax, max_ax]
    ax.loglog(xy_line, xy_line, '--k', alpha=0.2, linewidth=1)

    ax.set_xlabel(f'{label} conc\nin validation (mM)')
    ax.set_ylabel(f'{label} conc\nin model (mM)')

def save_fig(filename):
    file = os.path.join(OUTPUT_DIR, filename)
    plt.tight_layout()
    plt.savefig(file)
    plt.close('all')
    print(f'Saved to {file}')


if __name__ == '__main__':
    validation = load_validation()
    model = load_model()

    # Compare validation data AA increase to level of inhibition reduction in the model
    plot_ki_prediction(validation, model)
    save_fig('aa-ki-prediction.pdf')

    # Compare validation data to range of KI values produced in the model
    plot_ki_range(validation, model)
    save_fig('aa-ki-range.pdf')

    # Validation bar plot to reproduce Fig 1B from paper
    plot_bars([validation], [get_validation])
    save_fig('validation-bar.pdf')

    # Comparable bar plot from model output
    plot_bars([model], [get_model], log=True)
    save_fig('model-bar.pdf')

    # Side by side bar plot with validation and model data
    plot_bars([validation, model], [get_validation, get_model], log=True)
    save_fig('side-by-side-bar.pdf')

    # Side by side bar plot with validation and model data
    plot_bars([validation, model], [get_validation, get_model], normalize=True)
    save_fig('normalized-side-by-side-bar.pdf')

    # Scatter plot between validation and model
    plt.figure()
    plot_scatter(plt.gca(), validation, model)
    save_fig('scatter.pdf')

    # Scatter for each amino acid
    plot_sub_scatter(validation, model)
    save_fig('sub-scatter.pdf')
