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
    ('ILE', 'ilvA'),
    ('PRO', 'proB'),
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
    return validation[aa]['WT'], validation[aa][enz]

def get_model(model, aa, enz):
    aa_idx = AA_IDX[aa]
    control = model[0][aa_idx]
    conc = model[1][AA_IDX[aa], ENZ_IDX[enz], :]
    return control, conc


if __name__ == '__main__':
    validation = load_validation()
    model = load_model()

    _, axes = plt.subplots(len(COMPARISONS), 2, figsize=(8, 20))
    for i, (aa, enz) in enumerate(COMPARISONS):
        val = get_validation(validation, aa, enz)
        wcm = get_model(model, aa, enz)

        ax = axes[i, 0]
        data = [val[1] / val[0]] + list(wcm[1] / wcm[0])
        x = np.arange(len(data))
        ax.bar(x, data)
        ax.set_yscale('log')
        ax.set_ylabel(f'Normalized {aa} conc to WT')
        ax.set_xticks(x)
        ax.set_xticklabels(['Val'] + KI_FACTORS, rotation=45, fontsize=6)

        ax = axes[i, 1]
        data = list(val) + [wcm[0]] + list(wcm[1])
        x = np.arange(len(data))
        ax.bar(x, data)
        ax.set_yscale('log')
        ax.set_ylabel(f'{aa} conc')
        ax.set_xticks(x)
        ax.set_xticklabels(['Val WT', 'Val mutant', 'WCM WT'] + KI_FACTORS, rotation=45, fontsize=6)

    plt.tight_layout()
    plt.savefig(OUTPUT_FILE)
    print(f'Saved to {OUTPUT_FILE}')
