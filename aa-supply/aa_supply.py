#! /usr/bin/env python

"""
Script to test changes in supply for each amino acid when adding a single
amino acid to media and assuming an internal concentration increase of 2x.

Must supply out directory that will contain sim_data as first arg (eg. out/test)
"""

import os
import pickle
import sys

import numpy as np

from wholecell.utils import units


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))


# Check arg input
if len(sys.argv) != 2:
    print('Must provide out directory for sim_data load (eg. out/test)')
    sys.exit(1)

# Load sim_data
with open(os.path.join(sys.argv[1], 'kb', 'simData.cPickle'), 'rb') as f:
    sim_data = pickle.load(f)
metabolism = sim_data.process.metabolism
aa_ids = sim_data.molecule_groups.amino_acids

# Example data from a sim
enzymes = np.load(os.path.join(FILE_LOCATION, 'enzymes.npy'))
aa_conc = units.mol / units.L * np.load(os.path.join(FILE_LOCATION, 'aas.npy'))
dry_mass = units.fg * 351.2
aa_in_media = np.zeros(21, dtype=bool)
aa_in_media[-2] = True
time_step = 2

# Calculate original values
synthesis, enzyme_counts, saturation = metabolism.amino_acid_synthesis(enzymes, aa_conc)
aa_supply = time_step * (synthesis + metabolism.amino_acid_import(aa_in_media, dry_mass))

# Check change in supply when adding each amino acid to the media (assuming only its internal
# concentration has increased)
for aa_idx in range(21):
    if aa_idx == 19:
        continue

    aa_conc[aa_idx] *= 2
    aa_in_media[aa_idx] = True

    new_synthesis, new_enzymes, new_saturation = metabolism.amino_acid_synthesis(enzymes, aa_conc)
    new_aa_supply = time_step * (new_synthesis + metabolism.amino_acid_import(aa_in_media, dry_mass))

    aa_conc[aa_idx] /=2
    aa_in_media[aa_idx] = False

    print(f'Added {aa_ids[aa_idx]}:')
    fraction = (new_synthesis - synthesis) / synthesis
    for aa, frac in zip(aa_ids, fraction):
        if frac == 0:
            continue

        print(f'\t{aa}: {frac:.3}')
