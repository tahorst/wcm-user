#! /usr/bin/env python

"""
Load and explore regulation data.

TODO:
	tRNA control
	TF control
	AA regulation in pathways
"""

import csv
import os
import re
from typing import Dict, Tuple


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

AMINO_ACIDS = [
	'L-alanine',
	'L-arginine',
	'L-asparagine',
	'L-aspartate',
	'L-cysteine',
	'L-glutamine',
	'L-glutamate',
	'glycine',
	'L-histidine',
	'L-isoleucine',
	'L-leucine',
	'L-lysine',
	'L-methionine',
	'L-phenylalanine',
	'L-proline',
	'L-serine',
	'L-threonine',
	'L-tryptophan',
	'L-tyrosine',
	'L-valine',
	]


def load_regulation() -> Dict[str, Dict[str, Tuple[str, str]]]:
	"""
	Load regulation data.

	Returns:
		regulation: dictionary of regulatory pair mappings
			{regulator: {regulatee: (type, direction)}}
			type options:
				'Allosteric-Regulation-of-RNAP'
				'Compound-Mediated-Translation-Regulation'
				'Protein-Mediated-Attenuation'
				'Protein-Mediated-Translation-Regulation'
				'RNA-Mediated-Translation-Regulation'
				'Regulation'
				'Regulation-of-Enzyme-Activity'
				'Regulation-of-Translation'
				'Rho-Blocking-Antitermination'
				'Ribosome-Mediated-Attenuation'
				'Sigma-Factor'
				'Small-Molecule-Mediated-Attenuation'
				'Transcription-Factor-Binding'
				'Transcriptional-Attenuation'
			direction options:
				'+'
				'-'
	"""

	filename = os.path.join(FILE_LOCATION, 'regevents.tsv')
	with open(filename) as f:
		reader = csv.reader(f, delimiter='\t')

		regulation = {}
		for line in reader:
			if line[0].startswith('#'):
				continue

			reg_type = line[2]
			direction = line[3]
			regulator, regulatee = re.findall('(.*) [-+]*> (.*)', line[4])[0]

			reg = regulation.get(regulator, {})
			reg[regulatee] = (reg_type, direction)
			regulation[regulator] = reg

	return regulation

def summarize_aa_control(regulation: Dict[str, Dict[str, Tuple[str, str]]]) -> None:
	"""
	Show genes that have enzymatic activity regulated by amino acids.

	Args:
		regulation: dictionary of regulatory pair mappings, see load_regulation()
	"""

	print('Amino acid regulation of enzyme activity:')
	for aa in AMINO_ACIDS:
		regulated = [
			regulatee
			for regulatee, control in sorted(regulation[aa].items())
			if control[0] == 'Regulation-of-Enzyme-Activity' and control[1] == '-'
			]
		print(f'\t{aa}: {", ".join(regulated)}')


if __name__ == '__main__':
	regulation = load_regulation()
	summarize_aa_control(regulation)
