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

AMINO_ACIDS = {
	'L-alanine': 'Ala',
	'L-arginine': 'Arg',
	'L-asparagine': 'Asn',
	'L-aspartate': 'Asp',
	'L-cysteine': 'Cys',
	'L-glutamine': 'Gln',
	'L-glutamate': 'Glt',
	'glycine': 'Gly',
	'L-histidine': 'His',
	'L-isoleucine': 'Ile',
	'L-leucine': 'Leu',
	'L-lysine': 'Lys',
	'L-methionine': 'Met',
	'L-phenylalanine': 'Phe',
	'L-proline': 'Pro',
	'L-serine': 'Ser',
	'L-threonine': 'Thr',
	'L-tryptophan': 'Trp',
	'L-tyrosine': 'Tyr',
	'L-valine': 'Val',
	}
PATHWAYS = {
	aa: [
		f'{gene[:3].lower()}{letter}'
		for gene in regulated
		for letter in (gene[3:] if gene[3:] else [''])
		]
	for aa, regulated in {
		'L-alanine': ['IlvE', 'AvtA', 'AlaCA', 'IscS', 'SufS'],
		'L-arginine': ['ArgABCDEFGHI', 'AstC', 'GabT', 'CarAB'],
		'L-asparagine': ['AsnAB'],
		'L-aspartate': ['AspC'],
		'L-cysteine': ['CysEMK', 'GrxA'],
		'L-glutamine': ['GlnA'],
		'L-glutamate': ['GltBD', 'GdhA'],
		'glycine': ['GlyA'],
		'L-histidine': ['HisGIAHFBCD'],
		'L-isoleucine': ['IlvAGMIHCDE'],
		'L-leucine': ['LeuACDB, DmlA, TyrB, IlvE'],
		'L-lysine': ['LysCA', 'Asd', 'DapABDEF', 'SerC', 'ArgD'],
		'L-methionine': ['MetABCEH', 'MalY'],
		'L-phenylalanine': ['PheA', 'TyrAB', 'IlvE', 'AspC'],
		'L-proline': ['ProBAC'],
		'L-serine': ['SerACB'],
		'L-threonine': ['ThrBC'],
		'L-tryptophan': ['TrpEDCAB'],
		'L-tyrosine': ['TyrAB', 'PheA', 'AspC'],
		'L-valine': ['IlvIHGMBNCDE'],
		}.items()
	}


def load_regulation() -> Tuple[Dict[str, Dict[str, Tuple[str, str]]], Dict[str, Dict[str, Tuple[str, str]]]]:
	"""
	Load regulation data.

	Returns:
		regulators: dictionary of regulatory pair mappings
			{regulator: {regulatee: (type, direction)}}
		regulatees: dictionary of regulatory pair mappings
			{regulatee: {regulator: (type, direction)}}

	Notes:
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

		regulators = {}
		regulatees = {}
		for line in reader:
			if line[0].startswith('#'):
				continue

			reg_type = line[2]
			direction = line[3]
			regulator, regulatee = re.findall('(.*) [-+]*> (.*)', line[4])[0]

			# Map regulators to regulatees
			reg = regulators.get(regulator, {})
			reg[regulatee] = (reg_type, direction)
			regulators[regulator] = reg

			# Map regulatees to regulators
			reg = regulatees.get(regulatee, {})
			reg[regulator] = (reg_type, direction)
			regulatees[regulatee] = reg

	return regulators, regulatees

def summarize_aa_control(regulation: Dict[str, Dict[str, Tuple[str, str]]]) -> None:
	"""
	Show genes that have enzymatic activity regulated by amino acids.

	Args:
		regulation: dictionary of regulator to regulatee pair mappings,
			see load_regulation()
	"""

	print('Amino acid regulation of enzyme activity:')
	for aa in AMINO_ACIDS:
		regulated = [
			regulatee
			for regulatee, control in sorted(regulation[aa].items())
			if control[0] == 'Regulation-of-Enzyme-Activity' and control[1] == '-'
			]
		print(f'\t{aa}: {", ".join(regulated)}')

def summarize_pathway_control(regulation: Dict[str, Dict[str, Tuple[str, str]]]) -> None:
	"""
	Show amino acids that regulate enzymatic activity in synthesis pathways for
	other amino acids.

	Args:
		regulation: dictionary of regulatee to regulator pair mappings,
			see load_regulation()
	"""

	print('Pathway activity regulated by other amino acids:')
	for aa, enzymes in PATHWAYS.items():
		regulated = {
			AMINO_ACIDS[regulator]
			for enzyme in enzymes
			for regulator, control in sorted(regulation.get(enzyme, {}).items())
			if control[0] == 'Regulation-of-Enzyme-Activity' and control[1] == '-' and regulator in PATHWAYS
			}
		print(f'\t{aa}: {", ".join(sorted(regulated))}')


if __name__ == '__main__':
	regulators, regulatees = load_regulation()
	summarize_aa_control(regulators)
	summarize_pathway_control(regulatees)
