#! /usr/bin/env python

"""
Load and explore regulation data.
"""

import csv
import os
import re
from typing import Dict, Tuple


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

# Mapping for easier to read output and tRNA matching
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
# Mapping to all genes in each AA synthesis pathway curated from EcoCyc
# Some genes form complexes with each other or catalyze reactions in parallel
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
# Mapping to transcription factors that an amino acid can bind to curated from regulonDB
TRANSCRIPTION_FACTORS = {
	'alaS': ['L-alanine'],
	'argP': ['L-arginine', 'L-lysine'],  # TODO: lysine bound is inactive
	'argR': ['L-arginine'],
	'asnC': ['L-asparagine'],  # TODO: bound is inactive
	'gcvA': ['glycine'],  # TODO: bound is inactive
	'lrp': ['L-leucine'],
	'tyrR': ['L-phenylalanine', 'L-tryptophan', 'L-tyrosine'],
	'trpR': ['L-tryptophan'],
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

def summarize_trna_control(regulation: Dict[str, Dict[str, Tuple[str, str]]]) -> None:
	"""
	Show tRNA that regulate enzymatic expression in synthesis pathways for
	other amino acids.

	Args:
		regulation: dictionary of regulatee to regulator pair mappings,
			see load_regulation()
	"""

	print('Pathway expression regulated by other tRNA:')
	for aa, enzymes in PATHWAYS.items():
		regulated = {
			regulator[-4:-1]
			for enzyme in enzymes
			for regulator, control in sorted(regulation.get(enzyme, {}).items())
			if control[0] == 'Ribosome-Mediated-Attenuation' and control[1] == '-' and 'tRNA' in regulator
			}
		print(f'\t{aa}: {", ".join(sorted(regulated))}')

def summarize_tf_control(regulation: Dict[str, Dict[str, Tuple[str, str]]]) -> None:
	"""
	Show amino acids that regulate enzymatic expression in synthesis pathways
	for other amino acids through transcription factor binding.

	Args:
		regulation: dictionary of regulatee to regulator pair mappings,
			see load_regulation()
	"""

	print('Pathway expression regulated by other amino acids through TFs:')
	for aa, enzymes in PATHWAYS.items():
		regulated = {
			AMINO_ACIDS[ligand]
			for enzyme in enzymes
			for regulator, control in sorted(regulation.get(enzyme, {}).items())
			for ligand in TRANSCRIPTION_FACTORS.get(regulator, [])
			}
		print(f'\t{aa}: {", ".join(sorted(regulated))}')

def summarize_rnap_control(regulation: Dict[str, Dict[str, Tuple[str, str]]]) -> None:
	"""
	Show RNAP regulation (ppGpp or DksA) of enzymatic expression in synthesis
	pathways for amino acids.

	Args:
		regulation: dictionary of regulatee to regulator pair mappings,
			see load_regulation()
	"""

	print('Pathway expression regulated by RNAP allosteric effects:')
	for aa, enzymes in PATHWAYS.items():
		regulated = {
			regulator
			for enzyme in enzymes
			for regulator, control in sorted(regulation.get(enzyme, {}).items())
			if control[0] == 'Allosteric-Regulation-of-RNAP'
			}
		print(f'\t{aa}: {", ".join(sorted(regulated))}')

def summarize_other_control(regulation: Dict[str, Dict[str, Tuple[str, str]]]) -> None:
	"""
	Show other regulation of enzymes in synthesis pathways for amino acids
	that is not captured in other functions above.

	Args:
		regulation: dictionary of regulatee to regulator pair mappings,
			see load_regulation()
	"""

	print('Pathway regulated by other means:')
	for aa, enzymes in PATHWAYS.items():
		regulated = {
			f'{control[1]} {regulator} -> {enzyme} ({control[0]})'
			for enzyme in enzymes
			for regulator, control in sorted(regulation.get(enzyme, {}).items())
			if not (control[0] == 'Regulation-of-Enzyme-Activity' and control[1] == '-' and regulator in PATHWAYS)
			and not (control[0] == 'Ribosome-Mediated-Attenuation' and control[1] == '-' and 'tRNA' in regulator)
			and regulator not in TRANSCRIPTION_FACTORS
			and control[0] != 'Allosteric-Regulation-of-RNAP'
			}
		sep = '\n\t\t'
		print(f'\t{aa}:{sep}{sep.join(sorted(regulated))}')


if __name__ == '__main__':
	regulators, regulatees = load_regulation()
	summarize_aa_control(regulators)
	summarize_pathway_control(regulatees)
	summarize_trna_control(regulatees)
	summarize_tf_control(regulatees)
	summarize_rnap_control(regulatees)
	summarize_other_control(regulatees)
