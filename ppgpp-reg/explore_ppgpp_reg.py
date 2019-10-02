#! /usr/bin/env python

"""
Explore data regarding ppGpp regulation curated from EcoCyc.

Regulation file is not currently loaded into sim_data so it comes from a separate
file but should be integrated in the future.
"""

from __future__ import division

import cPickle
import csv
import os

import matplotlib.pyplot as plt
import numpy as np


# Directories
FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
OUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUT_DIR):
	os.mkdir(OUT_DIR)

# File paths
PPGPP_REG_FILE = os.path.join(DATA_DIR, 'ppgpp_regulation.tsv')
SIM_DATA_FILE = os.path.join(DATA_DIR, 'sim_data.cp')
SYNONYMS_FILE = os.path.join(DATA_DIR, 'gene_ids.tsv')


def load_regulation():
	with open(PPGPP_REG_FILE) as f:
		reader = csv.reader(f, delimiter='\t')
		header = reader.next()
		while header[0].startswith('#'):
			header = reader.next()
		data = np.array(list(reader))
		data[data == ''] = '0'

	genes = data[:, 0]
	ppgpp_reg = data[:, 1].astype(int)
	ppgpp_dksa_reg = data[:, 2].astype(int)

	return genes, ppgpp_reg, ppgpp_dksa_reg

def load_sim_data():
	with open(SIM_DATA_FILE) as f:
		sim_data = cPickle.load(f)

	return sim_data

def load_synonyms():
	synonyms = {}
	with open(SYNONYMS_FILE) as f:
		reader = csv.reader(f, delimiter='\t')
		headers = reader.next()
		for line in reader:
			if line[1]:
				synonyms.update({s: line[0] for s in line[1].strip('("")').split('" "')})

	return synonyms

def print_is_fraction(rna_data, key, idx):
	total = np.sum(rna_data[key])
	regulated = np.sum(rna_data[key][idx])
	print('\t{}: {:.1f}% ({}/{})'.format(key, 100 * regulated / total, regulated, total))

def plot_expression(expression, regulation, genes):
	print('Plotting expression in {} ...'.format(OUT_DIR))
	for exp, reg, gene in zip(expression, regulation, genes):
		plt.figure()

		plt.plot(exp, 'o')
		x = plt.xlim()
		if reg < 0:
			x = x[::-1]
		plt.plot(x, plt.ylim(), 'k--')
		plt.title(gene)

		plt.savefig(os.path.join(OUT_DIR, '{}.png'.format(gene)))
		plt.close('all')


if __name__ == '__main__':
	genes, ppgpp_reg, ppgpp_dksa_reg = load_regulation()
	sim_data = load_sim_data()
	synonyms = load_synonyms()

	replication = sim_data.process.replication
	transcription = sim_data.process.transcription
	complexation = sim_data.process.complexation

	gene_data = replication.geneData
	rna_expression = transcription.rnaExpression
	rna_data = transcription.rnaData

	# Duplicate entries
	unique_genes, counts = np.unique(genes, return_counts=True)
	duplicate_genes = unique_genes[counts > 1]
	print('{}/{} unique genes/total entries'.format(len(unique_genes), len(genes)))
	print('Duplicate genes:')
	for g in duplicate_genes:
		print('\t{}'.format(g))

	# Contradictory directions
	for g in duplicate_genes:
		idx = np.where(genes == g)
		ppgpp = ppgpp_reg[idx]
		ppgpp_dksa = ppgpp_dksa_reg[idx]
		reg = np.hstack((ppgpp, ppgpp_dksa))
		if -1 in reg and 1 in reg:
			print('Opposing regulation for {}'.format(g))
			print('\tppGpp reg: {}'.format(ppgpp))
			print('\tppGpp-DksA reg: {}'.format(ppgpp_dksa))

	# Check gene symbols in WCM
	wcm_genes = set(gene_data['symbol'])
	genes_not_found = [g for g in genes if g not in wcm_genes]
	print('\nGene symbol not found in wcm:')
	for g in genes_not_found:
		has_synonym = ' but in synonyms' if g in synonyms else ''
		print('\t{}{}'.format(g, has_synonym))

	# Trim genes not in wcm
	mask = np.array([g not in genes_not_found for g in genes])
	genes = genes[mask]
	ppgpp_reg = ppgpp_reg[mask]
	ppgpp_dksa_reg = ppgpp_dksa_reg[mask]
	total_reg = ppgpp_reg + ppgpp_dksa_reg
	total_reg = total_reg / np.fmax(1, np.abs(total_reg))
	n_genes = len(genes)

	# Convert to wcm framework
	symbol_to_id = {g['symbol']: g['name'] for g in gene_data}
	gene_id_to_rna_idx = {r['geneId']: i for i, r in enumerate(rna_data)}
	rna_id_to_monomer_id = {g['rnaId']: g['monomerId'] for g in gene_data}
	gene_ids = np.array([symbol_to_id[g] for g in genes])
	rna_data_idx = np.array([gene_id_to_rna_idx[g] for g in gene_ids])

	# Counts of RNAP, ribosomes, tRNA, synthetases
	synthetase_monomers = np.hstack([
		complexation.getMonomers(s)['subunitIds']
		if s in complexation.complexNames
		else [s]
		for s in transcription.synthetase_names
	])
	synthetase_monomers = np.array([m[:-3] for m in synthetase_monomers])
	synthetase_data = {
		'isSynthetase': np.array([rna_id_to_monomer_id[r[:-3]] in synthetase_monomers for r in rna_data['id']])
	}
	print('\nFraction regulated for specific groups:')
	print_is_fraction(rna_data, 'isRRna5S', rna_data_idx)
	print_is_fraction(rna_data, 'isRRna16S', rna_data_idx)
	print_is_fraction(rna_data, 'isRRna23S', rna_data_idx)
	print_is_fraction(rna_data, 'isRRna', rna_data_idx)
	print_is_fraction(rna_data, 'isTRna', rna_data_idx)
	print_is_fraction(rna_data, 'isRProtein', rna_data_idx)
	print_is_fraction(rna_data, 'isRnap', rna_data_idx)
	print_is_fraction(synthetase_data, 'isSynthetase', rna_data_idx)
	print('rProtein not regulated:')
	for i, r in enumerate(rna_data):
		if rna_data['isRProtein'][i] and i not in rna_data_idx:
			print('\t{}'.format(rna_id_to_monomer_id[r['id'][:-3]]))
	print('tRNA regulated:')
	for i, r in enumerate(rna_data):
		if rna_data['isTRna'][i] and i in rna_data_idx:
			print('\t{}'.format(r['id'][:-3]))
	print('Synthetases regulated:')
	for i, r in enumerate(rna_data):
		if synthetase_data['isSynthetase'][i] and i in rna_data_idx:
			print('\t{}'.format(rna_id_to_monomer_id[r['id'][:-3]]))

	# Expression in different conditions
	conditions = ['with_aa', 'basal', 'no_oxygen']
	expression = np.vstack([rna_expression[c] for c in conditions])[:, rna_data_idx]
	aa_consistent = np.sum(np.sign(expression[1, :] - expression[0, :]) == total_reg)
	anaerobic_consistent = np.sum(np.sign(expression[2, :] - expression[1, :]) == total_reg)
	print('\nRelative expression consistent for {:.1f}% of genes ({}/{}) in basal to with_aa'
		  .format(100 * aa_consistent / n_genes, aa_consistent, n_genes))
	print('Relative expression consistent for {:.1f}% of genes ({}/{}) in basal to no_oxygen'
		  .format(100 * anaerobic_consistent / n_genes, anaerobic_consistent, n_genes))
	plot_expression(expression.T, total_reg, genes)
