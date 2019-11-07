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
EXP_OUT_DIR = os.path.join(OUT_DIR, 'expression')
if not os.path.exists(OUT_DIR):
	os.mkdir(OUT_DIR)
if not os.path.exists(EXP_OUT_DIR):
	os.mkdir(EXP_OUT_DIR)

# File paths
PPGPP_REG_FILE = os.path.join(DATA_DIR, 'ppgpp_regulation.tsv')
SIM_DATA_FILE = os.path.join(DATA_DIR, 'sim_data.cp')
SYNONYMS_FILE = os.path.join(DATA_DIR, 'gene_ids.tsv')
FOLD_CHANGE_FILE = os.path.join(DATA_DIR, 'ppgpp_fc.tsv')


def load_regulation():
	with open(PPGPP_REG_FILE) as f:
		reader = csv.reader(f, delimiter='\t')
		header = reader.next()
		while header[0].startswith('#'):
			header = reader.next()
		data = np.array(list(reader))
		data[data == ''] = '0'

	genes = data[:, header.index('Gene')]
	ppgpp_reg = data[:, header.index('ppGpp')].astype(int)
	ppgpp_dksa_reg = data[:, header.index('DksA-ppGpp')].astype(int)
	curated = data[:, header.index('Curated Gene')]
	original_gene_mapping = {g: c for g, c in zip(genes, curated) if g != c}

	return genes, ppgpp_reg, ppgpp_dksa_reg, original_gene_mapping

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

def load_fc():
	valid_categories = {'A', 'B', 'C', 'D'}  # statistically significant change and not small RNA
	with open(FOLD_CHANGE_FILE) as f:
		reader = csv.reader(f, delimiter='\t')
		headers = reader.next()
		data = np.array(list(reader))

	data[data == 'N/A'] = '0'
	genes = data[:, headers.index('Gene')]
	early_expression = data[:, headers.index('1+2+ 5 min')].astype(float)
	late_expression = data[:, headers.index('1+2+ 10 min')].astype(float)
	discard_early = np.array([d not in valid_categories
		for d in data[:, headers.index('1+2+ 5 min Category')]], bool)
	discard_late = np.array([d not in valid_categories
		for d in data[:, headers.index('1+2+ 10 min Category')]])
	early_expression[discard_early] = 0
	late_expression[discard_late] = 0

	return genes, early_expression, late_expression

def print_is_fraction(rna_data, key, neg_idx, pos_idx):
	total = np.sum(rna_data[key])
	negative = np.sum(rna_data[key][neg_idx])
	positive = np.sum(rna_data[key][pos_idx])
	regulated = negative + positive
	print('\t{}: {:.1f}% ({}/{}), +:{}, -:{}'.format(key, 100 * regulated / total, regulated, total, positive, negative))

def plot_gene_comparison(gene, fc, exp):
	exp = exp / exp.max()  # normalize top to 1

	# Calculate fold change comparison data
	x_fc = [0, len(exp)-1]
	if fc <= 0:
		fc_low = exp.max()
	else:
		fc_low = exp.min()
	fc_high = fc_low * 2**fc
	y_fc = [fc_low, fc_high]

	# Plot data
	plt.figure()
	plt.plot(exp, 'o')
	plt.plot(x_fc, y_fc, 'x')

	# Format plot
	plt.legend(['WCM expression', 'FC data'])
	plt.xticks([0, 1, 2], ['AA', 'Basal', 'Anaerobic'])
	plt.ylabel('Normalized expression')

	# TODO: change directory if looping over all genes
	plt.savefig(os.path.join(OUT_DIR, '{}.png'.format(gene)))
	plt.close('all')

def plot_expression(expression, regulation, genes):
	print('\nPlotting expression in {} ...'.format(EXP_OUT_DIR))
	for exp, reg, gene in zip(expression, regulation, genes):
		plt.figure()

		plt.plot(exp, 'o')
		x = plt.xlim()
		if reg < 0:
			x = x[::-1]
		plt.plot(x, plt.ylim(), 'k--')
		plt.title(gene)

		plt.savefig(os.path.join(EXP_OUT_DIR, '{}.png'.format(gene)))
		plt.close('all')


if __name__ == '__main__':
	genes, ppgpp_reg, ppgpp_dksa_reg, original_gene_mapping = load_regulation()
	sim_data = load_sim_data()
	synonyms = load_synonyms()
	fc_genes, fc_early, fc_late = load_fc()

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
	total_reg = np.sign(ppgpp_reg + ppgpp_dksa_reg)
	n_genes = len(genes)

	# Convert to wcm framework
	symbol_to_id = {g['symbol']: g['name'] for g in gene_data}
	gene_id_to_rna_idx = {r['geneId']: i for i, r in enumerate(rna_data)}
	rna_id_to_monomer_id = {g['rnaId']: g['monomerId'] for g in gene_data}
	gene_ids = np.array([symbol_to_id[g] for g in genes])
	rna_data_idx = np.array([gene_id_to_rna_idx[g] for g in gene_ids])
	neg_rna_data_idx = np.array([gene_id_to_rna_idx[g] for g, r in zip(gene_ids, total_reg) if r < 0])
	pos_rna_data_idx = np.array([gene_id_to_rna_idx[g] for g, r in zip(gene_ids, total_reg) if r > 0])

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
	print_is_fraction(rna_data, 'isRRna5S', neg_rna_data_idx, pos_rna_data_idx)
	print_is_fraction(rna_data, 'isRRna16S', neg_rna_data_idx, pos_rna_data_idx)
	print_is_fraction(rna_data, 'isRRna23S', neg_rna_data_idx, pos_rna_data_idx)
	print_is_fraction(rna_data, 'isRRna', neg_rna_data_idx, pos_rna_data_idx)
	print_is_fraction(rna_data, 'isTRna', neg_rna_data_idx, pos_rna_data_idx)
	print_is_fraction(rna_data, 'isRProtein', neg_rna_data_idx, pos_rna_data_idx)
	print_is_fraction(rna_data, 'isRnap', neg_rna_data_idx, pos_rna_data_idx)
	print_is_fraction(synthetase_data, 'isSynthetase', neg_rna_data_idx, pos_rna_data_idx)
	print('Positive rProtein regulation:')
	for i, r in enumerate(rna_data):
		if rna_data['isRProtein'][i] and i in pos_rna_data_idx:
			print('\t{}'.format(rna_id_to_monomer_id[r['id'][:-3]]))
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

	# Compare to ppGpp FC data
	unique_genes = np.unique(genes)
	fc_gene_set = set(fc_genes)
	not_included = [g for g in unique_genes if g not in fc_gene_set and original_gene_mapping.get(g, '') not in fc_gene_set]
	n_total = len(unique_genes)
	n_included = n_total - len(not_included)
	print('\nFC data for {}/{} genes from EcoCyc that are in WCM'.format(n_included, n_total))
	print('Genes without FC data:')
	for g in not_included:
		print('\t{}'.format(g))

	# Get FC data for WCM genes
	fc_early_dict = {g: f for g, f in zip(fc_genes, fc_early)}
	fc_late_dict = {g: f for g, f in zip(fc_genes, fc_late)}
	fc_genes_in_wcm = np.array([
		g if g in fc_gene_set else original_gene_mapping.get(g)
		for g in genes
		if g in fc_gene_set or original_gene_mapping.get(g, '') in fc_gene_set
	])
	fc_mask = np.array([
		g in fc_gene_set or original_gene_mapping.get(g, '') in fc_gene_set
		for g in genes
	])
	reg_direction = np.sign(total_reg[fc_mask])

	## 5 min data (early)
	fc_in_wcm = np.array([fc_early_dict[g] for g in fc_genes_in_wcm])
	print('\nFC data for 5 min:')
	print('\t{} positive fold changes (mean FC: {:.2f})'
		  .format(np.sum(fc_in_wcm > 0), np.mean(fc_in_wcm[fc_in_wcm > 0])))
	print('\t{} negative fold changes (mean FC: {:.2f})'
		  .format(np.sum(fc_in_wcm < 0), np.mean(fc_in_wcm[fc_in_wcm < 0])))
	print('\t{} with no fold change'.format(np.sum(fc_in_wcm == 0)))
	fc_direction = np.sign(fc_in_wcm)
	fc_nonzero = fc_direction != 0
	consistent = reg_direction[fc_nonzero] == fc_direction[fc_nonzero]
	n_consistent = np.sum(consistent)
	print('\n\tFC consistent with curated regulation for {}/{} genes'.format(n_consistent, np.sum(fc_nonzero)))
	print('\tInconsistent genes with 5 and 10 min FC:')
	for g in fc_genes_in_wcm[fc_nonzero][~consistent]:
		print('\t\t{}: {}\t{}'.format(g, fc_early_dict[g], fc_late_dict[g]))

	## 10 min data (late)
	fc_in_wcm = np.array([fc_late_dict[g] for g in fc_genes_in_wcm])
	print('\nFC data for 10 min:')
	print('\t{} positive fold changes (mean FC: {:.2f})'
		  .format(np.sum(fc_in_wcm > 0), np.mean(fc_in_wcm[fc_in_wcm > 0])))
	print('\t{} negative fold changes (mean FC: {:.2f})'
		  .format(np.sum(fc_in_wcm < 0), np.mean(fc_in_wcm[fc_in_wcm < 0])))
	print('\t{} with no fold change'.format(np.sum(fc_in_wcm == 0)))
	fc_direction = np.sign(fc_in_wcm)
	fc_nonzero = fc_direction != 0
	consistent = reg_direction[fc_nonzero] == fc_direction[fc_nonzero]
	n_consistent = np.sum(consistent)
	print('\n\tFC consistent with curated regulation for {}/{} genes'.format(n_consistent, np.sum(fc_nonzero)))
	print('\tInconsistent genes with 5 and 10 min FC:')
	for g in fc_genes_in_wcm[fc_nonzero][~consistent]:
		print('\t\t{}: {}\t{}'.format(g, fc_early_dict[g], fc_late_dict[g]))

	# Plot SpoT expression comparison
	gene_of_interest = 'spoT'
	fc = fc_early_dict[gene_of_interest]
	exp = expression[:, genes==gene_of_interest]
	plot_gene_comparison(gene_of_interest, fc, exp)

	# Time intensive tasks
	plot_expression(expression.T, total_reg, genes)
