#! /usr/bin/env python

"""
Script to compare synthesis probabilities and expression calculated from old
parca methods to regulated by ppGpp.
"""

import argparse
import cPickle
import os

import matplotlib.pyplot as plt
import numpy as np

from wholecell.utils import units


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
OUT_DIR = os.path.join(FILE_LOCATION, 'out', 'ppgpp-regulation')
if not os.path.exists(OUT_DIR):
	os.mkdir(OUT_DIR)
SIM_DATA_FILE = os.path.join(DATA_DIR, 'sim_data_{}.cp')

CONDITIONS = ['with_aa', 'basal', 'no_oxygen']
MICROMOLAR = units.umol / units.L
DEFAULT_HASH = 'da3685f97'


def plot_ax(ax, old, new, highlighted, title, label, show_xlabel):
	"""
	Plot data on a given axis.
	"""

	for c, mask in highlighted.items():
		ax.plot(np.log10(old[mask]), np.log10(new[mask]), 'o', color=c, alpha=0.5)

	ax.set_title(title)

	x_min, x_max = ax.get_xlim()
	y_min, y_max = ax.get_ylim()
	ax_min = min(x_min, y_min)
	ax_max = max(x_max, y_max)
	ax.set_xlim([ax_min, ax_max])
	ax.set_ylim([ax_min, ax_max])
	if show_xlabel:
		ax.set_xlabel('Old {}'.format(label))
	ax.set_ylabel('ppGpp {}'.format(label))

def plot_synth_prob_comparison(sim_data, label, highlighted=None):
	"""
	Plot scatterplots between new and old methods for RNA synthesis probabilities.
	"""

	n_conditions = len(CONDITIONS)
	filename = os.path.join(OUT_DIR, '{}_synth_prob.png'.format(label))
	transcription = sim_data.process.transcription

	plt.figure(figsize=(5, 10))

	for i, condition in enumerate(CONDITIONS):
		old = transcription.rnaSynthProb[condition]
		ppgpp = sim_data.growthRateParameters.getppGppConc(sim_data.conditionToDoublingTime[condition])
		new = transcription.synth_prob_from_ppgpp(ppgpp)

		ax = plt.subplot(n_conditions, 1, i+1)
		if highlighted is None:
			highlighted = {'b': np.ones(len(new), bool)}
		show_xlabel = i == n_conditions - 1
		plot_ax(ax, old, new, highlighted, condition, 'synth prob', show_xlabel)

	plt.tight_layout()
	plt.savefig(filename)
	plt.close('all')

def plot_expression_comparison(sim_data, label, highlighted=None):
	"""
	Plot scatterplots between new and old methods for RNA expression.
	"""

	n_conditions = len(CONDITIONS)
	filename = os.path.join(OUT_DIR, '{}_expression.png'.format(label))
	transcription = sim_data.process.transcription
	km = transcription.ppgpp_km.asNumber(MICROMOLAR)

	plt.figure(figsize=(5, 10))

	for i, condition in enumerate(CONDITIONS):
		old = transcription.rnaExpression[condition]
		ppgpp = sim_data.growthRateParameters.getppGppConc(sim_data.conditionToDoublingTime[condition]).asNumber(MICROMOLAR)
		f = ppgpp**2 / (km**2 + ppgpp**2)
		new = transcription.exp_free * (1 - f) + transcription.exp_ppgpp * f

		ax = plt.subplot(n_conditions, 1, i+1)
		if highlighted is None:
			highlighted = {'b': np.ones(len(new), bool)}
		show_xlabel = i == n_conditions - 1
		plot_ax(ax, old, new, highlighted, condition, 'expression', show_xlabel)

	plt.tight_layout()
	plt.savefig(filename)
	plt.close('all')

def plot_split_expression_comparison(sim_data, label, highlighted=None):
	"""
	Plot scatterplots between new and old methods for RNA expression showing
	expression for free RNAP and ppGpp bound RNAP.
	"""

	n_conditions = len(CONDITIONS)
	filename = os.path.join(OUT_DIR, '{}_expression_split.png'.format(label))
	transcription = sim_data.process.transcription

	plt.figure(figsize=(8, 10))

	for i, condition in enumerate(CONDITIONS):
		old = transcription.rnaExpression[condition]

		for j, attr in enumerate(['exp_free', 'exp_ppgpp']):
			new = getattr(transcription, attr)
			ax = plt.subplot(n_conditions, 2, 2*i+j+1)
			if highlighted is None:
				highlighted = {'b': np.ones(len(new), bool)}
			show_xlabel = i == n_conditions - 1
			plot_ax(ax, old, new, highlighted, '{}: {}'.format(condition, attr), 'expression', show_xlabel)

	plt.tight_layout()
	plt.savefig(filename)
	plt.close('all')

def parse_args():
	# type: () -> argparse.Namespace
	"""
	Parses arguments from the command line.

	Returns:
		values of variables parsed from the command line
	"""

	parser = argparse.ArgumentParser()

	parser.add_argument('--hash',
		default=DEFAULT_HASH,
		help='Hash of sim_data object to use for analysis (default: {})'.format(DEFAULT_HASH))
	parser.add_argument('-l', '--label',
		default=None,
		help='Label to prepend to output files (default: selected hash)')

	return parser.parse_args()


if __name__ == '__main__':
	args = parse_args()

	with open(SIM_DATA_FILE.format(args.hash)) as f:
		sim_data = cPickle.load(f)
	transcription = sim_data.process.transcription
	rna_data = transcription.rnaData
	if args.label is None:
		label = '{}_'.format(args.hash)
	else:
		label = args.label

	# All genes
	plot_label = '{}all'.format(label)
	plot_synth_prob_comparison(sim_data, plot_label)
	plot_expression_comparison(sim_data, plot_label)
	plot_split_expression_comparison(sim_data, plot_label)

	# Highlight regulated genes
	rna_idx = {r[:-3]: i for i, r in enumerate(rna_data['id'])}
	neg_idx = np.array([rna_idx[r] for r, fc in zip(transcription.ppgpp_regulated_genes, transcription.ppgpp_fold_changes) if fc < 0])
	pos_idx = np.array([rna_idx[r] for r, fc in zip(transcription.ppgpp_regulated_genes, transcription.ppgpp_fold_changes) if fc > 0])
	neg_mask = np.zeros(len(rna_data), bool)
	neg_mask[neg_idx] = True
	pos_mask = np.zeros(len(rna_data), bool)
	pos_mask[pos_idx] = True
	highlighted = {
		'k': (~neg_mask) | (~pos_mask),
		'r': neg_mask,
		'g': pos_mask,
		}
	plot_label = '{}ppgpp_reg'.format(label)
	plot_synth_prob_comparison(sim_data, plot_label, highlighted=highlighted)
	plot_expression_comparison(sim_data, plot_label, highlighted=highlighted)
	plot_split_expression_comparison(sim_data, plot_label, highlighted=highlighted)

	# Highlight stable RNA
	stable_rna = rna_data['isTRna'] | rna_data['isRRna']
	highlighted = {
		'b': ~stable_rna,
		'r': stable_rna,
		}
	plot_label = '{}stable_rna'.format(label)
	plot_synth_prob_comparison(sim_data, plot_label, highlighted=highlighted)
	plot_expression_comparison(sim_data, plot_label, highlighted=highlighted)
	plot_split_expression_comparison(sim_data, plot_label, highlighted=highlighted)

	# Ribosome/RNAP related mRNA
	polymerizing_mrna = rna_data['isRProtein'] | rna_data['isRnap']
	highlighted = {
		'b': ~polymerizing_mrna,
		'r': polymerizing_mrna,
		}
	plot_label = '{}polymerizing_mrna'.format(label)
	plot_synth_prob_comparison(sim_data, plot_label, highlighted=highlighted)
	plot_expression_comparison(sim_data, plot_label, highlighted=highlighted)
	plot_split_expression_comparison(sim_data, plot_label, highlighted=highlighted)
