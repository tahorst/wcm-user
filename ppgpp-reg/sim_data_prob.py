#! /usr/bin/env python

"""
Script to compare synthesis probabilities and expression calculated from old
parca methods to regulated by ppGpp.
"""

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
SIM_DATA_FILE = os.path.join(DATA_DIR, 'sim_data_ppgpp_prob.cp')

CONDITIONS = ['with_aa', 'basal', 'no_oxygen']
MICROMOLAR = units.umol / units.L


def plot_ax(ax, old, new, highlighted, title, label, show_xlabel):
	"""
	Plot data on a given axis.
	"""

	for c, mask in highlighted.items():
		ax.plot(old[mask], new[mask], 'o', color=c, alpha=0.5)

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


if __name__ == '__main__':
	with open(SIM_DATA_FILE) as f:
		sim_data = cPickle.load(f)
	transcription = sim_data.process.transcription
	rna_data = transcription.rnaData

	# All genes
	plot_synth_prob_comparison(sim_data, 'all')
	plot_expression_comparison(sim_data, 'all')

	# Highlight regulated genes
	is_ppgpp_regulated = np.array([rna[:-3] in transcription.ppgpp_regulated_genes for rna in rna_data['id']])
	highlighted = {
		'b': ~is_ppgpp_regulated,
		'r': is_ppgpp_regulated,
		}
	plot_synth_prob_comparison(sim_data, 'ppgpp_reg', highlighted=highlighted)
	plot_expression_comparison(sim_data, 'ppgpp_reg', highlighted=highlighted)

	# Highlight stable RNA
	stable_rna = rna_data['isTRna'] | rna_data['isRRna']
	highlighted = {
		'b': ~stable_rna,
		'r': stable_rna,
		}
	plot_synth_prob_comparison(sim_data, 'stable_rna', highlighted=highlighted)
	plot_expression_comparison(sim_data, 'stable_rna', highlighted=highlighted)

	# Ribosome/RNAP related mRNA
	polymerizing_mrna = rna_data['isRProtein'] | rna_data['isRnap']
	highlighted = {
		'b': ~polymerizing_mrna,
		'r': polymerizing_mrna,
		}
	plot_synth_prob_comparison(sim_data, 'polymerizing_mrna', highlighted=highlighted)
	plot_expression_comparison(sim_data, 'polymerizing_mrna', highlighted=highlighted)
