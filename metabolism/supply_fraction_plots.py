#! /usr/bin/env python

"""
Plots to show the different components to the amino acid supply scaling function.
"""

from __future__ import division

import os

import matplotlib.pyplot as plt
import numpy as np


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUT_DIR):
	os.mkdir(OUT_DIR)


def plot_fraction(ax, x, y, title, c_basal, c_rich):
	# Plot
	ax.set_xscale('log')
	ax.plot(x, y)
	ax.axhline(1, linestyle='--', color='k', linewidth=0.5)

	# x-axis
	ax.set_xticks([], minor=True)
	ax.set_xticks([c_basal, c_rich])
	ax.set_xticklabels(['Basal\nTarget', 'Rich\nTarget'], fontsize=6)
	ax.set_xlabel('AA concentration', fontsize=7)

	# y-axis
	ax.set_ylim(0, 1.5)
	ax.set_yticks([0, 1])
	ax.set_yticklabels([0, 1], fontsize=6)
	ax.set_ylabel('{} fraction'.format(title), fontsize=7)

def plot_conditions(ax, c, fraction, c_target, c_other, labels):
	# Plot
	ax.set_xscale('log')
	ax.plot(c, fraction)
	ax.axhline(1, linestyle='--', color='k', linewidth=0.5)
	ax.axvline(c_target, linestyle='--', color='k', linewidth=0.5)

	# x-axis
	ax.set_xticks([], minor=True)
	ax.set_xticks([c_target, c_other])
	ax.set_xticklabels(labels, fontsize=6)
	ax.set_xlabel('AA concentration', fontsize=7)

	# y-axis
	ax.set_ylim(0, 2)
	ax.set_yticks([0, 1])
	ax.set_yticklabels([0, 1], fontsize=6)
	ax.set_ylabel('Total fraction', fontsize=7)


if __name__ == '__main__':
	c = np.logspace(-2, 4, 10001)
	c_basal = 5
	c_rich = 50
	f1 = 0.1
	f2 = 0.1
	KI = f1 * c_basal / (1 - f1)
	KM = (1 / f2 - 1) * c_rich

	supply = c / c - f1 + c_basal / (KM + c_basal)
	inhibited_synthesis = 1 / (1 + c / KI)
	import_rate = c / c - (supply + 1 / (1 + c_rich / KI) - f2)
	export_rate = c / (KM + c)
	basal_fraction = supply + inhibited_synthesis - export_rate
	rich_fraction = supply + inhibited_synthesis + import_rate - export_rate

	# Plot fractions
	plt.figure(figsize=(6, 6))

	## Base supply
	ax = plt.subplot(2, 2, 1)
	plot_fraction(ax, c, supply, 'Supply', c_basal, c_rich)

	## Inhibited synthesis
	ax = plt.subplot(2, 2, 2)
	plot_fraction(ax, c, inhibited_synthesis, 'Inhibited synthesis', c_basal, c_rich)

	## Import
	ax = plt.subplot(2, 2, 3)
	plot_fraction(ax, c, import_rate, 'Import', c_basal, c_rich)

	## Export
	ax = plt.subplot(2, 2, 4)
	plot_fraction(ax, c, export_rate, 'Export', c_basal, c_rich)

	## Save plot
	plt.tight_layout()
	plt.savefig(os.path.join(OUT_DIR, 'fraction_plot.pdf'))

	# Plot conditions
	plt.figure(figsize=(3, 6))

	## Basal
	ax = plt.subplot(2, 1, 1)
	plot_conditions(ax, c, basal_fraction, c_basal, c_rich, ['Basal\nTarget', 'Rich\nTarget'])

	## Rich
	ax = plt.subplot(2, 1, 2)
	plot_conditions(ax, c, rich_fraction, c_rich, c_basal, ['Rich\nTarget', 'Basal\nTarget'])

	## Save plot
	plt.tight_layout()
	plt.savefig(os.path.join(OUT_DIR, 'condition_fractions.pdf'))
