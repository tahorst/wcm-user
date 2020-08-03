#! /usr/bin/env python

"""
Plot outputs from NCA on a test data set with figure comparison from
https://www.eee.hku.hk/~cqchang/gNCA-fig.pdf.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np

import nca


BASE_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
NCA_DATA_DIR = os.path.join(DATA_DIR, 'nca')
NCA_TEST_TOPOLOGY_FILE = os.path.join(NCA_DATA_DIR, 'subnet1_top.tsv')
NCA_TEST_EXPRESSION_FILE = os.path.join(NCA_DATA_DIR, 'subnet1_data.tsv')

OUTPUT_DIR = os.path.join(BASE_DIR, 'out')
if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)


def test_nca(method: str) -> None:
    """
    Plots TFAs to match results from Chang et al. Bioinformatics. 2008. Fig 5.

    Args:
        method: NCA method to use (function must be defined in nca.py)
    """

    def plot_output(P, start, samples, label, method):
        tfs = np.array([1, 8, 9, 14, 15, 21, 30, 32, 33, 34, 35])
        data = P[tfs, start:start+samples]
        n_tfs = data.shape[0]

        plt.figure(figsize=(5, 15))
        for i, tf in enumerate(data):
            ax = plt.subplot(n_tfs, 1, i+1)
            ax.plot(range(len(tf)), tf)

        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f'{method}_{label}.png'))
        plt.close('all')

    with open(NCA_TEST_EXPRESSION_FILE) as f:
        reader = csv.reader(f, delimiter='\t')
        expression = np.array(list(reader), float)

    with open(NCA_TEST_TOPOLOGY_FILE) as f:
        reader = csv.reader(f, delimiter='\t')
        topology = np.array(list(reader), float)

    A, P = getattr(nca, method)(expression, topology)

    plot_output(P, 0, 14, 'elutriation', method)
    plot_output(P, 14, 18, 'alpha', method)
    plot_output(P, 32, 24, 'cdc', method)
    plot_output(P, 56, 13, 'cycle', method)


if __name__ == '__main__':
    for method in nca.METHODS:
        test_nca(method)
