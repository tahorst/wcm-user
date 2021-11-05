#!/usr/bin/env python
"""
Plot the expression included in the model based on output generated with
runscripts/reflect/model_inspection.py in the output file regulation.tsv.
"""

import matplotlib.pyplot as plt
import numpy as np


COLORS = ['black', 'orange', 'green', 'blue']
LABELS = ['Regulators', 'Positive', 'Positive and negative', 'Negative']


def plot_bar(data):
    x = np.arange(data.shape[1])
    width = 0.8 / data.shape[0]
    offset = 0.4 - width / 2
    for i, (d, color, label) in enumerate(zip(data, COLORS, LABELS)):
        plt.bar(x - offset + i * width, d, width=width, color=color, alpha=0.4, label=label)


if __name__ == '__main__':
    # Data rows line up with xlabels and columns line up with LABELS
    old_data = np.array([
        [0, 341, 0, 200],
        [22, 259, 49, 143],
        [22, 259, 49, 143],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        ]).T
    new_data = np.array([
        [0, 1070, 0, 879],
        [32, 599, 301, 468],
        [24, 525, 224, 352],
        [1, 193, 0, 201],
        [7, 0, 0, 31],
        ]).T
    xlabels = [
        'All regulatory\ninteractions',
        'All regulated\ngenes',
        'Genes regulated\nby transcription\nfactors',
        'Genes regulated\nby ppGpp\nregulation',
        'Genes regulated\nby transcriptional\nattenuation',
        ]

    x = np.arange(len(xlabels))

    plt.figure(figsize=(6, 6))

    plot_bar(new_data)
    plt.legend()
    plot_bar(old_data)

    plt.xticks(x, xlabels, fontsize=6)
    plt.ylabel('Count')

    plt.savefig('regulation.pdf')
