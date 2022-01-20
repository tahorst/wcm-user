#! /usr/bin/env python

"""
Example distributions based on summary statistic for # TODO: paper
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


# Data from Table S3
# (label, mean, std)
DOUBLING_TIME = [
    ('Rich', 22.50, 4.63),
    ('Glc+12AA', 26.61, 3.79),
    ('Glc+6AA', 30.14, 4.63),
    ('Glc', 37.66, 5.83),
    ]
# Data from # TODO
# (label, doublings/hr, CV)
GROWTH_RATE = [
    ('Rich', 2.4, 0.15),
    ('Glc', 0.9, 0.21),
    ]

FONT_SIZE = 9


def plot(x_min, x_max, normalization=1, filename='distributions'):
    plt.figure(figsize=(10, 5))
    ax = plt.gca()

    x = np.linspace(x_min, x_max, 1000)
    for label, mean, std in DOUBLING_TIME:
        new_mean = mean / normalization
        new_std = std / normalization
        pdf = stats.norm.pdf(x, new_mean, new_std)
        line = plt.plot(x, pdf, alpha=0.5, label=f'{label}: {new_mean:.1f} +/- {new_std:.2f}')[0]
        color = line.get_color()
        ax.axvline(new_mean, color=color, linestyle='--', alpha=0.5)

    ax.set_xlim([x_min, x_max])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Doubling time', fontsize=FONT_SIZE)
    ax.tick_params(labelsize=FONT_SIZE)
    ax.legend()

    plt.savefig(f'{filename}.pdf')


if __name__ == '__main__':
    plot(0, 80)
    plot(0, 4, normalization=DOUBLING_TIME[0][1], filename='distributions-normalized')
