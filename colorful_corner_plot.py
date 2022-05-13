import corner
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import os

CORNER_KWARGS = dict(
    smooth=0.9,
    label_kwargs=dict(fontsize=16),
    title_kwargs=dict(fontsize=16),
    quantiles=[0.05, 0.95],
    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
    plot_density=False,
    plot_datapoints=False,
    fill_contours=True,
    # show_titles=True,
    # max_n_ticks=3,
)


def overlaid_corner(samples_list, sample_labels, corner_range, param_labels, output_dir, filename):
    """Plots multiple corners on top of each other"""
    # get some constants
    n = len(samples_list)
    _, ndim = samples_list[0].shape
    max_len = max([len(s) for s in samples_list])
    cmap = plt.cm.get_cmap('jet', n)
    colors = [cmap(i) for i in range(n)]
    fig = corner.corner(
        samples_list[0],
        range=corner_range,
        color=colors[0],
        **CORNER_KWARGS
    )
    for idx in range(1, n):
        fig = corner.corner(
            samples_list[idx],
            fig=fig,
            weights=get_normalisation_weight(len(samples_list[idx]), max_len),
            labels=param_labels,

            color=colors[idx],
            **CORNER_KWARGS
        )
    plt.legend(
        handles=[
            mlines.Line2D([], [], color=colors[i], label=sample_labels[i])
            for i in range(n)
        ],
        fontsize=20, frameon=False,
        bbox_to_anchor=(1, ndim), loc="upper right"
    )
    plt.savefig(os.path.join(output_dir, filename+'_corner.png'))
    plt.close()


def get_normalisation_weight(len_current_samples, len_of_longest_samples):
    return np.ones(len_current_samples) * (len_of_longest_samples / len_current_samples)

