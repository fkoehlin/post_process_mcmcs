#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
.. module:: plot_1D_histograms
    :synopsis: plot 1D posterior histograms from an MCMC
.. moduleauthor:: Fabian Koehlinger <fabian.koehlinger@ipmu.jp>

Script for plotting all 1D posterior histograms (credibility intervals) from
a MontePython chain (from every sampler type).

Important: You must have translated your run into a default MontePython chain via

    python /path/to/montepython_public/montepython/MontePython.py info /path/to/your/MontePython/chain/{PC, NS, CH}

This script is self-consistent and can be called like this:

    python plot_1D_histograms /path/to/your/MontePython/chain/

    various other (mostly plotting-related) options can be set further below!
"""

import os
#import numpy as np
import math
import matplotlib.pyplot as plt
from plot_parameter_triangles import get_params_of_interest

def plot_histogram(path_in, path_out, key_params=[], plot_filetypes=['.pdf'], hist_kwargs={}):

    weights, points_cosmo, param_names, labels_chain = get_params_of_interest(path_in, key_params=key_params)

    # Find the appropriate number of columns and lines for the 1d posterior
    # plot
    num_columns = int(round(math.sqrt(len(param_names))))
    num_lines = int(math.ceil(len(param_names) * 1.0 / num_columns))

    # my solution:
    #num_columns = len(param_names) / 2
    #num_lines = 2 + len(param_names) % 2

    fig = plt.figure(figsize=(3 * num_columns, 3 * num_lines), tight_layout=True)

    for idx, param in enumerate(param_names):

        ax = fig.add_subplot(num_lines, num_columns, idx + 1, yticks=[])

        ax.set_title(r'$' + labels_chain[idx] + r'$')
        ax.hist(points_cosmo[:, idx], weights=weights, **hist_kwargs) #'ls' = hist_kwargs['ls'], 'normed' = hist_kwargs['normed'])

    for filetype in plot_filetypes:
        fname = os.path.join(path_out, 'histograms_1D' + filetype)
        fig.savefig(fname)
        print( 'Plot saved to: \n', fname)

    return

if __name__ == '__main__':

    import sys

    hist_kwargs = {'histtype': 'step',
                   'density': True,
                   'color': 'black',
                   'label': r'$\mathrm{fiducial}$',
                   'ls': '-',
                   'bins': 30
                  }


    path_in = sys.argv[1]

    # needs to be closed with '/' for glob.glob to work properly!
    if path_in[-1] != '/':
        path_in += '/'

    path_out = os.path.join(path_in, 'plots/')

    if not os.path.isdir(path_out):
        os.makedirs(path_out)

    plot_histogram(path_in, path_out, hist_kwargs=hist_kwargs)
