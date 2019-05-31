#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
.. module:: plot_2D_contours
    :synopsis: plot 2D posterior contours from a chain
.. moduleauthor:: Fabian Koehlinger <fabian.koehlinger@ipmu.jp>

Script for plotting adjustable 2D posterior contours (credibility intervals) from
a MontePython chain (from every sampler type).

Important: You must have translated your run into a default MontePython chain via

    python /path/to/montepython_public/montepython/MontePython.py info /path/to/your/MontePython/chain/{PC, NS, CH}

This script is self-consistent and can be called like this:

    python plot_2D_contours /path/to/your/MontePython/chain/

    various other (mostly plotting-related) options can be set further below!
"""

import os
import numpy as np
# for avoiding type 3 fonts:
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker
#from matplotlib import cm
#import mytriangle as triangle
#from scipy.ndimage import gaussian_filter
from corner import hist2d
#from utils import weighted_mean, minimum_credible_intervals
from plot_parameter_triangles import get_params_of_interest

def sigma(omega, b, alpha):

    #func = b*(omega/0.27)**(-alpha)
    # this is equivalent and avoids a negative exponent
    func = b*(0.3/omega)**(alpha)

    return func

def func2(alpha, omega, sigma8):

    func = sigma8*(omega/0.3)**alpha

    return func

def logsigma(omega, alpha, const):

    norm = 0.3
    #func = alpha*(np.log(norm)-np.log(omega))+const
    # the norm is arbitrary, so I can just use:
    func = -alpha*np.log(omega)+const

    return func

def make_2D_contour(path_out, paths_in=[], labels_ctr=[], colors=[], key_params=['Omega_m', 'sigma8'], ranges_2D_contour={}, fname_suffix='bla', levels=np.array([68.27, 95.45, 99.73]) / 100., plot_filetypes=['.pdf'], smooth=1., nbins=50):

    fig = plt.figure(tight_layout=True)
    ax = fig.add_subplot(111)

    min_x = np.zeros(len(paths_in))
    min_y = np.zeros_like(min_x)
    max_x = np.zeros_like(min_x)
    max_y = np.zeros_like(min_x)
    for idx, path_in in enumerate(paths_in):

        weights, points_cosmo, param_names, labels_TeX = get_params_of_interest(path_in, key_params=key_params)
        #print points_cosmo, points_cosmo.shape
        hist2d(points_cosmo[:, 0], points_cosmo[:, 1], bins=nbins, ax=ax, weights=weights, smooth=smooth, linewidths=1.5, color=colors[idx],
               plot_contours=True, fill_contours=True, plot_density=False, plot_datapoints=False, levels=levels, alpha=0.5) # labels_contour=labels_ctr[idx])

        min_x[idx] = points_cosmo[:, 0].min()
        max_x[idx] = points_cosmo[:, 0].max()
        min_y[idx] = points_cosmo[:, 1].min()
        max_y[idx] = points_cosmo[:, 1].max()

        # fake point for legends:
        ax.scatter(-1., -1., label=labels_ctr[idx])

    ax.set_xlabel(r'$' + labels_TeX[0] + '$') #, fontsize=fontsize_labels)
    ax.set_ylabel(r'$' + labels_TeX[1] + '$') #, fontsize=fontsize_labels)

    if len(ranges_2D_contour) != 0:
        min_x, max_x = ranges_2D_contour['x']
        min_y, max_y = ranges_2D_contour['y']
    else:
        min_x, max_x = min_x.min(), max_x.max()
        min_y, max_y = min_x.min(), max_y.max()

    ax.set_xlim([min_x, max_x])
    ax.set_ylim([min_y, max_y])

    ax.tick_params(axis='both', width=1) #labelsize=ticksize, width=1)
    #ax.legend(loc='upper right', scatterpoints=1, numpoints=1, frameon=False) #, fontsize=fontsize_legend)

    # some advanced legend:
    leg = ax.legend(loc='upper right', frameon=False) #, fontsize=fontsize_legend)
    # no lines in legend:
    for item in leg.legendHandles:
        item.set_visible(False)
    texts = [text for text in leg.get_texts()]
    for index, text in enumerate(texts):
        text.set_color(colors[index])

    for filetype in plot_filetypes:
        fname = os.path.join(path_out, fname_suffix + '_' + key_params[0] + '_vs_' + key_params[1] + filetype)
        fig.savefig(fname)
        print 'Plot saved to: \n', fname

    return

def make_color_2D_contour(chains=[], smooth=''):

    return

if __name__ == '__main__':

    import sys

    path_in = sys.argv[1]
    fname_suffix = sys.argv[2]

    path_out = os.path.join(path_in, 'plots/')

    if not os.path.isdir(path_out):
        os.makedirs(path_out)

    # supply min/max tuples per axis for single 2D contour plots:
    # x = Omega_m, y = sigma8
    ranges_2D_contour = {'x': (0.15, 0.60), 'y': (0.45, 1.05)}

    make_2D_contour(path_out, paths_in=[path_in], ranges_2D_contour=ranges_2D_contour, labels_ctr=['test'], colors=['blue'], key_params=['Omega_m', 'sigma8'], fname_suffix=fname_suffix, levels=np.array([68.27, 95.45]) / 100., plot_filetypes=['.pdf'], smooth=1., nbins=50)
