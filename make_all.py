#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
.. module:: make_all
    :synopsis: main script
.. moduleauthor:: Fabian Koehlinger <fabian.koehlinger@ipmu.jp>

Script for post-processing a MontePython chain (from every sampler type) fully
automatically and for plotting 1D histograms, 2D contours and 2D parameter
triangle plots.

Important: You must have translated your run into a default MontePython chain via

    python /path/to/montepython_public/montepython/MontePython.py info /path/to/your/MontePython/chain/{PC, NS, CH}

This is the main script, call it like this:

    python make_all.py /path/to/your/MontePython/chain/ model_name={'arbitrary string'} type_sampler={'MH', 'MN', 'NS', 'PC', 'CH"} chain_is={'2c', '2cosmo', '1c', '1cosmo'}

    various other (mostly plotting-related) options can be set further below!
"""

import os
import sys
import numpy as np
from write_parameter_table import write_table
from post_process_chain import post_process_chain_1cosmo, post_process_chain_2cosmos
from plot_parameter_triangles import plot_triangle_1cosmo, plot_triangles_2cosmos
from plot_1D_histograms import plot_histogram
from plot_2D_contours import make_2D_contour

if __name__ == '__main__':

    # TODO: use a proper parser?!
    path_to_chain = sys.argv[1]

    # needs to be closed with '/' for glob.glob to work properly!
    if path_to_chain[-1] != '/':
        path_to_chain += '/'

    model_name = sys.argv[2]

    path_out = path_to_chain
    path_plot = os.path.join(path_to_chain, 'plots/')

    # MH = Metropolis-Hastings, MN = MultiNest, CH = CosmoHammer, PC = PolyChord
    type_sampler = sys.argv[3]
    #print type_sampler

    # if you're dealing with a regular chain, just use any or no symbol
    # if you're dealing with a 2cosmos chain, use any of the following keywords:
    # '2c', '2cosmos', '2cosmo', '2COSMOS', '2COSMO', 'two_cosmos', 'two_cosmo'
    chain_is = sys.argv[4]

    if not os.path.isdir(path_plot):
        os.makedirs(path_plot)

    ### CUSTOMIZE ###

    # if you used Metropolis-Hastings as sampler, the first threshold per cent of the chain will be removed as burn-in:
    threshold = 0.3

    # specify filetype endings for plots:
    plot_filetypes = ['.pdf', '.png']

    # define some kwargs here:
    hist_kwargs = {'histtype': 'step',
                   'normed': True,
                   'color': 'black',
                   #'label': r'$\mathrm{fiducial}$',
                   'ls': '-'
                  }

    contour_kwargs = {'linestyles': '--'}

    legend_kwargs = {'fontsize': 14,
                     'bbox_to_anchor': (0.8, 0.8)
                    }

    levels = np.array([68.27, 95.45, 99.73]) / 100.
    # only use 68% and 95% contours:
    levels = levels[:2]

    # std for Gaussian filtering to smooth contours:
    # set to None for no smoothing!
    smooth = 1.

    # define key-parameters used in ...
    # must be consistent with CLASS and MontePython names!
    # if list is empty all parameters in chain will be used in plots!
    key_parameters = []
    #key_parameters = ['omega_cdm', 'ln10^{10}A_s', 'omega_b', 'n_s', 'h', 'c_min', 'A_IA', 'Omega_m', 'sigma8', 'S8']
    #key_parameters = ['Omega_m', 'sigma8', 'S8']

    # Options for 2D contour plots:
    # supply paths to several chains that are all plotted in one contour plot
    # this assumes that chains at those paths were also post-processed in the same way as here!
    #paths_to_chains = [path_to_chain, path_to_planck, path_to_kids450]
    #labels_ctr = ['KV450', 'Planck', 'KiDS-450']
    #colors_2D = ['blue', 'grey', 'green']
    paths_to_chains = [path_to_chain]
    #labels_2Dctrs = ['KV-450']
    labels_2Dctrs = [model_name]
    colors_2Dctrs = ['skyblue']
    #colors_2Dctrs = ['skyblue', 'limegreen', 'lightgrey']

    # supply min/max tuples per axis for single 2D contour plots:
    # x = Omega_m, y = sigma8
    ranges_2D_contour1 = {'x': (0.10, 0.50), 'y': (0.55, 1.30)}
    # x = Omega_m, y = S8
    ranges_2D_contour2 = {'x': (0.10, 0.50), 'y': (0., 1.)}

    ### CALLING THE SCRIPTS ###
    print '### Writing out parameter table. ###'
    write_table(path_to_chain, model_name=model_name, sampler=type_sampler, threshold=threshold)

    if chain_is in ['2c', '2cosmos', '2cosmo', '2COSMOS', '2COSMO', 'two_cosmos', 'two_cosmo']:
        print '### Converting chain to human-readable text and FITS. ###'
        post_process_chain_2cosmos(path_to_chain, model_name, sampler=type_sampler, threshold=threshold)
        print '### Creating triangle 2D parameter plot. ###'
        plot_triangles_2cosmos(path_to_chain, path_plot, fname_suffix=model_name, levels=levels, key_params=key_parameters, hist_kwargs=hist_kwargs, contour_kwargs=contour_kwargs, legend_kwargs=legend_kwargs, plot_filetypes=plot_filetypes, smooth=smooth)
        print '### Creating 1D parameter histograms. ###'
        plot_histogram(path_to_chain, path_plot, hist_kwargs=hist_kwargs, plot_filetypes=plot_filetypes)
        print '### Creating 2D contours: Omega_m vs. sigma8 ###'
        make_2D_contour(path_plot, paths_in=paths_to_chains, ranges_2D_contour=ranges_2D_contour1, labels_ctr=labels_2Dctrs, colors=colors_2Dctrs, key_params=['Omega_m_1', 'sigma8_1'], fname_suffix=model_name, levels=levels, plot_filetypes=plot_filetypes, smooth=smooth)
        print '### Creating 2D contours: Omega_m vs. S8 ###'
        make_2D_contour(path_plot, paths_in=paths_to_chains, ranges_2D_contour=ranges_2D_contour2,labels_ctr=labels_2Dctrs, colors=colors_2Dctrs, key_params=['Omega_m_1', 'S8_1'], fname_suffix=model_name, levels=levels, plot_filetypes=plot_filetypes, smooth=smooth)
    else:
        print '### Converting chain to human-readable text and FITS. ###'
        post_process_chain_1cosmo(path_to_chain, model_name, sampler=type_sampler, threshold=threshold)
        print '### Creating triangle 2D parameter plot. ###'
        plot_triangle_1cosmo(path_to_chain, path_plot, fname_suffix=model_name, levels=levels, key_params=key_parameters, hist_kwargs=hist_kwargs, contour_kwargs=contour_kwargs, legend_kwargs=legend_kwargs, plot_filetypes=plot_filetypes, smooth=smooth)
        print '### Creating 1D parameter histograms. ###'
        plot_histogram(path_to_chain, path_plot, hist_kwargs=hist_kwargs, plot_filetypes=plot_filetypes)
        print '### Creating 2D contours: Omega_m vs. sigma8 ###'
        make_2D_contour(path_plot, paths_in=paths_to_chains, ranges_2D_contour=ranges_2D_contour1, labels_ctr=labels_2Dctrs, colors=colors_2Dctrs, key_params=['Omega_m', 'sigma8'], fname_suffix=model_name, levels=levels, plot_filetypes=plot_filetypes, smooth=smooth)
        print '### Creating 2D contours: Omega_m vs. S8 ###'
        make_2D_contour(path_plot, paths_in=paths_to_chains, ranges_2D_contour=ranges_2D_contour2,labels_ctr=labels_2Dctrs, colors=colors_2Dctrs, key_params=['Omega_m', 'S8'], fname_suffix=model_name, levels=levels, plot_filetypes=plot_filetypes, smooth=smooth)
