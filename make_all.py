#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 16:48:45 2018

@author: fkoehlin
"""

import os
import sys
import numpy as np
from write_parameter_table import write_table
from post_process_chain import post_process_chain_1cosmo, post_process_chain_2cosmos
from plot_parameter_triangles import plot_triangle_1cosmo, plot_figures_2cosmos
from plot_1D_histograms import plot_histogram

if __name__ == '__main__':
    
    # will be made an option soon...
    chain_is = '1cosmo'
    
    # TODO: use a proper parser?!
    path_to_chain = sys.argv[1]
    model_name = sys.argv[2]
    
    path_out = path_to_chain
    path_plot = os.path.join(path_to_chain, 'plots/')
    
    if not os.path.isdir(path_plot):
        os.makedirs(path_plot)
        
    ### CUSTOMIZE ###
    # define some kwargs here:
    hist_kwargs = {'histtype': 'step',
                   'density': True,
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
    
    # define key-parameters used in ...
    # must be consistent with CLASS and MontePython names!
    # if list is empty all parameters in chain will be used in plots!
    key_parameters = []
    #key_parameters = ['Omega_m', 'sigma8', 'S8']
    
    
    ### CALLING THE SCRIPTS ###
    print '### Writing out parameter table. ###'
    write_table(path_to_chain)
    
    if chain_is in ['2cosmos', '2cosmo', '2COSMOS', '2COSMO', 'two_cosmos', 'two_cosmo']:
        post_process_chain_2cosmos(path_to_chain, model_name)
        plot_figures_2cosmos()
    else:
        print '### Converting chain to human-readable text and FITS. ###'
        post_process_chain_1cosmo(path_to_chain, model_name)
        print '### Creating triangle 2D parameter plot. ###'
        plot_triangle_1cosmo(path_to_chain, path_plot, fname_suffix=model_name, levels=levels, key_params=key_parameters, hist_kwargs=hist_kwargs, contour_kwargs=contour_kwargs, legend_kwargs=legend_kwargs)
        print '### Creating 1D parameter histograms. ###'
        plot_histogram(path_to_chain, path_plot, hist_kwargs=hist_kwargs)