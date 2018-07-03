#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 16:45:45 2018

@author: fkoehlin
"""

import os
#import numpy as np
import matplotlib.pyplot as plt
from plot_parameter_triangles import get_params_of_interest

def plot_histogram(path_in, path_out, key_params=[], hist_kwargs={}):
    
    weights, points_cosmo, param_names, labels_chain = get_params_of_interest(path_in, key_params=key_params)
    
    #print weights.shape
    #print points_cosmo.shape
    
    dim_x = len(param_names) / 2 
    dim_y = 2 + len(param_names) % 2
    
    #print dim_x, dim_y, len(param_names)
    #print param_names, len(param_names)
    
    fig, ax = plt.subplots(dim_y, dim_x, tight_layout=True)
    
    idx1 = 0
    idx2 = 0
    for idx, param in enumerate(param_names):
        
        ax[idx1, idx2].set_title(r'$' + labels_chain[idx] + r'$')
        ax[idx1, idx2].hist(points_cosmo[:, idx], weights=weights, **hist_kwargs) #'ls' = hist_kwargs['ls'], 'normed' = hist_kwargs['normed'])
    
        if idx2 == dim_x - 1:
            idx2 = 0
            idx1 += 1
        else:
            idx2 += 1
        
        #print idx1, idx2, idx
    
    # remove unused panels:    
    if len(param_names) % 2 != 0: 
        ax[-1, -1].axis('off')
        ax[-1, -2].axis('off')
        #print idx, idx1, idx2
        
    fname = os.path.join(path_out, 'histograms_1D.pdf')
    fig.savefig(fname)
    print 'Plot saved to: \n', fname
    
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
    path_out = os.path.join(path_in, 'plots/')
    
    if not os.path.isdir(path_out):
        os.makedirs(path_out)
    
    plot_histogram(path_in, path_out, hist_kwargs=hist_kwargs)