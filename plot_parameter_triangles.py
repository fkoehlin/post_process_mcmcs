#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
.. module:: plot_parameter_triangles
    :synopsis: plot all possible 2D posterior combinations from chain
.. moduleauthor:: Fabian Koehlinger <fabian.koehlinger@ipmu.jp>

Script for plotting all possible 2D posterior combinations from  a MontePython
chain (from every sampler type).

Important: You must have translated your run into a default MontePython chain via

    python /path/to/montepython_public/montepython/MontePython.py info /path/to/your/MontePython/chain/{PC, NS, CH}

This script is self-consistent and can be called like this:

    python plot_parameter_triangles.py /path/to/MontePython/chain/ fname_suffix={'arbitrary string'} chain_is={'2c', '2cosmo', '1c', '1cosmo'}

    various other (mostly plotting-related) options can be set further below!
"""

import os
import glob
import numpy as np
# for avoiding type 3 fonts:
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import corner

def get_params_of_interest(path_to_chain, key_params=[]):

    #fname = path_to_chain + 'chain_NS__accepted.txt'
    fname_chain = glob.glob(path_to_chain + '*_HEADER.txt')[0]
    fname_names = glob.glob(path_to_chain + '*_HEADER.paramnames')[0]

    data = np.loadtxt(fname_chain)
    weights = data[:, 0]
    mloglkl = data[:, 1]

    vals = data[:, 2:]

    names_chain = np.loadtxt(fname_names, dtype=str, delimiter='\t')
    all_labels = names_chain[:, 1]
    # drop column with LaTeX names and derived parameters
    param_names = names_chain[:, 0]
    # remove spaces in strings:
    for idx, param in enumerate(param_names):
        if param[-1] == ' ':
            param_names[idx] = param[:-1]

    chain = dict(zip(param_names, vals.T))
    names = dict(zip(param_names, all_labels))

    if len(key_params) == 0:
        points_cosmo = []
        for param in param_names:
            points_cosmo += [chain[param]]
        points_cosmo = np.asarray(points_cosmo)
        #print param_names
        labels_chain = all_labels
    else:
        points_cosmo = []
        labels_chain = []
        for param in key_params:
            points_cosmo += [chain[param]]
            labels_chain += [names[param]]
        points_cosmo = np.asarray(points_cosmo)
        param_names = key_params

    #print param_names
    #print labels_chain
    #print points_cosmo
    #exit()

    return weights, points_cosmo.T, np.asarray(param_names), labels_chain

def get_params_of_interest_2cosmos(path_to_chain, key_params=[]):

    #fname = path_to_chain + 'chain_NS__accepted.txt'
    fname_chain = glob.glob(path_to_chain + '*_HEADER.txt')[0]
    fname_names = glob.glob(path_to_chain + '*_HEADER.paramnames')[0]

    data = np.loadtxt(fname_chain)
    weights = data[:, 0]
    mloglkl = data[:, 1]

    vals = data[:, 2:]

    names_chain = np.loadtxt(fname_names, dtype=str, delimiter='\t')
    all_labels = names_chain[:, 1]
    # drop column with LaTeX names and derived parameters
    param_names_tmp = names_chain[:, 0]
    # remove spaces in strings:
    for idx, param in enumerate(param_names_tmp):
        if param[-1] == ' ':
            param_names_tmp[idx] = param[:-1]

    chain = dict(zip(param_names_tmp, vals.T))
    names = dict(zip(param_names_tmp, all_labels))

    if len(key_params) == 0:
        points_cosmo1 = []
        points_cosmo2 = []
        param_names = []
        labels_chain = []
        for idx, param in enumerate(param_names_tmp):
            if param[-2:] == '_1':
                points_cosmo1 += [chain[param]]
                # TODO: here's some REGEXP magic needed to remove the cosmo indices from the TeX labels...
                labels_chain += [all_labels[idx]]
                param_names += [param[:-2]]
            else:
                points_cosmo2 += [chain[param]]
        points_cosmo1 = np.asarray(points_cosmo1)
        points_cosmo2 = np.asarray(points_cosmo2)

    else:
        points_cosmo1 = []
        points_cosmo2 = []
        labels_chain = []
        for param in key_params:
            points_cosmo1 += [chain[param + '_1']]
            points_cosmo2 += [chain[param + '_2']]
            labels_chain += [names[param + '_1'][:-2]]
        points_cosmo1 = np.asarray(points_cosmo1)
        points_cosmo2 = np.asarray(points_cosmo2)
        param_names = key_params

    #print param_names
    #print labels_chain
    #print points_cosmo
    #exit()

    return weights, points_cosmo1.T, points_cosmo2.T, np.asarray(param_names), labels_chain

def plot_triangle_1cosmo(path_in, path_out, fname_suffix='bla', levels=np.array([68.27, 95.45, 99.73]) / 100., key_params=[], hist_kwargs={}, contour_kwargs={}, legend_kwargs={}, label_kwargs={}, plot_filetypes=['.pdf'], smooth=0.5, tick_labelsize=12):

    if len(key_params) == 0:
        fname_out = path_out + fname_suffix + '_all_params'
    else:
        fname_out = path_out + fname_suffix + '_key_params'

    weights, points_cosmo, param_names, labels_TeX = get_params_of_interest(path_in, key_params=key_params)
    #= chain1.values(), chain1.keys()

    #print points_cosmo.shape
    #print points_cosmo[:, 0].min(), points_cosmo[:, 0].max()
    #print param_names[0]
    #exit()
    # exact prior ranges (except for S8)
    plot_ranges = []
    labels = []
    for idx in xrange(len(param_names)):
        plot_ranges += [(points_cosmo[:, idx].min(), points_cosmo[:, idx].max())]
        labels += [r'$' + labels_TeX[idx] + r'$']

    '''
    # adjust prior ranges manually (for KV450):
    if 'omega_cdm' in param_names:
        idx_omega_cdm = int(np.where(param_names == 'omega_cdm')[0])
        plot_ranges[idx_omega_cdm] = [0.01, 0.99]
    if 'ln10^{10}A_s' in param_names:
        idx_lnAs = int(np.where(param_names == 'ln10^{10}A_s')[0])
        plot_ranges[idx_lnAs] = [1.7, 5.]
    if 'omega_b' in param_names:
        idx_omega_b = int(np.where(param_names == 'omega_b')[0])
        plot_ranges[idx_omega_b] = [0.01875, 0.02625]
    if 'n_s' in param_names:
        idx_ns = int(np.where(param_names == 'n_s')[0])
        plot_ranges[idx_ns] = [0.7, 1.3]
    if 'h' in param_names:
        idx_h = int(np.where(param_names == 'h')[0])
        plot_ranges[idx_h] = [0.64, 0.82]
    if 'A_IA' in param_names:
        idx_A_IA = int(np.where(param_names == 'A_IA')[0])
        plot_ranges[idx_A_IA] = [-6.0, 6.0]
    if 'c_min' in param_names:
        idx_c_min = int(np.where(param_names == 'c_min')[0])
        plot_ranges[idx_c_min] = [2., 3.13]
    if 'dc' in param_names:
        idx_dc = int(np.where(param_names == 'dc')[0])
        plot_ranges[idx_dc] = [-0.0006, 0.0006]
    if 'Ac' in param_names:
        idx_Ac = int(np.where(param_names == 'Ac')[0])
        plot_ranges[idx_Ac] = [0.62, 1.40]
    if 'D_z1' in param_names:
        idx_Dz1 = int(np.where(param_names == 'D_z1')[0])
        plot_ranges[idx_Dz1] = [-0.117, 0.117]
    if 'D_z2' in param_names:
        idx_Dz2 = int(np.where(param_names == 'D_z2')[0])
        plot_ranges[idx_Dz2] = [-0.069, 0.069]
    if 'D_z3' in param_names:
        idx_Dz3 = int(np.where(param_names == 'D_z3')[0])
        plot_ranges[idx_Dz3] = [-0.078, 0.078]
    if 'D_z4' in param_names:
        idx_Dz4 = int(np.where(param_names == 'D_z4')[0])
        plot_ranges[idx_Dz4] = [-0.036, 0.036]
    if 'D_z5' in param_names:
        idx_Dz5 = int(np.where(param_names == 'D_z5')[0])
        plot_ranges[idx_Dz5] = [-0.033, 0.033]

    '''
    #'''
    # adjust ranges manually for derived parameters: Omega_m, sigma8 and S8
    if 'Omega_m' in param_names:
        idx_Omega_m = int(np.where(param_names == 'Omega_m')[0])
        plot_ranges[idx_Omega_m] = [0.05, 0.55]

    if 'sigma8' in param_names:
        idx_sigma8 = int(np.where(param_names == 'sigma8')[0])
        plot_ranges[idx_sigma8] = [0.4, 1.3]

    if 'S8' in param_names:
        idx_S8 = int(np.where(param_names == 'S8')[0])
        plot_ranges[idx_S8] = [0.55, 0.90]
    #'''

    fig = corner.corner(points_cosmo, weights=weights, labels=labels, smooth=smooth, range=plot_ranges, plot_contours=True, hist_kwargs=hist_kwargs, levels=levels, plot_datapoints=False, plot_density=False, fill_contours=True, label_kwargs=label_kwargs)
    plt.legend(frameon=False, bbox_transform=plt.gcf().transFigure, **legend_kwargs)

    # for control of labelsize of x,y-ticks:
    for ax in fig.get_axes():
        #ax.tick_params(axis='both', which='major', labelsize=14)
        #ax.tick_params(axis='both', which='minor', labelsize=12)
        ax.tick_params(axis='both', labelsize=tick_labelsize)

    for filetype in plot_filetypes:
        plt.savefig(fname_out + filetype)
        print 'Plot saved to: \n', fname_out + filetype

    return

# TODO: modify!
def plot_triangles_2cosmos(path_in, path_out, fname_suffix='bla', levels=np.array([68.27, 95.45, 99.73]) / 100., key_params=[], hist_kwargs={}, contour_kwargs={}, legend_kwargs={}, label_kwargs={}, plot_filetypes=['.pdf'], smooth=0.5, tick_labelsize=12):

    if len(key_params) == 0:
        fname_out1 = path_out + fname_suffix + '_all_params'
        fname_out2 = path_out + fname_suffix + '_all_params_DIFF'
    else:
        fname_out1 = path_out + fname_suffix + '_key_params'
        fname_out2 = path_out + fname_suffix + '_key_params_DIFF'

    weights, points_cosmo1, points_cosmo2, param_names, labels_TeX = get_params_of_interest_2cosmos(path_in, key_params=key_params)
    points_diff = points_cosmo1 - points_cosmo2

    # set plot ranges between min/max values in chain:
    plot_ranges1 = []
    labels1 = []
    plot_ranges2 = []
    labels2 = []
    for idx in xrange(len(param_names)):
        plot_ranges1 += [(min(points_cosmo1[:, idx].min(), points_cosmo2[:, idx].min()), max(points_cosmo1[:, idx].max(), points_cosmo2[:, idx].max()))]
        labels1 += [r'$' + labels_TeX[idx] + r'$']

        plot_ranges2 += [(points_diff[:, idx].min(), points_diff[:, idx].max())]
        labels2 += [r'$\Delta \ ' + labels_TeX[idx] + r'$']
    #'''
    # adjust ranges manually for derived parameters: Omega_m, sigma8 and S8
    if 'Omega_m' in param_names:
        idx_Omega_m = int(np.where(param_names == 'Omega_m')[0])
        plot_ranges1[idx_Omega_m] = [0.05, 0.55]

    if 'sigma8' in param_names:
        idx_sigma8 = int(np.where(param_names == 'sigma8')[0])
        plot_ranges1[idx_sigma8] = [0.4, 1.3]

    if 'S8' in param_names:
        idx_S8 = int(np.where(param_names == 'S8')[0])
        plot_ranges1[idx_S8] = [0.55, 0.90]
    #'''

    hist_kwargs2 = hist_kwargs.copy()
    hist_kwargs2['color'] = 'blue'
    hist_kwargs2['linestyle'] = '--'

    fig1 = corner.corner(points_cosmo1, weights=weights, labels=labels1, smooth=smooth, range=plot_ranges1, plot_contours=True, hist_kwargs=hist_kwargs, levels=levels, plot_datapoints=False, plot_density=False, fill_contours=True, label_kwargs=label_kwargs)
    corner.corner(points_cosmo2, weights=weights, fig=fig1, labels=labels1, smooth=smooth, range=plot_ranges1, plot_contours=True, hist_kwargs=hist_kwargs2, levels=levels, plot_datapoints=False, plot_density=False, fill_contours=True, label_kwargs=label_kwargs, color='blue')
    fig1.legend(frameon=False, bbox_transform=plt.gcf().transFigure, **legend_kwargs)

    # for control of labelsize of x,y-ticks:
    for ax in fig1.get_axes():
        #ax.tick_params(axis='both', which='major', labelsize=14)
        #ax.tick_params(axis='both', which='minor', labelsize=12)
        ax.tick_params(axis='both', labelsize=tick_labelsize)

    for filetype in plot_filetypes:
        fig1.savefig(fname_out1 + filetype)
        print 'Plot saved to: \n', fname_out1 + filetype

    fig2 = corner.corner(points_diff, weights=weights, truths=np.zeros(len(labels1)), labels=labels2, smooth=smooth, range=plot_ranges2, plot_contours=True, hist_kwargs=hist_kwargs, levels=levels, plot_datapoints=False, plot_density=False, fill_contours=True, label_kwargs=label_kwargs)
    fig2.legend(frameon=False, bbox_transform=plt.gcf().transFigure, **legend_kwargs)

        # for control of labelsize of x,y-ticks:
    for ax in fig2.get_axes():
        #ax.tick_params(axis='both', which='major', labelsize=14)
        #ax.tick_params(axis='both', which='minor', labelsize=12)
        ax.tick_params(axis='both', labelsize=tick_labelsize)

    for filetype in plot_filetypes:
        fig2.savefig(fname_out2 + filetype)
        print 'Plot saved to: \n', fname_out2 + filetype

    return

if __name__ == '__main__':

    import sys

    # define some kwargs here:
    hist_kwargs = {'histtype': 'step',
                   'density': True,
                   'color': 'black',
                   'label': r'$\mathrm{fiducial}$',
                   'ls': '-'
                  }

    contour_kwargs = {'linestyles': '--'}

    legend_kwargs = {'fontsize': 14,
                     'bbox_to_anchor': (0.8, 0.8)
                    }

    # set fontsize of labels for all axes:
    label_kwargs = {'fontsize': 24}

    # set fontsize of x,y-ticks:
    tick_labelsize = 12

    # set confidence levels to plot:
    levels = np.array([68.27, 95.45, 99.73]) / 100.
    levels = levels[:2]

    path_in = sys.argv[1]

    # needs to be closed with '/' for glob.glob to work properly!
    if path_in[-1] != '/':
        path_in += '/'

    path_out = os.path.join(path_in, 'plots/')

    fname_suffix = sys.argv[2]

    chain_is = sys.argv[3]

    if not os.path.isdir(path_out):
        os.makedirs(path_out)

    #key_params = ['Omega_m', 'sigma8', 'S8']
    key_params = []

    if chain_is in ['2c', '2cosmos', '2cosmo', '2COSMOS', '2COSMO', 'two_cosmos', 'two_cosmo']:
        plot_triangles_2cosmos(path_in, path_out, fname_suffix=fname_suffix, levels=levels, key_params=key_params, label_kwargs=label_kwargs, tick_labelsize=tick_labelsize)
    else:
        plot_triangle_1cosmo(path_in, path_out, fname_suffix=fname_suffix, levels=levels, key_params=key_params, label_kwargs=label_kwargs, tick_labelsize=tick_labelsize)

    plt.show()
