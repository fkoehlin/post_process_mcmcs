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

def plot_triangle_1cosmo(path_in, path_out, fname_suffix='bla', levels=np.array([68.27, 95.45, 99.73]) / 100., key_params=[], hist_kwargs={}, contour_kwargs={}, legend_kwargs={}, plot_filetypes=['.pdf'], smooth=0.5):

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

    # adjust ranges for Omega_m, sigma8 and S8 manually:
    if 'Omega_m' in param_names:
        idx_Omega_m = int(np.where(param_names == 'Omega_m')[0])
        plot_ranges[idx_Omega_m] = [0.15, 0.60]

    if 'sigma8' in param_names:
        idx_sigma8 = int(np.where(param_names == 'sigma8')[0])
        plot_ranges[idx_sigma8] = [0.45, 1.05]

    if 'S8' in param_names:
        idx_S8 = int(np.where(param_names == 'S8')[0])
        plot_ranges[idx_S8] = [0.65, 0.90]

    corner.corner(points_cosmo, weights=weights, labels=labels, smooth=smooth, range=plot_ranges, plot_contours=True, hist_kwargs=hist_kwargs, levels=levels, plot_datapoints=False, plot_density=True)
    plt.legend(frameon=False, bbox_transform=plt.gcf().transFigure, **legend_kwargs)

    for filetype in plot_filetypes:
        plt.savefig(fname_out + filetype)
        print 'Plot saved to: \n', fname_out + filetype

    return

# TODO: modify!
def plot_figures_2cosmos(path_in1, path_in2, path_out, fname_suffix='bla', exclude_zbin=4, nzbins=5, levels=np.array([0.68, 0.95]), hist_kwargs1={}, hist_kwargs2={}, contour_kwargs1={}, contour_kwargs2={}, legend_kwargs={}, plot_filetypes=['.pdf'], smooth=0.5):

    labels1 = ''
    labels2 = ''
    for zbin1 in xrange(nzbins):
        for zbin2 in xrange(zbin1, nzbins):

            if zbin1 == exclude_zbin - 1 or zbin2 == exclude_zbin - 1:
                labels2 += '{:}{:}, '.format(zbin1 + 1, zbin2 + 1)
            else:
                labels1 += '{:}{:}, '.format(zbin1 + 1, zbin2 + 1)

    labels1 = r'$' + labels1[:-2] + r'$'
    labels2 = r'$' + labels2[:-2] + r'$'

    fname_out = path_out + fname_suffix + '_joint_chain'

    fname1 =  glob.glob(path_in1 + '*.fits')[0]
    weights_chain1, points_cosmo1_chain1, points_cosmo2_chain1 = get_params_of_interest(fname1)
    points_diff_chain1 = points_cosmo1_chain1[:, 0:3] - points_cosmo2_chain1[:, 0:3]

    '''
    if plot_noB:
        fname2 =  glob.glob(path_in2 + '*.fits')[0]
        weights_chain2, points_cosmo1_chain2, points_cosmo2_chain2 = get_params_of_interest(fname2)
        points_diff_chain2 = points_cosmo1_chain2[:, 0:3] - points_cosmo2_chain2[:, 0:3]
    '''
    # exact prior ranges (except for S8)
    rangePlot_1_vs_2 = [(-6,6), (0.4, 1.2), (0., 1.), (0.64, 0.82), (0.7, 1.3)]
    hist_kwargs1 = {}
    hist_kwargs1['histtype'] = 'step'
    hist_kwargs1['normed'] = True
    hist_kwargs1['color'] = 'black'
    hist_kwargs1['label'] = r'$z_i \times \, z_j = \{$' + labels1 + r'$\}$'

    hist_kwargs2 = {}
    hist_kwargs2['histtype'] = 'step'
    hist_kwargs2['normed'] = True
    hist_kwargs2['color'] = 'blue'
    hist_kwargs2['label'] = r'$z_i \times \, z_j = \{$' + labels2 + r'$\}$'

    labels = [r'$A_\mathrm{IA}$', r'$S_{8}$', r'$\Omega_{\rm m}$', r'$h$', r'$n_\mathrm{s}$']

    figure_1 = corner.corner(points_cosmo1_chain1, weights=weights_chain1, labels=labels, smooth=0.5, range=rangePlot_1_vs_2, plot_contours=True, hist_kwargs=hist_kwargs1, levels=levels, plot_datapoints=False, plot_density=True)
    figure_2 = corner.corner(points_cosmo2_chain1, weights=weights_chain1, fig=figure_1, smooth=0.5, range=rangePlot_1_vs_2, plot_contours=True, color='blue', hist_kwargs=hist_kwargs2, levels=levels, plot_datapoints=False, plot_density=True, ls='--')
    plt.legend(fontsize=fontsize_legend, frameon=False, bbox_to_anchor=(leg_x, leg_y), bbox_transform=plt.gcf().transFigure)

    for filetype in plot_filetypes:
        plt.savefig(fname_out + filetype)
        print 'Plot saved to: \n', fname_out + filetype

    #Look at differences (This accounts for the covariance!)
    fname_out = path_out + fname_suffix + '_diffs'

    rangePlotDiff = [(-5, 5), (-0.3, 0.3), (-0.5, 0.5)]
    hist_kwargs_diff1 = dict()
    hist_kwargs_diff1['histtype'] = 'step'
    hist_kwargs_diff1['normed'] = True
    hist_kwargs_diff1['color'] = 'black'
    hist_kwargs_diff1['label'] = r'$\mathrm{fiducial \ data}$'

    hist_kwargs_diff2 = dict()
    hist_kwargs_diff2['histtype'] = 'step'
    hist_kwargs_diff2['normed'] = True
    hist_kwargs_diff2['color'] = 'blue'
    hist_kwargs_diff2['label'] = r'$\mathrm{B-modes \ subtracted}$'
    hist_kwargs_diff2['linestyle'] = '--'

    contour_kwargs = {}
    contour_kwargs['linestyles'] = 'dashed'

    labels = [r'$\Delta A_\mathrm{IA}$',r'$\Delta S_8$', r'$\Delta \Omega_{\rm m}$']

    figure_diff1 = corner.corner(points_diff_chain1, weights=weights_chain1, labels=labels, smooth=smooth, range = rangePlotDiff, truths=[0, 0, 0], hist_kwargs=hist_kwargs_diff1, plot_contours=True, levels=levels, plot_datapoints=False, plot_density=False, color='black')
    '''
    if plot_noB:
        figure_diff2 = corner.corner(points_diff_chain2, weights=weights_chain2, fig=figure_diff1, labels=labels, smooth=0.5, range = rangePlotDiff, hist_kwargs=hist_kwargs_diff2, plot_contours=True, levels=levels, plot_datapoints=False, plot_density=False, color='blue', contour_kwargs=contour_kwargs)
        plt.legend(fontsize=fontsize_legend, frameon=False, bbox_to_anchor=(leg_x, leg_y), bbox_transform=plt.gcf().transFigure)
    '''

    for filetype in plot_filetypes:
        plt.savefig(fname_out + filetype)
        print 'Plot saved to: \n', fname_out + filetype

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

    levels = np.array([68.27, 95.45, 99.73]) / 100.
    levels = levels[:2]

    path_in = sys.argv[1]

    # needs to be closed with '/' for glob.glob to work properly!
    if path_in[-1] != '/':
        path_in += '/'

    path_out = os.path.join(path_in, 'plots/')

    fname_suffix = sys.argv[2]

    if not os.path.isdir(path_out):
        os.makedirs(path_out)

    #key_params = ['Omega_m', 'sigma8', 'S8']
    key_params = []

    plot_triangle_1cosmo(path_in, path_out, fname_suffix=fname_suffix, levels=levels, key_params=key_params)

    plt.show()