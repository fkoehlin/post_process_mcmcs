#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
.. module:: write_parameter_table
    :synopsis: write out a human-readable parameter table file
.. moduleauthor:: Fabian Koehlinger <fabian.koehlinger@ipmu.jp>

Script for writing out a human-readable parameter table file, 'parameter_table.txt', in the base folder
of the specified MontePython chain. This file contains useful info like posterior means, medians and
credibility intervals.

Important: You must have translated your run into a default MontePython chain via

    python /path/to/montepython_public/montepython/MontePython.py info /path/to/your/MontePython/chain/{PC, NS, CH}

This script is self-consistent and can be called like this:

    python write_parameter_table /path/to/MontePython/chain model_name={'arbitrary string'} sampler={'MH', 'NS', 'MN', 'PC', 'CH'}

    various other (mostly plotting-related) options can be set further below!
"""

import os
import sys
import glob
import numpy as np
from utils import minimum_credible_intervals, weighted_mean, quantile

def get_values_and_intervals(parameters, weights, labels, use_median=False):

    param_values = np.zeros((len(parameters), 7))
    confidence_values = np.zeros((len(parameters), 6))

    for idx, param in enumerate(parameters):

        print '-> Calculating histogram for {:}.'.format(labels[idx])

        if use_median:
            central_value = quantile(param, [0.5], weights=weights)[0]
        else:
            central_value = weighted_mean(param, weights=weights)

        # bounds returns [[-1sigma, +1sigma],[-2sigma, +2sigma], [-3sigma, +3sigma]]
        bounds = minimum_credible_intervals(param, central_value, weights, bins=50)

        param_values[idx, :] = np.concatenate(([central_value], bounds[:,0], bounds[:,1]))
        confidence_values[idx, :] = central_value + bounds.flatten()

    return param_values, confidence_values

def write_best_fit_to_file(fname, best_fit_params, fit_statistics, labels):
    """
    Store the bestfit parameters to a file
    """
    with open(fname, 'w') as bestfit_file:
        bestfit_file.write('# minimal \chi^2 = {:}, index in chain = {:.0f} \n'.format(fit_statistics[0], fit_statistics[3]))
        bestfit_file.write(
            '# %s\n' % ', '.join(['%16s' % label for label in labels]))
        # Removing scale factors in order to store true parameter values
        for idx in xrange(len(labels)):
            #bfvalue = chain[a[0], 2+i]*info.scales[i, i]
            bf_value = best_fit_params[idx]
            if bf_value > 0:
                bestfit_file.write(' %.6e\t' % bf_value)
            else:
                bestfit_file.write('%.6e\t' % bf_value)
        bestfit_file.write('\n')

    print 'File saved to: \n', fname

    return

def write_parameters_to_file(fname, best_fit_params, fit_statistics, param_values_mean, confidence_values_mean, param_values_median, confidence_values_median, labels, labels_tex):

    with open(fname, 'w') as f:
        f.write('# Best fitting values: \n')
        f.write('\chi^2 = {:.4f}, \chi^2_red = {:.4f} ({:} d.o.f.), index in chain = {:.0f} \n'.format(fit_statistics[0], fit_statistics[1], int(fit_statistics[2]), fit_statistics[3]))
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ')+'{:.4f} \n'.format(best_fit_params[index]))
        ### (weighted) MEAN ###
        f.write('\n'+'# parameter, MEAN, err_minus (68%), err_plus (68%), MEAN, err_minus (95%), err_plus (95%), MEAN, err_minus (99%), err_plus (99%) \n')
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ') + '{0:.4f} {1:.4f} +{2:.4f}, {0:.4f} {3:.4f} +{4:.4f}, {0:.4f} {5:.4f} +{6:.4f} \n'.format(param_values_mean[index, 0], param_values_mean[index, 1], param_values_mean[index, 4], param_values_mean[index, 2], param_values_mean[index, 5], param_values_mean[index, 3], param_values_mean[index, 6]))
        f.write('\n'+'# parameter, lower bound (68%), upper bound (68%), lower bound (95%), upper bound (95%), lower bound (99%), upper bound (99%) \n')
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ')+'1sigma >{:.4f}, 1sigma <{:.4f}, 2sigma >{:.4f}, 2sigma <{:.4f}, 3sigma >{:.4f}, 3sigma <{:.4f} \n'.format(confidence_values_mean[index, 0], confidence_values_mean[index, 1], confidence_values_mean[index, 2], confidence_values_mean[index, 3], confidence_values_mean[index, 4], confidence_values_mean[index, 5]))
        ### (weighted) MEDIAN ###
        f.write('\n'+'# parameter, MEDIAN, err_minus (68%), err_plus (68%), MEDIAN, err_minus (95%), err_plus (95%), MEDIAN, err_minus (99%), err_plus (99%) \n')
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ') + '{0:.4f} {1:.4f} +{2:.4f}, {0:.4f} {3:.4f} +{4:.4f}, {0:.4f} {5:.4f} +{6:.4f} \n'.format(param_values_median[index, 0], param_values_median[index, 1], param_values_median[index, 4], param_values_median[index, 2], param_values_median[index, 5], param_values_median[index, 3], param_values_median[index, 6]))
        f.write('\n'+'# parameter, lower bound (68%), upper bound (68%), lower bound (95%), upper bound (95%), lower bound (99%), upper bound (99%) \n')
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ')+'1sigma >{:.4f}, 1sigma <{:.4f}, 2sigma >{:.4f}, 2sigma <{:.4f}, 3sigma >{:.4f}, 3sigma <{:.4f} \n'.format(confidence_values_median[index, 0], confidence_values_median[index, 1], confidence_values_median[index, 2], confidence_values_median[index, 3], confidence_values_median[index, 4], confidence_values_median[index, 5]))
        ### (weighted) MEAN (TeX) ###
        f.write('\n'+'\n'+'\n'+'### TeX ###'+'\n'+'# parameter, MEAN, err_minus (68%), err_plus (68%), MEAN, err_minus (95%), err_plus (95%), MEAN, err_minus (99%), err_plus (99%) \n')
        for index, label in enumerate(labels_tex):
            name = label +':'
            f.write(name.ljust(20, ' ')+'{0:.2f}_{{{1:.2f}}}^{{+{2:.2f}}}, {0:.2f}_{{{3:.2f}}}^{{+{4:.2f}}}, {0:.2f}_{{{5:.2f}}}^{{+{6:.2f}}} \n'.format(param_values_mean[index, 0], param_values_mean[index, 1], param_values_mean[index, 4], param_values_mean[index, 2], param_values_mean[index, 5], param_values_mean[index, 3], param_values_mean[index, 6]))
        ### (weighted) MEDIAN (TeX) ###
        f.write('\n'+'\n'+'\n'+'### TeX ###'+'\n'+'# parameter, MEDIAN, err_minus (68%), err_plus (68%), MEDIAN, err_minus (95%), err_plus (95%), MEDIAN, err_minus (99%), err_plus (99%) \n')
        for index, label in enumerate(labels_tex):
            name = label +':'
            f.write(name.ljust(20, ' ')+'{0:.2f}_{{{1:.2f}}}^{{+{2:.2f}}}, {0:.2f}_{{{3:.2f}}}^{{+{4:.2f}}}, {0:.2f}_{{{5:.2f}}}^{{+{6:.2f}}} \n'.format(param_values_median[index, 0], param_values_median[index, 1], param_values_median[index, 4], param_values_median[index, 2], param_values_median[index, 5], param_values_median[index, 3], param_values_median[index, 6]))

    print 'File saved to: \n', fname

    return

def write_table(path_to_chain, model_name='bla', sampler='NS', threshold=0.3):

    if sampler == 'NS' or sampler == 'MN':
        fnames = [os.path.join(path_to_chain, 'chain_NS__accepted.txt')]
    elif sampler == 'MH':
        fnames = glob.glob(path_to_chain + '*.txt')
    elif sampler == 'CH':
        fnames = [os.path.join(path_to_chain, 'chain_CH__sampling.txt')]
    elif sampler == 'PC':
        fnames = [os.path.join(path_to_chain, 'chain_PC__accepted.txt')]
    else:
        print 'You must supply the type of sampler used for the MCMC (MH = Metropolis Hastings, MN = MultiNest, CH = CosmoHammer, PC = PolyChord).'

    # deal with multiple chains from MH run and combine them into one (also taking care of burn-in)
    counter = 0
    for fname in fnames:
        if fname not in glob.glob(path_to_chain + '*HEADER.txt') and fname != os.path.join(path_to_chain, 'parameter_table.txt'):
            data_tmp = np.loadtxt(fname)
            len_chain = data_tmp.shape[0]
            idx_gtr_threshold = int(threshold * len_chain)
            # remove first 30% of entries as burn-in from MH chain:
            # not necessary for NS and CH(?)!
            if sampler == 'MH':
                data_tmp = data_tmp[idx_gtr_threshold:, :]
            if counter == 0:
                data = data_tmp
            else:
                data = np.concatenate((data, data_tmp))
            counter += 1

    '''
    # remove first 30% of entries as burn-in from MH chain:
    # not necessary for NS and CH(?)!
    if sampler == 'MH':
        len_chain = data.shape[0]
        idx_gtr_threshold = int(threshold * len_chain)
        data = data[idx_gtr_threshold:, :]
    '''

    weights = data[:, 0]
    #print data.shape
    #print data[:, -1]
    # glob can expand names with *-operator!
    fname =  glob.glob(path_to_chain + '*_.paramnames')[0]
    #print fname
    names = np.loadtxt(fname, dtype=str, delimiter='\t')
    new_names = names.tolist() #[:-1, :] = names[:, :]

    #print np.shape(names)
    #print data.shape

    added_params = 0
    if 'Omega_m ' in names[:, 0] and 'sigma8 ' in names[:, 0]:
        idx_Om = np.where('Omega_m ' == names[:, 0])[0]
        idx_s8 = np.where('sigma8 ' == names[:, 0])[0]
        # +2 because of weights and mloglkl:
        S8 = data[:, idx_s8 + 2] * np.sqrt(data[:, idx_Om + 2] / 0.3)
        #print S8.mean()
        data = np.column_stack((data, S8))
        new_names.append(['S8', 'S_{8}'])

        added_params += 1

    elif 'Omega_m' in names[:, 0] and 'sigma8' in names[:, 0]:
        idx_Om = np.where('Omega_m' == names[:, 0])[0]
        idx_s8 = np.where('sigma8' == names[:, 0])[0]
        # +2 because of weights and mloglkl:
        S8 = data[:, idx_s8 + 2] * np.sqrt(data[:, idx_Om + 2] / 0.3)
        #print S8.mean()
        data = np.column_stack((data, S8))
        new_names.append(['S8', 'S_{8}'])

        added_params += 1
    #exit()

    for idx in xrange(2):

        if 'Omega_m_{:} '.format(idx + 1) in names[:, 0] and 'sigma8_{:} '.format(idx + 1) in names[:, 0]:
            idx_Om = np.where('Omega_m_{:} '.format(idx + 1) == names[:, 0])[0]
            idx_s8 = np.where('sigma8_{:} '.format(idx + 1) == names[:, 0])[0]
            # +2 because of weights and mloglkl:
            S8 = data[:, idx_s8 + 2] * np.sqrt(data[:, idx_Om + 2] / 0.3)
            #print S8.mean()
            data = np.column_stack((data, S8))
            new_names.append(['S8_{:}'.format(idx + 1), 'S_{{8, \, {:}}}'.format(idx + 1)])

            added_params += 1

        elif 'Omega_m_{:}'.format(idx + 1) in names[:, 0] and 'sigma8_{:}'.format(idx + 1) in names[:, 0]:
            idx_Om = np.where('Omega_m_{:}'.format(idx + 1) == names[:, 0])[0]
            idx_s8 = np.where('sigma8_{:}'.format(idx + 1) == names[:, 0])[0]
            # +2 because of weights and mloglkl:
            S8 = data[:, idx_s8 + 2] * np.sqrt(data[:, idx_Om + 2] / 0.3)
            #print S8.mean()
            data = np.column_stack((data, S8))
            new_names.append(['S8_{:}'.format(idx + 1), 'S_{{8, \, {:}}}'.format(idx + 1)])

            added_params += 1

    new_names = np.asarray(new_names, dtype=str)
    labels = new_names[:, 0]
    labels_tex = new_names[:, 1]

    for idx, label in enumerate(labels):
        if label[-1] == ' ':
            labels[idx] = label[:-1]
    #print new_names, new_names.shape

    #column_names = np.concatenate((np.asarray(['weights', 'mloglkl']), names[:, 0], np.asarray(['S8'])))
    chi2 = 2. * data[:, 1]
    min_chi2 = chi2.min()
    best_fit_index = np.where(data[:, 1] == data[:, 1].min())[0]
    #print best_fit_index
    #exit()
    best_fit_params = data[best_fit_index]
    #print data.shape
    #print best_fit_params, best_fit_params.shape
    #exit()
    fit_statistics = np.array([min_chi2, 0., 0., int(best_fit_index[0])])
    print 'Calculating histograms with central value = MEAN.'
    params_mean, conf_mean = get_values_and_intervals(data[:, 2:].T, weights, labels, use_median=False)
    print 'Calculating histograms with central value = MEDIAN.'
    params_median, conf_median = get_values_and_intervals(data[:, 2:].T, weights, labels, use_median=True)

    fname = os.path.join(path_to_chain, 'parameter_table.txt')
    write_parameters_to_file(fname, best_fit_params[0, 2:], fit_statistics, params_mean, conf_mean, params_median, conf_median, labels, labels_tex)

    fname = os.path.join(path_to_chain, model_name + '.bestfit')
    # remove S8 again so that this bestfit file can be used as MP's bestfit-file!!!
    write_best_fit_to_file(fname, best_fit_params[0, 2:best_fit_params.size - added_params], fit_statistics, labels[:len(labels) - added_params])

    return

if __name__ == '__main__':

    path_to_chain = sys.argv[1]

    model_name = sys.argv[2]

    # needs to be closed with '/' for glob.glob to work properly!
    if path_to_chain[-1] != '/':
        path_to_chain += '/'

    sampler = sys.argv[3]

    write_table(path_to_chain, model_name=model_name, sampler=sampler)
