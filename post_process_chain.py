#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
.. module:: post_process_chain
    :synopsis: turn the default Monte Python chain file into a (commented) text and FITS file
.. moduleauthor:: Fabian Koehlinger <fabian.koehlinger@ipmu.jp>

Script for translating the default Monte Python chain into a commented text and FITS file
in the base folder of the specified MontePython chain. The translated files will also
contain the parameter S8 = sigma8 * sqrt{Omega_m / 0.3} if sigma8 and Omega_m are among
the original parameters of the chain.

Important: You must have translated your run into a default MontePython chain via

This script is self-consistent and can be called like this:

    python post_process_chain.py /path/to/MontePython/chain/ model_name={'arbitrary string'} sampler={'MH', 'MN', 'NS', 'PC', 'CH'} chain_is={'2c', '2cosmo', '1c', '1cosmo'}
"""

import os
import sys
import glob
import numpy as np
import astropy.io.fits as fits
#import scipy.stats.distributions as dist

def post_process_chain_1cosmo(path_to_chain, model_name, sampler='NS', threshold=0.3):

    if sampler == 'NS' or sampler == 'MN':
        fnames = [os.path.join(path_to_chain, 'chain_NS__accepted.txt')]
    elif sampler == 'MH':
        fnames = glob.glob(path_to_chain + '*.txt')
    elif sampler == 'CH':
        fnames = [os.path.join(path_to_chain, 'chain_CH__sampling.txt')]
    elif sampler == 'PC':
        fnames = [os.path.join(path_to_chain, 'chain_PC__accepted.txt')]
    else:
        print( 'You must supply the type of sampler used for the MCMC (MH = Metropolis Hastings, MN = MultiNest, CH = CosmoHammer, PC = PolyChord).')

    # deal with multiple chains from MH run and combine them into one (also taking care of burn-in)
    counter = 0
    for fname in fnames:
        if fname not in glob.glob(path_to_chain + '*HEADER.txt') and fname != os.path.join(path_to_chain, 'parameter_table.txt'):
            data_tmp = np.loadtxt(fname)
            len_chain = data_tmp.shape[0]
            idx_gtr_threshold = int(threshold * len_chain)
            # remove first 30% as burn-in
            # not necessary for NS and CH(?)!
            if sampler == 'MH':
                data_tmp = data_tmp[idx_gtr_threshold:, :]
            if counter == 0:
                data = data_tmp
            else:
                data = np.concatenate((data, data_tmp))
            counter += 1

    '''
    # remove first 30% as burn-in
    # not necessary for NS and CH(?)!
    if sampler == 'MH':
        len_chain = data.shape[0]
        idx_gtr_threshold = int(threshold * len_chain)
        data = data[idx_gtr_threshold:, :]
    '''
    # glob can expand names with *-operator!
    fname =  glob.glob(path_to_chain + '*_.paramnames')[0]
    names = np.loadtxt(fname, dtype=str, delimiter='\t')

    # remove trailing spaces:
    for idx, name in enumerate(names[:, 0]):
        if name[-1] == ' ':
            names[idx, 0] = name[:-1]

    chain_dict = dict(zip(names[:, 0], data[:, 2:].T))
    new_names = names.tolist() #[:-1, :] = names[:, :]
    try:
        S8 = chain_dict['sigma8'] * np.sqrt(chain_dict['Omega_m'] / 0.3)
        data = np.column_stack((data, S8))
        new_names.append(['S8', 'S_{8}'])
    except:
        print( 'Could not calculate and append S8. \n Are Omega_m and sigma8 in the chain?')
    new_names = np.asarray(new_names, dtype=str)
    # better TeX:
    new_names[np.where(new_names[:, 0] == 'sigma8'), 1] = '\sigma_8'
    # \mathrm doesn't seem to work...
    new_names[np.where(new_names[:, 0] == 'Omega_m'), 1] = '\Omega_{m}'
    new_names[np.where(new_names[:, 0] == 'omega_cdm'), 1] = '\omega_{cdm}'
    new_names[np.where(new_names[:, 0] == 'omega_b'), 1] = '\omega_{b}'
    new_names[np.where(new_names[:, 0] == 'n_s'), 1] = 'n_{s}'
    new_names[np.where(new_names[:, 0] == 'ln10^{10}A_s'), 1] = '\ln 10^{10} A_s'

    column_names = np.concatenate((np.asarray(['weights', 'mloglkl']), names[:, 0], np.asarray(['S8'])))

    # some best-fitting statistics:
    '''
    chi2 = 2. * data[:, 1]
    min_chi2 = chi2.min()
    nparams = len(names) - nderived
    dof = ndata - nparams - 1
    min_chi2_red = min_chi2 / float(dof)

    print 'min(chi^2) = {:.3f}, min(chi^2)_red = {:.3f}, {:} d.o.f.'.format(min_chi2, min_chi2_red, int(dof))

    p_value = dist.chi2.sf(min_chi2, dof)

    print 'p-value = {:.6f}'.format(p_value)
    '''

    header = ''
    for elem in column_names:
        if elem[-1] == ' ':
            elem = elem[:-1]
        header += elem + ', '
    header = header[:-2]
    #print header
    #header =  header.tolist()
    fname = os.path.join(path_to_chain, model_name + '__HEADER.txt')
    np.savetxt(fname, data, header=header)
    print( 'Data saved to: \n', fname)
    if sampler == 'MH':
        print( 'ATTENTION: the above file does not contain the first {:}% of entries from the original chain (burn-in).'.format(threshold * 100))

    # for consistency with getdist:
    fname = os.path.join(path_to_chain, model_name + '__HEADER.paramnames')
    np.savetxt(fname, new_names, delimiter='\t', fmt='%s')

    # write out best-fitting statistics:
    '''
    header = 'min(chi^2), min(chi^2)_red, d.o.f., number of datapoints, number of free parameters, p-value'
    fname = path_to_chain + 'chi2_info.txt'
    savedata = np.column_stack((min_chi2, min_chi2_red, dof, ndata, nparams, p_value))
    np.savetxt(fname, savedata, header=header)
    print 'Data saved to: \n', fname
    '''

    # save chain also as FITS file:
    columns = []
    for index_col, column_name in enumerate(column_names):
        columns.append(fits.Column(name=column_name, format='D', array=data[:, index_col]))

    coldefs = fits.ColDefs(columns)
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    prihdu = fits.PrimaryHDU()
    thdulist = fits.HDUList([prihdu, tbhdu])
    fname = os.path.join(path_to_chain, model_name + '.fits')
    if not os.path.isfile(fname):
        thdulist.writeto(fname)
        print( 'Data saved to: \n', fname)
        if sampler == 'MH':
            print( 'ATTENTION: the above file does not contain the first {:}% of entries from the original chain (burn-in).'.format(threshold * 100))
    else:
        print( 'FITS file already exists. Proceeding without overwriting it!')

    return

def post_process_chain_2cosmos(path_to_chain, model_name, sampler='NS', threshold=0.3):

    if sampler == 'NS' or sampler == 'MN':
        fname = os.path.join(path_to_chain, 'chain_NS__accepted.txt')
    elif sampler == 'MH':
        fname = glob.glob(path_to_chain + '*.txt')[0]
    elif sampler == 'CH':
        fname = os.path.join(path_to_chain, 'chain_CH__sampling.txt')
    elif sampler == 'PC':
        fname = os.path.join(path_to_chain, 'chain_PC__accepted.txt')
    else:
        print( 'You must supply the type of sampler used for the MCMC (MH = Metropolis Hastings, MN = MultiNest, CH = CosmoHammer, PC = PolyChord).')

    data = np.loadtxt(fname)

    # remove first 30% of samples as burn-in
    # not necessary for NS and CH(?)!
    if sampler == 'MH':
        len_chain = data.shape[0]
        idx_gtr_threshold = int(threshold * len_chain)
        data = data[idx_gtr_threshold:, :]

    fname =  glob.glob(path_to_chain + '*_.paramnames')[0]
    names = np.loadtxt(fname, dtype=str, delimiter='\t')
    # remove trailing spaces:
    for idx, name in enumerate(names[:, 0]):
        if name[-1] == ' ':
            names[idx, 0] = name[:-1]

    new_names = names.tolist() #[:-1, :] = names[:, :]

    chain_dict = dict(zip(names[:, 0], data[:, 2:].T))

    try:
        S8_1 = chain_dict['sigma8_1'] * np.sqrt(chain_dict['Omega_m_1'] / 0.3)
        #print S8_1.mean()
        S8_2 = chain_dict['sigma8_2'] * np.sqrt(chain_dict['Omega_m_2'] / 0.3)
        #print S8_2.mean()
        data = np.column_stack((data, S8_1))
        data = np.column_stack((data, S8_2))
        new_names.append(['S8_1', 'S_{8, \, 1}'])
        new_names.append(['S8_2', 'S_{8, \, 2}'])
    except:
        print( 'Could not calculate and append S8. \n Are Omega_m and sigma8 in the chain?')
    new_names = np.asarray(new_names, dtype=str)
    #print new_names, new_names.shape

    column_names = np.concatenate((np.asarray(['weights', 'mloglkl']), names[:, 0], np.asarray(['S8_1', 'S8_2'])))
    chi2 = 2. * data[:, 1]
    min_chi2 = chi2.min()
    '''
    nparams = len(names) - nderived
    if run == 'E':
        ndata = nzbins * (nzbins + 1) / 2. * nells
        dof = ndata - nparams - 1
    elif run == 'B':
        ndata = nzbins * (nzbins + 1) / 2. * nells
        dof = ndata - nparams - 1
    elif run == 'EB':
        ndata = nzbins * (nzbins + 1) / 2. * nells + nzbins * (nzbins + 1) / 2. * nells
        dof = ndata - nparams - 1
    elif run == 'CF' or run == 'LS' or run == 'SS' or run == 'LSSS' or run == 'SSLS' or run == '2c':
        ndata = nzbins * (nzbins + 1) * nells
        dof = ndata - nparams - 1
    min_chi2_red = min_chi2 / float(dof)

    print 'min(chi^2) = {:.3f}, min(chi^2)_red = {:.3f}, {:} d.o.f.'.format(min_chi2, min_chi2_red, int(dof))

    p_value = dist.chi2.sf(min_chi2, dof)

    print 'p-value = {:.6f}'.format(p_value)
    '''

    header = ''
    for elem in column_names:
        if elem[-1] == ' ':
            elem = elem[:-1]
        header += elem + ', '
    header = header[:-2]
    #print header
    #header =  header.tolist()
    fname = os.path.join(path_to_chain, model_name + '__HEADER.txt')
    np.savetxt(fname, data, header=header)
    print( 'Data saved to: \n', fname)
    if sampler == 'MH':
        print( 'ATTENTION: the above file does not contain the first {:}% of entries from the original chain (burn-in).'.format(threshold * 100))

    # for consistency with getdist:
    fname = os.path.join(path_to_chain, model_name + '__HEADER.paramnames')
    np.savetxt(fname, new_names, delimiter='\t', fmt='%s')

    '''
    header = 'min(chi^2), min(chi^2)_red, d.o.f., number of datapoints, number of free parameters, p-value'
    fname = path_to_chain + 'chi2_info.txt'
    savedata = np.column_stack((min_chi2, min_chi2_red, dof, ndata, nparams, p_value))
    np.savetxt(fname, savedata, header=header)
    print 'Data saved to: \n', fname
    '''

    # save chain also as FITS file:
    columns = []
    for index_col, column_name in enumerate(column_names):
        columns.append(fits.Column(name=column_name, format='D', array=data[:, index_col]))

    coldefs = fits.ColDefs(columns)
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    prihdu = fits.PrimaryHDU()
    thdulist = fits.HDUList([prihdu, tbhdu])
    fname = os.path.join(path_to_chain, model_name + '.fits')

    if not os.path.isfile(fname):
        thdulist.writeto(fname)
        print( 'Data saved to: \n', fname)
        if sampler == 'MH':
            print( 'ATTENTION: the above file does not contain the first {:}% of entries from the original chain (burn-in).'.format(threshold * 100))
    else:
        print( 'FITS file already exists. Proceeding without overwriting it!')

    return

if __name__ == '__main__':

    path_to_chain = sys.argv[1]

    # needs to be closed with '/' for glob.glob to work properly!
    if path_to_chain[-1] != '/':
        path_to_chain += '/'

    model_name = sys.argv[2]
    sampler = sys.argv[3]
    chain_is = sys.argv[4]

    if chain_is in ['2c', '2cosmos', '2cosmo', '2COSMOS', '2COSMO', 'two_cosmos', 'two_cosmo']:
        post_process_chain_2cosmos(path_to_chain, model_name, sampler=sampler)
    else:
        post_process_chain_1cosmo(path_to_chain, model_name, sampler=sampler)
