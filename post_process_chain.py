#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:16:20 2017

@author: fkoehlin
"""

import os
import sys
import glob
import numpy as np
import astropy.io.fits as fits
#import scipy.stats.distributions as dist

def post_process_chain_1cosmo(path_to_chain, model_name):
    
    fname = os.path.join(path_to_chain, 'chain_NS__accepted.txt')
    
    data = np.loadtxt(fname)
    S8 = data[:, -1] * np.sqrt(data[:, -2] / 0.3)
    
    data = np.column_stack((data, S8))
    # glob can expand names with *-operator!
    fname =  glob.glob(path_to_chain + '*_.paramnames')[0]
    names = np.loadtxt(fname, dtype=str, delimiter='\t')
    new_names = names.tolist() #[:-1, :] = names[:, :]
    new_names.append(['S8', 'S_{8}'])
    new_names = np.asarray(new_names, dtype=str)
    # better TeX:
    new_names[np.where(new_names[:, 0] == 'sigma8 '), 1] = '\sigma_8 '
    # \mathrm doesn't seem to work...
    #new_names[np.where(new_names[:, 0] == 'Omega_m '), 1] = '\Omega_\mathrm{m} '
    #new_names[np.where(new_names[:, 0] == 'omega_cdm '), 1] = '\omega_\mathrm{cdm} ' 
    #new_names[np.where(new_names[:, 0] == 'omega_b '), 1] = '\omega_\mathrm{b} ' 
    #new_names[np.where(new_names[:, 0] == 'n_s '), 1] = 'n_\mathrm{s} ' 
    new_names[np.where(new_names[:, 0] == 'ln10^{10}A_s '), 1] = '\ln 10^{10} A_s ' 
    
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
    fname = os.path.join(path_to_chain, model_name + '_chain_NS__accepted_HEADER.txt')
    np.savetxt(fname, data, header=header)
    print 'Data saved to: \n', fname
    
    # for consistency with getdist:
    fname = os.path.join(path_to_chain, model_name + '_chain_NS__accepted_HEADER.paramnames')
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
        print 'Data saved to: \n', fname
    else:
        print 'FITS file already exists. Proceeding without overwriting it!'
    

    return 

def post_process_chain_2cosmos(path_to_chain, model_name):
                
    fname = path_to_chain + 'chain_NS__accepted.txt'
    data = np.loadtxt(fname)
    #print data
    #print data[:, -1]
    S8_1 = data[:, -3] * np.sqrt(data[:, -4] / 0.3)
    #print S8_1.mean()
    S8_2 = data[:, -1] * np.sqrt(data[:, -2] / 0.3)
    #print S8_2.mean()
    data = np.column_stack((data, S8_1))
    data = np.column_stack((data, S8_2))
    # glob can expand names with *-operator!
    fname =  glob.glob(path_to_chain + '*_.paramnames')[0]
    #print fname
    names = np.loadtxt(fname, dtype=str, delimiter='\t')
    #print names
    names = np.loadtxt(fname, dtype=str, delimiter='\t')
    #print names, names.shape
    new_names = names.tolist() #[:-1, :] = names[:, :]
    new_names.append(['S8_1', 'S_{8, \, 1}'])
    new_names.append(['S8_2', 'S_{8, \, 2}'])
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
    fname = os.path.join(path_to_chain, model_name + '_chain_NS__accepted_HEADER.txt')
    np.savetxt(fname, data, header=header)
    print 'Data saved to: \n', fname
    
    # for consistency with getdist:
    fname = os.path.join(path_to_chain, model_name + '_chain_NS__accepted_HEADER.paramnames')
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
        print 'Data saved to: \n', fname
    else:
        print 'FITS file already exists. Proceeding without overwriting it!'
        
    return

if __name__ == '__main__':    
    
    path_to_chain = sys.argv[1]
    
    # needs to be closed with '/' for glob.glob to work properly!
    if path_to_chain[-1] != '/':
        path_to_chain += '/'

    model_name = sys.argv[2]
    chain_is = sys.argv[3]

    # number of derived parameters included in chain (e.g. Omega_m, sigma8):
    #nderived = sys.argv[3]
    
    # number of independent data points in data-vector:
    #ndata = sys.argv[4]
    
    if chain_is in ['2cosmos', '2cosmo', '2COSMOS', '2COSMO', 'two_cosmos', 'two_cosmo']:
        post_process_chain_2cosmos(path_to_chain, model_name)
    else:
        post_process_chain_1cosmo(path_to_chain, model_name)