#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 16:12:20 2018

@author: fkoehlin
"""

import numpy as np

def cm2inch(value):

    return value / 2.54

# Bayesian way of defining confidence intervals:
def minimum_credible_intervals(values, central_value, weights, bins=20):
    """
    Extract minimum credible intervals (method from Jan Haman) FIXME

    copy & paste from Monte Python (2.1.2) with own modifications
    --> checked that this function returns same output as Monte Python; modifications are all okay!!!

    """
    #histogram = info.hist
    #bincenters = info.bincenters
    #levels = info.levels

    histogram, bin_edges = np.histogram(values, bins=bins, weights=weights, density=False)
    bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])

    # Defining the sigma contours (1, 2 and 3-sigma)
    levels = np.array([68.27, 95.45, 99.73])/100.

    bounds = np.zeros((len(levels), 2))
    j = 0
    delta = bincenters[1]-bincenters[0]
    left_edge = max(histogram[0] - 0.5*(histogram[1]-histogram[0]), 0.)
    right_edge = max(histogram[-1] + 0.5*(histogram[-1]-histogram[-2]), 0.)
    failed = False
    for level in levels:
        norm = float(
            (np.sum(histogram)-0.5*(histogram[0]+histogram[-1]))*delta)
        norm += 0.25*(left_edge+histogram[0])*delta
        norm += 0.25*(right_edge+histogram[-1])*delta
        water_level_up = np.max(histogram)*1.0
        water_level_down = np.min(histogram)*1.0
        top = 0.

        iterations = 0
        while (abs((top/norm)-level) > 0.0001) and not failed:
            top = 0.
            water_level = (water_level_up + water_level_down)/2.
            #ontop = [elem for elem in histogram if elem > water_level]
            indices = [i for i in range(len(histogram))
                       if histogram[i] > water_level]
            # check for multimodal posteriors
            '''
            if ((indices[-1]-indices[0]+1) != len(indices)):
                print'Could not derive minimum credible intervals for this multimodal posterior!' + \
                     'Please try running longer chains or reducing the number of bins (default: 40)'
                failed = True
                break
            '''
            top = (np.sum(histogram[indices]) -
                   0.5*(histogram[indices[0]]+histogram[indices[-1]]))*(delta)

            # left
            if indices[0] > 0:
                top += (0.5*(water_level+histogram[indices[0]]) *
                        delta*(histogram[indices[0]]-water_level) /
                        (histogram[indices[0]]-histogram[indices[0]-1]))
            else:
                if (left_edge > water_level):
                    top += 0.25*(left_edge+histogram[indices[0]])*delta
                else:
                    top += (0.25*(water_level + histogram[indices[0]]) *
                            delta*(histogram[indices[0]]-water_level) /
                            (histogram[indices[0]]-left_edge))

            # right
            if indices[-1] < (len(histogram)-1):
                top += (0.5*(water_level + histogram[indices[-1]]) *
                        delta*(histogram[indices[-1]]-water_level) /
                        (histogram[indices[-1]]-histogram[indices[-1]+1]))
            else:
                if (right_edge > water_level):
                    top += 0.25*(right_edge+histogram[indices[-1]])*delta
                else:
                    top += (0.25*(water_level + histogram[indices[-1]]) *
                            delta * (histogram[indices[-1]]-water_level) /
                            (histogram[indices[-1]]-right_edge))

            if top/norm >= level:
                water_level_down = water_level
            else:
                water_level_up = water_level
            # safeguard, just in case
            iterations += 1
            if (iterations > 1e4):
                print('The loop to check for sigma deviations was taking too long to converge.')
                failed = True
                break

        # min
        if failed:
            bounds[j][0] = np.nan
        elif indices[0] > 0:
            bounds[j][0] = bincenters[indices[0]] - delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-histogram[indices[0]-1])
        else:
            if (left_edge > water_level):
                bounds[j][0] = bincenters[0]-0.5*delta
            else:
                bounds[j][0] = bincenters[indices[0]] - 0.5*delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-left_edge)

        # max
        if failed:
            bounds[j][1] = np.nan
        elif indices[-1] < (len(histogram)-1):
            bounds[j][1] = bincenters[indices[-1]] + delta*(histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-histogram[indices[-1]+1])
        else:
            if (right_edge > water_level):
                bounds[j][1] = bincenters[-1]+0.5*delta
            else:
                bounds[j][1] = bincenters[indices[-1]] + \
                    0.5*delta*(histogram[indices[-1]]-water_level) / \
                    (histogram[indices[-1]]-right_edge)

        j += 1

    for elem in bounds:
        for j in (0, 1):
            elem[j] -= central_value

    return bounds

def weighted_mean(values, weights=None):

    if weights is None:
        weights = np.ones_like(values)

    return np.sum(weights*values)/np.sum(weights)

def quantile(x, q, weights=None):
    """
    Like numpy.percentile, but:
    * Values of q are quantiles [0., 1.] rather than percentiles [0., 100.]
    * scalar q not supported (q must be iterable)
    * optional weights on x
    """
    if weights is None:
        return np.percentile(x, [100. * qi for qi in q])
    else:
        idx = np.argsort(x)
        xsorted = x[idx]
        cdf = np.add.accumulate(weights[idx])
        cdf /= cdf[-1]

        return np.interp(q, cdf, xsorted).tolist()
