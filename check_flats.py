#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import warnings

from fits import get_fits_header, get_fits_data


def detect_flats(pattern):
    for filename in pattern:
        header = get_fits_header(filename)
        data = get_fits_data(filename)
        NBINS = 1000
        histogram = plt.hist(data.flatten(), NBINS)
        data_split = np.array_split(data, 100)
        data_split_avg = [np.mean(arr) for arr in data_split]
        for n in data_split_avg: 
            goodflat_condition = 10000<n<55000
            if goodflat_condition == True:
                pass
            elif n >= 65536:
                warnings.warn('SATURATED FLAT ' + filename)
            else: 
                warnings.warn('NONLINEAR FLAT REGIME ' + filename)


##        plt.show()
                   
