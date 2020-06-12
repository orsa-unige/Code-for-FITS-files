#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import warnings

from fits import get_fits_data


def detect_flats(filenames):
    for filename in filenames:
        data = get_fits_data(filename)
        NBINS = 1000
        histogram = plt.hist(data.flatten(), NBINS)
        data_split = np.array_split(data, 100)
        data_split_avg = [np.mean(arr) for arr in data_split]
        for n in data_split_avg: 
            goodflat_condition = 10000<n<55000
            if goodflat_condition == True:
                is_it_a_good_flat = True
            elif n >= 65536:
                is_it_a_good_flat = False
                warnings.warn('SATURATED FLAT ' + filename)
            else:
                 is_it_a_good_flat = False
                 warnings.warn('NONLINEAR FLAT REGIME ' + filename)
        print(is_it_a_good_flat, filename)

##        plt.show()
                   
