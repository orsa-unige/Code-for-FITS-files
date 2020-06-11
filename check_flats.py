#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import warnings

from fits import get_fits_data


def detect_flats(pattern):
    for filename in pattern:
        data = get_fits_data(filename)
        NBINS = 1000
        histogram = plt.hist(data.flatten(), NBINS)
        data_split = np.array_split(data, 100)
        data_split_avg = [np.mean(arr) for arr in data_split]
        for n in data_split_avg: 
           if goodflat_condition == True:
                answer = True
            elif n >= 65536:
                answer = False
                warnings.warn('SATURATED FLAT ' + filename)
            else:
                 answer = False
                 warnings.warn('NONLINEAR FLAT REGIME ' + filename)
        print(answer, filename)

##        plt.show()
                   
