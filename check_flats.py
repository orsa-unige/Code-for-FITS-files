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
        windowsize_r = 10
        windowsize_c = 10
        for r in range(0,data.shape[0] - windowsize_r, windowsize_r):
            for c in range(0,data.shape[1] - windowsize_c, windowsize_c):
                window = data[r:r+windowsize_r,c:c+windowsize_c]
                for w in window:
                    mean = w.mean()
        if mean > 10000 and mean < 55000:
            print('FLAT ok')
        else:
            warnings.warn('NONLINEAR FLAT REGIME ') + filename)
##        plt.show()
                   
