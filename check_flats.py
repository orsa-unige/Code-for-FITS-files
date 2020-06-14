#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import warnings

def detect_flats(data, size, saturated = 65536, min_val = 10000, max_val = 55000):
    lenr, lenc = int(data.shape[0]/size), int(data.shape[1]/size)
    NBINS = 1000
    histogram = plt.hist(data.flatten(), NBINS)
    data_split = np.array([data[i*size:(i+1)*size,j*size:(j+1)*size]
          for (i,j) in np.ndindex(lenr,lenc)]).reshape(lenr,lenc,size,size)
    data_split_avg = [np.mean(arr) for arr in data_split]
    for avg in data_split_avg:
        goodflat_condition = min_val<avg<max_val
        if goodflat_condition == True:
            is_it_a_good_flat = True
            print(is_it_a_good_flat)
            continue
        elif avg >= saturated:
            is_it_a_good_flat = False
            print(is_it_a_good_flat)
            warnings.warn('SATURATED FLAT ')
            break
        else:
            is_it_a_good_flat = False
            print(is_it_a_good_flat)
            warnings.warn('NONLINEAR FLAT REGIME ')
            break
##    plt.show()
    return(is_it_a_good_flat)
        
