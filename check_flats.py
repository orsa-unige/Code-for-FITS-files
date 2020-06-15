#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
#import matplotlib.pyplot as plt
#import warnings
from astropy import log # Questo è uguale a warnings però lo uso già dappertutto

def detect_flats(data, size, min_val=10000, max_val=55000, saturated=65536):
    # lenr, lenc = int(data.shape[0]/size), int(data.shape[1]/size)
    # NBINS = 1000
    # histogram = plt.hist(data.flatten(), NBINS)
    # data_split = np.array([data[i*size:(i+1)*size,j*size:(j+1)*size]
    #       for (i,j) in np.ndindex(lenr,lenc)]).reshape(lenr,lenc,size,size)
    data_split = np.array_split(data, size) # Questo evita lenr lenc, gestisce lui la divisione.
    data_split_avg = [np.mean(arr) for arr in data_split]
    for avg in data_split_avg:
        #goodflat_condition = min_val<avg<max_val # Alla fine lo usi solo una volta
        if min_val < avg < max_val: # Se è True, non serve esplicitarlo
            is_it_a_good_flat = True
            log.info(is_it_a_good_flat) # Uso log
            continue # Brava!!!
        else: # Va bene come l'hai fatto tu. Questa è un'alternativa:
            is_it_a_good_flat = False
            if avg >= saturated: log.warning('SATURATED FLAT ')
            else: log.warning('NONLINEAR FLAT REGIME ')
            break
##    plt.show()
    return(is_it_a_good_flat)


# Se vogliamo scorrere la lista una sola volta senza stare a cercare
# in particolare se è saturato o non lineare, allora conviene solo
# qualcosa tipo:

def check_counts(data, size, min_val=10000, max_val=55000):
    data_split = np.array_split(data, size)
    data_split_avg = [ np.mean(arr) for arr in data_split ]    
    if all([ min_val < avg < max_val for avg in data_split_avg ]):
        log.info("Ok")
        return True
    else:
        log.warn("Some chunk is is saturated or non linear")
        return False


