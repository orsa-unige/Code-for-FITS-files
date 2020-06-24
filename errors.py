#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
from astropy import log
from astropy.wcs import WCS
from astropy.io import fits

# our functions
from fits import get_fits_data, get_fits_header
from photometry import load_catalog, set_apertures, do_photometry                             


def signal_noise(filenames):
    for filename in filenames:
        catalog = load_catalog(filename)
        apers = set_apertures(catalog)
        data = get_fits_data(filename)
        header = get_fits_header(filename)
        wcs = WCS(header)
        phot_table = do_photometry(data, apers=apers, wcs=wcs)
        pixar = apers[0].to_pixel(wcs)
        pixan = apers[1].to_pixel(wcs)
        bkg_mean = phot_table['aperture_sum_1'] / pixan.area
        R_star = phot_table['residual_aperture_sum']     
        R_sky = bkg_mean * pixar.area
        n_pix_aperture = pixan.area
        
        with open('instruments.json') as jfile:
            instrument_dict = json.load(jfile)

        instrument_id = header['INSTRUME']
        
        if instrument_id == 'Mexman':
            key = instrument_dict['Mexman']
            gain = key['gain']
        elif instrument_id == 'DFOSC_FASU':
            log.warning('gain == null ')
        elif instrument_id == 'STL-11000 3 CCD Camera':
            key = instrument_dict['STL-11000 3 CCD Camera']
            gain = key['gain']
        else:
            log.warning('INSTRUMENT NOT FOUND ')
            
        
    return(R_star, R_sky, n_pix_aperture, gain)
         
