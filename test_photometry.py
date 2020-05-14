#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.vizier import Vizier
from photutils import CircularAperture, aperture_photometry, CircularAnnulus
import matplotlib.pyplot as plt
import astropy.units as u
from astroquery.mast import Catalogs
from astropy.table import Table

# our functions
from reduction import get_fits_header, get_fits_data

def img_coords(pattern):
    for filename in pattern:
        header = get_fits_header(filename)
        w = WCS(header)
        data = get_fits_data(filename)
        mean, median, std = sigma_clipped_stats(data, sigma=3)
        daofind = DAOStarFinder(fwhm=3, threshold=5.*std)
        sources = daofind(data - median)
        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output
       # Pixel coordinates of the sources
        c = dict();
        c['x1'] = np.array(sources['xcentroid'])
        c['y1'] = np.array(sources['ycentroid'])
        c['loc'] = (c['x1'],c['y1'])
        return(c)


def circular_aperture_img(pattern):
    for filename in pattern:
        header = get_fits_header(filename)
        wcs = WCS(header)
        c = img_coords(pattern)
        loc = c['loc']
        e = dict();
        e['data'] = get_fits_data(filename)
        e['positions'] = np.transpose(loc)
        e['apertures'] = CircularAperture(e['positions'], r=6.8)
        e['annulus_aperture'] = CircularAnnulus(e['positions'], r_in=10, r_out=15)
        phot_table = aperture_photometry(e['data'], e['apertures'])
        for col in phot_table.colnames:
            phot_table[col].info.format = '%.8g'
        bkg_mean = phot_table['aperture_sum'] / e['annulus_aperture'].area
        bkg_sum = bkg_mean * e['apertures'].area
        final_sum = phot_table['aperture_sum'] - bkg_sum
        phot_table['residual_aperture_sum'] = final_sum
        print(phot_table) 
##        print(phot_table[29]) #GJ3470 location

        return(e)

def circular_aperture_catalog(pattern):
    for filename in pattern:
        header = get_fits_header(filename)
        catalog = Catalogs.query_region(f'{header["CRVAL1"]} {header["CRVAL2"]}',
                                        radius = '0.2 deg',
                                        catalog = 'TIC')
        dattab = Table(catalog)
        radec = (f['catalog']['ra','dec','Bmag'])
        mask = radec['Bmag']<19.0
        mag_radec = radec[mask]
        positions_cat = SkyCoord(mag_radec['ra'],
                             mag_radec['dec'],
                             frame='icrs',
                             unit=(u.deg,u.deg))
        aperture_cat = SkyCircularAperture(positions_cat,
                                           r=4.5*u.arcsec) 
        annulus_cat = SkyCircularAnnulus(positions_cat,
                                         r_in=5*u.arcsec,
                                         r_out=8*u.arcsec)
        apers = [aperture_cat, annulus_cat]
        return(apers)

def plot_img(pattern):
    e = circular_aperture_img(pattern)
    plt.figure()
    plt.imshow(e['data'], origin='lower',vmin=700,vmax=10000)
    e['apertures'].plot(color='blue', lw=1.5, alpha=0.5)
    e['annulus_aperture'].plot(color='white',lw=2)
    plt.colorbar()

    return(plt.show())
    
    

    

def main():
    pattern = sys.argv[1:]
    circular_aperture_img(pattern)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 :
        print(" Usage:  "+sys.argv[0]+" fits file or list of fits files")
        sys.exit()
    main()     
