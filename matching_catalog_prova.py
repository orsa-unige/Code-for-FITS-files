#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from photutils import detect_threshold
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
from photutils import deblend_sources
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
import numpy as np
from astropy.visualization import SqrtStretch
from astroquery.vizier import Vizier
from astropy.coordinates import match_coordinates_sky
from scipy.ndimage import rotate
import math

def choose_hdu(filename):
    finfo = fits.info(filename, output=False) # Returns a list of tuples.
    finfo_list = [item for item in finfo if 'COMPRESSED_IMAGE' in item]
    if not finfo_list:
        return finfo[0][0] # 0 if not compressed
    else:
        return finfo_list[0][0]

def get_data(pattern): # Gets file data, rotates the image and extracts sources
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)

        # Gets file data and rotates the image
        data = fits.getdata(filename, ext=0)
        rot = rotate(data,-90,reshape=False)

        # extracts the light sources from the image (basing on sigma, FWHM, thrsholds...) 
        threshold = detect_threshold(rot, nsigma=2.)
        sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
        kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
        kernel.normalize()
        mean, median, std = sigma_clipped_stats(rot, sigma=3)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
        sources = daofind(rot - median)
        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output

       # Pixel coordinates of the sources
        x1 = np.array(sources['xcentroid'])
        y1 = np.array(sources['ycentroid'])

def get_target_coord(pattern):
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
        catalog_object_coord = SkyCoord.from_name(header['OBJECT'], parse=True)
        xobj = catalog_object_coord.ra.degree
        yobj = catalog_object_coord.dec.degree
        print(xobj,yobj)        
    
        
def main():
    pattern = sys.argv[1:]
    get_data(pattern)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 :
        print(" Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()
    main()
