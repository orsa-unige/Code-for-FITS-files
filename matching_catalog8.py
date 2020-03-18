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


def rotate_img(pattern): # Gets file data, rotates the image and extracts sources
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)

        # Gets file data and rotates the image
        data = fits.getdata(filename, ext=0)
        rot = rotate(data,-90,reshape=False)
        return(rot)

def get_target_coord(pattern):
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
        catalog_object_coord = SkyCoord.from_name(header['OBJECT'], parse=True)
        xobj = catalog_object_coord.ra.degree
        yobj = catalog_object_coord.dec.degree
        return(xobj,yobj)


def detect_sources(pattern):# extracts the light sources from the image (basing on sigma, FWHM, thrsholds...) 
    rot = rotate_img(pattern)
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
    return(x1,y1)

def pix_to_deg(pattern):
    coords = detect_sources(pattern)
    x1 = coords[0]
    y1 = coords[1]
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
        
    coords_obj = get_target_coord(pattern)
    xobj = coords_obj[0]
    yobj = coords_obj[1]
    
    Bin = float(header['CCDSUM'][0])
    Res = 0.25 #arcsec/px
    f_arcsec = Bin*Res #arcsec/px
    f = f_arcsec/3600

    for n in x1:
        dx = x1 - x1[30]
        x_deg = xobj + (f*dx)

    for n in y1:
        dy = y1 - y1[30]
        y_deg = yobj + (f*dy)
    return(x_deg,y_deg)

    
def get_catalog():
    coo = SkyCoord.from_name('GJ3470')
    rad = 40*u.arcmin
    cat_id = 'I/284/out'
    v=Vizier(catalog=cat_id,columns=["RAJ2000","DEJ2000","Plx","RAJ2000","DEJ2000"],column_filters={'Imag':'>13'})
    v.ROW_LIMIT = -1
    tab = v.query_region(coo, radius=rad, catalog = cat_id)[0]
    x2 = tab['RAJ2000']
    y2 = tab['DEJ2000']
    print(coo.galactic)
    return(x2,y2)

def plot_imag_and_catalog(pattern):
    fig1,ax = plt.subplots()
    coords_deg = pix_to_deg(pattern)
    x_deg = coords_deg[0]
    y_deg = coords_deg[1]
    coords_obj = get_target_coord(pattern)
    x_obj = coords_obj[0]
    y_obj = coords_obj[1]
    ax.plot(x_deg,y_deg, 'ok', ms=5)
    ax.plot(x_obj, y_obj, 'oy', mfc='none',ms=10)
    fig2, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlim(119.65,119.85)
    ax.set_ylim(15.3,15.5)
    coords_cat = get_catalog()
    x2 = coords_cat[0]
    y2 = coords_cat[1]
    ax.plot(x2,y2,'or',mfc='none',ms=7)
    ax.plot(x_deg,y_deg,'ob',ms=4)
    ax.plot(x_obj,y_obj,'oy',ms=4)
    plt.show()
    return plt.show()

def main():
    pattern = sys.argv[1:]
    plot_imag_and_catalog(pattern)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 :
        print(" Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()
    main()
