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
        data = fits.getdata(filename, ext=0)
        rot = rotate(data,-90,reshape=False)
        return(rot)

def get_target_coord(pattern):
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
        catalog_object_coord = SkyCoord.from_name(header['OBJECT'], parse=True)
        c = dict();
        c['ra_obj'] = catalog_object_coord.ra.degree
        c['dec_obj'] = catalog_object_coord.dec.degree
        return(c)


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
    d = dict();
    d['x1'] = np.array(sources['xcentroid'])
    d['y1'] = np.array(sources['ycentroid'])
    return(d)

def pix_to_deg(pattern):
    c = get_target_coord(pattern)
    d = detect_sources(pattern)
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
    
    Bin = float(header['CCDSUM'][0])
    Res = 0.25 #arcsec/px
    f_arcsec = Bin*Res #arcsec/px
    f = f_arcsec/3600
    p = dict();

    for n in d['x1']:
        dx = d['x1'] - d['x1'][30]
        p['ra_img_deg'] = c['ra_obj'] + (f*dx)

    for n in d['y1']:
        dy = d['y1'] - d['y1'][30]
        p['dec_img_deg'] = c['dec_obj'] + (f*dy)
    return(p)

def get_catalog():
    coo = SkyCoord.from_name('GJ3470')
    rad = 40*u.arcmin
    cat_id = 'I/284/out'
    v=Vizier(catalog=cat_id,columns=["RAJ2000","DEJ2000","Plx","RAJ2000","DEJ2000"],column_filters={'Imag':'>13'})
    v.ROW_LIMIT = -1
    tab = v.query_region(coo, radius=rad, catalog = cat_id)[0]
    g = dict();
    g['ra_cat'] = tab['RAJ2000']
    g['dec_cat'] = tab['DEJ2000']
    print(coo.galactic)
    return(g)

def plot_imag_and_catalog(pattern):
    c = get_target_coord(pattern)
    d = detect_sources(pattern)
    p = pix_to_deg(pattern)
    g = get_catalog()
    fig1,ax = plt.subplots()
    ax.plot(p['ra_img_deg'],p['dec_img_deg'], 'ok', ms=5)
    ax.plot(c['ra_obj'], c['dec_obj'], 'oy', mfc='none',ms=10)
    
    fig2, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlim(119.65,119.85)
    ax.set_ylim(15.3,15.5)
    ax.plot(g['ra_cat'],g['dec_cat'],'or',mfc='none',ms=7)
    ax.plot(p['ra_img_deg'],p['dec_img_deg'],'ob',ms=4)
    ax.plot(c['ra_obj'], c['dec_obj'],'oy',ms=4)
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
