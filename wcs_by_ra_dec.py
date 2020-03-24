#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
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
import math
from astropy.coordinates import FK5

def choose_hdu(filename):
    finfo = fits.info(filename, output=False) # Returns a list of tuples.
    finfo_list = [item for item in finfo if 'COMPRESSED_IMAGE' in item]
    if not finfo_list:
        return finfo[0][0] # 0 if not compressed
    else:
        return finfo_list[0][0]

def get_header(pattern):
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
        return(header)

def get_target_coords(pattern):
    header = get_header(pattern)
    d = dict();
    d['target_coords'] = SkyCoord.from_name(header['OBJECT'], parse=True)
    d['ra_target'] = d['target_coords'].ra.degree
    d['dec_target'] = d['target_coords'].dec.degree
    return(d)

def get_wcs(pattern):
    for filename in pattern:
        header = get_header(pattern)
        c = SkyCoord(header['RA'], header['DEC'], unit=(u.hourangle, u.deg))
        w = WCS(filename)
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.crpix = [header['NAXIS1']/2, header['NAXIS2']/2]
        w.wcs.crval = [c.ra.degree, c.dec.degree]
        w.wcs.cd = [[0, 0.00015], [0.00015, 0]]
        return(w)
    
def get_sources_coords(pattern):
    for filename in pattern:
        data = fits.getdata(filename, ext=0)
        threshold = detect_threshold(data, nsigma=2.)
        sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
        kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
        kernel.normalize()
        mean, median, std = sigma_clipped_stats(data, sigma=3)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
        sources = daofind(data - median)
        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output
       # Pixel coordinates of the sources
        w = get_wcs(pattern)
        c = dict();
        x1 = np.array(sources['xcentroid'])
        y1 = np.array(sources['ycentroid'])
        data1 = np.array([x1,y1])
        world = w.wcs_pix2world(x1,y1,0)
        c['ra_img'] = world[0]
        c['dec_img'] = world[1]
        return(c)
    
def get_catalog(pattern):
    d =  get_target_coords(pattern)
    cooEQ =  d['target_coords'].fk5.transform_to(FK5(equinox="J2019.9"))
    rad = 40*u.arcmin
    cat_id = 'I/284/out'
    v=Vizier(catalog=cat_id,columns=["RAJ2000","DEJ2000","Plx","RAJ2000","DEJ2000"],column_filters={'Imag':'>14'})
    v.ROW_LIMIT = -1
    tab = v.query_region(cooEQ, radius=rad, catalog = cat_id)[0]
    g = dict();
    g['ra_cat'] = tab['RAJ2000']
    g['dec_cat'] = tab['DEJ2000']
    return(g)

def plot_img(pattern):
    fig1, ax = plt.subplots()
    c = get_sources_coords(pattern)
    ax.plot(c['ra_img'],c['dec_img'], 'ok', ms=5)

    fig2, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlim(120,120.20)
    ax.set_ylim(15.25,15.45)
    g = get_catalog(pattern)
    ax.plot(g['ra_cat'],g['dec_cat'],'or',mfc='none',ms=7)
    ax.plot(c['ra_img'],c['dec_img'],'ob',ms=4)
    plt.show()
    return plt.show()

def main():
    pattern = sys.argv[1:]
    plot_img(pattern)
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 :
        print(" Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()
    main()
