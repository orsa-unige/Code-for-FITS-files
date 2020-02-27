#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
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


# the code verify if the file (or list of files) is or is not compressed and gets the headers

def choose_hdu(filename):
    finfo = fits.info(filename, output=False) # Returns a list of tuples.
    finfo_list = [item for item in finfo if 'COMPRESSED_IMAGE' in item]
    if not finfo_list:
        return finfo[0][0] # 0 if not compressed
    else:
        return finfo_list[0][0]

# finds the stars in the image by a given FWHM,
# gets the stars coordinates, transform them from pixel to deg

def get_wcs(pattern):
    for filename in pattern:

        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
        
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
       
        x = np.array(sources['xcentroid'])
        y = np.array(sources['ycentroid'])
        OB = header['OBJECT']
        cobj = SkyCoord.from_name(OB, parse=True)
        cobjx = cobj.ra.degree
        cobjy = cobj.dec.degree
        refx = x[1]
        refy = y[1]
        
        Bin = float(header['CCDSUM'][0])
        Res = 0.25 #arcsec/px
        f_arcsec = Bin*Res #arcsec/px
        f = f_arcsec/3600
        
        for n in x:
            dx = x - refx
            xscal = cobjx + (f*dx)

        for n in y:
            dy = y - refy
            yscal = cobjy + (f*dy)
            
        fig1,ax=plt.subplots()
        c1 = SkyCoord(ra=xscal*u.degree, dec=yscal*u.degree) # ra and dec of each source
        ax.plot(xscal,yscal,'ok',ms=5)
        
        result = Vizier.query_object(OB)
        interesting_table = result['I/305/out']
        x1 = np.array(interesting_table['RAJ2000'])
        y1 = np.array(interesting_table['DEJ2000'])
        c2 = SkyCoord(ra=x1*u.degree, dec = y1*u.degree)
        fig2, ax = plt.subplots()
        ax.plot(x1,y1,'or',mfc='none',ms=10)

        idx, d2d, d3d = c1.match_to_catalog_sky(c2)
        fig3,ax=plt.subplots()
        ax.plot(xscal,yscal,'ok',ms=5)
        ax.plot(x1,y1,'or',mfc='none',ms=10)
        plt.show()

def main():
    pattern = sys.argv[1:]
    get_wcs( pattern)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 :
        print(" Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()
    main()
