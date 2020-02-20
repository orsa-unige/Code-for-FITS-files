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
from astropy.visualization.mpl_normalize import ImageNormalize

def choose_hdu(filename):
    finfo = fits.info(filename, output=False) # Returns a list of tuples.
    finfo_list = [item for item in finfo if 'COMPRESSED_IMAGE' in item]
    if not finfo_list:
        return finfo[0][0] # 0 if not compressed
    else:
        return finfo_list[0][0]

def get_wcs(pattern):
    for filename in pattern:

        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
        print(header)
        w = WCS(header)
        print(w)

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

        x_unscaled = np.array(sources['xcentroid'])
        y_unscaled = np.array(sources['ycentroid'])
        x_min = min(x_unscaled)
        y_min = min(y_unscaled)
        x_max = max(x_unscaled)
        y_max = max(y_unscaled)
        x1 = ((x_unscaled - x_min)/x_max)*90
        y1 = ((y_unscaled - y_min)/y_max)*90
        fig,ax=plt.subplots()
        ax.plot(x1,y1,'ok',ms=5)
            
        plt.show()
        c1 = SkyCoord(ra=x1*u.degree, dec=y1*u.degree)
        print(c1)
        coo = SkyCoord.from_name('GJ3470')
    

        
def main():
    pattern = sys.argv[1:]
    get_wcs( pattern)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 :
        print(" Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()
    main()

