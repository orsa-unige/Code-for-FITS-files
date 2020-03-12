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

       # Gets the coordinates (deg) of the target  
        catalog_object_coor = SkyCoord.from_name(header['OBJECT'], parse=True)
        xobj = catalog_object_coor.ra.degree
        yobj = catalog_object_coor.dec.degree

        # Transforms pixels into degs
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

        # plots the file's coordinates
        fig1,ax=plt.subplots()
        ax.plot(x_deg,y_deg, 'ok', ms=5)
        ax.plot(xobj, yobj, 'oy', mfc='none',ms=10)

        # Gets a catalog by name and its objects coordinates
        coo = SkyCoord.from_name('GJ3470')
        rad = 40*u.arcmin
        cat_id = 'I/284/out'
        v=Vizier(catalog=cat_id,columns=["RAJ2000","DEJ2000","Plx","RAJ2000","DEJ2000"],column_filters={'RAJ2000':'!null','DEJ2000':'!null'}) 
        v.ROW_LIMIT = -1
        tab = v.query_region(coo, radius=rad, catalog = cat_id)[0]
        x2 = tab['RAJ2000']
        y2 = tab['DEJ2000']

        # plots the file's coordinates (blu), the target (yeallow) and the catalog (red & black) 
        fig2, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.plot(x2,y2,'.k',ms=1)
        pmdd = .5
        cond=((x2-xobj)**2+(y2-yobj)**2)<pmdd**2
        ax.plot(x2[cond],y2[cond],'.r',ms=1)
        ax.plot(x_deg,y_deg,'ob',ms=2)
        ax.plot(xobj,yobj,'oy',ms=3)
        plt.show()
        print(coo.galactic)

        
def main():
    pattern = sys.argv[1:]
    get_wcs( pattern)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 :
        print(" Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()
    main()
