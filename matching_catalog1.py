#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.wcs import WCS

def get_wcs(pattern):
    for filename in pattern:
        with fits.open(filename, 'update') as hdul:
            hdr = hdul[0].header
            key = hdr['OBJECT']
            if key == 'wasp16':
                w = WCS(filename)
                print(filename)
                print(key)
                print(w)

                image_data = fits.getdata(filename, ext=0)
                from photutils import detect_threshold
                threshold = detect_threshold(image_data, 2)

                from astropy.convolution import Gaussian2DKernel
                from astropy.stats import gaussian_fwhm_to_sigma
                from photutils import detect_sources
                sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
                kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
                kernel.normalize()
                segm = detect_sources(image_data, threshold, npixels=5, filter_kernel=kernel)

                import numpy as np
                import matplotlib.pyplot as plt
                from astropy.visualization import SqrtStretch
                from astropy.visualization.mpl_normalize import ImageNormalize
                norm = ImageNormalize(stretch=SqrtStretch())
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
                ax1.imshow(image_data, origin='lower', cmap='Greys_r', norm=norm)
                ax1.set_title('Data')
                cmap = segm.cmap(random_state=12345)
                ax2.imshow(segm, origin='lower', cmap=cmap)
                ax2.set_title('Segmentation Image')
                plt.show()

##                from astropy.coordinates import SkyCoord
##                from astropy import units as u
##                RA = hdr['RA']
##                DEC = hdr ['DEC']
##                print(RA,DEC)
##                c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
##                catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  Come ottengo RA e DEC del catalogo?
##                idx, d2d, d3d = c.match_to_catalog_sky(catalog)
##                matches = catalog[idx]
##                (matches.separation_3d(c) == drd).all()
##                dra = (mathces.ra - c.ra).arcmin
##                ddec = (matches.dec - c.dec).arcmin
##                
                    
                
def main():
    pattern = sys.argv[1:]
    get_wcs( pattern)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 :
        print(" Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()

    main()
