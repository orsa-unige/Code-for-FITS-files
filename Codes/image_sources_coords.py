#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def get_sources_coords(target, pattern):
    ''' Takes objects coords from .fits files and transform them to RA and DEC
        using the WCS keywords'''
    from astropy.io import fits
    from astropy.stats import sigma_clipped_stats
    from photutils import DAOStarFinder
    import numpy as np
    from manipulate_header import wcs, get_header
    w = wcs(target, pattern)
    header = get_header(pattern)
    
    for filename in pattern:
        data = fits.getdata(filename, ext=0)
        mean, median, std = sigma_clipped_stats(data, sigma=3)
        daofind = DAOStarFinder(fwhm=3, threshold=5.*std)
        sources = daofind(data - median)
        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output
       # Pixel coordinates of the sources
        x1 = np.array(sources['xcentroid'])
        y1 = np.array(sources['ycentroid'])
        world = w.wcs_pix2world(x1,y1,0)
        c = dict();
        world_x = world[0]
        world_y = world[1]
        c['RA_img'] = world_x
        c['DEC_img'] = world_y
        return(c)

def main():
    target = sys.argv[1]
    pattern = sys.argv[2:]
    get_sources_coords(target, pattern)
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3 :
        print(" Usage:  "+sys.argv[0]+" target name, catalog name and fits file or list of fits files")
        sys.exit()
    main()     
    

