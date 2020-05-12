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


def choose_hdu(filename):
    ''' Verifies if fits file is compressed or not. only .fits or .fits.fz
        format are accepted'''
    finfo = fits.info(filename, output=False) # Returns a list of tuples.
    finfo_list = [item for item in finfo if 'COMPRESSED_IMAGE' in item]
    if not finfo_list:
        return finfo[0][0] # 0 if not compressed
    else:
        return finfo_list[0][0]
    

def get_header(pattern):
    '''Returns the HEADER parameter, that will be used to obtain coordinates
        and catalog stuff'''
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
        return(header)

def wcs(pattern):
    header = get_header(pattern)
    w = WCS(header)
    return(w)


def img_coords(pattern):
    w = wcs(pattern)
    header = get_header(pattern)
    for filename in pattern:
        data = fits.getdata(filename, ext=0)
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
        world = w.wcs_pix2world(c['x1'],c['y1'],0)
        world_x = world[0]
        world_y = world[1]
        c['RA_img'] = world_x
        c['DEC_img'] = world_y
        return(c)

def get_target(pattern):
    header = get_header(pattern)
    obj = header['OBJECT']
    d = dict();
    d['target_coords'] = SkyCoord.from_name(obj)
    d['ra_target'] = d['target_coords'].ra.degree
    d['dec_target'] = d['target_coords'].dec.degree
    return(d)

def circular_aperture(pattern):
    for filename in pattern:
        c = img_coords(pattern)
        data = fits.getdata(filename)
        loc = c['loc']
        e = dict();
        e['positions'] = np.transpose(loc)
        apertures = CircularAperture(e['positions'], r=6.8)
        header = get_header(pattern)
        wcs = WCS(header)
        plt.figure()
        plt.imshow(data, origin='lower',vmin=700,vmax=10000)
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        phot_table = aperture_photometry(data, apertures)
        annulus_aperture = CircularAnnulus(e['positions'], r_in=10, r_out=15)
        annulus_aperture.plot(color='white',lw=2)
        plt.colorbar()

        annulus_masks = annulus_aperture.to_mask(method='center')
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        phot = aperture_photometry(data, apertures)
        phot['annulus_median'] = bkg_median
        phot['aper_bkg'] = bkg_median * apertures.area
        phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
        for col in phot.colnames:
            phot[col].info.format = '%.8g'
        print(phot)
        print(phot[29])

        plt.show()
        return(e)

def plot(pattern):
    header = get_header(pattern)
    w = WCS(header)
    for filename in pattern:
        data = fits.getdata(filename, ext=0)
        mean, median, std = sigma_clipped_stats(data, sigma=3)
        daofind = DAOStarFinder(fwhm=3, threshold=5.*std)
        sources = daofind(data - median)
        x = np.array(sources['xcentroid'])
        y = np.array(sources['ycentroid'])
    world = w.wcs_pix2world(x,y,0)
    positions = SkyCoord(ra=world[0]*u.deg,dec=world[1]*u.deg).to_string('hmsdms')
    catalogData = Catalogs.query_region('119.77433151 15.39145556', radius =
                                        '0.7 deg',catalog = 'TIC')
    apertures = CircularAperture(world, r=0.5)
    dattab = Table(catalogData)
    radec = (catalogData['ra','dec','Bmag'])
    mask = radec['Bmag']<19.0
    mag_radec = radec[mask]
    
    fig = plt.figure()
    ax = plt.subplot(projection=w)
    im = ax.imshow(data, origin = 'lower', clim = (700,10000))
    fig.colorbar(im)
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    ax.scatter(mag_radec['ra']*u.deg, mag_radec['dec']*u.deg,
           facecolors='none', edgecolors='r', linewidths=0.5,
          transform=ax.get_transform('icrs'))

    plt.show()
    
    

    

def main():
    pattern = sys.argv[1:]
    circular_aperture(pattern)
    plot(pattern)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2 :
        print(" Usage:  "+sys.argv[0]+" fits file or list of fits files")
        sys.exit()
    main()     
