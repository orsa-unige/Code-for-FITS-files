#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.wcs import utils
from astropy.wcs import WCS
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits

## L normale
xy1 = (np.array([11, 11, 11, 11, 12, 13]),np.array([14, 13, 12, 11, 11, 11]))
coord = SkyCoord(ra=[0.11, 0.11, 0.11, 0.11, 0.12, 0.13],
                 dec=[0.14, 0.13, 0.12, 0.11, 0.11, 0.11],
                 unit=(u.deg,u.deg))
ra = coord.ra.deg
dec = coord.dec.deg

center = SkyCoord(ra=0.11,dec=0.11,unit=(u.deg,u.deg))

a = 11
b = 11

x1 = xy1[0]
y1 = xy1[1]
    
def rotate():
    d = dict();
    angle = (0)
    angle_rad = np.deg2rad(angle)
    d['x2'] = a + np.cos(angle_rad) * (x1 - a) - np.sin(angle_rad) * (y1 - b)
    d['y2'] = b + np.sin(angle_rad) * (x1 - a) + np.cos(angle_rad) * (y1 - b)
    print(np.cos(angle_rad),np.sin(angle_rad))
    return d

def wcs():
    d = rotate()
    xyrot = (np.array(d['x2']),np.array(d['y2']))
    wcs = utils.fit_wcs_from_points(xyrot, coord, center)
    print(wcs)
    header = fits.Header()
    header.extend(wcs.to_header(), update=True)
    w = WCS(header)
    e = dict();
    world = w.wcs_pix2world(d['x2'],d['y2'],0)
    e['xr'] = world[0] 
    e['yr'] = world[1]
    return(e)


def plot():
    d = rotate()
    e = wcs()

    fig1,ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlim(4,15)
    ax.set_ylim(4,15)
    ax.plot(d['x2'],d['y2'],'ok',ms=5)

    fig2,ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlim(0.075,0.150)
    ax.set_ylim(0.075,0.150)
    ax.plot(e['xr'],e['yr'],'ok',ms=5)
    ax.plot(ra,dec,'or',ms=1)
    
    plt.show()



def main():
    plot()
    
if __name__ == '__main__':

    main()     

