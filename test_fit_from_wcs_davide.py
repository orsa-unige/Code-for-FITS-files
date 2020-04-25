#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.wcs import utils
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord

## L normale


def rotate(angle):
    angle_rad = np.deg2rad(float(angle))

    # Center of "pixel L"
    a = 1500
    b = 1000

    # # Coordinates of "pixel L"
    # x1,y1 = ( np.array([110, 110, 110, 110, 120, 130]),
    #           np.array([140, 130, 120, 110, 110, 110]) )

    x1,y1 = (np.array([2810.156, 2810.156,  650.236, 1820.927, 3425.779, 2750.369,
                    212.422, 1146.91 , 27.055, 2100.888,  648.149,   22.212,
                    2003.314,  727.098,  248.91 ,  409.998, 1986.931,  128.925,
                    1106.654, 1502.67 ]),
             np.array([1670.347, 1670.347,  360.325,  165.663,  900.922,  700.148,
                    1416.235, 1372.364,  398.823,  580.316,  317.952,  733.984,
                    339.024,  234.29 , 1241.608,  293.545, 1794.522, 1365.706,
                    583.135,   25.306]))

    # Coordinates of "pixel L", rotated by angle_rad and centered in a,b
    xyrot  = ( a + np.cos(angle_rad) * (x1 - a) - np.sin(angle_rad) * (y1 - b),
                   b + np.sin(angle_rad) * (x1 - a) + np.cos(angle_rad) * (y1 - b) )

    return xyrot


def wcs(xyrot):

    # Center of "sky L"
    #center = SkyCoord(ra=0.11, dec=0.11, unit=(u.deg,u.deg))

    # Coordinates of "sky L" expressed in ra,dec
    # ra,dec = ([0.11, 0.11, 0.11, 0.11, 0.12, 0.13],
    #           [0.14, 0.13, 0.12, 0.11, 0.11, 0.11] )

    center = SkyCoord(246.7368408, 43.480712949, frame = 'icrs', unit = (u.deg,u.deg))

    ra,dec = (np.array([246.75001315, 246.75001315, 246.72033646, 246.72303144,
                   246.74164072, 246.73540614, 246.73379121, 246.73761455,
                   246.7179495 , 246.73051123, 246.71970072, 246.7228646 ,
                   246.72647213, 246.7188386 , 246.7314031 , 246.71821002,
                   246.74785534, 246.73265223, 246.72579817, 246.71943263]),
             np.array([43.48690547,  43.48690547,  43.46792989,  43.48075238,
                   43.49560501,  43.48903538,  43.46045875,  43.47030776,
                   43.46132376,  43.48252763,  43.46802566,  43.46035331,
                   43.48218262,  43.46908299,  43.46131665,  43.46560591,
                   43.47791234,  43.45973025,  43.47208325,  43.47779988]))

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg,u.deg))

    # Fitting "pixel L" and "sky L"
    wcs = utils.fit_wcs_from_points(xyrot, coord, center, sip_degree=2)

    # Coordinates of "sky L" expressed in pixel
    x,y = wcs.all_pix2world(ra,dec, 0)

    cd = wcs.wcs.cd # rotation matrix
    # norm = np.linalg.norm(cd) # norm
    # cd_norm = cd / norm

    ###
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(xyrot[0],xyrot[1],'ob') # pixel L
    ax.plot(x,y,'xr')  # sky L in pixel
    #ax.plot(ra,dec,'xg')  # sky L in radec
    plt.savefig("filname.ps") # scusa, non mi va il terminale tk...
    plt.close('all')

    return cd


def main():
    '''
    Main function
    '''
    angle = sys.argv[1]
    xyrot = rotate(angle)
    wcs(xyrot)


if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 2 :    # argv[0] is the filename.
        print("Usage:  "+sys.argv[0]+" <angle>")
        sys.exit()

    main()
