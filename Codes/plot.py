#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def plot_img(target,pattern):
    ''' Plots RA and DEC of the files.
        In the second plot, the catalog will be plotted with red circles,
        the image objects with blu spheres, and the target coords with
        the yellow dot'''
    import matplotlib.pyplot as plt
    from astropy.coordinates import FK5
    from cat_and_target_coords import get_target_coords, get_catalog
    from image_sources_coords import get_sources_coords
    fig1, ax = plt.subplots()
    c = get_sources_coords(target, pattern)
    ax.plot(c['RA_img'],c['DEC_img'], 'ok', ms=5)

    fig2, ax = plt.subplots()
    ax.set_aspect('equal')
    d = get_target_coords(target, pattern)
    g = get_catalog(target, pattern)
    ax.plot(g['RA_catalog'],g['DEC_catalog'],'or',mfc='none',ms=7)
    ax.plot(c['RA_img'],c['DEC_img'],'ob',ms=4)
    ax.plot(d['RA_target'],d['DEC_target'],'oy',ms=4)
    plt.show()
    return plt.show()

def main():
    target = sys.argv[1]
    pattern = sys.argv[2:]
    plot_img(target, pattern)
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print(" Usage:  "+sys.argv[0]+" target name, catalog name and fits file or list of fits files")
        sys.exit()
    main()   
