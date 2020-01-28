#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.time import Time

def update_keyword(key, valore, pattern):

    for filename in pattern:
        with fits.open(filename,'update') as hdul:
            hdr = hdul[0].header
            if key in hdr.keys():
                print(filename+":"+key+" "+hdr[key]+" --> "+valore)
            else:
                print(filename+": "+key+" undefined --> "+valore)
            hdr[key] = valore
            data_osservazione = hdr['DATE-OBS*']
            dateandtime = data_osservazione[0]
            dataeora = dateandtime[0:22]
            t = Time(dataeora, format='isot', scale='utc')
            mjd = t.mjd
            if "MJD-OBS" in hdr.keys():
                print("Already existing")
            else:
                print("Creating now")
                hdr['MJD-OBS'] = mjd
            print(mjd)


def main():
    key = sys.argv[1]
    valore = sys.argv[2]
    pattern = sys.argv[3:]

    update_keyword(key, valore, pattern)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 4:
        print("Usage:  "+sys.argv[0]+"KEYWORD value <list of FITS files>")
        sys.exit()

    main()
    

