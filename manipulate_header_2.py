#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits

def update_keyword(key, valore, pattern):
    for filename in pattern:
        with fits.open(filename,'update') as hdul:
            hdr = hdul[0].header
            if key in hdr:
                print(filename+":"+key+" "+hdr[key]+" --> "+valore)
            else:
                print(filename+": "+key+" undefined --> "+valore)


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
