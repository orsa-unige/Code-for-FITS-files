#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits

def update_keyword(key, valore, pattern):
    '''Creates or updates a keyword in a list of fits files.
    Provide: KEYWORD value [list of FITS files]
    '''
    
    for filename in pattern:  # Così non c'è bisogno di glob! Vedi (b).
        with fits.open(filename,'update') as hdulist:
            hdr = hdulist[0].header
            if key.upper() in hdr.keys(): # hdr ha un sacco di metodi, vedi (c).
                print(filename+": "+key+" "+hdr[key]+" --> "+valore)
            else:
                print(filename+": "+key+" undefined --> "+valore)
            hdr[key] = valore


def main():        
    key = sys.argv[1] # Keyword in maiuscolo. Vedi (a).
    valore = sys.argv[2]      # Valore
    pattern = sys.argv[3:]    # File(s). "3:" significa "dal 3 in poi".

    update_keyword(key, valore, pattern)

    
if __name__ == '__main__':
    import sys

    if len(sys.argv) < 4 :    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" KEYWORD value <list of FITS files>")
        sys.exit()

    main()

# vedi:

# (a) https://stackoverflow.com/questions/9257094/how-to-change-a-string-into-uppercase
# (b) https://stackoverflow.com/questions/29992414/unexpected-output-with-glob-glob
# (c) https://docs.astropy.org/en/stable/io/fits/api/headers.html#header
