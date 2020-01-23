#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from astropy.io import fits

def main():

    if len(sys.argv) < 4 :    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" KEYWORD value <list of FITS files>")
        sys.exit()
        
    key = sys.argv[1].upper() # Keyword in maiuscolo. Vedi (a).
    valore = sys.argv[2]      # Valore
    pattern = sys.argv[3:]    # File(s). "3:" significa "dal 3 in poi".
    
    for filename in pattern:  # Così non c'è bisogno di glob! Vedi (b).
        with fits.open(filename,'update') as hdulist:
            hdr = hdulist[0].header
            if key in hdr.keys(): # hdr ha un sacco di metodi, vedi (c).
                print(filename+": "+key+" "+hdr[key]+" --> "+valore)
            else:
                print(filename+": "+key+" undefined --> "+valore)
            hdr[key] = valore


if __name__ == "__main__":
    main()

# vedi:

# (a) https://stackoverflow.com/questions/9257094/how-to-change-a-string-into-uppercase
# (b) https://stackoverflow.com/questions/29992414/unexpected-output-with-glob-glob
# (c) https://docs.astropy.org/en/stable/io/fits/api/headers.html#header
