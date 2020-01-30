#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits

def read_header(pattern):
    for filename in pattern:
        with fits.open(filename,'update') as hdul:
            hdr = hdul[0].header
    return hdr

    x = read_header(pattern)



def update_keyword(pattern, key, valore):
    x
    print(hdr)
    if key in hdr:
        print(filename+":"+key+" "+hdr[key]+" --> "+valore)
    else:
        print(filename+": "+key+" undefined --> "+valore)
        hdr[key] = valore

   

def main():
    pattern = sys.argv[1:]

    read_header(pattern)
    if __name__ == '__main__':
        import sys
    if len(sys.argv) < 2:
        print("Usage:  "+sys.argv[0]+"<list of FITS files>")
        sys.exit()
    pattern = sys.argv[1:]
    key = sys.argv[2]  
    valore = sys.argv[3]

    update_keyword(pattern, key, valore)

    if __name__ == '__main__':
        import sys
        if len(sys.argv) < 4 :    # C'è anche lo [0] che è il nome del file :)
            print("Usage:  "+sys.argv[0]+" KEYWORD value <list of FITS files>")
        sys.exit()

    main()     
