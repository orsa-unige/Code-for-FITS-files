#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.time import Time # Lo metto qui e non in MJD(), sennò lo importa ogni volta che esegue la funzione.

def read_func(func, key, valore, pattern): # Serve anche la funzione da chiamare

    for filename in pattern:

        with fits.open(filename, 'update') as hdul:
            hdr = hdul[0].header
            
            if func == "update_keyword":                
                pass # ho gia- tutto come mi serve

            elif func == "update_mjd":
                data_osservazione=hdr['DATE-OBS'] # Ho tolto l'asterisco di DATE-OBS*, non so perché fosse lì.
                # dateandtime = data_osservazione # [0] Ah, forse con l'asterisco fa un array.
                dataeora= data_osservazione[0:22]
                t = Time(dataeora, format='isot', scale='utc')
                key = 'MJD-OBS' # Ecco la chiave su cui chiamare la nostra funzione update_keyword
                x = t.mjd # Perfetto. Ora ho key, valore e filename
                valore = str(t.mjd)

            if key in hdr:
                print(filename, key,":", hdr[key], "--->", valore)
            else:
                print(filename, key,":", "undefined --> ", valore)
                hdr[key] = valore # Indentazione! voglio che lo sostituisca in ogni caso, e non solo nel caso "else"
                # print(hdr)  # Non voglio vedere tutto l'header!
                # Non mi faccio dare indietro niente. 
                 
def main():

    func = sys.argv[1] # Attenzione, aggiungo un parametro per decidere che funzione usare!
    key = sys.argv[2] # cambio l-ordine dei parametri 
    valore = sys.argv[3]  # nel caso di mjd non serve... ci metto ciccio
    pattern = sys.argv[4:]
    
    read_func(func, key, valore, pattern)
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 1 :  # A te una gestione migliore del numero di argomenti!
        print("Usage:  "+sys.argv[0]+" update_keyword KEYWORD value <list of FITS files>")
        print("Usage:  "+sys.argv[0]+" update_mjd chiave_a_caso altra_chiave_a_caso <list of FITS files>")
        sys.exit()

    main()
