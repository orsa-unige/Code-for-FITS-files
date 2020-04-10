#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules
from astropy import units as u
from astropy.coordinates import SkyCoord, FK5
from astroquery.vizier import Vizier

# Our modules
from manipulate_header import get_obj, get_header, EQUINOX


def get_target_coords(target,pattern):
    ''' Get target coords by header 'OBJECT' or by input and gives his RA
        and DEC'''
    o = get_obj(target, pattern)
    header = get_header(pattern)
    eq = EQUINOX(pattern)
    d = dict();
    d['target_coords'] = SkyCoord.from_name(o['obj']).fk5.transform_to(FK5(equinox=eq))        
    d['RA_target'] = d['target_coords'].ra.degree
    d['DEC_target'] = d['target_coords'].dec.degree
    return(d)


def get_catalog(target, pattern):
    '''Gets a catalog from target's name and gives his objects coordinates
        in RA and DEC'''
    d = get_target_coords(target, pattern)
    cooEQ = d['target_coords']
    rad = 40*u.arcmin
    cat_id = 'I/284/out'
    v=Vizier(catalog=cat_id,columns=["RAJ2000","DEJ2000","Plx","RAJ2000","DEJ2000"],column_filters={'Imag':'>14'})
    v.ROW_LIMIT = -1
    tab = v.query_region(cooEQ, radius=rad, catalog = cat_id)[0]
    g = dict();
    g['RA_catalog'] = tab['RAJ2000']
    g['DEC_catalog'] = tab['DEJ2000']
    return(g)


def main():
    target = sys.argv[1]
    pattern = sys.argv[2:]
    get_catalog(target, pattern)

    
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print(" Usage:  "+sys.argv[0]+" target name, catalog name and fits file or list of fits files")
        sys.exit()
    main()     
    

