#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def choose_hdu(filename):
    ''' Verifies if fits file is compressed or not. only .fits or .fits.fz
        format are accepted'''
    from astropy.io import fits
    finfo = fits.info(filename, output=False) # Returns a list of tuples.
    finfo_list = [item for item in finfo if 'COMPRESSED_IMAGE' in item]
    if not finfo_list:
        return finfo[0][0] # 0 if not compressed
    else:
        return finfo_list[0][0]
    

def get_header(pattern):
    '''Returns the HEADER parameter, that will be used to obtain coordinates
        and catalog stuff'''
    from astropy.io import fits
    for filename in pattern:
        which_hdu = choose_hdu(filename)
        header = fits.getheader(filename, which_hdu)
        return(header)

def MJD(pattern):
    '''Verifies if 'MJD' header is present and, if not,
        creates it by the oservation date & time, provided by headers
        'DATE-OBS','MJD-OBS' or 'MJD' '''
    from astropy.time import Time
    header = get_header(pattern)
    if 'MJD-OBS' in header or 'MJD' in header:
        pass
    elif 'JD' in header:
        jd = header['JD']
        t = Time (jd, format='jd')
        mjd = t.mjd
    elif 'DATE-OBS' in header:
        DO = header['DATE-OBS']
        t = Time(DO, format='isot')
        mjd = t.mjd
    header['MJD'] = mjd
    return(mjd)

def EQUINOX(pattern):
    ''' Header 'EQUINOX' is used to provide a catalog with the same equinox
        of the file. This function is used to verify if the header EQUINOX (if
        present), has the right format (JXXXX.X -> jyear_str). If not, it will
        be changed to be readable'''
    header = get_header(pattern)
    from astropy.time import Time
    if 'EQUINOX' in header: #check format! Must contain "J"
        equinox = str(header['EQUINOX'])
        condition = equinox.startswith('J')
        if condition is True:
            pass
        else:
            t = Time(equinox,format='jyear')
    elif 'JD' in header: 
        jd = header['JD']
        t = Time(jd,format='jd')
    elif 'MJD-OBS' in header or 'MJD' in header: 
        mjd = header['MJD-OBS'][0:22] if 'MJD-OBS' in header else header['MJD'][0:22] 
        t = Time(mjd,format='mjd')
    elif 'DATE-OBS' in header:
        DO = header['DATE-OBS']
        t = Time(DO,format='isot')

    eq = t.jyear_str
    header['EQUINOX'] = eq
    return(eq)

def get_obj(target, filename):
    ''' WCS keywords need target's name. Since it won't always be written in
        'OBJECT' header, it must be provided as second
        element of the input line'''
    header = get_header(filename)                                                                                                                                  
    o = dict();                                                                                                                                                
 #   o['x_target_pix'] = x_target   #in pix                                                                                                                                        
 #   o['y_target_pix'] = y_target   #in pix                                                                                                                                                     
    o['obj'] = header['OBJECT'] if 'OBJECT' in header else target                                                                                                                                                                                                                                                                                                                                                                                                      
    return(o)

def which_instrument(pattern):
    ''' This is used to build different rotational matrices basing on different
        instruments used for the observations

        Parameters: an angle (rad) is required to build the matrix. Default
        is pi/2'''
    import numpy as np
    header = get_header(pattern)
    if 'CCDXBIN' in header:
        plate = (0.25 * header['CCDXBIN'])/3600
    else:
        binning = float(header['CCDSUM'][0])
        plate = (0.25*binning)/3600
        
    angle = np.pi/2
    flip = 1 if header['TELESCOP'] == '0.84m' and header['INSTRUME'] == 'Mexman' else -1
    cd = np.array([[plate*np.cos(angle), plate*np.sin(angle)*flip],
                    [plate*np.sin(angle), plate*np.cos(angle)]])
    return(cd)
    


def wcs(target, pattern):
    ''' Provides WCS keywords to convert pixel coordinates of the files
        to sky coordinates. It uses the rotational matrix obtained in previous
        function (which_instrument)'''
    from astropy.wcs import WCS
    import numpy as np
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    header = get_header(pattern)
    cd = which_instrument(pattern)

    for filename in pattern:
        w = WCS(header)
        
        if 'CTYPE' in header and 'CRVAL' in header and 'CRPIX' in header and 'CD1_1 CD1_2' in header and 'CD2_1 CD2_2' in header and 'NAXIS' in header: #cerco WCS in header 
            pass
        else:
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            
            w.wcs.cd = cd

            if 'RA' in header and 'DEC' in header:
                c = SkyCoord(header['RA'], header['DEC'], unit=(u.hourangle, u.deg))
                w.wcs.crval = [c.ra.degree, c.dec.degree]
                w.wcs.crpix = [header['NAXIS1']/2, header['NAXIS2']/2]
            else:
                o = get_obj(target, pattern)
                c = SkyCoord.from_name(o['obj'])
                w.wcs.crval = [c.ra.degree, c.dec.degree]
                w.wcs.crpix = [header['NAXIS1']/2, header['NAXIS2']/2] #o, in alternativa, x_target e y_target date in input

        header.extend(w.to_header(), update=True)
        return(w)

    
def main():
    target = sys.argv[1]
    pattern = sys.argv[2:]
    wcs(target, pattern)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3 :
        print(" Usage:  "+sys.argv[0]+" target name and fits file or list of fits files")
        sys.exit()
    main()     
