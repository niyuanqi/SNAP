#################################################################
# Name:     Vizier.py                                           #
# Author:   Yuan Qi Ni                                          #
# Version:  February 9, 2017                                    #
# Function: Program contains routines for obtaining differential#
#           photometric catalogs by querying Vizier database.   #
#           Currently supports AAVSO-APASS-DR9 and USNO-B1.     #
#################################################################

import numpy as np

#Query AAVSO APASS DR9 catalog
def aavso(radeg,decdeg,fovam,band,out=False):

    import urllib as url
    
    str1 = 'http://webviz.u-strasbg.fr/viz-bin/asu-tsv/?-source=II/336' 
    str2 = '&-c.ra={:4.6f}&-c.dec={:4.6f}&-c.bm={:4.7f}/{:4.7f}&-out.max=unlimited'.format(radeg,decdeg,fovam,fovam)
     
    # Make sure str2 does not have any spaces or carriage returns/line feeds when you # cut and paste into your code  
    sr = str1+str2 
    f = url.urlopen(sr) 
    
    # Read from the object, storing the page's contents in 's'. 
    s = f.read() 
    f.close()
    
    #write text file
    if out:
        #outname = 'AAVSO-APASS-DR9_ra{:4.6f}dec{:4.6f}bm{:3.1f}x{:3.1f}.cat'.format(radeg,decdeg,fovam,fovam)
        outname = out
        f = open(outname, 'w')
        f.write(s)
        f.close()
        
    return aavso_static(s, band)

#read static Vizier AAVSO APASS DR9 catalog
def aavso_static(lines, band):
    RAcol = 1
    DEcol = 2
    bands = {'B-V':8, 'B': 12, 'V': 10, 'i': 18, 'r': 16}
    banderrs = {'B-V':9, 'B': 13, 'V': 11, 'i': 19, 'r': 17}
    null = '      '

    #parse text
    sl = lines.splitlines()
    sl = sl[53:-1] # get rid of header
    name = np.array([]) 
    rad = np.array([]) # RA in degrees 
    ded = np.array([]) # DEC in degrees 
    rmag = np.array([]) # rmag
    rmagerr = np.array([]) # rmag error

    for k in sl:
        kw = k.split('\t')
        if kw[0] != '':  
            name = np.append(name,kw[0])
            rad = np.append(rad,float(kw[1])) 
            ded = np.append(ded,float(kw[2]))
            # deal with case where no mag is reported
            if (kw[bands[band]] != null) and (kw[banderrs[band]]) != null:
                rmag = np.append(rmag,float(kw[bands[band]]))
                rmagerr = np.append(rmagerr,float(kw[banderrs[band]]))
            else:
                rmag = np.append(rmag,np.nan)
                rmagerr = np.append(rmagerr,np.nan)
        
    return name,rad,ded,rmag,rmagerr,lines

#query Vizier USNO-B1 catalog
def usnoB(radeg,decdeg,fovam,band,out=False): # RA/Dec in decimal degrees/J2000.0 FOV in arc min. 

    import urllib as url
    
    str1 = 'http://webviz.u-strasbg.fr/viz-bin/asu-tsv/?-source=USNO-B1' 
    str2 = '&-c.ra={:4.6f}&-c.dec={:4.6f}&-c.bm={:4.7f}/{:4.7f}&-out.max=unlimited'.format(radeg,decdeg,fovam,fovam)

    RAcol = 1
    DEcol = 2
    bands = {'B': 12}
    #banderrs = {'B': 13}
    null = '     '
    empt = ''
     
    # Make sure str2 does not have any spaces or carriage returns/line feeds when you # cut and paste into your code  
    sr = str1+str2 
    f = url.urlopen(sr)
    
    # Read from the object, storing the page's contents in 's'. 
    s = f.read() 
    f.close()

    #write text file
    if out:
        outname = 'USNO-B1_ra{:4.6f}dec{:4.6f}bm{:4.1f}x{:4.1f}.cat'.format(radeg,decdeg,fovam,fovam)
        f = open(outname, 'w')
        f.write(s)
        f.close()
    
    sl = s.splitlines() 
    sl = sl[46:] # get rid of header 
    name = np.array([]) 
    rad = np.array([]) # RA in degrees 
    ded = np.array([]) # DEC in degrees 
    rmag = np.array([]) # rmage 
    
    for k in sl:
        kw = k.split('\t')
        if kw[0] != '':  
            name = np.append(name,kw[0])
            rad = np.append(rad,float(kw[1])) 
            ded = np.append(ded,float(kw[2]))
            # deal with case where no mag is reported
            #if (kw[bands[band]] != null) and (kw[banderrs[band]]) != null:
            if (kw[bands[band]] != null) and (kw[bands[band]] != empt):
                print "Mag", kw[bands[band]], "len", len(kw[bands[band]])
                rmag = np.append(rmag,float(kw[bands[band]]))
                #rmagerr = np.append(rmag,float(kw[banderrs[band]]))
            else:
                rmag = np.append(rmag,np.nan)
                #rmagerr = np.append(rmagerr,np.nan)
        
    return name,rad,ded,rmag,s
