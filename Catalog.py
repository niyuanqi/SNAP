#################################################################
# Name:     Catalog.py                                          #
# Author:   Yuan Qi Ni                                          #
# Version:  May, 27, 2016                                       #
# Function: Program contains functions for loading photometric  #
#           reference star catalog data for differential        #
#           photometry.                                         #
#################################################################

#essential modules
import numpy as np

#function: load AAVSO catalog from online
def catAAVSO(radeg,decdeg,fovam,band,out=False):
    import os
    from Vizier import aavso, aavso_static

    #output file supplied?
    if out:
        #supplied output filename
        #check if file already present
        if os.path.isfile(out):
            #file already present
            f = open(out, 'r')
            s = f.read()
            f.close()
            ID, RA, DEC, catM, catMerr, catLines = aavso_static(s, band) 
        else:
            #fetch file from online
            ID, RA, DEC, catM, catMerr, catLines =  aavso(radeg,decdeg,fovam,band,out)
    else:
        #output file not supplied
        #fetch file from online
        ID, RA, DEC, catM, catMerr, catLines =  aavso(radeg,decdeg,fovam,band,out)
            
    #filter out nans
    index = np.invert(np.logical_or(np.isnan(catM), np.isnan(catMerr)))
    return ID[index], RA[index], DEC[index], catM[index], catMerr[index]

#function: load stable star location
def catDiff(catname, band=False):
    #load stable reference star location
    ID,RA,DEC,B,Berr,V,Verr,i,ierr = np.loadtxt(catname, unpack=True, comments=';')
    if not hasattr(ID, '__iter__'):
        ID,RA,DEC,B,Berr,V,Verr,i,ierr = np.array([ID]),np.array([RA]),np.array([DEC]),np.array([B]),np.array([Berr]),np.array([V]),np.array([Verr]),np.array([i]),np.array([ierr])

    #standard filters
    bands = {'V':0,'B':1,'I':2}
    #catalog magnitudes
    M = [V,B,i]
    Merr = [Verr,Berr,ierr]
    if band:
        #choose magnitude using band argument
        catM = [M[bands[band]]]
        catMerr = [Merr[bands[band]]]
    else:
        #give back all magnitudes
        catM = M
        catMerr = Merr

    #filter out bad values
    j=0
    while (j<len(ID)):
        if any([(val[j]=='NA' or val[j]=='-0') for val in list(catM)+list(catMerr)]):
            ID,RA,DEC = np.delete([ID,RA,DEC],j,1)
            catM,catMerr = np.delete([catM,catMerr],j,2)
        else:
            j = j+1
    RA,DEC = RA.astype(float),DEC.astype(float)
    #return catalog magnitudes
    if band:
        catM =  np.array(catM).squeeze().astype(float)
        catMerr =  np.array(catMerr).squeeze().astype(float)
        if catM.shape == ():
            catM = np.array([catM])
            catMerr = np.array([catMerr])
        return ID, RA, DEC, catM, catMerr
    else:
        catM = [cat.astype(float) for cat in catM]
        catMerr = [cat.astype(float) for cat in catMerr]
        return ID, RA, DEC, catM, catMerr

#function: load DPRS catalog (give flag 'dprs' to MagCalc)
def catDPRS(catname, band=False):
    #load values with DPRS conventions
    ID,RA,DEC,V,Verr,B,Berr,g,gerr,r,rerr,i,ierr = np.loadtxt(catname, unpack=True, skiprows=2)

    #standard DPRS filters
    bands = {'V':0,'B':1,'G':2,'R':3,'I':4}
    #catalog magnitudes
    M = [V,B,g,r,i]
    Merr = [Verr,Berr,gerr,rerr,ierr]
    if band:
        #choose magnitude using band argument
        catM = M[bands[band]]
        catMerr = Merr[bands[band]]
    else:
        #give back all magnitudes
        catM = M
        catMerr = Merr
    
    #return catalog magnitudes
    return ID, RA, DEC, catM, catMerr

#function: load phot catalog (give flag 'phot' to MagCalc)
def catPhot(catname,band=False):
    #load values with phot conventions
    RA,RAerr,DEC,DECerr,Nobs,V,Verr,B,Berr,g,gerr,r,rerr,i,ierr = np.loadtxt(catname, unpack=True, skiprows=2, delimiter=",", dtype='string')
    ID = range(len(RA))

    #standard phot filters
    bands = {'V':0,'B':1,'G':2,'R':3,'I':4}
    #catalog magnitudes
    M = [V,B,g,r,i]
    Merr = [Verr,Berr,gerr,rerr,ierr]
    if band:
        #choose magnitude using band argument
        catM = [M[bands[band]]]
        catMerr = [Merr[bands[band]]]
    else:
        #give back all magnitudes
        catM = M
        catMerr = Merr

    #filter out bad values
    j=0
    while (j<len(ID)):
        if any([(val[j]=='NA' or val[j]=='-0') for val in list(catM)+list(catMerr)]):
            ID,RA,RAerr,DEC,DECerr,Nobs = np.delete([ID,RA,RAerr,DEC,DECerr,Nobs],j,1)
            catM,catMerr = np.delete([catM,catMerr],j,2)
        else:
            j = j+1
    RA,RAerr,DEC,DECerr,Nobs = RA.astype(float),RAerr.astype(float),DEC.astype(float),DECerr.astype(float),Nobs.astype(int)
    #return catalog magnitudes
    if band:
        catM =  np.array(catM).squeeze().astype(float)
        catMerr =  np.array(catMerr).squeeze().astype(float)
        return ID, RA, DEC, catM, catMerr
    else:
        catM = [cat.astype(float) for cat in catM]
        catMerr = [cat.astype(float) for cat in catMerr]
        return ID, RA, DEC, catM, catMerr
