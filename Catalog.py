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
