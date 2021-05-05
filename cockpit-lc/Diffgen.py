#################################################################
# Name:     DiffFits.py                                         #
# Author:   Yuan Qi Ni                                          #
# Date:     July 14, 2017                                       #
# Function: Program uses DiffIm routine to subtract images.     #
#           Update /raw files and ObjData.py before running.    #
#################################################################

#essential modules
import numpy as np
import os
from glob import glob
import math
from astropy.io import fits

#essential files from SNAP
from SNAP.DiffIm import make_diff_image
from SNAP.Analysis.LCRoutines import*
from SNAP.MagCalc import*
from SNAP.Catalog import*
from SNAP.Photometry import*
from SNAP.PSFlib import*
#essential imports
from ContextManager import cd
#essential data
from ObjData import *

#reference files
bands = ['B','V','I']
bindex = {'B':0, 'V':1, 'I':2}
refs = ['../ref/'+Brefname, '../ref/'+Vrefname, '../ref/'+Irefname]

#current working directory
wd = os.getcwd()
#make directory for diff images
with cd(wd+"/../"):
    if not os.path.isdir("diff"): os.mkdir('diff')
    if not os.path.isdir("conv"): os.mkdir('conv')

#for each band
for i in range(len(bands)):
    #get all band files
    files = sorted(glob('../raw/'+prefix+bands[i]+'*.fits'))
    N = len(files)
    for n, filename in enumerate(files):
        print ""
        print "Performing image subtraction on file "+str(n+1)+"/"+str(N)
        print filename
        #output filename
        diffname = '.'.join(filename.split('.')[:-1])+".diff.fits"
        diffname = '../diff/'+'/'.join(diffname.split('/')[2:])
        convname = '.'.join(filename.split('.')[:-1])+".conv.fits"
        convname = '../conv/'+'/'.join(convname.split('/')[2:])
        #other parameters
        fo = filename.split('/')[2]
        fo = '.'.join(fo.split('.')[2:5])
        band = fo[0]

        #subtract if not already subtracted
        if os.path.exists(diffname) and os.path.exists(convname):
            print "Already subtracted "+filename
        else:
            print "Subtracting "+filename
            #retrieve parameters from image
            Mtest = True
            try: #try to load image
                image, to, wcs, hdr = loadFits(filename, year=year, getwcs=True, gethdr=True, verbosity=0)
                
            except FitsError:
                #image critically failed to load
                Mtest = False
                so = "FITS_ERROR"
                to = 0
                print "Critical error loading image!"
            if Mtest:
                try:
                    print "Extracting psf"
                    PSF, PSFerr, Med, Noise = magnitude(image, image, wcs, cattype, catname, (ra,dec), radius=size, psf=1, name=name, band=band, fwhm=5.0, limsnr=SNRnoise, satmag=satlvl, refmag=rellvl, fitsky=True, satpix=satpix, verbosity=0, diagnosis=True)
                    #print image fwhm
                    fwhm = np.mean(E2moff_toFWHM(*PSF[:-1]))
                    print "Image fwhm", fwhm
                    if fwhm == 0:
                        raise PSFError('Unable to perform photometry on reference stars.')
                    print ""
                    print "Performing subtraction, generating files"
                    print diffname, convname
                    make_diff_image(filename, refs[i], diffname, convname,
                                    tmp_fwhm=ref_fwhms[i], src_fwhm=fwhm,
                                    tmpdir="DITemp"+str(n))
                
                except PSFError:
                    Mtest = False
                    print "PSF can't be extracted!"
                    print "Not performing subtraction."
                except: #General catastrophic failure
                    Mtest = False
                    print "Unknown catastrophic failure!"
                    print "Not performing subtraction."
                

