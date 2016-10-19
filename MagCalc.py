#################################################################
# Name:     MagCalc.py                                          #
# Author:   Yuan Qi Ni                                          #
# Version:  August 25, 2016                                     #
# Function: Developed from LimMag.py routine. Program contains  # 
#           essential functions for computing magnitude         #
#           magnitude at position of specified source in image  #
#           file using provided differential photometric        #
#           reference star catalog.                             #
#################################################################

#run sample
#python MagCalc.py -c phot -o N300-1.Q0.SN -b 'B' -p 14.263303:-37.039900 -r 2000 -fwhm 5 -vvv -n 3.0 N300-1.Q0.B.151010_1604.A.033278.005604N3646.0060.nh.crop.fits N300_1_Q0_SN.csv
#python MagCalc.py -c diff -o KSP-N300-Nova -b 'B' -p 13.789218:-37.704572 -r 1000 -fwhm 5 -n 3 -s 14.0 -vv N300-1.Q2.B.151009_0015.S.000859.005606N3754.0060.nh.fits N300-1.Q2.diff.cat

#essential modules
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse

#essential files
from Catalog import*
from Photometry import*
from Astrometry import*

#class: exception to detect invalid images
class ImageError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

#function: compute magnitude of object in image with catalog
def magnitude(filename, cat, catname, (RA,DEC), radius=500, name='object', band='V', fwhm=5.0, limsnr=0.0, satmag=14.0, verbosity=0):
    #load HDU image
    if verbosity > 0:
        print "loading hdu"
    hdulist = fits.open(filename)
    info = hdulist.info()
    image = hdulist[0].data
    header = hdulist[0].header
    #close HDU image
    hdulist.close()
    
    #print hdulist header
    if verbosity > 2:
        print "\n info \n"
        print header
    #calculate observation time in day of year
    try:
        time = isot_day(Time(header['DATE-OBS']))
    except KeyError:
        time = 0

    #construct world coordinate system
    if verbosity > 0:
        print "loading world coordinate system"
    wcs = WCS(filename)
    
    #convert position of source to pixel 
    print RA, DEC
    X, Y = wcs.all_world2pix(RA, DEC, 0)
    X, Y  = int(X), int(Y)
    if verbosity > 0:
        print "Source located at: " + str(X) + ", " + str(Y)
    
    #load photometric reference stars catalog
    if verbosity > 0:
        print "loading catalog"
    if cat == 'phot':
        ID, RA, DEC, catM, catMerr = catPhot(catname,band=band)
    elif cat == 'dprs':
        ID, RA, DEC, catM, catMerr = catDPRS(catname,band=band)
    elif cat == 'diff':
        ID, RA, DEC, catM, catMerr = catDiff(catname,band=band)
        print catM

    #convert position of catalog stars to world coordinates
    catX, catY = wcs.all_world2pix(RA, DEC, 0)
    catX, catY = catX.astype(float), catY.astype(float)
    #select catalog stars within some radius of object
    index = dist(catX,catY,X,Y) < radius
    ID, catX, catY, catM, catMerr = ID[index], catX[index], catY[index], catM[index], catMerr[index]
    #select catalog stars within edges
    index = np.logical_and(catX > 15, image.shape[1]-catX > 15)
    ID, catX, catY, catM, catMerr = ID[index], catX[index], catY[index], catM[index], catMerr[index]
    index = np.logical_and(catY > 15, image.shape[0]-catY > 15)
    ID, catX, catY, catM, catMerr = ID[index], catX[index], catY[index], catM[index], catMerr[index]
    #select unsaturated catalog stars
    index = catM > satmag
    ID, catX, catY, catM, catMerr = ID[index], catX[index], catY[index], catM[index], catMerr[index]
    #select bright enough catalog stars
    index = catM < 19.0
    ID, catX, catY, catM, catMerr = ID[index], catX[index], catY[index], catM[index], catMerr[index]
    if verbosity > 0:
        print "Selected catalog star IDs:"
        print ID
        for i in range(len(ID)):
            print RA[int(i)], DEC[int(i)]

    #aperture photometry on catalog stars
    n = len(catX)
    catI = np.zeros(n) #intensity list
    catSN = np.zeros(n) #signal to noise list
    catpopt = [] #fits to catalog stars
    catperr = [] #fit errors
    if verbosity > 0:
        print "Calculating intensity for "+str(n)+" catalog stars."

    #calculate photometry for each reference star
    for i in range(n):
        if verbosity > 0:
            print "Computing "+str(i+1)+"/"+str(n)
        #position of star in catalog
        x0, y0 = catX[i], catY[i]
        #calculate intensity and SN ratio with reduced verbosity
        PSFpopt, PSFperr, X2dof, skyN = PSFextract(image, x0, y0, fwhm=fwhm, verbosity=verbosity-1)
        I, SN = photometry(image, x0, y0, PSFpopt, PSFperr, skyN, verbosity=verbosity-1)
        #check if reference stars are valid
        if I == 0 or SN == 0 or skyN == 0:
            raise ImageError('Unable to perform photometry on reference stars.')
        #save intensity and SN ratio
        catI[i] = I
        catSN[i] = SN
        #save catalog star fits
        catpopt.append(PSFpopt)
        catperr.append(PSFperr)
    catpopt = np.array(catpopt)
    catperr = np.array(catperr)
    #calculate average psf among reference stars
    w = 1/np.square(catperr)
    catpopt = (catpopt*w).sum(0)/w.sum(0)
    catperr = np.sqrt(1/w.sum(0))
    if verbosity > 0:
        print "Average PSF [A,a,b,X0,Y0,B] = "+str(catpopt)
        print "parameter errors = "+str(catperr)
    
    #calculate photometry for source object
    if verbosity > 0:
        print "Computing magnitude of source "+name
    PSFpopt, PSFperr, X2dof, skyNo = PSFfit(image, catpopt, catperr, X, Y, verbosity=verbosity)
    Io, SNo = photometry(image, X, Y, PSFpopt, PSFperr, skyNo, verbosity=verbosity)
    #check if source is valid
    mo = float('NaN')
    mo_err = float('NaN')
    if Io != 0 and SNo != 0 and skyN != 0:
        #convert position to world coordinates
        Xo, Yo = PSFpopt[3], PSFpopt[4]
        RAo, DECo = wcs.all_pix2world(Xo, Yo, 0)
    
        #calculate magnitude of object wrt each reference star
        mi = catM - 2.512*np.log10(Io/catI)
        mi_err = np.sqrt(np.square((2.512/np.log(10))*(1/catSN))+np.square(catMerr))
        #calculate weighted mean
        w = 1/np.square(mi_err)
        mo = np.sum(mi*w)/np.sum(w)
        mo_rand = np.sqrt(1/np.sum(w))
        mo_err = np.sqrt(np.square((2.512/np.log(10))*(1/SNo)) + mo_rand**2)
    else:
        print "No source candidate detected"
        RAo, DECo = float('NaN'), float('NaN')

    if limsnr != 0 and skyNo != 0:
        #sky noise properly estimated, calculate limiting magnitude
        ru = 10.0 #recursion seed
        rl = 0.1 #recursion seed
        mlim, SNlim, expu, expd = limitingM(ru, rl, limsnr, catpopt, np.mean(catSN), skyNo, catM, catMerr, catSN, catI, verbosity)
        #return calculated magnitude, magnitude errors, and limiting magnitude
        return time, RAo, DECo, mo, mo_err, mlim
    elif limsnr != 0:
        #no sky noise estimate
        return time, RAo, DECo, mo, mo_err, float('NaN')
    else:
        #return calculated magnitude and magnitude errors
        return time, RAo, DECo, mo, mo_err

#function: recursively calculates limiting magnitude by scaling PSF to SN3.0
def limitingM(ru, rl, limsnr, popt, sno, skyN, catM, catMerr, catSN, catI, verbosity=0, level=0):
    if len(popt) == 5:
        #PSF is moffat
        A,a,b,X0,Y0 = popt[0], popt[1], popt[2], popt[3], popt[4]
        FWHM = moff_toFWHM(a,b)
        #estimate upper bound for factor shift
        f = sno/(limsnr)
        #calculate snr for many synthetic sources near limiting snr
        n = 10 #resolution of parameter space
        Aest = A/f
        A_trials = np.linspace(rl*Aest,ru*Aest,n)
        I_trials = np.zeros(len(A_trials))
        SN_trials = np.zeros(len(A_trials))
        #compute optimal aperture radius (90% source light)
        frac = 0.9
        opt_r = a*np.sqrt(np.power(1 - frac,1/(1-b)) - 1)/FWHM
        #check if wings are too large to be sensical
        opt_r = min(opt_r, 3.0)
        #do photometry over synthetic sources for each A
        if verbosity > 0:
            print "Computing "+str(n)+" synthetic sources to find mlim"
            print "PSF parameters: " + str(popt)
            print "A bound: " + str(Aest)
        for j in range(len(A_trials)):
            if verbosity > 2:
                print "Computing "+str(j+1)+"/"+str(n)
            popt_trial = [A_trials[j],a,b,X0,Y0]
            #get apertures around synthetic star
            aperture = ap_synth(D2moff, popt_trial, opt_r*FWHM)
            #photometry over synthetic aperture
            I = np.sum(aperture)
            sigma = np.sqrt(I + (skyN**2)*aperture.size)
            SN = I/sigma
            #append to trials
            I_trials[j] = I
            SN_trials[j] = SN
        #calculate I for source closest to limiting snr
        Ilim = I_trials[np.argmin(np.square(SN_trials - limsnr))]
        SNlim = SN_trials[np.argmin(np.square(SN_trials - limsnr))]
        SNUbound = max(SN_trials)
        SNLbound = min(SN_trials)
        #calculate magnitude wrt each reference star
        mlim = catM - 2.512*np.log10(Ilim/np.array(catI))
        mlim_err = np.sqrt(np.square((2.512/np.log(10))*(1/catSN))+np.square(catMerr))
        w = 1/np.square(mlim_err)
        mlim = np.sum(mlim*w)/np.sum(w)
        #output comments
        if verbosity > 0:
            print "maximum SN="+str(SNUbound)
            print "minimum SN="+str(SNLbound)
            print "mlim calculated at SN="+str(SNlim)
            print "mlim calculated at R="+str(opt_r)+"FWHM"
        #check convergence
        if abs(SNlim - limsnr) > 0.1:
            #prevent stack overflow
            if level+1 > 10:
                if verbosity > 0:
                    print "Convergence terminated to prevent stack overflow."
                return float('NaN'), float('NaN'), float('NaN'), float('NaN')
            if SNUbound > limsnr and SNLbound < limsnr:
                #take closest points and converge more
                if verbosity > 0:
                    print "mlim Monte Carlo inconvergent, refine.\n"
                ru = A_trials[SN_trials>limsnr][0]/Aest
                rl = A_trials[SN_trials<limsnr][-1]/Aest
                return limitingM(ru, rl, limsnr, popt, sno, skyN, catM, catMerr, catSN, catI, verbosity, level+1)
            elif SNUbound < limsnr:
                #need to rise more to converge
                if verbosity > 0:
                    print "mlim Monte Carlo inconvergent, upstep.\n"
                return limitingM(ru*10.0, rl, limsnr, popt, sno, skyN, catM, catMerr, catSN, catI, verbosity, level+1)
            elif SNLbound > limsnr:
                #need to drop more to converge
                if verbosity > 0:
                    print "mlim Monte Carlo inconvergent, downstep.\n"
                return limitingM(ru, rl/10.0, limsnr, popt, sno, skyN, catM, catMerr, catSN, catI, verbosity, level+1)
        else:
            #convergent
            return mlim, SNlim, ru, rl

        
#function: main function, command line execution
def main():
    #command line arguments
    parser = argparse.ArgumentParser(description="Find Photometric Magnitude")
    parser.add_argument("filename", type=str, help="fits image containing source")
    parser.add_argument("-c", "--catalog", type=str, default='phot', help="reference stars catalog convention")
    parser.add_argument("catname", type=str, help="tab separated reference stars catalog file")
    parser.add_argument("-r", "--radius", type=float, default=1000.0, help="radius in which to take reference stars")
    parser.add_argument("-o", "--source", type=str, default='object', help="target source name")
    parser.add_argument("-b", "--band", type=str, default='V', help="image filter band")
    parser.add_argument("-p", "--position", type=str, help="RA:DEC as deg:deg")
    parser.add_argument("-fwhm", type=float, default=5.0, help="image fwhm upper bound")
    parser.add_argument("-n", "--noiseSNR", type=float, default=0.0, help="limiting signal to noise threshold")
    parser.add_argument("-s", "--satMag", type=float, default=14.0, help="CCD saturation, reference star magnitude upper bound")
    parser.add_argument("-v", "--verbosity", action="count", default=0)
    args = parser.parse_args()
    
    #extract RA, DEC from position argument
    RA, DEC = [float(coord) for coord in args.position.split(':')]
    #compute position, magnitude and error
    if args.noiseSNR != 0:
        time, RA, DEC, M, Merr, Mlim = magnitude(args.filename, args.catalog, args.catname, (RA,DEC), args.radius, args.source, args.band, args.fwhm, args.noiseSNR, args.satMag, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, M, Merr, Mlim
    else:
        time, RA, DEC, M, Merr = magnitude(args.filename, args.catalog, args.catname, (RA,DEC), args.radius, args.source, args.band, args.fwhm, args.noiseSNR, args.satMag, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, M, Merr   
        
#command line execution
if __name__ == "__main__":
    #call main function
    main()
