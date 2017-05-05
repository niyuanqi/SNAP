#################################################################
# Name:     MagCalc.py                                          #
# Author:   Yuan Qi Ni                                          #
# Version:  April 28, 2016                                      #
# Function: Program contains essential functions for computing  #
#           magnitude at position of specified source in image  #
#           file using provided differential photometric        #
#           reference star catalog.                             #
#################################################################

#run sample
#python MagCalc.py -c phot -o N300-1.Q0.SN -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 N300-1.Q0.B.151010_1604.A.033278.005604N3646.0060.nh.crop.fits N300_1_Q0_SN.csv
#python MagCalc.py -c diff -o KSP-N300-Nova -b 'B' -p 13.789218:-37.704572 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv N300-1.Q2.B.151009_0015.S.000859.005606N3754.0060.nh.fits N300-1.Q2.diff.cat
#python MagCalc.py -c dprs -o KSP-OT-1 -b 'B' -p 140.92247:-21.969278 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv N2784-7.Q1.B.150402_2125.S.015081.092205N2208.0060.nh.fits N2784-7.Q1.DPRS.cat

#essential modules
import numpy as np

#essential files
from Catalog import*
from Photometry import*
from Astrometry import*

#band definitions
bands = {'U':0, 'B':1, 'V':2, 'R':3, 'I':4}
fluxes = [1810, 4260, 3640, 3080, 2550] #Jansky

#class: exception to clarify cause of crash as inability to extract psf on image
class PSFError(Exception):
    def __init__(self, value):
        #value is error message
        self.value = value
    def __str__(self):
        #set error message as value
        return repr(self.value)

#class: exception to clarify cause of crash as invalid fits file
class FitsError(Exception):
    def __init__(self, value):
        #value is error message
        self.value = value
    def __str__(self):
        #set error message as value
        return repr(self.value)

def loadFits(filename, verbosity=0):
    """
    #################################################################
    # Desc: Load fits file for MagCalc.                             #
    # ------------------------------------------------------------- #
    # Imports: astropy.io.fits, astropy.time.Time, astropy.wcs.WCS  #
    # ------------------------------------------------------------- #
    # Input                                                         #
    # ------------------------------------------------------------- #
    #  filename: str fits filename to be opened                     #
    # verbosity; int counts verbosity level                         #
    # ------------------------------------------------------------- #
    # Output                                                        #
    # ------------------------------------------------------------- #
    # image: numpy array containing image data                      #
    #  time: float time in days since start of year YYYY            #
    #   wcs: astropy wcs object, world coordinate system on image   #
    #################################################################
    """

    from astropy.io import fits
    from astropy.time import Time
    from astropy.wcs import WCS
    
    try: #try to load image data
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
        #if verbosity > 2:
        #    print "\n info \n"
        #    print header
    except:
        raise FitsError('Unable to load fits data.')

    try: #calculate observation time in day of year
        time = isot_day(Time(header['DATE-OBS']))
    except KeyError:
        time = 0

    try: #try to load WCS
        if verbosity > 0:
            print "loading world coordinate system"
        wcs = WCS(filename)
    except:
        raise FitsError('Unable to load wcs data.')

    return image, time, wcs

def magnitude(image, wcs, cat, catname, (RAo,DECo), radius=500, name='object', band='V', fwhm=5.0, limsnr=0.0, satmag=14.0, verbosity=0):
    """
    #####################################################################
    # Desc: Compute magnitude of object in image using ref catalog.     #
    # ----------------------------------------------------------------- #
    # Imports:                                                          #
    # ----------------------------------------------------------------- #
    # Input                                                             #
    # ----------------------------------------------------------------- #
    #     image: numpy array containing image data.                     #
    #       wcs: astropy wcs object, world coordinate system on image.  #
    #       cat: str catalog type (phot, dprs, or diff from Catalog.py) #
    #   catname: str catalog name.                                      #
    # RAo, DECo: float equatorial coordinate of object in degrees.      #
    #    radius; float radius around object in which to take ref stars. #
    #      name; str name of object.                                    #
    #      band; char observational filter of data.                     #
    #      fwhm; float estimate of FWHM on image.                       #
    #    limsnr; float signal to noise ratio defining detection limit,  #
    #            if 0.0, then no detection limits are calculated.       #
    #    satmag; float magnitude below which reference stars are        #
    #            considered to be saturated and hence not used.         #
    # verbosity; int counts verbosity level.                            #
    # ----------------------------------------------------------------- #
    # Output                                                            #
    # ----------------------------------------------------------------- #
    #  RAo, DECo: float measured equatorial coordinate of source.       #
    #    Io, SNo: float measured intensity and SNR of source.           #
    # mo, mo_err: float calibrated magnitude and error of source.       #
    #       mlim; float detection limit at source position.             #
    #####################################################################
    """
    
    #convert position of source to pixel 
    Xo, Yo = wcs.all_world2pix(RAo, DECo, 0)
    Xo, Yo  = int(Xo), int(Yo)
    if verbosity > 0:
        print "Source located at: " + str(Xo) + ", " + str(Yo)

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
    index = dist(catX,catY,Xo,Yo) < radius
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
    if len(ID) == 0:
        raise ImageError('No reference stars in image.')
    if verbosity > 0:
        #output selected catalog stars
        print "Selected catalog star IDs:"
        print ID
        for i in range(len(ID)):
            print RA[int(i)], DEC[int(i)]

    if verbosity > 3:
        #plot image of catalog positions
        plt.imshow(image, cmap='Greys', vmax=0.0001*np.amax(image), vmin=0)
        plt.scatter(catX, catY)
        plt.scatter(Xo,Yo,c='r')
        plt.show()
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
            print "Magnitude "+str(catM[i])
        #position of star in catalog
        x0, y0 = catX[i], catY[i]
        #calculate intensity and SN ratio with reduced verbosity
        PSFpopt, PSFperr, X2dof, skypopt, skyN = PSFextract(image, x0, y0, fwhm=fwhm, verbosity=verbosity-1)
        I, SN = photometry(image, x0, y0, PSFpopt, skypopt, skyN, verbosity=verbosity-1)
        #check if reference stars are valid
        if I == 0 or SN == 0 or skyN == 0:
            raise PSFError('Unable to perform photometry on reference stars.')
        #save intensity and SN ratio
        catI[i] = I
        catSN[i] = SN
        #save catalog star fits
        catpopt.append(PSFpopt)
        catperr.append(PSFperr)
    catpopt = np.array(catpopt)
    catperr = np.array(catperr)
    """
    #diagnostic for SNR, N, I calculation routine
    #checks for wrong correlation between intensity and noise
    import matplotlib.pyplot as plt
    plt.title("SNR calculation for various reference stars")
    plt.scatter(catM, np.log(catSNr), c='b', label="Calculated using fit residuals")
    plt.scatter(catM, np.log(catSN), c='r', label="Calculated using intensity")
    plt.ylabel("log SNR")
    plt.xlabel("Mag (~log I)")
    plt.legend()
    plt.show()
    plt.title("Noise under various reference stars")
    plt.scatter(catM, np.log(catI/catSNr), c='b', label="Calculated using fit residuals")
    plt.scatter(catM, np.log(catI/catSN), c='r', label="Calculated using intensity")
    plt.ylabel("log Noise")
    plt.xlabel("Mag (~log I)")
    plt.legend()
    plt.show()
    plt.title("Intensity of various reference stars")
    plt.scatter(catM, np.log(catI), c='b')
    plt.ylabel("log I")
    plt.xlabel("Mag (catalog)")
    plt.legend()
    plt.show()
    """

    #calculate average psf among reference stars
    w = 1/np.square(catperr)
    catpopt = (catpopt*w).sum(0)/w.sum(0)
    catperr = np.sqrt(1/w.sum(0))
    if verbosity > 0:
        print "Average PSF [A,a,b,X0,Y0,B] = "+str(catpopt)
        print "parameter errors = "+str(catperr)
        print "\nMean SN of reference stars:",np.mean(catSN)
        print ""

    #calculate photometry for source object
    if verbosity > 0:
        print "Computing magnitude of source "+name
    PSFpopt, PSFperr, X2dof, skypopto, skyNo = PSFfit(image, catpopt, catperr, Xo, Yo, verbosity=verbosity)
    Io, SNo = photometry(image, Xo, Yo, PSFpopt, skypopto, skyNo, verbosity=verbosity)
    #check if source is valid
    if Io != 0 and SNo != 0 and skyNo != 0:
        #calculate relative flux of object wrt each reference star
        Ir = fluxes[bands[band]]*np.power(10,-catM/2.512)*Io/catI
        print Ir
        Io_err = Io/SNo
        catI_err = catI/catSN
        Ir_err = Ir*np.sqrt(np.square(Io_err/Io)+np.square(catI_err/catI)+np.square(np.log(10)*catMerr/2.512))
        print Ir_err

        #calculate weighted mean
        w = 1/np.square(Ir_err)
        I = np.sum(Ir*w)/np.sum(w)
        I_rand = np.sqrt(np.sum(w))
        I_err = np.sqrt(Io_err**2 + I_rand**2)
        SN = I/I_err

        print I, I_err, SN
    else:
        #no valid source
        Io, SNo = float('NaN'), float('NaN')
    #try to compute magnitude if source is present
    if Io != float('NaN') and Io > 0 and skyNo != 0:
        #convert position to world coordinates
        Xp, Yp = PSFpopt[3], PSFpopt[4]
        RAo, DECo = wcs.all_pix2world(Xp, Yp, 0)
    
        #calculate magnitude of object wrt each reference star
        mi = catM - 2.512*np.log10(Io/catI)
        mi_err = np.sqrt(np.square((2.512/np.log(10))*(1/catSN))+np.square(catMerr))
        #calculate weighted mean
        w = 1/np.square(mi_err)
        mo = np.sum(mi*w)/np.sum(w)
        mo_rand = np.sqrt(1/np.sum(w))
        mo_err = np.sqrt(np.square((2.512/np.log(10))*(1/SNo)) + mo_rand**2)
    else:
        mo, mo_err = float('NaN'), float('NaN')
        RAo, DECo = wcs.all_pix2world(Xo, Yo, 0)

    if limsnr != 0 and skyNo != 0 and PSFverify(catpopt, catpopt[-2], catpopt[-1]):
        #sky noise properly estimated, calculate limiting magnitude
        ru = 10.0 #recursion seed
        rl = 0.1 #recursion seed
        mlim, SNlim, expu, expd = limitingM(ru, rl, limsnr, catpopt, np.mean(catSN), skyNo, catM, catMerr, catSN, catI, verbosity)
        #return calculated magnitude, magnitude errors, and limiting magnitude
        return RAo, DECo, I, SN, mo, mo_err, mlim
    elif limsnr != 0:
        #no sky noise estimate
        return RAo, DECo, I, SN, mo, mo_err, float('NaN')
    else:
        #return calculated magnitude and magnitude errors
        return RAo, DECo, I, SN, mo, mo_err

#function: recursively calculates limiting magnitude by scaling PSF to SN3.0
def limitingM(ru, rl, limsnr, popt, sno, skyN, catM, catMerr, catSN, catI, verbosity=0, level=0):
    """
    ##########################################################################
    # Desc: Recursively calculates limiting magnitude by scaling PSF to SNR. #
    # ---------------------------------------------------------------------- #
    # Imports:                                                               #
    # ---------------------------------------------------------------------- #
    # Input                                                                  #
    # ---------------------------------------------------------------------- #
    #    ru, rl: Upper and lower ratio of fixed scaling estimate.            #
    #            Estimate A_est obtained from popt, and Monte Carlo          #
    #            values taken between ru*A_est and rl*A_est.                 #
    #    limsnr: float signal to noise ratio defining detection limit.       #
    #      popt: iterable floats (len 5) containing PSF on image.            #
    #       sno: average signal to noise of ref stars used to get popt.      #
    #      skyN: sky noise in annulus at source position.                    #
    #      catM: list of magnitudes of reference stars.                      #
    #   catMerr: list of magnitude errors of reference stars.                #
    #     catSN: list of SNRs of reference stars.                            #
    #      catI: list of intensities of reference stars.                     #
    #     level; recursion level of computation, min 0, max 10.              #
    # verbosity; int counts verbosity level.                                 #
    # ---------------------------------------------------------------------- #
    # Output                                                                 #
    # ---------------------------------------------------------------------- #
    #  RAo, DECo: float measured equatorial coordinate of source.            #
    #    Io, SNo: float measured intensity and SNR of source.                #
    # mo, mo_err: float calibrated magnitude and error of source.            #
    #       mlim; float detection limit at source position.                  #
    ##########################################################################
    """
    
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
        #check if psf is too small
        opt_r = max(opt_r, 1.0/FWHM)
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
            #sigma = np.sqrt(I + (skyN**2)*aperture.size)
            #at SN <= 5, noise dominated, I/(skyN**2)*aperture.size < 0.1
            sigma = np.sqrt(aperture.size*skyN**2)
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
        
#command line execution
if __name__ == "__main__":

    import argparse
    
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
    parser.add_argument("-n", "--noiseSNR", type=float, default=0.0, help="signal to noise at detection limit")
    parser.add_argument("-s", "--satMag", type=float, default=14.0, help="CCD saturation, reference star magnitude upper bound")
    parser.add_argument("-v", "--verbosity", action="count", default=0)
    args = parser.parse_args()
    
    #extract RA, DEC from position argument
    RA, DEC = [float(coord) for coord in args.position.split(':')]

    #load fits file, get relevant data
    image, time, wcs = loadFits(args.filename, args.verbosity)
    
    #compute position, magnitude and error
    if args.noiseSNR != 0:
        RA, DEC, I, SN, M, Merr, Mlim = magnitude(image, wcs, args.catalog, args.catname, (RA,DEC), args.radius, args.source, args.band, args.fwhm, args.noiseSNR, args.satMag, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, I, SN, M, Merr, Mlim
    else:
        RA, DEC, I, SN, M, Merr = magnitude(image, wcs, args.catalog, args.catname, (RA,DEC), args.radius, args.source, args.band, args.fwhm, args.noiseSNR, args.satMag, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, I, SN, M, Merr
