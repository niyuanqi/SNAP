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
#python MagCalc.py -c phot -o N300-1.Q0.SN -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 -f 16.0 --fit_sky N300-1.Q0.B.151010_1604.A.033278.005604N3646.0060.nh.crop.fits N300_1_Q0_SN.csv
#python MagCalc.py -c diff -o KSP-N300-Nova -b 'B' -p 13.789218:-37.704572 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv --fit_sky N300-1.Q2.B.151009_0015.S.000859.005606N3754.0060.nh.fits N300-1.Q2.diff.cat
#python MagCalc.py -c dprs -o KSP-OT-1 -b 'B' -p 140.92247:-21.969278 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv --fit_sky N2784-7.Q1.B.150402_2125.S.015081.092205N2208.0060.nh.fits N2784-7.Q1.DPRS.cat

#essential modules
import numpy as np

#band definitions
bands = {'U':0, 'B':1, 'V':2, 'R':3, 'I':4}
fluxes = [1810, 4260, 3640, 3080, 2550] #Jansky
fluxes = [flux*1e6 for flux in fluxes]

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

def loadFits(filename, year=2016, getwcs=False, verbosity=0):
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
    #    getwcs; boolean whether to return wcs                      #
    #      year; int year to measure time to                        #
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

    from Astrometry import*
    
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
        time = isot_day(Time(header['DATE-OBS']), int(year))
    except KeyError:
        time = 0

    if getwcs:
        try: #try to load WCS
            if verbosity > 0:
                print "loading world coordinate system"
            wcs = WCS(filename)
        except:
            raise FitsError('Unable to load wcs data.')

        return image, time, wcs
    else:
        return image, time

def magnitude(image, catimage, wcs, cat, catname, (RAo,DECo), radius=500, aperture=None, name='object', band='V', fwhm=5.0, limsnr=0.0, satmag=14.0, refmag=19.0, fitsky=True, verbosity=0):
    """
    #####################################################################
    # Desc: Compute magnitude of object in image using ref catalog.     #
    # ----------------------------------------------------------------- #
    # Imports:                                                          #
    # ----------------------------------------------------------------- #
    # Input                                                             #
    # ----------------------------------------------------------------- #
    #     image: numpy array containing image data on which to measure  #
    #            source photometry.                                     #
    #  catimage: numpy array containing image data on which to measure  #
    #            reference star photometry.                             #
    #       wcs: astropy wcs object, world coordinate system on image.  #
    #       cat: str catalog type (phot, dprs, or diff from Catalog.py) #
    #   catname: str catalog name.                                      #
    # RAo, DECo: float equatorial coordinate of object in degrees.      #
    #    radius; float radius around object in which to take ref stars. #
    #  aperture; float aperture around object in which to integrate     #
    #            light. If None; PSF is integrated. If not a positive   #
    #            number; then PSF is used to get Kron aperture.         #
    #      name; str name of object.                                    #
    #      band; char observational filter of data.                     #
    #      fwhm; float estimate of FWHM on image.                       #
    #    limsnr; float signal to noise ratio defining detection limit,  #
    #            if 0.0, then no detection limits are calculated.       #
    #    satmag; float magnitude below which reference stars are        #
    #            considered to be saturated and hence not used.         #
    #    fitsky; boolean, if True; fit for planar sky around source to  #
    #            be subtracted from image before fitting/integrating.   #
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

    #essential imports
    import matplotlib.pyplot as plt
    
    #essential files
    from Catalog import*
    from Photometry import*
    from Astrometry import*
    
    #convert position of source to pixel 
    Xo, Yo = wcs.all_world2pix(RAo, DECo, 0)
    Xo, Yo = int(Xo), int(Yo)
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
    elif cat == 'aavso':
        fovam = 2.0*radius*0.4/60.0 #arcmin radius in KMT scaling
        ID, RA, DEC, catM, catMerr = catAAVSO(RAo,DECo,fovam,band,out=catname)

    #convert position of catalog stars to world coordinates
    catX, catY = wcs.all_world2pix(RA, DEC, 0)
    catX, catY = catX.astype(float), catY.astype(float)
    #select catalog stars within some radius of object
    index = dist(catX,catY,Xo,Yo) < radius
    #select catalog stars within edges
    index = np.logical_and(index, np.logical_and(catX > 15, catimage.shape[1]-catX > 15))
    index = np.logical_and(index, np.logical_and(catY > 15, catimage.shape[0]-catY > 15))
    #select unsaturated catalog stars
    index = np.logical_and(index, catM > satmag)
    #select bright enough catalog stars
    index = np.logical_and(index, catM < refmag)
    #crop values to mask
    ID, catX, catY, catRA, catDEC, catM, catMerr = ID[index], catX[index], catY[index], RA[index], DEC[index], catM[index], catMerr[index]
    if len(ID) == 0:
        raise PSFError('No reference stars in image.')
    if verbosity > 0:
        #output selected catalog stars
        print "Selected catalog star IDs:"
        for i in range(len(ID)):
            print ID[int(i)], catX[int(i)], catY[int(i)]
            print RA[int(i)], DEC[int(i)], catM[int(i)], catMerr[int(i)]

    if verbosity > 3:
        #plot image of catalog positions
        plt.imshow(catimage, cmap='Greys', vmax=0.0001*np.amax(catimage), vmin=0)
        plt.scatter(catX, catY)
        plt.scatter(Xo,Yo,c='r')
        plt.show()
    #aperture photometry on catalog stars
    n = len(catX)
    catI = np.zeros(n) #intensity list
    catSN = np.zeros(n) #signal to noise list
    catPSF = [] #PSF fits to catalog stars
    catPSFerr = [] #fit errors
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
        PSFpopt, PSFperr, X2dof, skypopt, skyN = PSFextract(catimage, x0, y0, fwhm=fwhm, fitsky=fitsky, verbosity=verbosity-1)
        PSF, PSFerr = [PSFpopt[1], PSFpopt[2]], [PSFperr[1], PSFperr[2]]
        #check preferred intensity calculation method
        if aperture is None:
            #integrate PSF directly
            I, SN = PSF_photometry(catimage, x0, y0, PSFpopt, skypopt, skyN, verbosity=verbosity-1)
        else:
            #perform aperture photometry
            if aperture <= 0:
                #use FWHM of catPSF to define Kron aperture
                I, SN = Ap_photometry(catimage, x0, y0, skypopt, skyN, PSF=PSF, fitsky=fitsky, verbosity=verbosity-1)
            else:
                #use aperture given directly
                I, SN = Ap_photometry(catimage, x0, y0, skypopt, skyN, radius=aperture, fitsky=fitsky, verbosity=verbosity-1)
        #check if reference stars are valid
        if I == 0 or SN == 0 or skyN == 0:
            raise PSFError('Unable to perform photometry on reference stars.')
        #save intensity and SN ratio
        catI[i] = I
        catSN[i] = SN
        #save catalog star fits
        catPSF.append(PSF)
        catPSFerr.append(PSFerr)
    catPSF = np.array(catPSF)
    catPSFerr = np.array(catPSFerr)

    """
    if verbosity > 1:
        #diagnostic for SNR, N, I calculation routine
        #checks for wrong correlation between intensity and noise
        import matplotlib.pyplot as plt
        plt.title("SNR calculation for various reference stars")
        plt.scatter(catM, np.log(catSN), c='r')
        plt.plot(catM, max(np.log(catSN))-(catM-min(catM))/2)
        plt.ylabel("log SNR")
        plt.xlabel("Mag (~log I)")
        plt.legend()
        plt.show()
        plt.title("Noise under various reference stars")
        plt.scatter(catM, np.log(catI/catSN), c='r')
        plt.plot(catM, max(np.log(catI/catSN))-(catM-min(catM))/2)
        plt.ylabel("log Noise")
        plt.xlabel("Mag (~log I)")
        plt.legend()
        plt.show()
        plt.title("Intensity of various reference stars")
        plt.scatter(catM, np.log(catI), c='r')
        plt.plot(catM, max(np.log(catI))-(catM-min(catM)))
        plt.ylabel("log I")
        plt.xlabel("Mag (catalog)")
        plt.legend()
        plt.show()
    """
        
    #calculate average psf among reference stars
    w = 1/np.square(catPSF)
    catPSF = (catPSF*w).sum(0)/w.sum(0)
    catPSFerr = np.sqrt(1/w.sum(0))
    if verbosity > 0:
        print "Average PSF [a,b] = "+str(catPSF)
        print "Average FWHM = "+str(moff_toFWHM(*catPSF))
        print "parameter errors = "+str(catPSFerr)
        print "\nMean SN of reference stars:",np.mean(catSN)
        print ""

    #calculate photometry for source object
    if verbosity > 0:
        print "Computing magnitude of source "+name

    #extract PSF to as great a degree as possible from source
    PSFpopt, PSFperr, X2dof, skypopto, skyNo = PSFfit(image, catPSF, catPSFerr, Xo, Yo, fitsky=fitsky, verbosity=verbosity)
    #check preferred intensity calculation method
    if aperture is None:
        #integrate PSF directly
        Io, SNo = PSF_photometry(image, Xo, Yo, PSFpopt, skypopto, skyNo, verbosity=verbosity)
    else:
        #perform aperture photometry
        if aperture <= 0:
            #use FWHM of catPSF to define Kron aperture
            Io, SNo = Ap_photometry(image, Xo, Yo, skypopto, skyNo, PSF=catPSF, fitsky=fitsky, verbosity=verbosity)
        else:
            #use aperture given directly
            Io, SNo = Ap_photometry(image, Xo, Yo, skypopto, skyNo, radius=aperture, fitsky=fitsky, verbosity=verbosity)
        
    #check if source is valid
    if Io != 0 and SNo != 0 and skyNo != 0:
        #calculate relative flux of object wrt each reference star
        Ir = fluxes[bands[band]]*np.power(10,-catM/2.512)*Io/catI
        Io_err = Io/SNo
        catI_err = catI/catSN
        Ir_err = Ir*np.sqrt(np.square(1/SNo)+np.square(1/catSN)+np.square(np.log(10)*catMerr/2.512))

        #calculate weighted mean
        w = 1/np.square(Ir_err)
        I = np.sum(Ir*w)/np.sum(w)
        I_rand = np.sqrt(1/np.sum(w))
        I_err = np.sqrt((I*Io_err/Io)**2 + I_rand**2)

        if verbosity > 0:
            print "Contribution of intrinsic error:", Io_err/Io
            print "Contribution of ref star scatter:", I_rand/I
            
        #SN = I/I_err
    else:
        #no valid source
        I, SNo = float('NaN'), float('NaN')
    #try to compute magnitude if source is present
    if I != float('NaN') and I > 0 and skyNo != 0:
        #convert position to world coordinates
        Xp, Yp = PSFpopt[3], PSFpopt[4]
        RAo, DECo = wcs.all_pix2world(Xp, Yp, 0)
    
        #calculate magnitude of object wrt each reference star
        #mi = catM - 2.512*np.log10(Io/catI)
        #mi_err = np.sqrt(np.square((2.512/np.log(10))*(1/catSN))+np.square(catMerr))
        #calculate weighted mean
        #w = 1/np.square(mi_err)
        #mo = np.sum(mi*w)/np.sum(w)
        #mo_rand = np.sqrt(1/np.sum(w))
        #mo_err = np.sqrt(np.square((2.512/np.log(10))*(1/SNo)) + mo_rand**2)

        mo = -2.512*np.log10(I/fluxes[bands[band]])
        mo_err = (2.512/np.log(10))*(I_err/I)
    else:
        mo, mo_err = float('NaN'), float('NaN')
        RAo, DECo = wcs.all_pix2world(Xo, Yo, 0)

    if limsnr != 0 and skyNo != 0 and PSFverify(catPSF):
        #sky noise properly estimated, calculate limiting magnitude
        ru = 10.0 #recursion seed
        rl = 0.1 #recursion seed
        mlim, SNlim, expu, expd = limitingM(ru, rl, limsnr, catPSF, np.mean(catSN), skyNo, catM, catMerr, catSN, catI, verbosity)
        #return calculated magnitude, magnitude errors, and limiting magnitude
        return RAo, DECo, I, SNo, mo, mo_err, mlim
    elif limsnr != 0:
        #no sky noise estimate
        return RAo, DECo, I, SNo, mo, mo_err, float('NaN')
    else:
        #return calculated magnitude and magnitude errors
        return RAo, DECo, I, SNo, mo, mo_err

#function: recursively calculates limiting magnitude by scaling PSF to SN3.0
def limitingM(ru, rl, limsnr, PSF, sno, skyN, catM, catMerr, catSN, catI, verbosity=0, level=0):
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
    #       PSF: iterable floats (len 2) containing PSF on image.            #
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

    from Photometry import*
    
    if len(PSF) == 2:
        #PSF is moffat
        a,b = PSF[0], PSF[1]
        A = np.mean(catSN)*(b-1)/np.pi/np.square(a)
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
            print "PSF parameters: " + str(PSF)
            print "A bound: " + str(Aest)
        for j in range(len(A_trials)):
            if verbosity > 2:
                print "Computing "+str(j+1)+"/"+str(n)
            popt_trial = [A_trials[j],a,b,0,0]
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
                return limitingM(ru, rl, limsnr, PSF, sno, skyN, catM, catMerr, catSN, catI, verbosity, level+1)
            elif SNUbound < limsnr:
                #need to rise more to converge
                if verbosity > 0:
                    print "mlim Monte Carlo inconvergent, upstep.\n"
                return limitingM(ru*10.0, rl, limsnr, PSF, sno, skyN, catM, catMerr, catSN, catI, verbosity, level+1)
            elif SNLbound > limsnr:
                #need to drop more to converge
                if verbosity > 0:
                    print "mlim Monte Carlo inconvergent, downstep.\n"
                return limitingM(ru, rl/10.0, limsnr, PSF, sno, skyN, catM, catMerr, catSN, catI, verbosity, level+1)
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
    parser.add_argument("-r", "--radius", type=float, default=1000.0, help="pixel radius in which to take reference stars")
    parser.add_argument("-a", "--aperture", type=float, default=None, help="Aperture on which to perform aperture photometry. If not given, default is PSF photometry. If not a positive number, will use PSF to define Kron radius as aperture.")
    parser.add_argument("-o", "--source", type=str, default='object', help="target source name")
    parser.add_argument("-y", "--year", type=int, default=2016, help="year in which source was observed")
    parser.add_argument("-b", "--band", type=str, default='V', help="image filter band")
    parser.add_argument("-p", "--position", type=str, help="RA:DEC as deg:deg")
    parser.add_argument("-fwhm", type=float, default=5.0, help="image fwhm upper bound")
    parser.add_argument("-n", "--noiseSNR", type=float, default=0.0, help="signal to noise at detection limit")
    parser.add_argument("-s", "--satMag", type=float, default=14.0, help="CCD saturation, reference star magnitude upper bound")
    parser.add_argument("-f", "--refMag", type=float, default=19.0, help="Reliable lower bound for reference star brightness")
    parser.add_argument("-d", "--diffIm", type=str, default=None, help="Difference fits image containing source, which if given will be used instead to perform source photometry. Original image will be used for reference star photometry. Difference image wcs must match original image.")
    parser.add_argument("--fit_sky", action='store_const', const=True, default=False, help="Give this flag if it is desirable to fit for and subtract planar sky around the source.")
    parser.add_argument("-v", "--verbosity", action="count", default=0)
    args = parser.parse_args()
    
    #extract RA, DEC from position argument
    RA, DEC = [float(coord) for coord in args.position.split(':')]
    
    #load fits file, get relevant data
    catimage, time, wcs = loadFits(args.filename, getwcs=True, verbosity=args.verbosity)
    if args.diffIm is not None:
        #load difference image for source photometry
        image, t = loadFits(args.diffIm, verbosity=args.verbosity)
    else:
        #use original image for source photometry
        image = catimage
    
    #compute position, magnitude and error
    if args.noiseSNR != 0:
        RA, DEC, I, SN, M, Merr, Mlim = magnitude(image, catimage, wcs, args.catalog, args.catname, (RA,DEC), args.radius, args.aperture, args.source, args.band, args.fwhm, args.noiseSNR, args.satMag, args.refMag, args.fit_sky, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, I, SN, M, Merr, Mlim
    else:
        RA, DEC, I, SN, M, Merr = magnitude(image, catimage, wcs, args.catalog, args.catname, (RA,DEC), args.radius, args.aperture, args.source, args.band, args.fwhm, args.noiseSNR, args.satMag, args.refMag, args.fit_sky, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, I, SN, M, Merr
