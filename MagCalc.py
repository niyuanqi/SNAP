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
#python -m SNAP.MagCalc -c phot -o N300-1.Q0.SN -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 -f 16.0 --fit_sky conv.fits N300_1_Q0_SN.csv -d diff.fits
#python MagCalc.py -c phot -o N300-1.Q0.SN -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 -f 16.0 --fit_sky N300-1.Q0.B.151010_1604.A.033278.005604N3646.0060.nh.crop.fits N300_1_Q0_SN.csv
#python MagCalc.py -c diff -o KSP-N300-Nova -b 'B' -p 13.789218:-37.704572 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv --fit_sky N300-1.Q2.B.151009_0015.S.000859.005606N3754.0060.nh.fits N300-1.Q2.diff.cat
#python MagCalc.py -c dprs -o KSP-OT-1 -b 'B' -p 140.92247:-21.969278 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv --fit_sky N2784-7.Q1.B.150402_2125.S.015081.092205N2208.0060.nh.fits N2784-7.Q1.DPRS.cat

#essential modules
import numpy as np

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

def loadFits(filename, year=2016, getwcs=False, gethdr=False, verbosity=0):
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
    #    gethdr; boolean whether to return header                   #
    #      year; int year to measure time to                        #
    # ------------------------------------------------------------- #
    # Output                                                        #
    # ------------------------------------------------------------- #
    # image: numpy array containing image data                      #
    #  time: float time in days since start of year YYYY            #
    #   wcs: astropy wcs object, world coordinate system on image   #
    #   hdr: astropy header object                                  #
    #################################################################
    """

    from astropy.io import fits
    from astropy.time import Time
    from astropy.wcs import WCS

    from Astrometry import isot_day
    
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

    retlist = [image, time]

    if getwcs:
        try: #try to load WCS
            if verbosity > 0:
                print "loading world coordinate system"
            wcs = WCS(filename)
            retlist += [wcs]
        except:
            raise FitsError('Unable to load wcs data.')
    if gethdr:
        retlist += [header]
    return retlist

def magnitude(image, catimage, wcs, cat, catname, (RAo,DECo), radius=500, over_intens=None, aperture=None, psf='1', name='object', band='V', fwhm=5.0, limsnr=3.0, satmag=14.0, refmag=19.0, fitsky=True, satpix=40000.0, verbosity=0, diagnosis=False):
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
    #   over_in; float intensity calculation override. If given,        #
    #            this value will be used for the source intensity.      #
    #            is performed on the source.                            #
    #  aperture; float aperture around object in which to integrate     #
    #            light. If None; PSF is integrated. If not a positive   #
    #            number; then PSF is used to get Kron aperture.         #
    #       psf; string selects what parameters to fit the source with. #
    #            (e.g., '1' fits PSF central position only,             #
    #                   '2' fits PSF central position and height,       #
    #                   '3' fits PSF central position, height, shape.)  #
    #      name; str name of object.                                    #
    #      band; char observational filter of data.                     #
    #      fwhm; float estimate of FWHM on image.                       #
    #    limsnr; float signal to noise ratio defining detection limit,  #
    #            if 0.0, then no detection limits are calculated.       #
    #    satmag; float magnitude below which reference stars are        #
    #            considered to be saturated and hence not used.         #
    #    refmag; float magnitude above which reference stars are        #
    #            considered to be reliable, and therefore used.         #
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
    
    #essential functions
    import Catalog as ctlg
    import PSFlib as plib
    import Photometry as pht
    from Analysis.Cosmology import bands, flux_0

    #Single object? Generalize to multiple object. (Compatibility)
    if not hasattr(name, '__iter__'):
        #Make all listable objects into list
        name = [name]
        RAo = [RAo]
        DECo = [DECo]
        psf = [psf]
        fitsky = [fitsky]
        if over_intens is not None:
            over_intens = [over_intens]
    #number of sources to perform photometry on
    Nobj = len(name)
    
    #convert position of source to pixel
    Xo, Yo = np.zeros(Nobj), np.zeros(Nobj)
    for i in range(Nobj):
        Xo[i], Yo[i] = wcs.all_world2pix(RAo[i], DECo[i], 0)
        if verbosity > 0:
            print "Source "+str(i+1)+" located at: "+str(Xo[i])+", "+str(Yo[i])
        
    #load photometric reference stars catalog
    if verbosity > 0:
        print "loading catalog"
    if cat == 'phot':
        ID, RA, DEC, catM, catMerr = ctlg.catPhot(catname,band=band)
    elif cat == 'dprs':
        ID, RA, DEC, catM, catMerr = ctlg.catDPRS(catname,band=band)
    elif cat == 'diff':
        ID, RA, DEC, catM, catMerr = ctlg.catDiff(catname,band=band)
    elif cat == 'aavso':
        fovam = 2.0*radius*0.4/60.0 #arcmin radius in KMT scaling
        if band == 'I':
            if verbosity > 0:
                print "Performing AAVSO i -> I band conversion (Jodri 2006)"
            IDi, RAi, DECi, catMi, catMierr = ctlg.catAAVSO(RAo[0],DECo[0],fovam,'i',out=catname)
            IDr, RAr, DECr, catMr, catMrerr = ctlg.catAAVSO(RAo[0],DECo[0],fovam,'r',out=catname)
            ID, RA, DEC, catM, catMerr = [], [], [], [], []
            for i in range(len(IDi)):
                #for each ID in i band
                if IDi[i] in IDr:
                    #if also in r band
                    j = list(IDr).index(IDi[i]) #here it is
                    #get I band from i and r
                    ID.append(IDi[i])
                    RA.append(RAi[i])
                    DEC.append(DECi[i])
                    #Jodri 2006 general stars transform
                    catMI = 1.083*catMi[i] - 0.083*catMr[j] - 0.376
                    catMIerr = np.sqrt(((catMr[j]-catMi[i])*0.006)**2 + (0.004)**2 + (1.083*catMierr[i])**2 + (0.083*catMrerr[j])**2)
                    catM.append(catMI)
                    catMerr.append(catMIerr)
            ID, RA, DEC, catM, catMerr = np.array(ID), np.array(RA), np.array(DEC), np.array(catM), np.array(catMerr)
            ID, RA, DEC, catM, catMerr = ctlg.catAAVSO(RAo[0],DECo[0],fovam,'i',out=catname)
        else:
            ID, RA, DEC, catM, catMerr = ctlg.catAAVSO(RAo[0],DECo[0],fovam,band,out=catname)
    
    #convert position of catalog stars to pixels
    catX, catY = wcs.all_world2pix(RA, DEC, 0)
    catX, catY = catX.astype(float), catY.astype(float)
    #select catalog stars within some radius of object
    index = pht.dist(catX,catY,np.mean(Xo),np.mean(Yo)) < radius
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
            print catRA[int(i)], catDEC[int(i)], catM[int(i)], catMerr[int(i)]
    #number of selected catalog stars
    Ncat = len(ID)

    if satpix == 0:
        #measure saturation level: works if there is saturated star 
        sat = satpix(image)
    if verbosity > 3:
        #essential extra import
        import matplotlib.pyplot as plt
        #plot image of catalog positions
        plt.imshow(catimage, cmap='Greys', vmax=0.001*np.amax(catimage), vmin=0)
        plt.scatter(catX, catY)
        plt.scatter(Xo, Yo, c='r')
        for i in range(Ncat):
            plt.text(catX[i], catY[i], ID[i])
        for i in range(Nobj):
            plt.text(Xo[i], Yo[i], name[i])
        plt.show()
    #photometry on catalog stars
    catMags = [] #magnitude list
    catMagerrs = []
    catPSFs = [] #PSF fits to catalog stars
    catPSFerrs = [] #fit errors
    catSkyMs = [] #background average count
    catSkyNs = [] #background noise
    catXs = []
    catYs = []
    catX2dofs = []
    catIDs = []
    if verbosity > 0:
        print "Extracting PSF of "+str(Ncat)+" catalog stars."
    
    #calculate PSF for each reference star
    for i in range(Ncat):
        if verbosity > 0:
            print "\nComputing PSF of "+str(i+1)+"/"+str(Ncat)
        #position of reference star
        x0, y0 = catX[i], catY[i]
        #calculate intensity and SN ratio
        #verbosity is reduced for catalog stars
        try:
            PSFpopt, PSFperr, X2dof, skypopt, skyN = pht.PSFextract(catimage, x0, y0, fwhm=fwhm, fitsky=fitsky[0], sat=satpix, verbosity=verbosity-1)
            PSF, PSFerr = PSFpopt[1:5], PSFperr[1:5]
        except:
            PSFpopt = [0]*7
        
        #Take only reference stars whose fits are sane 
        if plib.E2moff_verify(PSFpopt, x0, y0):
            PSF[3] = PSF[3] % 180.0 #principle angle
            #break x,y degeneracy in theta
            if PSF[0] > PSF[1]:
                #switched up x,y axes
                PSF[0], PSF[1] = PSF[1], PSF[0]
                if PSF[3]<90:
                    PSF[3] += 90
                if PSF[3]>90:
                    PSF[3] -= 90
                if verbosity > 2:
                    print "Swapping theta to enforce ax<ay"
            #save magnitude of catalog star
            catMags.append(catM[i])
            catMagerrs.append(catMerr[i])
            #save catalog star fit
            catPSFs.append(PSF)
            catPSFerrs.append(PSFerr)
            #save sky details
            catSkyMs.append(skypopt[0]*PSFpopt[5]+skypopt[1]*PSFpopt[6]+skypopt[2])
            catSkyNs.append(skyN)
            #save ref star position
            catXs.append(PSFpopt[5])
            catYs.append(PSFpopt[6])
            #save fit X2/dof
            catX2dofs.append(X2dof)
            catIDs.append(ID[i])
        else:
            if verbosity > 0:
                #say something about fit being bad for this particular star
                print "\nReference star ID"+str(ID[i])+" fit unacceptable"
                print "Criminal located at position "+str([x0,y0])+".\n"

    ncat = len(catMags) #number of good stars
    if verbosity > 0:
        print "\nNumber of reference stars used: "+str(ncat)+"/"+str(Ncat)
    if float(ncat)/Ncat < 0.5:
        #over half reference stars are invalid... how??
        raise PSFError('Unable to perform photometry on reference stars.')
    catMags = np.array(catMags)
    catMagerrs = np.array(catMagerrs)
    catPSFs = np.array(catPSFs)
    catPSFerrs = np.array(catPSFerrs)
    catSkyMs = np.array(catSkyMs)
    catSkyNs = np.array(catSkyNs)
    catXs = np.array(catXs)
    catYs = np.array(catYs)
    catX2dofs = np.array(catX2dofs)
    catIDs = np.array(catIDs)
    #print catPSFs.T[3]
    
    #calculate average psf among reference stars
    w = 1/np.square(catPSFerrs)
    catPSF = (catPSFs*w).sum(0)/w.sum(0)
    catPSFerr = np.sqrt(1/w.sum(0))
    #calculate average sky parameters
    skyval = catSkyMs.mean() #mean constant background
    noise = catSkyNs.mean() #mean std from background
    if verbosity > 0:
        print "Average PSF [ax, ay, b, theta] =",str(catPSF)
        print "parameter errors =",str(catPSFerr)
        print "Average FWHMx,FWHMy =",str(plib.E2moff_toFWHM(*catPSF[:-1]))
        print "Average background sky count =",str(skyval)
        print "Average noise in background =",str(noise)
        print ""

    if diagnosis:
        if verbosity > 0:
            print "Returning image data for diagnosis"
        return catPSF, catPSFerr, skyval, noise

    #Integration using common PSF
    catIs = np.zeros(ncat)
    catSNs = np.zeros(ncat)
    skyNs = np.zeros(ncat)
    for i in range(ncat):
        if verbosity > 0:
            print "Computing intensity of "+str(i+1)+"/"+str(ncat)
        #position of star in catalog
        x0, y0 = catXs[i], catYs[i]
        #calculate intensity and SN ratio with reduced verbosity
        PSFpopt, PSFperr, X2dof, skypopt, skyN = pht.PSFscale(catimage, catPSF, catPSFerr, x0, y0, fitsky=fitsky[0], sat=satpix, verbosity=verbosity-1)
        #check preferred intensity calculation method
        if aperture is None:
            #integrate PSF directly
            I, SN = pht.PSF_photometry(catimage, x0, y0, PSFpopt, PSFperr, '1', skypopt, skyN, verbosity=verbosity-1)
        elif aperture <= 0:
            #perform aperture photometry
            #use FWHM of catPSF to define Kron aperture
            I, SN = pht.Ap_photometry(catimage, x0, y0, skypopt, skyN, PSF=catPSF, fitsky=True, verbosity=verbosity-1)
        else:
            #use aperture given directly
            I, SN = pht.Ap_photometry(catimage, x0, y0, skypopt, skyN, radius=aperture, fitsky=True, verbosity=verbosity-1)
        #check if reference stars are valid
        if I == 0 or SN == 0 or skyN == 0:
            raise PSFError('Unable to perform photometry on reference stars.')
        #save intensity and SN ratio
        catIs[i]=I
        catSNs[i]=SN
        skyNs[i]=skyN

    if verbosity > 0:
        print "Mean SN of reference stars:",np.mean(catSNs)
        print "Mean background noise:", np.mean(skyNs), np.std(skyNs)
        print ""
    if verbosity > 1:
        #essential extra import
        import MagPlot as magplt
        
        #diagnostic for SNR, N, I calculation routine
        #instrumental magnitudes
        insMags = -2.5*np.log10(catIs)
        insMagerrs = (2.5/np.log(10))*(1/catSNs)
        #checks for right correlation between intensity and SNR
        magplt.sn_corr_plot(insMags, catSNs)
        #checks for right correlation between intensity and noise
        catNs = catIs/catSNs
        magplt.noise_corr_plot(insMags, catNs)
        #photometric solution between instrumental magnitudes vs catalog
        magplt.phot_sol(insMags, insMagerrs, catMags, catMagerrs)
        if band == 'B':
            magplt.Bcol_corr(cat, catname, catIDs, RAo, DECo, radius, insMags, insMagerrs, catMags, catMagerrs)
        elif band == 'I':
            magplt.Icol_corr(cat, catname, catIDs, RAo, DECo, radius, insMags, insMagerrs, catMags, catMagerrs)
        #check reference star fit qualities
        magplt.X2_hist(catX2dofs)
        print ""
        
    #calculate photometry for source object
    #extract PSF to as great a degree as needed from source
    if Nobj == 1 or aperture is not None:
        PSFpopt, PSFperr, X2dof, skypopto, skyNo = [0]*7, [0]*7, 0, [0]*3, 0
        if verbosity > 0:
            print "Computing photometry of source "+name[0]

        if aperture is not None:
            PSFpopt, PSFperr, X2dof, skypopto, skyNo = pht.PSFextract(image, Xo[0], Yo[0], fwhm, fitsky=fitsky[0], sat=satpix, verbosity=verbosity)
        elif psf[0] == '1':
            PSFpopt, PSFperr, X2dof, skypopto, skyNo = pht.PSFscale(image, catPSF, catPSFerr, Xo[0], Yo[0], fitsky=fitsky[0], sat=satpix, verbosity=verbosity)
        elif psf[0] == '2':
            PSFpopt, PSFperr, X2dof, skypopto, skyNo = pht.PSFfit(image, catPSF, catPSFerr, Xo[0], Yo[0], fitsky=fitsky[0], sat=satpix, verbosity=verbosity)
        elif psf[0] == '3':
            PSFpopt, PSFperr, X2dof, skypopto, skyNo = pht.PSFextract(image, Xo[0], Yo[0], fwhm, fitsky=fitsky[0], sat=satpix, verbosity=verbosity)
        else:
            PSFpopt, PSFperr, X2dof, skypopto, skyNo = pht.PSFscale(image, catPSF, catPSFerr, Xo[0], Yo[0], fitsky=fitsky[0], sat=satpix, verbosity=verbosity)
        
        PSFpopt, PSFperr = [PSFpopt], [PSFperr]
        #check preferred intensity calculation method
        if aperture is None:
            #integrate PSF directly
            Io, SNo = pht.PSF_photometry(image, Xo[0], Yo[0], PSFpopt[0], PSFperr[0], psf[0], skypopto, skyNo, verbosity=verbosity)
        else:
            #perform aperture photometry
            if aperture <= 0:
                #use FWHM of catPSF to define Kron aperture
                Io, SNo = pht.Ap_photometry(image, Xo[0], Yo[0], skypopto, skyNo, PSF=catPSF, fitsky=fitsky[0], verbosity=verbosity)
            else:
                #use aperture given directly
                Io, SNo = pht.Ap_photometry(image, Xo[0], Yo[0], skypopto, skyNo, radius=aperture, fitsky=fitsky[0], verbosity=verbosity)
        Io, SNo = [Io], [SNo]
    #deal with mult-object photometry
    elif Nobj > 0:
        PSFpopt, PSFperr, X2dof, skypopto, skyNo = pht.PSFmulti(image, catPSF, catPSFerr, psf, Xo, Yo, fitsky=fitsky, sat=satpix, verbosity=verbosity)

        #check preferred intensity calculation method
        if aperture is None:
            #integrate PSF directly
            Io, SNo = np.zeros(Nobj), np.zeros(Nobj)
            for i in range(Nobj):
                if verbosity > 0:
                    print "Computing photometry of source "+name[i]
                Io[i], SNo[i] = pht.PSF_photometry(image, Xo[i], Yo[i], PSFpopt[i], PSFperr[i], psf[i], skypopto, skyNo, verbosity=verbosity)
        else:
            #perform aperture photometry
            print "Multi-object aperture photometry under construction!"
            print "Try again later, when it is finished."
            Io, SNo = [0]*Nobj, [0]*Nobj
    else:
        #Invalid fit selected, don't fit source
        if verbosity > 0:
            print "No source selected."
            Io, SNo, skyNo = [0]*Nobj, [0]*Nobj, 0
    #Override intensity value using other calculation
    if over_intens is not None:
        Io = over_intens
    
    #Process each source
    I,RAo,DECo,mo,mo_err = np.zeros(Nobj),np.zeros(Nobj),np.zeros(Nobj),np.zeros(Nobj),np.zeros(Nobj)
    for i in range(Nobj):
        #check if source is valid
        if Io[i] != 0 and SNo[i] != 0 and skyNo != 0:
            #calculate relative flux of object wrt each reference star
            Ir = flux_0[bands[band]]*1e6*np.power(10,-catMags/2.5)*Io[i]/catIs
            Io_err = Io[i]/SNo[i]
            catI_err = catIs/catSNs
            Ir_err = Ir*np.sqrt(np.square(1/catSNs)+np.square(np.log(10)*catMagerrs/2.5))

            #calculate weighted mean
            w = 1/np.square(Ir_err)
            I[i] = np.sum(Ir*w)/np.sum(w)
            I_rand = np.sqrt(1/np.sum(w))
            I_err = np.sqrt((I[i]*Io_err/Io[i])**2 + I_rand**2)
            JN = I[i]/Io[i] #jansky - data number conversion, uJy/#

            #if verbosity > 0:
                #print "Contribution of intrinsic error:", Io_err/Io[i]
                #print "Contribution of ref star scatter:", I_rand/I[i]
                
            #SN = I/I_err
        else:
            #no valid source
            I[i], SNo[i] = float('NaN'), float('NaN')
        #try to compute magnitude if source is present
        if I[i] != float('NaN') and I[i] > 0 and skyNo != 0:
            #convert position to world coordinates
            if psf[i][0] != 's':
                Xp, Yp = PSFpopt[i][5], PSFpopt[i][6]
            else:
                if psf[i][1] == 'n':
                    Xp, Yp = PSFpopt[i][3], PSFpopt[i][4]
                else:
                    Xp, Yp = PSFpopt[i][2], PSFpopt[i][3]
                    
            RAo[i], DECo[i] = wcs.all_pix2world(Xp, Yp, 0)
            #calculate magnitude from flux
            mo[i] = -2.5*np.log10(I[i]/(flux_0[bands[band]]*1e6))
            mo_err[i] = (2.5/np.log(10))*(I_err/I[i])
        else:
            #bad source
            mo[i], mo_err[i] = float('NaN'), float('NaN')
            RAo[i], DECo[i] = wcs.all_pix2world(Xo[i], Yo[i], 0)

    if limsnr != 0 and skyNo != 0:
        #sky noise properly estimated, calculate limiting magnitude
        ru = 10.0 #recursion seed
        rl = 0.1 #recursion seed
        mlim, SNlim, expu, expd = limitingM(ru, rl, limsnr, catPSF, catPSFerr, skyNo, catMags, catMagerrs, catSNs, catIs, verbosity)
        #return calculated magnitude, magnitude errors, and limiting magnitude
        return RAo, DECo, I, SNo, mo, mo_err, mlim
    elif limsnr != 0:
        #no sky noise estimate
        return RAo, DECo, I, SNo, mo, mo_err, float('NaN')
    else:
        #return calculated magnitude and magnitude errors
        return RAo, DECo, I, SNo, mo, mo_err

#function: recursively calculates limiting magnitude by scaling PSF to SN3.0
def limitingM(ru, rl, limsnr, PSF, PSFerr, skyN, catM, catMerr, catSN, catI, verbosity=0, level=0):
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
    
    import PSFlib as plib
    
    if len(PSF) == 4:
        #extract values from PSF
        ax, ay = abs(PSF[0]), abs(PSF[1])
        b, theta = PSF[2], PSF[3]
        axerr, ayerr = PSFerr[0], PSFerr[1]
        thetaerr, berr = PSFerr[2], PSFerr[3]
        FWHMx, FWHMy = plib.E2moff_toFWHM(ax,ay,b)
        #estimate height of average catalog PSF
        A = np.mean(catSN)*(b-1)/np.pi/ax/ay
        #estimate upper bound for factor to apply to A to get to source height
        f = np.mean(catSN)/(limsnr)
        #calculate snr for many synthetic sources near limiting snr
        n = 10 #resolution of parameter space
        Aest = A/f
        A_trials = np.linspace(rl*Aest,ru*Aest,n)
        I_trials = np.zeros(len(A_trials))
        SN_trials = np.zeros(len(A_trials))
        #compute optimal aperture radius (90% source light)
        frac = 0.9
        #do photometry over synthetic sources for each A
        if verbosity > 0:
            print "Computing "+str(n)+" synthetic sources to find mlim"
            print "PSF parameters: " + str(PSF)
            print "A bound: " + str(Aest)
        for j in range(len(A_trials)):
            if verbosity > 2:
                print "Computing "+str(j+1)+"/"+str(n)
            #integrate synthetic PSF
            I = plib.E2moff_integrate(A_trials[j],ax,ay,b,frac)
            #get aperture size around synthetic star
            ap_size = plib.E2moff_apsize(ax,ay,b,frac)
            #at SN <= 5, noise dominated, I/(skyN**2*aperture.size) < 0.1
            sigma = np.sqrt(np.absolute(I) + ap_size*skyN**2)
            SN = I/sigma
            #append to trials
            I_trials[j] = I
            SN_trials[j] = SN
        #calculate I for source closest to limiting snr
        idx = np.argmin(np.square(SN_trials - limsnr))
        Ilim, SNlim = I_trials[idx], SN_trials[idx]
        SNUbound = max(SN_trials)
        SNLbound = min(SN_trials)
        #calculate magnitude wrt each reference star
        mlim = catM - 2.5*np.log10(Ilim/np.array(catI))
        mlim_err = np.sqrt(np.square((2.5/np.log(10))*(1/catSN))+np.square(catMerr))
        w = 1/np.square(mlim_err)
        mlim = np.sum(mlim*w)/np.sum(w)
        #output comments
        if verbosity > 0:
            print "maximum SN="+str(SNUbound)
            print "minimum SN="+str(SNLbound)
            print "mlim calculated at SN="+str(SNlim)
        #check convergence
        if abs(SNlim - limsnr) > 0.1:
            #prevent stack overflow after 20 recursions
            if level+1 > 20:
                if verbosity > 0:
                    print "Convergence terminated to prevent stack overflow."
                return float('NaN'), float('NaN'), float('NaN'), float('NaN')
            #not yet convergent, narrow down SN range
            if SNUbound > limsnr and SNLbound < limsnr:
                #take closest points and converge more
                if verbosity > 0:
                    print "mlim Monte Carlo inconvergent, refine.\n"
                ru = A_trials[SN_trials>limsnr][0]/Aest
                rl = A_trials[SN_trials<limsnr][-1]/Aest
                return limitingM(ru, rl, limsnr, PSF, PSFerr, skyN, catM, catMerr, catSN, catI, verbosity, level+1)
            elif SNUbound < limsnr:
                #need to rise more to converge
                if verbosity > 0:
                    print "mlim Monte Carlo inconvergent, upstep.\n"
                return limitingM(ru*10.0, rl, limsnr, PSF, PSFerr, skyN, catM, catMerr, catSN, catI, verbosity, level+1)
            elif SNLbound > limsnr:
                #need to drop more to converge
                if verbosity > 0:
                    print "mlim Monte Carlo inconvergent, downstep.\n"
                return limitingM(ru, rl/10.0, limsnr, PSF, PSFerr, skyN, catM, catMerr, catSN, catI, verbosity, level+1)
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
    parser.add_argument("-psf", type=str, default='1', help="Degrees of freedom to fit PSF of source. All reference stars are fit will full moffat fit to get PSF on image. 1=only relax height of PSF for source. 2=relax height and centroid of PSF. 3=relax height, centroid, and shape. sn=Sersic galaxy fit with general n. s1.324=Sersic galaxy fit with n=1.324 or whatev.")
    parser.add_argument("-o", "--source", type=str, default='object', help="target source name")
    parser.add_argument("-y", "--year", type=int, default=2016, help="year in which source was observed")
    parser.add_argument("-b", "--band", type=str, default='V', help="image filter band")
    parser.add_argument("-p", "--position", type=str, help="RA:DEC as deg:deg")
    parser.add_argument("-fwhm", type=float, default=5.0, help="image fwhm upper bound")
    parser.add_argument("-n", "--noiseSNR", type=float, default=0.0, help="signal to noise at detection limit")
    parser.add_argument("-s", "--satMag", type=float, default=14.0, help="CCD saturation, reference star magnitude upper bound")
    parser.add_argument("-sp", "--satpix", type=float, default=40000.0, help="CCD upper valid pixel count. If given value is 0, code can determine satpix (only for images containing at least one star which has saturated CCD full well capacity).")
    parser.add_argument("-f", "--refMag", type=float, default=19.0, help="Reliable lower bound for reference star brightness")
    parser.add_argument("-d", "--diffIm", type=str, default=None, help="Difference fits image containing source, which if given will be used instead to perform source photometry. Original image will be used for reference star photometry. Difference image wcs and psf must match original image. If not, matching is required in preprocessing.")
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
        RA, DEC, I, SN, M, Merr, Mlim = magnitude(image, catimage, wcs, args.catalog, args.catname, (RA,DEC), args.radius, args.aperture, args.psf, args.source, args.band, args.fwhm, args.noiseSNR, args.satMag, args.refMag, args.fit_sky, args.satpix, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, I, SN, M, Merr, Mlim
    else:
        RA, DEC, I, SN, M, Merr = magnitude(image, catimage, wcs, args.catalog, args.catname, (RA,DEC), args.radius, args.aperture, args.psf, args.source, args.band, args.fwhm, args.noiseSNR, args.satMag, args.refMag, args.fit_sky, args.satpix, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, I, SN, M, Merr
