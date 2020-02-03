#################################################################
# Name:     ClickMag.py                                         #
# Author:   Yuan Qi Ni                                          #
# Version:  December 5, 2019                                    #
# Function: Program contains functions for doing photometry by  #
#           clicking on an image.                               #
#################################################################

#Example:
#python -m SNAP.ClickMag -c phot -b 'I' -r 3000 -fwhm 5 -v -n 3.0 -s 15.0 -f 16.0 Ideepestcrop.fits N300_1_Q0_SN.csv
#Enter aperture '0' for Kron aperture photometry

#function: click on an image to receive magnitudes of sources
def ClickMag(image, wcs, cat, catname, radius=500, band='V', fwhm=5.0, limsnr=3.0, satmag=14.0, refmag=19.0, satpix=40000.0, verbosity=0, diagnosis=False):
    """
    #####################################################################
    # Desc: Compute magnitude of objects in image using ref catalog.    #
    #       Select objects by clicking on them.                         #
    # ----------------------------------------------------------------- #
    # Imports:                                                          #
    # ----------------------------------------------------------------- #
    # Input                                                             #
    # ----------------------------------------------------------------- #
    #     image: numpy array containing image data on which to measure  #
    #            source and reference star photometry.                  #
    #       wcs: astropy wcs object, world coordinate system on image.  #
    #       cat: str catalog type (phot, dprs, or diff from Catalog.py) #
    #   catname: str catalog name.                                      #
    #    radius; float radius around object in which to take ref stars. #
    #      band; char observational filter of data.                     #
    #      fwhm; float estimate of FWHM on image.                       #
    #    limsnr; float signal to noise ratio defining detection limit,  #
    #            if 0.0, then no detection limits are calculated.       #
    #    satmag; float magnitude below which reference stars are        #
    #            considered to be saturated and hence not used.         #
    #    refmag; float magnitude above which reference stars are        #
    #            considered to be reliable, and therefore used.         #
    # verbosity; int counts verbosity level.                            #
    # ----------------------------------------------------------------- #
    # Output                                                            #
    # ----------------------------------------------------------------- #
    #  RAo, DECo: float measured equatorial coordinates of sources.     #
    #    Io, SNo: float measured intensities and SNRs of sources.       #
    # mo, mo_err: float calibrated magnitudes and errors of sources.    #
    #       mlim; float detection limits at source positions.           #
    #####################################################################
    """
    #essential imports
    import matplotlib.pyplot as plt
    import numpy as np

    #essential functions
    import Catalog as ctlg
    import PSFlib as plib
    import Photometry as pht
    from MagCalc import PSFError, magnitude
    from Analysis.Cosmology import bands, flux_0

    #missing
    #(RAo,DECo) CHECK
    #aperture=None
    #psf='1'
    #name='object' CHECK
    #fitsky=True
    
    #steps
    
    #load all reference stars on image DONE
    #plot image DONE
    #plot reference stars in blue DONE
    #centroid clicked sources -> replot DONE
    #after each click, ask for photometry method
    #evaluate ra and dec coordinates DONE
    #evaluate magnitudes of each using MagCalc DONE

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
        fovam = 4.0*radius*0.4/60.0 #arcmin radius in KMT scaling
        RAo, DECo = [wcs.wcs.crval[0]], [wcs.wcs.crval[1]]
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
    print max(catM)
    #select catalog stars within edges
    index = np.logical_and(catX > 80, image.shape[1]-catX > 80)
    index = np.logical_and(index, np.logical_and(catY > 80, image.shape[0]-catY > 80))
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
        print "Showing catalog star IDs:"
        for i in range(len(ID)):
            print ID[int(i)], catX[int(i)], catY[int(i)]
            print catRA[int(i)], catDEC[int(i)], catM[int(i)], catMerr[int(i)]
    #number of selected catalog stars
    Ncat = len(ID)
    print Ncat
    
    #plot image of field
    fig = plt.figure()
    plt.imshow(image, cmap='Greys', vmax=0.0001*np.amax(image), vmin=0)
    plt.scatter(catX, catY, c='b', marker='+')
    plt.tight_layout()
    
    x_cens, y_cens = [], []
    psfs = []
    apers = []
    fit_skys = []
    #function: record click position on matplotlib image
    def onclick(event):
        #output click position
        if verbosity > 0:
            print 'Clicked pixel: x=%d, y=%d'%(event.xdata, event.ydata)
        #record clicked position
        x_click, y_click = event.xdata, event.ydata
        #obtain an aperture around clicked star
        intens, x, y = pht.ap_get(image, x_click, y_click, 0, 1*fwhm)
        #computer centroid of clicked star
        x_cen = np.sum(intens*x)/intens.sum()
        y_cen = np.sum(intens*y)/intens.sum()
        print 'Centroid pixel: x=%d, y=%d'%(x_cen, y_cen)
        #Should we keep the star
        keep = raw_input("Keep? (y/n)")
        if keep == 'y' or keep == "":
            #record centroid position
            x_cens.append(x_cen)
            y_cens.append(y_cen)
            #plot centroid position
            plt.scatter(x_cen, y_cen, marker='+', c='r', s=80)
            fig.canvas.draw()
            
            if len(x_cens) == 1:
                #First source -> prompt parameters
                use_prev = 'n'
            else:
                #Prompt for previous parameter selection
                use_prev = raw_input("Use previous photometry parameters? (y/n)")
            
            if use_prev == 'y' or use_prev == "":
                #Use previous parameters
                apers.append(apers[-1])
                psfs.append(psfs[-1])
                fit_skys.append(fit_skys[-1])
            else:
                #Prompt for aperture size
                aper = raw_input("Aperture=? (float / empty for PSF)")
                if aper == "":
                    #PSF photometry parameters
                    psf = raw_input("PSF type=? ('1', '2', '3', 's<sersic index>')")
                    apers.append(None)
                    psfs.append(psf)
                else:
                    #Aperture photometry parameters
                    aper = float(aper)
                    apers.append(aper)
                    psfs.append('1')
                #Prompt for planar sky fitting
                fit_sky = raw_input("Fit planar background? (y/n)")
                if fit_sky == 'y' or fit_sky == "":
                    fit_skys.append(True)
                else:
                    fit_skys.append(False)
        else:
            #Drop clicked position and keep listening
            print "Dropping"
        print ""
        print "Click another source"
    #plot image and listen for a click
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    if verbosity > 0:
        print ""
        print "Click on the sources"
        print '---------------------------------------'
    plt.show()
    #convert source centroid positions to ra, dec degrees
    x_cens = np.array(x_cens)
    y_cens = np.array(y_cens)
    ra_cens, dec_cens = wcs.all_pix2world(x_cens, y_cens, 0)

    #compute magnitude of each source
    RAs, DECs = np.zeros(len(x_cens)),np.zeros(len(x_cens))
    Is, SNs = np.zeros(len(x_cens)),np.zeros(len(x_cens))
    Ms, Merrs = np.zeros(len(x_cens)),np.zeros(len(x_cens))
    if limsnr != 0:
        #also compute limiting magnitudes
        Mlims = np.zeros(len(x_cens))
        for i in range(len(ra_cens)):
            RA, DEC, I, SN, M, Merr, Mlim = magnitude(image, image, wcs, cat, catname, (ra_cens[i],dec_cens[i]), radius, apers[i], psfs[i], "Source"+str(i), band, fwhm, limsnr, satmag, refmag, fit_skys[i], satpix, verbosity)
            #output position, magnitude
            print "Output for Source "+str(i)
            print RA, DEC, I, SN, M, Merr, Mlim
            RAs[i], DECs[i], Is[i], SNs[i], Ms[i], Merrs[i], Mlims[i] = RA[0], DEC[0], I[0], SN[0], M[0], Merr[0], Mlim
        retlist = [RAs, DECs, Is, SNs, Ms, Merrs, Mlims]
    else:
        #don't compute limiting magnitudes
        for i in range(len(ra_cens)):
            RA, DEC, I, SN, M, Merr = magnitude(image, image, wcs, cat, catname, (ra_cens[i],dec_cens[i]), radius, apers[i], psfs[i], "Source"+str(i), band, fwhm, limsnr, satmag, refmag, fit_skys[i], satpix, verbosity)
            #output position, magnitude
            print "Output for Source "+str(i)
            print RA, DEC, I, SN, M, Merr
            RAs[i], DECs[i], Is[i], SNs[i], Ms[i], Merrs[i] = RA[0], DEC[0], I[0], SN[0], M[0], Merr[0]
        retlist = [RAs, DECs, Is, SNs, Ms, Merrs]

    #plot diagnostic
    if verbosity > 0:
        #plot image of field
        fig = plt.figure()
        plt.imshow(image, cmap='Greys', vmax=0.0001*np.amax(image), vmin=0)
        plt.scatter(catX, catY, c='b', marker='+', label="catalog")
        plt.scatter(x_cens, y_cens, marker='+', c='r', s=80, label="centroids")
        x_meas, y_meas = wcs.all_world2pix(RAs, DECs, 0)
        plt.scatter(x_meas, y_meas, marker='+', c='g', s=80, label="measured")
        print len(x_cens), len(x_meas)
        plt.tight_layout()
        plt.show()

    #return source measurements
    return retlist

#command line execution
if __name__ == "__main__":

    import argparse
    import numpy as np
    
    #command line arguments
    parser = argparse.ArgumentParser(description="Find Photometric Magnitude")
    parser.add_argument("filename", type=str, help="fits image containing source")
    parser.add_argument("-c", "--catalog", type=str, default='phot', help="reference stars catalog convention")
    parser.add_argument("catname", type=str, help="tab separated reference stars catalog file")
    parser.add_argument("-r", "--radius", type=float, default=1000.0, help="pixel radius in which to take reference stars")
    parser.add_argument("-b", "--band", type=str, default='V', help="image filter band")
    parser.add_argument("-fwhm", type=float, default=5.0, help="image fwhm upper bound")
    parser.add_argument("-n", "--noiseSNR", type=float, default=0.0, help="signal to noise at detection limit")
    parser.add_argument("-s", "--satMag", type=float, default=14.0, help="CCD saturation, reference star magnitude upper bound")
    parser.add_argument("-sp", "--satpix", type=float, default=40000.0, help="CCD upper valid pixel count. If given value is 0, code can determine satpix (only for images containing at least one star which has saturated CCD full well capacity).")
    parser.add_argument("-f", "--refMag", type=float, default=19.0, help="Reliable lower bound for reference star brightness")
    parser.add_argument("-v", "--verbosity", action="count", default=0)
    args = parser.parse_args()

    from MagCalc import loadFits
    
    #load fits file, get relevant data
    image, time, wcs = loadFits(args.filename, getwcs=True, verbosity=args.verbosity)
    
    #compute position, magnitude and error
    if args.noiseSNR != 0:
        RA, DEC, I, SN, M, Merr, Mlim = ClickMag(image, wcs, args.catalog, args.catname, args.radius, args.band, args.fwhm, args.noiseSNR, args.satMag, args.refMag, args.satpix, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, I, SN, M, Merr, Mlim
        np.savetxt("source_list.cat", np.array([RA, DEC, I, SN, M, Merr, Mlim]).T)
    else:
        RA, DEC, I, SN, M, Merr = ClickMag(image, wcs, args.catalog, args.catname, args.radius, args.band, args.fwhm, args.noiseSNR, args.satMag, args.refMag, args.satpix, args.verbosity)
        #output position, magnitude
        print time, RA, DEC, I, SN, M, Merr
        np.savetxt("source_list.cat", np.array([RA, DEC, I, SN, M, Merr]).T)
