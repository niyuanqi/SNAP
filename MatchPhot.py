#################################################################
# Name:     MatchPhot.py                                        #
# Author:   DECuan Qi Ni                                        #
# Version:  February 9, 2017                                    #
# Function: Program contains routine for matching photometry    #
#           between given catalog file and standard AAVSO-APASS #
#           reference star catalog.                             #
#################################################################

#essential modules
import numpy as np 

#essential imports
from Astrometry import*
from Vizier import*

#linear photometric solution
def linMatch(x, a, b, x_err=0, a_err=0, b_err=0):
    m = a*x + b
    m_err = np.sqrt(b_err**2 + np.square(a_err*x)+np.square(x_err*a))
    return m, m_err

#linear function
def lin(x, a, b):
    return a*x + b

#distance metric on image
def dist(x1, y1, x2, y2):
    return np.sqrt(np.square(x1-x2)+np.square(y1-y2))

#match photometry between catalog and reference
def matchSExPhot(catname, band, refname=None, fitsname=None, pos=None, verbosity=0):

    from astropy.io import fits
    from astropy.wcs import WCS
    from scipy.optimize import curve_fit
    
    #get centroid positions in image from SExtractor Catalog
    if verbosity > 1:
        print "Loading SExtractor Catalog"
    N, F, F_err, RA, DEC = np.loadtxt(catname, unpack=True)
    #fitting will be in log axes
    F_err = F_err/(F*np.log(10.0))
    F = np.log10(F)

    fovam = 60.0 # size of KMTNet square search field in arc min

    if refname is not None:
        if verbosity > 1:
            print "Reference Catalog Provided"
        f = open(refname, 'r')
        s = f.read()
        f.close()
        #parse AAVSO catalog
        name,rad,ded,rmag = aavso_static(s,band)
    elif fitsname is not None:
        if verbosity > 1:
            print "Fits image field coordinates provided, searching AAVSO"
        #import fits image corresponding to catalog 
        s1 = fits.open(fitsname)
        # Read central pixel from the FITS file
        cnx = s1[0].header['CRPIX1']
        cny = s1[0].header['CRPIX2']
        s1.close()
        #get wcs from FITS file
        w = WCS(fitsname)
        radeg, decdeg = w.all_pix2world(cnx, cny, 0)
        radeg, decdeg = float(radeg), float(decdeg)
        #Query AAVSO catalog
        name,rad,ded,rmag,s = aavso(radeg,decdeg,fovam,band,out=True)
    elif pos is not None:
        if verbosity > 1:
            print "Field coordinates provided, searching AAVSO"
        #Query AAVSO catalog
        name,rad,ded,rmag,s = aavso(pos[0],pos[1],fovam,band,out=True)
    else:
        print "Insufficient information to retrieve AAVSO catalog, try -h flag and give one of three optional arguments."
        return

    if verbosity > 1:
            print "Matching star fields"
    wh = np.where(rmag < 16.0)[0] # select only bright stars r < 16 mag.
    rad, ded, rmag = rad[wh], ded[wh], rmag[wh]
    wh = np.where(rmag > 14.0)[0] # select only unsaturated stars r > 14 mag.
    rad, ded, rmag = rad[wh], ded[wh], rmag[wh]
    #match approximately the magnitudes from SExtractor
    ram, dem, RAm, DECm, magm, Fm, Fm_err = [],[],[],[],[],[],[]
    #for each SExtractor object
    for i in range(len(N)):
        #check if object in aavso within 1arcsec
        for j in range(len(rmag)):
            if not np.isnan(rmag)[j]:
                #use small angle approximation for close stars
                if smallAngle([rad[j], ded[j]],[RA[i], DEC[i]])*13600 < 2:
                    #add object to matched list
                    ram.append(rad[j])
                    dem.append(ded[j])
                    RAm.append(RA[i])
                    DECm.append(DEC[i])
                    magm.append(rmag[j])
                    Fm.append(F[i])
                    Fm_err.append(F_err[i])
    Fm = np.array(Fm)
    Fm_err = np.array(Fm_err)

    if verbosity > 1:
        print "Fitting AAVSO magnitude to SExtractor flux"
    #fit photometric solution
    popt, pcov = curve_fit(lin, magm, Fm, sigma = Fm_err)
    perr = np.sqrt(np.diag(pcov))

    #invert linear function
    popt = [1/popt[0], -popt[1]/popt[0]]
    perr = [perr[0]/popt[0]**2, (-popt[1]/popt[0])*np.sqrt((perr[1]/popt[1])**2+(perr[0]/popt[0])**2)]

    Ffit = Fm
    Mfit, Mfit_err = linMatch(Ffit, popt[0], popt[1], Fm_err, perr[0], perr[1])
    X2dof = np.sum(np.square((Mfit-magm)/Mfit_err))/(len(Mfit)-2.0)

    if verbosity > 1:
        print "Linear Photmetric Solution: "+str(popt)
        print "Errors: "+str(perr)
        print "Fit X2/dof: "+str(X2dof)
        print "Writing photometric fit parameters"
    #get base of catalog name for output
    base = catname[:-4]
    
    #output fit parameters in txt file
    s = "# Photometric linear fit: maps flux to "+band+" band magnitude."
    s += "\n# catalog stars from aavso taken from:\n"
    if refname is not None:
        s += '# '+refname
    else:
        s += '# AAVSO-APASS-DR9_ra{:4.6f}dec{:4.6f}bm{:4.1f}x{:4.1f}.cat'.format(radeg,decdeg,fovam,fovam)
    s += '\n# fit X2/dof = '+str(X2dof)[:4]
    s += '\n# magnitude = a*flux + b'
    s += '\n# row 1: a\tb'
    s += '\n# row 2: errors'
    s += '\n# ----------------------'
    s += "\n  "+'\t'.join("%.4f" % x for x in popt)+'\n  '+'\t'.join("%.4f" % x for x in perr)
    f = open(base+'.phot','w')
    f.write(s)
    f.close()

    if verbosity > 0:
        if verbosity > 1:
            print "Writing Catalog with corrected magnitudes"

        #compute corrected magnitudes for all objects in SExtractor catalog
        M, M_err = linMatch(F, popt[0], popt[1], F_err, perr[0], perr[1])
        #output fit parameters in txt file
        s = "Photometric magnitudes in "+band+" band."
        s += '\ncalculated from SExtractor fluxes.'
        s += "\nused photometric solution from:"
        s += '\n'+base+'.phot'
        s += '\nmagnitude = a*flux + b'
        s += '\nrow: M\tM_err\tRA\tDEC'
        s += '\n---------------------------------'
        #output catalog file with corrected magnitudes
        np.savetxt(base+".txt", np.array([M, M_err, RA, DEC]).T, fmt='%.3f', delimiter='\t', header=s)
        
    if verbosity > 1:

        import matplotlib.pyplot as plt 
        
        #plot matching star field
        plt.title("AAVSO matched to SExtractor 14<m<16")
        plt.plot(RAm, DECm, 'b.', label='source')
        plt.plot(ram, dem, 'r+', label='AAVSO') 
        plt.locator_params(axis='x',nbins=4) 
        plt.locator_params(axis='y',nbins=4) 
        plt.tick_params('x',pad=10) 
        plt.xlabel('x pixel') 
        plt.ylabel('y pixel') 
        plt.ticklabel_format(useOffset=False) 
        plt.axis('scaled') 
        plt.legend()
        #plt.ylim(0,9200)
        #plt.xlim(0,9200)
        plt.show()
        
        #plot fit residuals
        plt.plot(Fm, magm, 'r+')
        plt.errorbar(Ffit, Mfit, yerr=Mfit_err, label="fit\nm="+str(popt[0])[:7]+"pm"+str(perr[0])[:6]+"\nb="+str(popt[1])[:6]+"pm"+str(perr[1])[:5]+"\nX2dof="+str(X2dof)[:4])
        plt.xlabel("SExtractor")
        plt.ylabel("AAVSO-APASS-DR9")
        plt.axis('equal')
        plt.legend()
        plt.show()

    #return photometric fit parameters
    return popt, perr
        
#command line execution
if __name__ == "__main__":
    
    import argparse
    
    #command line arguments
    parser = argparse.ArgumentParser(description="Find Photometric Magnitude")
    parser.add_argument("catname", type=str, help="catalog whose photometry needs to be matched to AAVSO")
    parser.add_argument("band", type=str, help="catalog band to be matched (B, V or I)")
    parser.add_argument("-fits", "--fitsname", type=str, help="Program searches online using FITS header location to find suitable catalog.")
    parser.add_argument("-ref", "--refname", type=str, help="AAVSO reference catalog file from Vizier.")
    parser.add_argument("-p", "--position", type=str, help="RA:DEC as deg:deg")
    parser.add_argument("-v", "--verbosity", action="count", default=0, help="Default outputs photometric fit parameters into .phot file. -v writes SExtractor catalog magnitudes into .txt file. -vv outputs noisy diagnostic messages and plots.")
    args = parser.parse_args()

    if args.refname is not None:
        #compute photometric solution
        popt, perr = matchPhot(args.catname, args.band, refname=args.refname, verbosity=args.verbosity)
        print "Linear Photmetric Solution: "+str(popt)
        print "Errors: "+str(perr)
    elif args.fitsname is not None:
        #compute photometric solution
        popt, perr = matchPhot(args.catname, args.band, fitsname=args.fitsname, verbosity=args.verbosity)
        print "Linear Photmetric Solution: "+str(popt)
        print "Errors: "+str(perr)
    elif args.position is not None:
        #extract RA, DEC from position argument
        RA, DEC = [float(coord) for coord in args.position.split(':')]
        pos = [RA, DEC]
        #compute photometric solution
        popt, perr = matchPhot(args.catname, args.band, pos=pos, verbosity=args.verbosity)
        print "Linear Photmetric Solution: "+str(popt)
        print "Errors: "+str(perr)
    else:
        print "Insufficient information to retrieve AAVSO catalog, try -h flag and give one of three optional arguments."
