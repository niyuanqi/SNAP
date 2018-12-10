#################################################################
# Name:     MagPlot.py                                          #
# Author:   Yuan Qi Ni                                          #
# Version:  April 28, 2016                                      #
# Function: Program contains essential functions for plotting   #
#           instrumental magnitude relationships.               #
#################################################################

#essential imports
import numpy as np
import matplotlib.pyplot as plt

#function: plot intensity vs SNs
def sn_corr_plot(insMags, catSNs):
    plt.title("Measured SNR of reference stars")
    plt.scatter(insMags, np.log(catSNs), c='r')
    fit = np.polyfit(insMags, np.log(catSNs), 1)
    plt.plot(insMags, np.polyval(fit, insMags), zorder=2)
    plt.ylabel("log SNR")
    plt.xlabel("Mag (~log I)")
    plt.show()

#function: plot intensity vs noise
def noise_corr_plot(insMags, catNs):
    plt.title("Measured noise under reference stars")
    plt.scatter(insMags, np.log(catNs), c='r')
    fit = np.polyfit(insMags, np.log(catNs), 1)
    plt.plot(insMags, np.polyval(fit, insMags), zorder=2)
    plt.ylabel("log Noise")
    plt.xlabel("Mag (~log I)")
    plt.show()

#function: plot magnitude solution
def phot_sol(insMags, insMagerrs, catMags, catMagerrs):
    #essential extra imports
    from scipy.optimize import curve_fit
    from SNAP.Analysis.LCFitting import linfunc
    plt.title("Photometric solution of reference stars")
    plt.errorbar(catMags, insMags, xerr=catMagerrs, yerr=insMagerrs, fmt='r+', zorder=1)
    #fit photometric solution
    popt, pcov = curve_fit(linfunc,catMags,insMags,p0=[1,-29.1])
    perr = np.sqrt(np.diag(pcov))
    photsol = linfunc(catMags, *popt)
    print "Photometric solution:",popt, perr
    plt.plot(catMags, photsol, zorder=2)
    plt.ylabel("-2.5 log I")
    plt.xlabel("Mag (catalog)")
    plt.show()

#function: plot B band color correlation
def col_corr(cat, catname, catIDs, insMags, insMagerrs, catMags, catMagerrs):
    #essential extra imports
    from scipy.optimize import curve_fit
    from SNAP.Analysis.LCFitting import linfunc
    import Catalog as ctlg
    #load V band data
    if cat == 'phot':
        IDV, RAV, DECV, catMV, catMerrV = ctlg.catPhot(catname,band='V')
    elif cat == 'dprs':
        IDV, RAV, DECV, catMV, catMerrV = ctlg.catDPRS(catname,band='V')
    elif cat == 'diff':
        IDV, RAV, DECV, catMV, catMerrV = ctlg.catDiff(catname,band='V')
    elif cat == 'aavso':
        fovam = 2.0*radius*0.4/60.0 #arcmin radius in KMT scaling
        IDV, RAV, DECV, catMV, catMerrV = ctlg.catAAVSO(RAo[0],DECo[0],fovam,'V',out=catname)
    #photometric solution color dependence
    plt.title("B band dependence on B-V")
    dI = catMags - insMags
    dI_err = np.sqrt(catMagerrs**2 + insMagerrs**2)
    V , Verr = [], []
    for i in range(len(catMV)):
        if IDV[i] in catIDs:
            V.append(catMV[i])
            Verr.append(catMerrV[i])
    BV = catMags - V
    BV_err = np.sqrt(catMagerrs**2 + np.square(Verr))
    plt.errorbar(BV, dI, xerr=BV_err, yerr=dI_err, fmt='r+', zorder=1)
    #fit color dependence
    popt, pcov = curve_fit(linfunc,BV,dI,p0=[0.27,28.1])
    perr = np.sqrt(np.diag(pcov))
    colsol = linfunc(BV, *popt)
    print "Color correlation:",popt, perr
    plt.plot(BV, colsol, zorder=2)
    plt.ylabel("B - b")
    plt.xlabel("B - V")
    plt.show()

#function: plot Chi2 histogram
def X2_hist(catX2dofs):
    plt.title("Reference star PSF fit qualities")
    nums, bins, patches = plt.hist(catX2dofs, bins=30)
    plt.xlabel("X2/dof")
    plt.ylabel("Counts")
    plt.title("Reference star fit qualities")
    plt.show()
