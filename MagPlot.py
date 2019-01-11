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
def col_corr(cat, catname, catIDs, RAo, DECo, radius, insMags, insMagerrs, catMags, catMagerrs):
    #essential extra imports
    from scipy.optimize import curve_fit
    from SNAP.Analysis.LCFitting import linfunc
    import Catalog as ctlg
    #load V band data
    if cat == 'aavso':
        fovam = 2.0*radius*0.4/60.0 #arcmin radius in KMT scaling
        IDBV, RABV, DECBV, catBV, catBVerr = ctlg.catAAVSO(RAo[0],DECo[0],fovam,'B-V',out=catname)
        B, Berr = [], []
        KB, KBerr = [], []
        BV, BV_err = [], []
        for i in range(len(catBV)):
            if IDBV[i] in catIDs:
                Bid = list(catIDs).index(IDBV[i])
                B.append(catMags[Bid])
                Berr.append(catMagerrs[Bid])
                KB.append(insMags[Bid])
                KBerr.append(insMagerrs[Bid])
                BV.append(catBV[i])
                BV_err.append(catBVerr[i])
        B, Berr = np.array(B), np.array(Berr)
        KB, KBerr = np.array(KB), np.array(KBerr)
        BV, BV_err = np.array(BV), np.array(BV_err)
    else:
        #fetch V band magnitudes
        if cat == 'phot':
            IDV, RAV, DECV, catMV, catMerrV = ctlg.catPhot(catname,band='V')
        elif cat == 'dprs':
            IDV, RAV, DECV, catMV, catMerrV = ctlg.catDPRS(catname,band='V')
        elif cat == 'diff':
            IDV, RAV, DECV, catMV, catMerrV = ctlg.catDiff(catname,band='V')
        #compute B-V
        B, Berr = [], []
        KB, KBerr = [], []
        V , Verr = [], []
        for i in range(len(catMV)):
            if IDV[i] in catIDs:
                Bid = list(catIDs).index(IDV[i])
                B.append(catMags[Bid])
                Berr.append(catMagerrs[Bid])
                KB.append(insMags[Bid])
                KBerr.append(insMagerrs[Bid])
                V.append(catMV[i])
                Verr.append(catMerrV[i])
        B, Berr = np.array(B), np.array(Berr)
        KB, KBerr = np.array(KB), np.array(KBerr)
        V, Verr = np.array(V), np.array(Verr)
        BV = B-V
        BV_err = np.sqrt(np.square(Berr) + np.square(Verr))
    #photometric solution color dependence
    plt.title("B band dependence on B-V")
    dI = B - KB
    dI_err = np.sqrt(Berr**2 + KBerr**2)
    plt.errorbar(BV, dI, xerr=BV_err, yerr=dI_err, fmt='r+', zorder=1)
    #average B-V color
    BV_err = [BV_err[i] if BV_err[i] > 0 else 0.0005 for i in range(len(BV))]
    w = 1/np.square(BV_err)
    BV_mean = np.sum(BV*w)/np.sum(w)
    BV_merr = np.sqrt(1/np.sum(w))
    print "Average color (B-V):", BV_mean, "+/-", BV_merr 
    #fit color dependence
    popt, pcov = curve_fit(linfunc,BV,dI,sigma=dI_err,p0=[0.27,28.1],absolute_sigma=True)
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
