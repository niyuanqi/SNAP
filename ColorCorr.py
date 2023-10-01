#################################################################
# Name:     ColorCorr.py                                        #
# Author:   Yuan Qi Ni                                          #
# Version:  July 13, 2020                                       #
# Function: Program contains functions which demonstrates the   #
#           color - intrumental magnitude correlation.          #
#           See Park et al. 2017 for details.                   #
#################################################################

#essential imports
import numpy as np

#function: plot B band color correlation
def Bcol_corr(cat, catname, catIDs, RAo, DECo, radius, insMags, insMagerrs, catMags, catMagerrs, plot=True):
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
    dI = B - KB
    dI_err = np.sqrt(Berr**2 + KBerr**2)
    #average B-V color
    BV_err = [BV_err[i] if BV_err[i] > 0 else 0.0005 for i in range(len(BV))]
    w = 1/np.square(BV_err)
    BV_mean = np.sum(BV*w)/np.sum(w)
    BV_merr = np.sqrt(1/np.sum(w))
    print "Average color (B-V):", BV_mean, "+/-", BV_merr 
    
    #make color correlation plot
    if plot:
        #essential additional import
        import matplotlib.pyplot as plt

        #fit color dependence
        plt.title("B band dependence on B-V")
        plt.errorbar(BV, dI, xerr=BV_err, yerr=dI_err, fmt='k+', zorder=1)
        popt, pcov = curve_fit(linfunc,BV,dI,p0=[0.27,27.8],
                               sigma=dI_err,absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        colsol = linfunc(BV, *popt)
        #mask out 3sig deviators 
        mask = np.absolute(dI-colsol) < 3*np.std(dI-colsol)
        plt.scatter(BV[mask], dI[mask], c='r')
        popt, pcov = curve_fit(linfunc,BV[mask],dI[mask],p0=[0.27,27.8],
                               sigma=dI_err[mask],absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        colsol = linfunc(BV[mask], *popt)
        print "Color correlation:",popt, perr
        print "Nstar:",len(BV[mask])
        print "Pearson:",np.corrcoef(BV[mask],dI[mask])
        plt.plot(BV[mask], colsol, zorder=2)
        plt.ylabel("B - inst")
        plt.xlabel("B - V")
        plt.show()
        
    #return mean color of reference stars
    return BV_mean, BV_merr 
    

#function: plot B band color correlation
def Icol_corr(cat, catname, catIDs, RAo, DECo, radius, insMags, insMagerrs, catMags, catMagerrs, plot=True):
    #essential extra imports
    from scipy.optimize import curve_fit
    from SNAP.Analysis.LCFitting import linfunc
    import Catalog as ctlg
    #fetch V band magnitudes
    if cat == 'phot':
        IDV, RAV, DECV, catMV, catMerrV = ctlg.catPhot(catname,band='V')
    elif cat == 'dprs':
        IDV, RAV, DECV, catMV, catMerrV = ctlg.catDPRS(catname,band='V')
    elif cat == 'diff':
        IDV, RAV, DECV, catMV, catMerrV = ctlg.catDiff(catname,band='V')
    elif cat == 'aavso':
        fovam = 2.0*radius*0.4/60.0 #arcmin radius in KMT scaling
        IDV, RAV, DECV, catMV, catMerrV = ctlg.catAAVSO(RAo[0],DECo[0],fovam,'V',out=catname)
    #compute V-I
    I, Ierr = [], []
    KI, KIerr = [], []
    V , Verr = [], []
    for i in range(len(catMV)):
        if IDV[i] in catIDs:
            Iid = list(catIDs).index(IDV[i])
            I.append(catMags[Iid])
            Ierr.append(catMagerrs[Iid])
            KI.append(insMags[Iid])
            KIerr.append(insMagerrs[Iid])
            V.append(catMV[i])
            Verr.append(catMerrV[i])
    I, Ierr = np.array(I), np.array(Ierr)
    KI, KIerr = np.array(KI), np.array(KIerr)
    V, Verr = np.array(V), np.array(Verr)
    VI = V-I
    VI_err = np.sqrt(np.square(Verr) + np.square(Ierr))
    #photometric solution color dependence
    dI = I - KI
    dI_err = np.sqrt(Ierr**2 + KIerr**2)
    #average B-V color
    VI_err = [VI_err[i] if VI_err[i] > 0 else 0.0005 for i in range(len(VI))]
    w = 1/np.square(VI_err)
    VI_mean = np.sum(VI*w)/np.sum(w)
    VI_merr = np.sqrt(1/np.sum(w))
    print "Average color (V-I):", VI_mean, "+/-", VI_merr 

    #make color correlation plot
    if plot:
        #essential additional import
        import matplotlib.pyplot as plt
        
        #fit color dependence
        plt.title("I band dependence on V-I")
        plt.errorbar(VI, dI, xerr=VI_err, yerr=dI_err, fmt='k+', zorder=1)
        popt, pcov = curve_fit(linfunc,VI,dI,p0=[0.27,27.8],
                               sigma=dI_err,absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        colsol = linfunc(VI, *popt)
        #mask out 3sig deviators 
        mask = np.absolute(dI-colsol) < 3*np.std(dI-colsol)
        plt.scatter(VI[mask], dI[mask], c='r')
        popt, pcov = curve_fit(linfunc,VI[mask],dI[mask],p0=[0.27,27.8],
                               sigma=dI_err[mask],absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        colsol = linfunc(VI[mask], *popt)
        print "Color correlation:",popt, perr
        print "Nstar:",len(VI[mask])
        print "Pearson:",np.corrcoef(VI[mask],dI[mask])
        plt.plot(VI[mask], colsol, zorder=2)
        plt.ylabel("i - inst")
        plt.xlabel("V - i")
        plt.show()

    #return mean color of reference stars
    return VI_mean, VI_merr 
