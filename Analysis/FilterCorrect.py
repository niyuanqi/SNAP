#################################################################
# Name:     Scorrect.py                                         #
# Author:   Yuan Qi Ni                                          #
# Version:  Feb. 6, 2019                                        #
# Function: Program contains various routines for manipulating  #
#           multi-band light curve files and arrays.            #
#           Primarily for correcting filter discrepancy related #
#           magnitude - color dependencies using linear fit or  #
#           spectral correction methods.                        #
#################################################################

#essential modules
import numpy as np

def neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.log_likelihood(y)

#function: return corrected B band magnitude based Spectral Corrections
def SBcorrectMag(ts, mags, errs, tcorr, Scorr, tdiv=0, interp='GP',
                 Bcol=0, Vcol=1, SBVega=0, mBVr=0, mBVrerr=0,
                 Sinterp='GP', Scorr_err=None):
    '''
    #######################################################################
    # Input                                                               #
    # ------------------------------------------------------------------- #
    #   mags: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitudes in float.              #
    #                                                                     #
    #   errs: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitude errors in float.        #
    #                                                                     #
    #     ts: list of time arrays (eg. [tB, tV, tI]) where each is an     #
    #         array of time (in float) corresponding to the light curve.  #
    #                                                                     #
    #   Bcol: index of B band column                                      #
    #   Vcol: index of V band column                                      #
    #                                                                     #
    #   tdiv: dividing time, before which apply spectral information is   #
    #         insufficient, and we should apply simple BV correction.     #
    #  tcorr: epochs at which S corrections were measured                 #
    #  Scorr: spectral correction measured from spectra                   #
    # SBVega: spectral correction for vega
    #                                                                     #
    #   mBVr: mean color of reference stars used mBVr=<B-V>_r             #
    #mBVrerr: error in above color                                        #
    # ------------------------------------------------------------------- #
    # Output                                                              #
    # ------------------------------------------------------------------- #
    #    tcs: list of time arrays.                                        #
    #                                                                     #
    #  magcs: list of corrected light curves.                             #
    #                                                                     #
    #  errcs: list of corrected errors.                                   #
    #######################################################################
    '''
    import copy
    tcs, magcs, errcs = copy.deepcopy(ts), copy.deepcopy(mags), copy.deepcopy(errs)
    #B band correlation with B-V
    c = 0.27
    #B band to be corrected
    Bin, Bin_err = mags[Bcol], errs[Bcol]
    Bout, Bout_err = magcs[Bcol], errcs[Bcol]
    #mask times over which Scorrs are valid
    mask = [ts[Bcol]>=tdiv]

    if Sinterp == 'GP':
        from scipy.optimize import minimize
        import george
        from george import kernels

        # matern kernel in time-band space
        rt = 3.0
        mu, sigma = np.mean(Scorr), np.sqrt(np.var(Scorr))
        kernel = sigma*kernels.Matern32Kernel(metric=[rt], ndim=1)
        #initialize gaussian process
        gp = george.GP(kernel, mean=mu)
        gp.compute(tcorr, Scorr_err)  # You always need to call compute once.
        initial_params = gp.get_parameter_vector()
        bounds = gp.get_parameter_bounds()
        #train gaussian process
        r = minimize(neg_log_like, initial_params, method="L-BFGS-B",
                     bounds=bounds, args=(Scorr, gp))
        gp.set_parameter_vector(r.x)
        gp.get_parameter_dict()
        print r.x
        #predict using gaussian process
        scorr, scorr_var = gp.predict(Scorr, tcs[Bcol][mask])
    else:
        #correct B band using Bout = Bin + Scorr
        scorr = np.interp(tcs[Bcol][mask],tcorr,Scorr)
    #correct B band using S correction
    Bout[mask] = Bin[mask] + scorr - SBVega - c*mBVr
    Bout_err[mask] = np.sqrt(Bin_err[mask]**2 + (c*mBVrerr)**2)
    
    #mask times over which Scorrs are invalid
    mask = [ts[Bcol]<tdiv]
    #correct B band using Bout = (Bout-Vin)*c + Bin
    if interp == 'GP':
        from SEDAnalysis import SEDinterp
        #Construct V band Gaussian Process interpolator
        gp = SEDinterp(ts[Vcol][0], ['V'], [ts[Vcol]],
                       [mags[Vcol]], [errs[Vcol]], retGP=True)[0]
        Vin, Vin_var = gp.predict(mags[Vcol], ts[Bcol][mask])
        Vin_err = np.sqrt(np.diag(Vin_var))
    else:
        #Interpolate linearly
        Vin = np.interp(ts[Bcol][mask], ts[Vcol], mags[Vcol])
        Vin_err = np.interp(ts[Bcol][mask], ts[Vcol], errs[Vcol])
    Bout[mask] = (Bin[mask] - c*Vin - c*mBVr)/(1.-c)
    Bout_err[mask] = np.sqrt(np.square(c*Vin_err)+np.square(Bin_err[mask])
                             +np.square(c*mBVrerr))/(1.-c)

    #return corrected magnitudes
    magcs[Bcol] = Bout
    errcs[Bcol] = Bout_err
    #return corrected light curves
    return tcs, magcs, errcs

#function: return corrected B band magnitude based on V band correlation
def BVcorrectMag(ts, mags, errs, interp='GP', Bcol=0, Vcol=1, mBVr=0, mBVrerr=0):
    '''
    #######################################################################
    # Input                                                               #
    # ------------------------------------------------------------------- #
    #   mags: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitudes in float.              #
    #                                                                     #
    #   errs: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitude errors in float.        #
    #                                                                     #
    #     ts: list of time arrays (eg. [tB, tV, tI]) where each is an     #
    #         array of time (in float) corresponding to the light curve.  #
    #                                                                     #
    #   Bcol: index of B band column                                      #
    #   Vcol: index of V band column                                      #
    #                                                                     #
    #   mBVr: mean color of reference stars used mBVr=<B-V>_r             #
    #mBVrerr: error in above color                                        #
    # ------------------------------------------------------------------- #
    # Output                                                              #
    # ------------------------------------------------------------------- #
    #    tcs: list of time arrays.                                        #
    #                                                                     #
    #  magcs: list of corrected light curves.                             #
    #                                                                     #
    #  errcs: list of corrected errors.                                   #
    #######################################################################
    '''
    import copy
    tcs, magcs, errcs = copy.deepcopy(ts), copy.deepcopy(mags), copy.deepcopy(errs)
    #B band correlation with B-V
    c = 0.27
    #correct B band using Bout = (Bout-Vin)*c + Bin
    Bin, Bin_err = mags[Bcol], errs[Bcol]
    if Sinterp == 'GP':
        from SEDAnalysis import SEDinterp
        #Construct V band Gaussian Process interpolator
        gp = SEDinterp(ts[Vcol][0], ['V'], [ts[Vcol]],
                       [mags[Vcol]], [errs[Vcol]], retGP=True)[0]
        Vin, Vin_var = gp.predict(mags[Vcol], ts[Bcol][mask])
        Vin_err = np.sqrt(np.diag(Vin_var))
    else:
        #Interpolate linearly
        Vin = np.interp(ts[Bcol][mask], ts[Vcol], mags[Vcol])
        Vin_err = np.interp(ts[Bcol][mask], ts[Vcol], errs[Vcol])    
    Bout = (Bin - c*Vin - c*mBVr)/(1.-c)
    Bout_err = np.sqrt(np.square(c*Vin_err)+np.square(Bin_err)
                       +np.square(c*mBVrerr))/(1.-c)
    magcs[Bcol] = Bout
    errcs[Bcol] = Bout_err
    #return corrected light curves
    return tcs, magcs, errcs

#function: return corrected B band limiting magnitude based reference star B-V
def BVcorrectLim(ts, lims, Bcol=0, mBVr=0):
    '''
    #######################################################################
    # Input                                                               #
    # ------------------------------------------------------------------- #
    #   lims: list of det limits (eg. in different bands [B, V, I])       #
    #         where each is an array of magnitudes in float.              #
    #                                                                     #
    #     ts: list of time arrays (eg. [tB, tV, tI]) where each is an     #
    #         array of time (in float) corresponding to the light curve.  #
    #                                                                     #
    #   Bcol: index of B band column                                      #
    #                                                                     #
    #   mBVr: mean color of reference stars used mBVr=<B-V>_r             #
    # ------------------------------------------------------------------- #
    # Output                                                              #
    # ------------------------------------------------------------------- #
    #    tcs: list of time arrays.                                        #
    #                                                                     #
    #  limcs: list of corrected det limits.                               #
    #######################################################################
    '''
    import copy
    tcs, limcs = copy.deepcopy(ts), copy.deepcopy(lims)
    #B band correlation with B-V
    c = 0.27
    #correct B band using Bout = (B-Vin)*c + Bin
    Bin = lims[Bcol]
    Bout = Bin - c*mBVr
    limcs[Bcol] = Bout
    #return corrected light curves
    return tcs, limcs

#function: return corrected B band magnitude based on V band correlation
#only to be applied over short times, where B-V doesn't vary much
def BVcorrectFlux(ts, mags, errs, te, Fe, SNe, plot=True):
    '''
    #######################################################################
    # Input                                                               #
    # ------------------------------------------------------------------- #
    #   mags: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitudes in float.              #
    #                                                                     #
    #   errs: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitude errors in float.        #
    #                                                                     #
    #     ts: list of time arrays (eg. [tB, tV, tI]) where each is an     #
    #         array of time (in float) corresponding to the light curve.  #
    # ------------------------------------------------------------------- #
    # Output                                                              #
    # ------------------------------------------------------------------- #
    #    tcs: list of time arrays.                                        #
    #                                                                     #
    #  magcs: list of corrected light curves.                             #
    #                                                                     #
    #  errcs: list of corrected errors.                                   #
    #######################################################################
    '''
    import copy
    from scipy.optimize import curve_fit
    from SNAP.Analysis.Cosmology import flux_0
    from SNAP.Analysis.LCFitting import linfunc
    from SNAP.Analysis.LCRoutines import LCcolors
    
    tce, Fce, SNce = copy.deepcopy(te), copy.deepcopy(Fe), copy.deepcopy(SNe)
    #B band correlation with B-V
    c = 0.27
    #correct B band
    Bin, Bin_err = Fe[0], Fe[0]/SNe[0]
    #compute B-V color
    tdiff, C, C_err = LCcolors(ts, mags, errs)
    tBV, BV, BV_err = tdiff[0], C[0], C_err[0]
    #take only relevant interval
    mask = np.logical_and(tBV>te[0][0], tBV<te[0][-1])
    tBV, BV, BV_err = tBV[mask], BV[mask], BV_err[mask]
    #interpolate
    popt, pcov = curve_fit(linfunc,tBV,BV,p0=[0.0,0.0],sigma=BV_err,absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    #color at flux epochs
    BVe = popt[0]*te[0] + popt[1]
    print popt, perr
    BVe_err = np.sqrt(perr[1]**2 + (te[0]*perr[0])**2)
    #plot fit
    if plot:
        import matplotlib.pyplot as plt
        print "Plotting color interpolation"
        plt.title("B-V fit interpolation")
        plt.errorbar(tBV, BV, yerr=BV_err, fmt='k-')
        plt.errorbar(te[0], BVe, yerr=BVe_err, fmt='g-')
        plt.xlabel("Time")
        plt.ylabel("B-V")
        plt.show()
    #correct flux
    Bout = Bin*np.power(10., (-c/(1.-c))*(BVe/2.5))
    Bout_err = Bout*np.sqrt((Bin_err/Bin)**2+(c*np.log(10)*BVe_err/(2.5*(1.-c)))**2)
    Fce[0] = Bout
    SNce[0] = Bout/Bout_err
    #return corrected light curves
    return tce, Fce, SNce

#function: return corrected B band magnitude based Spectral Corrections
def SIcorrectMag(ts, mags, errs, tcorr, Scorr, tdiv=0, interp='GP',
                 Icol=2, Vcol=1, SIVega=0, mVIr=0, mVIrerr=0,
                 Sinterp='GP', Scorr_err=None):
    '''
    #######################################################################
    # Input                                                               #
    # ------------------------------------------------------------------- #
    #   mags: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitudes in float.              #
    #                                                                     #
    #   errs: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitude errors in float.        #
    #                                                                     #
    #     ts: list of time arrays (eg. [tB, tV, tI]) where each is an     #
    #         array of time (in float) corresponding to the light curve.  #
    #                                                                     #
    #   Vcol: index of V band column                                      #
    #   Icol: index of I band column                                      #
    #                                                                     #
    #   tdiv: dividing time, before which apply spectral information is   #
    #         insufficient, and we should apply simple BV correction.     #
    #  tcorr: epochs at which S corrections were measured                 #
    #  Scorr: spectral correction measured from spectra                   #
    # SIVega: spectral correction for vega
    #                                                                     #
    #   mVIr: mean color of reference stars used mBVr=<B-V>_r             #
    #mVIrerr: error in above color                                        #
    # ------------------------------------------------------------------- #
    # Output                                                              #
    # ------------------------------------------------------------------- #
    #    tcs: list of time arrays.                                        #
    #                                                                     #
    #  magcs: list of corrected light curves.                             #
    #                                                                     #
    #  errcs: list of corrected errors.                                   #
    #######################################################################
    '''
    import copy
    tcs, magcs, errcs = copy.deepcopy(ts), copy.deepcopy(mags), copy.deepcopy(errs)
    #B band correlation with B-V
    c = 0.0
    #B band to be corrected
    Iin, Iin_err = mags[Icol], errs[Icol]
    Iout, Iout_err = magcs[Icol], errcs[Icol]
    #mask times over which Scorrs are valid
    mask = [ts[Icol]>=tdiv]
    if Sinterp == 'GP':
        from scipy.optimize import minimize
        import george
        from george import kernels

        # matern kernel in time-band space
        rt = 3.0
        mu, sigma = np.mean(Scorr), np.sqrt(np.var(Scorr))
        kernel = sigma*kernels.Matern32Kernel(metric=[rt], ndim=1)
        #initialize gaussian process
        gp = george.GP(kernel, mean=mu)
        gp.compute(tcorr, Scorr_err)  # You always need to call compute once.
        initial_params = gp.get_parameter_vector()
        bounds = gp.get_parameter_bounds()
        #train gaussian process
        r = minimize(neg_log_like, initial_params, method="L-BFGS-B",
                     bounds=bounds, args=(Scorr, gp))
        gp.set_parameter_vector(r.x)
        gp.get_parameter_dict()
        print r.x
        #predict using gaussian process
        scorr, scorr_var = gp.predict(Scorr, tcs[Icol][mask])
    else:
        #correct B band using Bout = Bin + Scorr
        scorr = np.interp(tcs[Icol][mask],tcorr,Scorr)
    #correct I band using S correction
    Iout[mask] = Iin[mask] + scorr - SIVega - c*mVIr
    Iout_err[mask] = np.sqrt(Iin_err[mask]**2 + (c*mVIrerr)**2)
    
    #mask times over which Scorrs are invalid
    mask = [ts[Icol]<tdiv]
    #correct B band using Bout = (Bout-Vin)*c + Bin
    if interp == 'GP':
        from SEDAnalysis import SEDinterp
        #Construct V band Gaussian Process interpolator
        gp = SEDinterp(ts[Vcol][0], ['V'], [ts[Vcol]],
                       [mags[Vcol]], [errs[Vcol]], retGP=True)[0]
        Vin, Vin_var = gp.predict(mags[Vcol], ts[Icol][mask])
        Vin_err = np.sqrt(np.diag(Vin_var))
    else:
        #Interpolate linearly
        Vin = np.interp(ts[Icol][mask], ts[Vcol], mags[Vcol])
        Vin_err = np.interp(ts[Icol][mask], ts[Vcol], errs[Vcol])
    Iout[mask] = (Iin[mask] + c*Vin - c*mVIr)/(1.+c)
    Iout_err[mask] = np.sqrt(np.square(c*Vin_err)+np.square(Iin_err[mask])
                             +np.square(c*mVIrerr))/(1.-c)

    #return corrected magnitudes
    magcs[Icol] = Iout
    errcs[Icol] = Iout_err
    #return corrected light curves
    return tcs, magcs, errcs

#function: return corrected B band limiting magnitude based reference star B-V
def VIcorrectLim(ts, lims, Icol=2, mVIr=0):
    '''
    #######################################################################
    # Input                                                               #
    # ------------------------------------------------------------------- #
    #   lims: list of det limits (eg. in different bands [B, V, I])       #
    #         where each is an array of magnitudes in float.              #
    #                                                                     #
    #     ts: list of time arrays (eg. [tB, tV, tI]) where each is an     #
    #         array of time (in float) corresponding to the light curve.  #
    #                                                                     #
    #   Icol: index of I band column                                      #
    #                                                                     #
    #   mVIr: mean color of reference stars used mVIr=<V-I>_r             #
    # ------------------------------------------------------------------- #
    # Output                                                              #
    # ------------------------------------------------------------------- #
    #    tcs: list of time arrays.                                        #
    #                                                                     #
    #  limcs: list of corrected det limits.                               #
    #######################################################################
    '''
    import copy
    tcs, limcs = copy.deepcopy(ts), copy.deepcopy(lims)
    #I band correlation with V-I
    c = 0.0
    #correct I band using Iout = (Vin-Iout)*c + Iin
    Iin = lims[Icol]
    Iout = Iin - c*mVIr
    limcs[Icol] = Iout
    #return corrected light curves
    return tcs, limcs
