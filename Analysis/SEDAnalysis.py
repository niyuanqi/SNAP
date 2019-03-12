#################################################################
# Name:     SEDAnalysis.py                                      #
# Author:   Yuan Qi Ni                                          #
# Version:  Mar. 12, 2019                                       #
# Function: Program contains various models to analyse and fit  #
#           SED constructed from photometry.                    #
#################################################################

#essential modules
import numpy as np

#################################################################
# Constructing Photometric SED                                  #
#################################################################

#function: mask a list
def lmask(lst, msk):
    return [lst[i] for i in xrange(len(msk)) if msk[i]]

#function: SED gaussian interpolator
def SEDinterp(t, bands, SED_ts, SED_lcs, SED_errs=None,
              bandmask=None, interp='Gauss'):

    from scipy.stats import norm
    
    #check if number of light curves given is same as number of wavelengths
    Nw = len(bands)
    if len(SED_ts)!= Nw or len(SED_lcs) != Nw:
        print "Must give one light curve for each wavelength specified."
        return
    if bandmask is None:
        bandmask = np.ones(Nw).astype(bool)
    #check coverage
    for i in range(Nw):
        if t > max(SED_ts[i]) or t < min(SED_ts[i]):
            print "Incomplete coverage in lc "+str(i)+", masking..."
            bandmask[i] = False
    #check if some bands are masked
    SED_ts = lmask(SED_ts, bandmask)
    SED_lcs = lmask(SED_lcs, bandmask)
    if SED_errs is not None:
        SED_errs = lmask(SED_errs, bandmask)
    Nm = int(sum(bandmask))
    
    #Create 1D array over wavelength
    fluxes = np.zeros(Nm)
    for i in range(Nm):
        if interp == 'Gauss':
            #blur by gaussian function with std=1day
            gauss_win = norm.pdf(SED_ts[i], t, 1.0)
            fluxes[i] = np.sum(gauss_win*SED_lcs[i])/np.sum(gauss_win)
        elif interp == 'lin':
            fluxes[i] = np.interp(t, SED_ts[i], SED_lcs[i])
    if SED_errs is not None:
        flux_errs = np.zeros(Nm)
        for i in range(Nm):
            if interp == 'Gauss':
                #blur by gaussian function with std=1day
                gauss_win = norm.pdf(SED_ts[i], t, 1.0)
                flux_errs[i] = np.sum(gauss_win*SED_errs[i])/np.sum(gauss_win)
            elif interp == 'lin':
                flux_errs[i] = np.interp(t, SED_ts[i], SED_errs[i])

    #Return SED
    if SED_errs is None:
        return bandmask, fluxes
    else:
        return bandmask, fluxes, flux_errs

#################################################################
# Calibrating Spectroscopy Using Photometry                     #
#################################################################

#function: Chi squared difference between spectrum and photometric SED
def err_specSED(spec, SED_wave, SED_flux, SED_flux_err):
    spec_flux = spec(SED_wave, flux_unit='jy')
    return (spec_flux - SED_flux)/SED_flux_err

#function: calibrate spectrum using photometric SED
def SEDcalib(p0, spec, SED_filts, SED_flux, SED_flux_err):

    from scipy.optimize import leastsq
    from SNAP.Analysis.SpecAnalysis import filter_flux
    
    #get flux from spectrum in each filter
    spec_flux = np.zeros(len(SED_filts))
    for i in range(len(SED_filts)):
        flux = filter_flux(SED_filts[i], spec)
        spec_flux[i] = flux
    errfunc = lambda a: (a*spec_flux - SED_flux)/SED_flux_err
    #fit using scipy least-squares minimizer
    popt, ier = leastsq(errfunc, p0, full_output=0, maxfev=100000)
    return popt[0]
