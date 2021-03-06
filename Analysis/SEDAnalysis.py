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

#Gaussian process negative log likelihood minimizer
def neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.log_likelihood(y)

#function: SED gaussian interpolator
def SEDinterp(t, bands, SED_ts, SED_lcs, SED_errs=None,
              bandmask=None, interp='GP', retGP=False):

    from scipy.stats import norm
    from scipy.optimize import minimize
    if interp == 'GP':
        import george
        from george import kernels
    
    #check if number of light curves given is same as number of wavelengths
    Nw = len(bands)
    if len(SED_ts)!= Nw or len(SED_lcs) != Nw:
        print "Must give one light curve for each wavelength specified."
        return
    if bandmask is None:
        bandmask = np.ones(Nw).astype(bool)
        #check coverage
    for i in range(Nw):
        if bandmask[i] and (t > max(SED_ts[i]) or t < min(SED_ts[i])):
            print "Incomplete coverage in lc "+str(i)+", masking..."
            bandmask[i] = False
    #check if some bands are masked
    SED_ts = lmask(SED_ts, bandmask)
    SED_lcs = lmask(SED_lcs, bandmask)
    if SED_errs is not None:
        SED_errs = lmask(SED_errs, bandmask)
    Nm = int(sum(bandmask))

    #interpolate photometric SED
    if interp == 'GP':
        fluxes = np.zeros(Nm)
        flux_errs = np.zeros(Nm)
        #Gaussian Process requires errors to be given
        if SED_errs is None:
            print "GP option requires input errorbars."
            return
        gps = []
        for i in range(Nm):
            # matern kernel in time-band space
            rt = 3.0
            mu, sigma = np.mean(SED_lcs[i]), np.sqrt(np.var(SED_lcs[i]))
            kernel = sigma*kernels.Matern32Kernel(metric=[rt], ndim=1)
            #initialize gaussian process
            gp = george.GP(kernel, mean=mu)
            # You always need to call compute once.
            gp.compute(SED_ts[i], SED_errs[i])
            print("Initial log likelihood: {0}".format(
                gp.log_likelihood(SED_lcs[i])))
            initial_params = gp.get_parameter_vector()
            bounds = gp.get_parameter_bounds()
            #train gaussian process
            r = minimize(neg_log_like, initial_params, method="L-BFGS-B",
                         bounds=bounds, args=(SED_lcs[i], gp))
            gp.set_parameter_vector(r.x)
            gp.get_parameter_dict()
            print("GP trained parameters: {0}".format(r.x))
            gps.append(gp)
        if retGP:
            #return gaussian process
            return gps
        else:
            #predict using gaussian process
            for i in range(Nm):
                flux, flux_var = gps[i].predict(SED_lcs[i], t)
                flux_err = np.sqrt(np.diag(flux_var))
                fluxes[i] = flux
                flux_errs[i] = flux_err
    else:    
        #Create 1D array over wavelength
        fluxes = np.zeros(Nm)
        for i in range(Nm):
            if interp == 'linear':
                fluxes[i] = np.interp(t, SED_ts[i], SED_lcs[i])
            elif interp == 'nearest':
                i_near = np.argmin(np.square(SED_ts[i]-t))
                fluxes[i] = SED_lcs[i][i_near]
        if SED_errs is not None:
            flux_errs = np.zeros(Nm)
            for i in range(Nm):
                if interp == 'linear':
                    flux_errs[i] = np.interp(t, SED_ts[i], SED_errs[i])
                elif interp == 'nearest':
                    i_near = np.argmin(np.square(SED_ts[i]-t))
                    flux_errs[i] = SED_errs[i][i_near]
                
    #Return SED
    if SED_errs is None:
        return bandmask, fluxes
    else:
        return bandmask, fluxes, flux_errs

#function: integrate SED using trapezoidal rule
def SEDtrap(wave, flux, fluxerr=None, N=100):
    #convert A to Hz
    freq = 3.0e8/(wave*1e-10)
    if fluxerr is None:
        #simple trapezoidal rule
        return -1*np.trapz(flux, freq)
    else:
        #bootstrap SED fluxes
        SEDtrials = np.zeros((len(wave),N))
        #for each wave, draw N monte carlo values
        for i in range(len(wave)):
            #draw from gaussian centered at M, with sigma M_err
            SEDtrials[i] = np.random.normal(flux[i],np.absolute(fluxerr[i]),N)
        SEDtrials = SEDtrials.T

        #initializations
        integs = np.zeros(N)
        #For each SED, integrate
        for j, SED in enumerate(SEDtrials):
            integs[j] = -1*np.trapz(SED, freq)
            
        #return average integral and error
        integ_mean = np.mean(integs)
        integ_err = np.std(integs)
        return integ_mean, integ_err

#function: generate SED using blackbody
def genBlackbod(wave, T, r, Terr=None, rerr=None, N=100):
    
    if Terr is None:
        return planck(wave,T)*r
    else:
        #bootstrap T and r
        Ttrials = np.random.normal(T,np.absolute(Terr),N)
        rtrials = np.random.normal(r,np.absolute(rerr),N)

        #initializations
        SEDtrials = np.zeros((N, len(wave)))
        #for each trial, generate SED
        for j in range(N):
            SEDtrials[j] = planck(wave,Ttrials[j])*rtrials[j]

        #return average SED and error
        SED_mean = np.mean(SEDtrials, 0)
        SED_err = np.std(SEDtrials, 0)
        return SED_mean, SED_err

#function: generate Rayleigh-Jeans tail based off of given data.
def genRJtail(wave, phot_wave, flux, flux_err=None):
    
    #get best fit RJtail
    RJtail = lambda x, a: a/x**2
    a = flux*phot_wave**2
    if flux_err is not None:
        aerr = flux_err*phot_wave**2
    #evaluate RJtail
    tail_flux = RJtail(wave, a)
    tail_err = RJtail(wave, aerr)
    return tail_flux, tail_err

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
    popt, ier = leastsq(errfunc, p0, full_output=0, maxfev=1000000)
    return popt[0]

#################################################################
# Blackbody SED template                                        #
#################################################################

#function: jy to flam (erg/s/cm2/AA)
def Jy_toFLAM(x, flux_jy):
    #constants
    wave = x*1e-8 #angstrom to cm
    c = 2.99792458e10 #cm/s
    flux_flam = flux_jy*1e-23*c/wave**2/1e8
    return flux_flam

#function: flam to jy (1e-23 erg/s/cm2/Hz)
def FLAM_toJy(x, flux_flam):
    #constants
    wave = x*1e-8 #angstrom to cm
    c = 2.99792458e10 #cm/s
    flux_jy = flux_flam*1e8*wave**2/c/1e-23
    return flux_jy

#function: photlam to flam (erg/s/cm2/AA)
def PHOTLAM_toFLAM(x, flux_phot):
    #constants
    wave = x*1e-8 #angstrom to cm
    h = 6.6260755e-27 #erg*s
    c = 2.99792458e10 #cm/s
    flux_flam = flux_phot*h*c/wave
    return flux_flam
    
#function: stefan Boltzmann's law
def SBlaw(T):
    sb = 5.67051e-5 #erg/s/cm2/K4
    #black body total flux
    integ = sb*np.power(T,4) #ergs/s/cm2
    return integ

#function: black body distribution (wavelength)
def blackbod(x, T):
    #constants
    wave = x*1e-8 #angstrom to cm
    h = 6.6260755e-27 #erg*s
    c = 2.99792458e10 #cm/s
    k = 1.380658e-16 #erg/K
    freq = c/wave #Hz
    #planck distribution
    p_rad = (2*h*freq**3/c**2)/(np.exp(h*freq/(k*T))-1.0) #erg/s/cm2/rad2/Hz power per area per solid angle per frequency
    #integrated planck over solid angle
    p_int = np.pi*p_rad #erg/s/cm2/Hz #power per area per frequency [fnu]
    return p_int

#function: normalized planck distribution (wavelength)
def planck(x, T):
    #black body total flux
    integ = SBlaw(T) #erg/s/cm2
    #blackbody distribution
    p_int = blackbod(x,T) #erg/s/cm2/Hz #power per area per frequency [fnu]
    #normalized planck distribution
    return (p_int/integ) #1/Hz, luminosity density

#function: fit planck's law for black body temperature, received fraction
def fitBlackbod(waves, fluxes, fluxerrs=None, plot=False, ptitle=""):

    from scipy.optimize import curve_fit

    #blackbody flux function
    BBflux = lambda x, T, r : planck(x,T)*r
    
    #estimate temperature
    est = [10000.0, 1e14]
    est = [10000.0, 1e20]
    #fit blackbody temperature
    if fluxerrs is not None:
        popt, pcov = curve_fit(BBflux, waves, fluxes, sigma=fluxerrs, p0=est, absolute_sigma=True)
    else:
        popt, pcov = curve_fit(BBflux, waves, fluxes, p0=est, maxfev=10000000)
    perr = np.sqrt(np.diag(pcov))
    T, Terr = popt[0], perr[0] #K
    r, rerr = popt[1], perr[1] #dimensionless
    #plot fit if given
    if plot:
        import matplotlib.pyplot as plt

        print "Temperature [K]:", T, Terr
        print "Received/Emitted:", r, rerr
        if fluxerrs is not None:
            plt.errorbar(waves, fluxes, yerr=fluxerrs, fmt='g+')
        else:
            plt.plot(waves, fluxes, color='g')
        w = np.linspace(min(waves), max(waves), 100)
        plt.plot(w, BBflux(w, T, r), c='b',
                 label="T = {:.0f} ({:.0f}) K\nr = {:.3f} ({:.3f})".format(
                     T, Terr, r, rerr))
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux")
        plt.title(ptitle)
        plt.legend(loc='lower right')
        plt.tight_layout()
        plt.show()
    w = np.linspace(min(waves), max(waves), 100)
    #return blackbody temperature
    return T, Terr, r, rerr, w, BBflux(w, T, r)

#function: fit Rayleigh-Jeans tail
def fitRJtail(waves, fluxes, fluxerrs):

    from scipy.optimize import curve_fit

    #estimate constant in F = a(lambda)^-2
    a, aerr = fluxes*waves**2, fluxerrs*waves**2
    #take a weighted sum
    w = 1/np.square(aerr)
    a_mean = np.sum(w*a)/np.sum(w)
    a_err = np.sqrt(1/np.sum(w))
    return a_mean, a_err

#function: L, T fluxes derived using blackbody spectrum
def BBflux(Lc,Teff,wave,z,DM):
    #give wave in observer frame
    
    #luminosity distance [pc -> cm]
    dl = 10*np.power(10, DM/5.0)*3.086*10**18
    Area = 4.0*np.pi*np.square(dl) #cm^2
    #kasen model in observer band
    Lc_wave = planck(wave/(1.0+z),Teff)*Lc/Area
    #Lc_angle_wave = np.nan_to_num(Lc_angle_wave)
    #ergs/s/Hz/cm^2, luminosity density in observer frame
    #return in uJy, 10**29 uJy = 1 ergs/s/Hz/cm^2
    return Lc_wave*10**29

#function: color temperature
def colorT(color, band1, band2, z=0):

    from SNAP.Analysis.Cosmology import wave_0, bands, Flux_toMag

    #trial temperatures
    N = 10000
    Ts = np.logspace(3, 5, N) #K
    Ls = np.ones(N)
    #1-2 color for each temperature
    Tflux1 = blackbod(wave_0[bands[band1]]/(1.0+z),Ts)
    Tflux2 = blackbod(wave_0[bands[band2]]/(1.0+z),Ts)
    Tcolors = -2.5*np.log10(Tflux1/Tflux2)
    #Tflux1 = BBflux(Ls, Ts, wave_0[bands[band1]], z, 1) #uJy
    #Tflux2 = BBflux(Ls, Ts, wave_0[bands[band2]], z, 1) #uJy
    #Tmag1 = Flux_toMag(band1, Tflux1*1e-6) #mag
    #Tmag2 = Flux_toMag(band2, Tflux2*1e-6) #mag
    #Tcolors = Tmag1-Tmag2
    #closest color
    i_min = np.argmin(np.abs(Tcolors-color))
    return Ts[i_min]

#function: color temperature and radius
def colorTR(mag1, mag2, band1, band2, z, DM):
    
    from SNAP.Analysis.Cosmology import wave_0, bands, Mag_toFlux
    
    #luminosity distance [pc -> cm]
    dl = 10*np.power(10, DM/5.0)*3.086*10**18
    Area = 4.0*np.pi*np.square(dl) #cm^2
    #calculate luminosity density in band1
    flux_den1 = Mag_toFlux(band1, mag1) #Jy <-> 1e-23 ergs/s/Hz/cm^2
    lum_den1 = flux_den1*Area*1e-23 #ergs/s/Hz
    flux_den2 = Mag_toFlux(band2, mag2) #Jy <-> 1e-23 ergs/s/Hz/cm^2
    lum_den2 = flux_den2*Area*1e-23 #ergs/s/Hz

    
    #calculate color temperature
    color = mag1 - mag2
    T = colorT(color, band1, band2, z)
    #predicted flux density in band1
    flux_den1 = blackbod(wave_0[bands[band1]]/(1.+z),T) #erg/s/cm2/Hz 
    flux_den2 = blackbod(wave_0[bands[band2]]/(1.+z),T) #erg/s/cm2/Hz 
    #calculate radius
    area1 = lum_den1/flux_den1 #cm^2
    area2 = lum_den2/flux_den2 #cm^2
    area = 0.5*(area1+area2)
    #area = area2
    R = np.sqrt(area/(4*np.pi)) #cm

    #return color temperature and radius
    return T, R
#function: color radius
def colorR(T, mag1, mag2, band1, band2, z, DM):
    
    from SNAP.Analysis.Cosmology import wave_0, bands, Mag_toFlux
    
    #luminosity distance [pc -> cm]
    dl = 10*np.power(10, DM/5.0)*3.086*10**18
    Area = 4.0*np.pi*np.square(dl) #cm^2
    #calculate luminosity density in band1
    flux_den1 = Mag_toFlux(band1, mag1) #Jy <-> 1e-23 ergs/s/Hz/cm^2
    lum_den1 = flux_den1*Area*1e-23 #ergs/s/Hz
    flux_den2 = Mag_toFlux(band2, mag2) #Jy <-> 1e-23 ergs/s/Hz/cm^2
    lum_den2 = flux_den2*Area*1e-23 #ergs/s/Hz
    
    #predicted flux density in band1
    flux_den1 = blackbod(wave_0[bands[band1]]/(1.+z),T) #erg/s/cm2/Hz 
    flux_den2 = blackbod(wave_0[bands[band2]]/(1.+z),T) #erg/s/cm2/Hz 
    #calculate radius
    area1 = lum_den1/flux_den1 #cm^2
    area2 = lum_den2/flux_den2 #cm^2
    area = 0.5*(area1+area2)
    #area = area2
    R = np.sqrt(area/(4*np.pi)) #cm

    #return color temperature and radius
    return R
#function: get blackbody function from temperature and radius
def TRblackbod(T, R, z, DM, EBV=0):
    
    from SNAP.Analysis.Cosmology import wave_0, bands, Mag_toFlux
    import extinction as extn
    import astropy.units as u

    #luminosity distance [pc -> cm]
    dl = 10*np.power(10, DM/5.0)*3.086*10**18
    Area = 4.0*np.pi*np.square(dl) #cm^2
    area = 4.0*np.pi*R**2
    #predicted flux density at wave
    flux_obs = lambda wave: extn.apply(extn.fm07(wave*u.AA, EBV*3.1), 
                                       blackbod(wave/(1.+z),T)*1e23*area/Area) #Jy
    #flux_obs = lambda wave: blackbod(wave/(1.+z),T)*1e23*area/Area #Jy
    return flux_obs

#function: fit blackbody with extinction
def fitExtBlackbod(waves, fluxes, fluxerrs=None, EBV=None, plot=False, ptitle=""):

    import extinction as extn
    from scipy.optimize import curve_fit
    import astropy.units as u

    if EBV is None:
        #blackbody+extinction(Rv=3.1) flux function
        exBBflux = lambda x, T, r, EBV: extn.apply(extn.fm07(x*u.AA, EBV*3.1), planck(x,T)*r)
        #estimate temperature
        est = [10000.0, 1e20, 10]
    else:
        exBBflux = lambda x, T, r: extn.apply(extn.fm07(x*u.AA, EBV*3.1), planck(x,T)*r)
        est = [10000.0, 1e25]
        if EBV < 0:
            est = [1000.0, 1e18]
        else:
            est = [10000.0, 1e20]

    #fit blackbody temperature
    if fluxerrs is not None:
        popt, pcov = curve_fit(exBBflux, waves, fluxes, sigma=fluxerrs, p0=est, absolute_sigma=True, maxfev=1000000)
    else:
        popt, pcov = curve_fit(exBBflux, waves, fluxes, p0=est, maxfev=10000000)
    perr = np.sqrt(np.diag(pcov))
    T, Terr = popt[0], perr[0] #K
    r, rerr = popt[1], perr[1] #dimensionless
    if EBV is None:
        EBV, EBVerr = popt[2], perr[2] #mag
        print "E(B-V) for Rv=3.1:", EBV, EBVerr
    else:
        EBVerr = 0
    #plot fit if given
    if plot:
        import matplotlib.pyplot as plt

        print "Temperature [K]:", T, Terr
        print "Received/Emitted:", r, rerr
        if fluxerrs is not None:
            plt.errorbar(waves, fluxes, yerr=fluxerrs, fmt='g+')
        else:
            plt.plot(waves, fluxes, color='g')
        w = np.linspace(min(waves), max(waves), 100)
        plt.plot(w, exBBflux(w, *popt), c='r',
                 label="T = {:.0f} ({:.0f}) K\nr = {:.3f} ({:.3f})\nEBV = {:.3f} ({:.3f})".format(
                     T, Terr, r, rerr, EBV, EBVerr))
        #plt.plot(w, extn.remove(extn.fm07(w*u.AA, EBV*3.1), exBBflux(w, *popt)), c='b',
        #         label="dereddened")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux")
        plt.title(ptitle)
        plt.legend(loc='lower right')
        plt.tight_layout()
        plt.show()
    w = np.linspace(min(waves), max(waves), 100)
    #return blackbody temperature
    return T, Terr, r, rerr, EBV, EBVerr, w, exBBflux(w, *popt)
