#################################################################
# Name:     SpecAnalysis.py                                     #
# Author:   Yuan Qi Ni                                          #
# Version:  Feb. 19, 2019                                       #
# Function: Program contains various models to analyse and fit  #
#           supernova spectra.                                  #
#################################################################

#essential modules
import numpy as np

#################################################################
# Processing Spectra using Synphot                              #
#################################################################

#function: loading spectra
def spec_fits(filename, get_err=False, get_meta=False):
    from specutils.io import read_fits
    import astropy.units as u
    from astropy.io import fits
    from synphot import SourceSpectrum
    from synphot.models import Empirical1D

    #load fits file
    spec_fits = read_fits.read_fits_spectrum1d(filename)
    #read header info
    spec_meta = fits.getheader(filename)
    #check for multispec
    if isinstance(spec_fits, list):
        spec = spec_fits[1]
        if get_err:
            spec_err = spec_fits[3]
        print spec_meta['SITEID']
    else:
        spec = spec_fits
    #spectrum flux data
    flux = spec.data
    if spec_meta['BUNIT']=='erg/cm2/s/A':
        flux = flux * u.erg/u.cm**2/u.s/u.AA
        if get_err:
            flux_err = spec_err.data
            flux_err = flux_err * u.erg/u.cm**2/u.s/u.AA
    wave = spec.dispersion
    if isinstance(wave, u.Quantity):
        wave = wave.value * u.AA
    #create synphot spectrum object
    spec = SourceSpectrum(Empirical1D, points=wave,
                          lookup_table=flux, keep_neg=True)
    if get_err:
        spec_err = SourceSpectrum(Empirical1D, points=wave,
                          lookup_table=flux_err, keep_neg=True)
    #check if metadata is needed
    if get_meta:
        from astropy.coordinates import SkyCoord
        
        #calculate time observed
        if 'T' in spec_meta['DATE-OBS']:
            isot_obs = spec_meta['DATE-OBS']
        elif 'UT' in spec_meta:
            isot_obs = spec_meta['DATE-OBS']+'T'+spec_meta['UT']
        else:
            isot_obs = spec_meta['DATE-OBS']+'T'+spec_meta['UT-TIME']
        #exposure time
        exp_obs = spec_meta['EXPTIME']
        #RA, DEC observed
        RA, DEC = spec_meta['RA'], spec_meta['DEC']
        coord = SkyCoord(RA, DEC, unit=(u.hourangle, u.deg))
        RA, DEC = coord.ra.degree, coord.dec.degree
        #return spectrum and metadata
        #return spectrum
        if get_err:
            return spec, spec_err, isot_obs, exp_obs, (RA,DEC)
        else:
            return spec, isot_obs, exp_obs, (RA,DEC)
    else:
        #return spectrum
        if get_err:
            return spec, spec_err
        else:
            return spec

#function: calculate zero point of filter system
def filter_vega_zp(filterfile, filterzero):
    '''
    ######################################################
    # Input                                              #
    # -------------------------------------------------- #
    # filterfile: file containing filter function        #
    # filterzero: flux corresponding to mag=0            #
    # -------------------------------------------------- #
    # Output                                             #
    # -------------------------------------------------- #
    # mag_0: magnitude of Vega in filter system.         #
    ######################################################
    '''
    from synphot import SourceSpectrum, SpectralElement, Observation
    
    #load Vega spectrum
    spec = SourceSpectrum.from_vega()
    #Load ascii filter function
    if filterfile.split('/')[-1][:6] == 'Bessel':
        filt = SpectralElement.from_file(filterfile, wave_unit='nm')
    else:
        filt = SpectralElement.from_file(filterfile, wave_unit='AA')
    wave = filt.waveset
    #Synthetic observation
    obs = Observation(spec, filt)
    flux = obs.effstim(flux_unit='jy', waverange=(wave[0],wave[-1]))
    #Calibrate to zero point (zp)
    mag_0 = -2.512*np.log10(flux/filterzero)
    return mag_0

#function: evaluates flux of spectrum in filter system
def filter_flux(filterfile, syn_spec, syn_err=None, wrange=None):
    '''
    ######################################################
    # Input                                              #
    # -------------------------------------------------- #
    # filterfile: file containing filter function        #
    #   syn_spec: synphot spectrum object of source      #
    #     wrange: (start,end) observed wavelength range  #
    # -------------------------------------------------- #
    # Output                                             #
    # -------------------------------------------------- #
    # flux: flux of source in filter system.             #
    ######################################################
    '''
    from synphot import SpectralElement, Observation, Empirical1D, SourceSpectrum
    import astropy.units as u
    
    #Load ascii filter function
    if filterfile.split('/')[-1][:6] == 'Bessel':
        filt = SpectralElement.from_file(filterfile, wave_unit='nm')
    else:
        filt = SpectralElement.from_file(filterfile, wave_unit='AA')
    if wrange is None:
        fwave = filt.waveset
        swave = syn_spec.waveset
        wrange = (max(fwave[0],swave[0]),min(fwave[-1],swave[-1]))
        wlengths = fwave[np.logical_and(fwave>wrange[0], fwave<wrange[1])]
    else:
        wlengths = fwave[np.logical_and(fwave>wrange[0], fwave<wrange[1])]
    #Synthetic observation
    obs = Observation(syn_spec, filt, force='extrap')
    #flux = obs.effstim(flux_unit='jy', waverange=wrange).value
    flux = obs.effstim(flux_unit='jy', wavelengths=wlengths).value
    #Synthetic observation of error spectrum
    if syn_err is not None:
        #square filter and error spectrum
        filt2 = SpectralElement(Empirical1D, points=fwave,
                          lookup_table=np.square(filt(fwave)))
        pseudo_flux = np.square(syn_err(swave,flux_unit='jy')).value
        syn_err2 = SourceSpectrum(Empirical1D, points=swave,
                                   lookup_table=pseudo_flux)
        #sum errors in quadrature
        obs = Observation(syn_err2, filt2, force='extrap')
        flux_err = np.sqrt(obs.effstim(waverange=wrange).value)
        return flux, flux_err
    else:
        return flux

#function: evaluates magnitudes of spectrum in filter system
def filter_mag(filterfile, filterzero, syn_spec, syn_err=None, wrange=None):
    '''
    ######################################################
    # Input                                              #
    # -------------------------------------------------- #
    # filterfile: file containing filter function        #
    # filterzero: flux [Jy] corresponding to mag=0       #
    #   syn_spec: synphot spectrum object of source      #
    #     wrange: (start,end) observed wavelength range  #
    # -------------------------------------------------- #
    # Output                                             #
    # -------------------------------------------------- #
    # mag: magnitude of source in filter system.         #
    ######################################################
    '''
    from synphot import SpectralElement, Observation, Empirical1D, SourceSpectrum
    import astropy.units as u
    
    #Load ascii filter function
    if filterfile.split('/')[-1][:6] == 'Bessel':
        filt = SpectralElement.from_file(filterfile, wave_unit='nm')
    else:
        filt = SpectralElement.from_file(filterfile, wave_unit='AA')
    if wrange is None:
        fwave = filt.waveset
        swave = syn_spec.waveset
        wrange = (max(fwave[0],swave[0]),min(fwave[-1],swave[-1]))
    #Synthetic observation
    obs = Observation(syn_spec, filt, force='extrap')
    flux = obs.effstim(flux_unit='jy', waverange=wrange).value
    #Calibrate magnitude with zero point
    mag = -2.512*np.log10(flux/filterzero)
    #Synthetic observation of error spectrum
    if syn_err is not None:
        #square filter and error spectrum
        filt2 = SpectralElement(Empirical1D, points=fwave,
                          lookup_table=np.square(filt(fwave)))
        pseudo_flux = np.square(syn_err(swave,flux_unit='jy')).value
        syn_err2 = SourceSpectrum(Empirical1D, points=swave,
                                   lookup_table=pseudo_flux)
        #sum errors in quadrature
        obs = Observation(syn_err2, filt2, force='extrap')
        flux_err = np.sqrt(obs.effstim(waverange=wrange).value)
        #convert to magnitude error
        mag_err = (2.5/np.log(10))*(flux_err/flux)
        return mag, mag_err
    else:
        return mag

#################################################################
# Spectral Filter Corrections                                   #
#################################################################

#function: evaluate S correction on Vega spectrum from filt1 to filt2
def Scorr_vega(filt1,filt2, zerof1,zerof2):
    from synphot import SourceSpectrum
    
    #load Vega spectrum
    spec = SourceSpectrum.from_vega()
    #Scorr = -2.5log(flux2/flux1) such that mag2 = mag1 + Scorr
    mag1 = filter_mag(filt1, zerof1, spec)
    mag2 = filter_mag(filt2, zerof2, spec)
    scorr = mag2 - mag1
    return scorr

#function: evaluate S correction on spec from filt1 to filt2
def Scorr(filt1,filt2, zerof1,zerof2, syn_spec, syn_err=None, wrange1=None,wrange2=None):
    #Scorr = -2.5log(flux2/flux1) such that mag2 = mag1 + Scorr
    if syn_err is not None:
        #if error spectrum is given
        mag1, err1 = filter_mag(filt1, zerof1, syn_spec, syn_err, wrange1)
        mag2, err2 = filter_mag(filt2, zerof2, syn_spec, syn_err, wrange2)
        scorr = mag2 - mag1
        scorr_err = np.sqrt(err1**2 + err2**2)
        return scorr, scorr_err
    else:
        #if error spectrum is not given
        mag1 = filter_mag(filt1, zerof1, syn_spec, wrange1)
        mag2 = filter_mag(filt2, zerof2, syn_spec, wrange2)
        scorr = mag2 - mag1
        return scorr

#################################################################
# Spectral Line Fitting                                         #
#################################################################

#function: bin spectrum to resolution
def bin_spec(wave, spec, R=400):
    spec_filt = np.zeros(len(spec))
    #bin resolved intervals
    for i in range(len(spec)):
        dw = wave[i]/R
        mask = np.logical_and(wave >= wave[i]-dw/2., wave < wave[i]+dw/2)
        spec_filt[i] = np.mean(spec[mask])
    return spec_filt

#function: Skewed Gaussian line profile
def lpNorm(x, A, mu, sig, skew, b):
    from scipy.stats import skewnorm, beta

    return A*skewnorm.pdf((x-mu)/sig, skew) + b
    
#function: fit the center of a line feature
def fitLine(center, width, spec, spec_err=None, plot=False):
    from scipy.optimize import curve_fit, minimize
    
    spec_wave = spec.waveset.value
    b_est = max(spec(spec_wave, flux_unit='flam').value*1.e14)
    mask = np.logical_and(spec_wave < center+width/2.,
                          spec_wave > center-width/2.)
    spec_wave = spec_wave[mask]
    spec_flux = spec(spec_wave, flux_unit='flam').value*1.e14

    #fit line center
    
    est = [-1.5, center, width/2.0, 3, b_est]
    est = [min(spec_flux)-b_est, center, width/2.0, 3, b_est]
    print est
    if spec_err is None:
        popt, pcov = curve_fit(lpNorm, spec_wave, spec_flux, p0=est,
                               maxfev=100000)
    else:
        spec_flux_err = spec_err(spec_wave, flux_unit='flam').value
        popt, pcov = curve_fit(lpNorm, spec_wave, spec_flux, p0=est,
                               sigma=spec_flux_err, maxfev=100000)
    perr = np.sqrt(np.diag(pcov))
    min_wave = minimize(lambda x: 1e20*lpNorm(x, *popt), [popt[1]])['x'][0]
    print popt
    #plot best fit
    if plot:
        import matplotlib.pyplot as plt
        
        plt.plot(spec_wave, spec_flux, c='g')
        plt.plot(spec_wave, lpNorm(spec_wave, *popt), c='r')
        plt.title("Skewed Gaussian Line Fit")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux [flam]")
        plt.show()

    #return line center
    return min_wave, perr[1]

def lin(x, a, b):
    return a*x + b

def voigt(x, y):
    # The Voigt function is also the real part of
    # w(z) = exp(-z^2) erfc(iz), the complex probability function,
    # which is also known as the Faddeeva function. Scipy has
    # implemented this function under the name wofz()
    from scipy.special import wofz
    z = x + 1j * y
    I = wofz(z).real
    return I


def Voigt(nu, alphaD, alphaL, nu_0, A):
    # The Voigt line shape in terms of its physical parameters
    # alphaD, alphaL half widths at half max for Doppler and Lorentz(not FWHM)
    # A - scaling factor
    f = np.sqrt(np.log(2))
    x = (nu - nu_0) / alphaD * f
    y = alphaL / alphaD * f
    V = A * f / (alphaD * np.sqrt(np.pi)) * voigt(x, y)
    return V

def DVoigt(nu, aD1, aL1, nu1, A1, aD2, aL2, nu2, A2, a, b):
    # The Voigt line shape in terms of its physical parameters
    # aD, aL are half widths at half max for Doppler and Lorentz(not FWHM)
    # A - scaling factor
    #first Voigt
    V1 = Voigt(nu, aD1, aL1, nu1, A1)
    #second Voigt
    V2 = Voigt(nu, aD2, aL2, nu2, A2)
    #linear background
    backg = b + a * nu
    V = V1+V2+backg
    return V

#function: measure equivalent width of Na I D lines
def fitNaID(spec, spec_err=None, r=15, z=0, plot=False, params=True):
    #r is annulus size around lines for fitting
    from scipy.optimize import curve_fit
    
    #load spectrum data
    spec_wave = spec.waveset.value
    spec_flux = spec(spec_wave, flux_unit='flam').value*1e14

    #Na I D locations
    line2 = 5890.0*(1.+z)
    line1 = 5896.0*(1.+z)

    #crop a fitting window around the Na I D
    mask = np.logical_and(spec_wave > line2-r, spec_wave < line1+r)
    #estimate background level
    b_est = 0.5*(spec_flux[mask][0] + spec_flux[mask][-1])

    #fitting doublet profile
    est = [0.4, 0.8, line2, -1.0,
           0.4, 0.8, line1, -0.5,
           0, b_est]
    print est
    if spec_err is None:
        popt, pcov = curve_fit(DVoigt, spec_wave[mask], spec_flux[mask],
                               p0=est, maxfev=1000000)
    else:
        spec_flux_err = spec_err(spec_wave, flux_unit='flam').value
        popt, pcov = curve_fit(DVoigt, spec_wave[mask], spec_flux[mask],
                               p0=est, sigma=spec_flux_err[mask],
                               maxfev=1000000)
    perr = np.sqrt(np.diag(pcov))
    print popt, perr

    #calculate equivalent widths
    wave2 = np.linspace(line2-r, line2+r, 1000)
    V2 = Voigt(wave2, popt[0], popt[1], popt[2], popt[3])
    backg2 = popt[9] + popt[8] * wave2
    EW2 = np.trapz(-V2/backg2, wave2)
    
    wave1= np.linspace(line1-r, line1+r, 1000)
    V1 = Voigt(wave1, popt[4], popt[5], popt[6], popt[7])
    backg1 = popt[9] + popt[8] * wave1
    EW1 = np.trapz(-V1/backg1, wave1)
    print "D2: EW =",EW2
    print "D1: EW =",EW1

    #calculate galactic extinction from Poznanski et al. 2012
    EBV2 = np.power(10, 2.16*EW2-1.91)
    err2 = EBV2*np.log(10)*0.15
    EBV1 = np.power(10, 2.47*EW1-1.76)
    err1 = EBV1*np.log(10)*0.17
    #SFD 1998 extinction
    print "D2: E(B-V) =",np.power(10, 2.16*EW2-1.91),err2
    print "D1: E(B-V) =",np.power(10, 2.47*EW1-1.76),err1

    #plotting the fit
    if plot:
        import matplotlib.pyplot as plt
        
        plt.plot(spec_wave, spec_flux, c='g')
        #plt.plot(spec_wave[mask], DVoigt(spec_wave[mask], *popt), c='r')
        plt.plot(wave2, V2+backg2, c='r')
        plt.plot(wave1, V1+backg1, c='b')
        plt.title("Doublet Voigt Line Fit")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux [flam]")
        plt.show()

    if params:
        return popt, perr

#function: Na I D upper limit
def limNaID(popt, spec, spec_err=None, r=15, sn=3, z=0, plot=False):
    #r is annulus size around lines for fitting
    from scipy.optimize import curve_fit
        
    #load spectrum data
    spec_wave = spec.waveset.value
    spec_flux = spec(spec_wave, flux_unit='flam').value*1e14

    #Na I D locations
    line2 = 5890.0*(1.+z)
    line1 = 5896.0*(1.+z)

    #crop a fitting window around the Na I D
    mask = np.logical_and(spec_wave > line2-r, spec_wave < line1+r)
    #estimate background level
    b_est = 0.5*(spec_flux[mask][0] + spec_flux[mask][-1])

    #fit a line to determine background level in spectrum
    plin, pcov = curve_fit(lin, spec_wave[mask], spec_flux[mask], p0=[0, b_est], maxfev=1000000)
    perr = np.sqrt(np.diag(pcov))
    lin_flux = lin(spec_wave[mask], *plin)
    #determine noise level per pixel in spectrum
    noise = np.std(spec_flux[mask] - lin_flux)

    #generate Voigt profiles until S/N = 3
    n = 500
    A = np.linspace(-0.1, -1, n)
    
    SNR1 = np.zeros(n)
    wave1= np.linspace(line1-r, line1+r, 1000)
    backg1 = lin(wave1, *plin)
    for i in range(n):
        V1 = Voigt(wave1, popt[0], popt[1], line1, A[i])
        EW1 = np.trapz(-V1/backg1, wave1)
        #signal to noise ratio
        EWmask1 = np.logical_and(wave1>=line1-0.5*EW1, wave1<=line1+0.5*EW1)
        SNR1[i] = np.mean(V1[EWmask1]/noise)
    i1 = np.argmin(np.absolute(SNR1 + sn))
    A1 = A[i1]
    SNR1 = SNR1[i1]
    V1 = Voigt(wave1, popt[0], popt[1], line1, A1)
    EW1 = np.trapz(-V1/backg1, wave1)
        
    SNR2 = np.zeros(n)
    wave2= np.linspace(line2-r, line2+r, 1000)
    backg2 = lin(wave2, *plin)
    for i in range(n):
        V2 = Voigt(wave2, popt[4], popt[5], line2, A[i])
        EW2 = np.trapz(-V2/backg2, wave2)
        #signal to noise ratio
        EWmask2 = np.logical_and(wave2>=line2-0.5*EW2, wave2<=line2+0.5*EW2)
        SNR2[i] = np.mean(V2[EWmask2]/noise)
    i2 = np.argmin(np.absolute(SNR2 + sn))
    A2 = A[i2]
    SNR2 = SNR2[i2]
    V2 = Voigt(wave1, popt[4], popt[5], line2, A2)
    EW2 = np.trapz(-V2/backg2, wave2)

    #convert EW to galactic extinction from Poznanski et al. 2012
    print "SNR2=", SNR2
    print "SNR1=", SNR1
    print "D2: EW =",EW2
    print "D1: EW =",EW1
    EBV2 = np.power(10, 2.16*EW2-1.91)
    err2 = EBV2*np.log(10)*0.15
    EBV1 = np.power(10, 2.47*EW1-1.76)
    err1 = EBV1*np.log(10)*0.17
    #SFD 1998 extinction
    print "D2: E(B-V) =",np.power(10, 2.16*EW2-1.91),err2
    print "D1: E(B-V) =",np.power(10, 2.47*EW1-1.76),err1

    #plotting the fit
    if plot:
        import matplotlib.pyplot as plt
        
        plt.plot(spec_wave, spec_flux, c='g')
        plt.plot(spec_wave[mask], lin_flux, c='k')
        plt.plot(wave2, V2+backg2, c='r')
        plt.plot(wave1, V1+backg1, c='b')
        plt.title("Doublet Voigt Line Fit")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux [flam]")
        plt.show()

#function: measure equivalent width of a line
def EWdirect(spec, lim1, lim2, r=20, spec_err=None, plot=False):
    #lim1 and lim2 are limits surrounding the line
    #r is annulus size for fitting
    
    #load spectrum data
    spec_wave = spec.waveset.value
    spec_flux = spec(spec_wave, flux_unit='flam').value
    
    #crop a window around the line
    lmask = np.logical_and(spec_wave > lim1, spec_wave < lim2)
    #crop an annulus around the line
    amask = np.logical_and(spec_wave > lim1-30, spec_wave < lim2+30)
    amask = np.logical_and(amask, np.logical_not(lmask))
    
    #fit continuum emission using 2nd order polynomial
    cont_poly = np.polyfit(spec_wave[amask], spec_flux[amask], 3)

    #measure equivalent width
    cont_flux = np.polyval(cont_poly, spec_wave[lmask])
    ratio = (cont_flux - spec_flux[lmask])/cont_flux
    EW = np.trapz(ratio, spec_wave[lmask])
    print "EW:", EW

    #plotting the fit
    if plot:
        import matplotlib.pyplot as plt

        plt.plot(spec_wave[lmask], ratio)
        plt.show()
        
        plt.plot(spec_wave, spec_flux, c='g')
        plt.plot(spec_wave[amask],
                 np.polyval(cont_poly, spec_wave[amask]),
                 c='r')
        plt.title("Direct Integration")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux [flam]")
        plt.show()
