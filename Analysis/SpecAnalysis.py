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
        else:
            isot_obs = spec_meta['DATE-OBS']+'T'+spec_meta['UT']
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
    #Synthetic observation
    obs = Observation(syn_spec, filt, force='extrap')
    flux = obs.effstim(flux_unit='jy', waverange=wrange).value
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

#function: Skewed Gaussian line profile
def lpNorm(x, A, mu, sig, skew, b):
    from scipy.stats import skewnorm, beta

    return A*skewnorm.pdf((x-mu)/sig, skew) + b
    
#function: fit the center of a line feature
def fitLine(center, width, spec, spec_err=None, plot=False):
    from scipy.optimize import curve_fit, minimize
    
    spec_wave = spec.waveset.value
    b_est = max(spec(spec_wave, flux_unit='flam').value)
    mask = np.logical_and(spec_wave < center+width/2.,
                          spec_wave > center-width/2.)
    spec_wave = spec_wave[mask]
    spec_flux = spec(spec_wave, flux_unit='flam').value

    #fit line center
    
    est = [-1.0, center, width/2.0, 0.0, b_est]
    if spec_err is None:
        popt, pcov = curve_fit(lpNorm, spec_wave, spec_flux, p0=est,
                               maxfev=100000)
    else:
        spec_flux_err = spec_err(spec_wave, flux_unit='flam').value
        popt, pcov = curve_fit(lpNorm, spec_wave, spec_flux, p0=est,
                               sigma=spec_flux_err, maxfev=100000)
    perr = np.sqrt(np.diag(pcov))
    min_wave = minimize(lambda x: 1e20*lpNorm(x, *popt), [popt[1]])['x'][0]
    
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


