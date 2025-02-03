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
        spec = spec_fits[0]
        if get_err:
            spec_err = spec_fits[3]
        #print spec_meta['SITEID']
    else:
        spec = spec_fits
    #check units
    if 'BUNIT' not in spec_meta or spec_meta['BUNIT']=='erg/cm2/s/A' or spec_meta['BUNIT']=='erg/cm2/s/Angstrom':
        unit = u.erg/u.cm**2/u.s/u.AA
    elif spec_meta['BUNIT']=='adu' or spec_meta['BUNIT']=='DU/PIXEL':
        unit = 1
    #spectrum flux data
    flux = spec.data*unit
    if get_err:
        flux_err = spec_err.data*unit
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
    if get_meta == 'mjd':
        from astropy.time import Time
        datestr = 'MJD-OBS'
        isot_obs = Time(spec_meta[datestr], format='mjd', scale='utc').isot
        return spec, isot_obs
    elif get_meta:
        from astropy.coordinates import SkyCoord
        #calculate time observed
        if 'DATE-OBS' in spec_meta:
            datestr = 'DATE-OBS'
            if 'T' in spec_meta[datestr]:
                isot_obs = spec_meta[datestr]
            elif 'UT' in spec_meta:
                isot_obs = spec_meta[datestr]+'T'+spec_meta['UT']
            elif 'UT-TIME' in spec_meta:
                isot_obs = spec_meta[datestr]+'T'+spec_meta['UT-TIME']
            elif 'TIME-OBS' in spec_meta:
                isot_obs = spec_meta[datestr]+'T'+spec_meta['TIME-OBS']
        #exposure time
        exp_obs = spec_meta['EXPTIME']
        #RA, DEC observed
        RA, DEC = spec_meta['RA'], spec_meta['DEC']
        coord = SkyCoord(RA, DEC, unit=(u.hourangle, u.deg))
        RA, DEC = coord.ra.degree, coord.dec.degree
        #return spectrum and metadata
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

#function: loading spectra
def spec_dat(filename):
    import astropy.units as u
    from synphot import SourceSpectrum
    from synphot.models import Empirical1D
    
    wave, flux = np.loadtxt(filename, unpack=True)
    #create synphot spectrum object
    spec = SourceSpectrum(Empirical1D, points=wave*u.AA,
                          lookup_table=flux*u.erg/u.cm**2/u.s/u.AA,
                          keep_neg=True)
    return spec

#function: loading spectra
def spec_txt(filename):
    import astropy.units as u
    from synphot import SourceSpectrum
    from synphot.models import Empirical1D
    
    wave, flux = np.loadtxt(filename, unpack=True, skiprows=1, delimiter=', ')
    #create synphot spectrum object
    spec = SourceSpectrum(Empirical1D, points=wave*u.AA,
                          lookup_table=flux*u.erg/u.cm**2/u.s/u.AA,
                          keep_neg=True)
    return spec

#################################################################
# Synthetic Observations                                        #
#################################################################

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
        fwave = filt.waveset
        wlengths = fwave[np.logical_and(fwave>wrange[0]*u.AA, fwave<wrange[1]*u.AA)]
    #Synthetic observation
    obs = Observation(syn_spec, filt, force='extrap')
    flux = obs.effstim(flux_unit='jy', waverange=wrange).value
    #flux = obs.effstim(flux_unit='jy', wavelengths=wlengths).value
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
def filter_mag(filterfile, filterzero, syn_spec, syn_err=None, z=0, wrange=None):
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
    waves = np.linspace(wrange[0], wrange[-1], 10000)
    filt = SpectralElement(Empirical1D, points=fwave*(1.+z),
                           lookup_table=filt(fwave))
        
    #Synthetic observation
    obs = Observation(syn_spec, filt, force='extrap')
    flux = obs.effstim(flux_unit='jy', waverange=wrange).value
    #flux = obs.effstim(flux_unit='jy', wavelengths=waves).value
    #Calibrate magnitude with zero point
    mag = -2.5*np.log10(flux/filterzero)
    #Synthetic observation of error spectrum
    if syn_err is not None:
        #square filter and error spectrum
        filt2 = SpectralElement(Empirical1D, points=fwave*(1.+z),
                                lookup_table=np.square(filt(fwave)))
        pseudo_flux = np.square(syn_err(swave,flux_unit='jy')).value
        syn_err2 = SourceSpectrum(Empirical1D, points=swave,
                                  lookup_table=pseudo_flux)
        #sum errors in quadrature
        obs = Observation(syn_err2, filt2, force='extrap')
        #flux_err = np.sqrt(obs.effstim(waverange=wrange).value)
        flux_err = np.sqrt(obs.effstim(wavelengths=waves).value)
        #convert to magnitude error
        mag_err = (2.5/np.log(10))*(flux_err/flux)
        return mag, mag_err
    else:
        return mag

#################################################################
# Spectral Filter Corrections                                   #
#################################################################

#function: evaluate S correction on Vega spectrum from filt1 to filt2
def Scorr_vega(filt1,filt2, zerof1,zerof2, wrange1=None,wrange2=None):
    from synphot import SourceSpectrum
    
    #load Vega spectrum
    spec = SourceSpectrum.from_vega()
    #Scorr = -2.5log(flux2/flux1) such that mag2 = mag1 + Scorr
    mag1 = filter_mag(filt1, zerof1, spec, wrange=wrange1)
    mag2 = filter_mag(filt2, zerof2, spec, wrange=wrange2)
    scorr = mag2 - mag1
    return scorr

#function: evaluate S correction on spec from filt1 to filt2
def Scorr(filt1,filt2, zerof1,zerof2, syn_spec, syn_err=None, wrange1=None,wrange2=None):
    #Scorr = -2.5log(flux2/flux1) such that mag2 = mag1 + Scorr
    if syn_err is not None:
        #if error spectrum is given
        mag1, err1 = filter_mag(filt1, zerof1, syn_spec, syn_err, wrange=wrange1)
        mag2, err2 = filter_mag(filt2, zerof2, syn_spec, syn_err, wrange=wrange2)
        scorr = mag2 - mag1
        scorr_err = np.sqrt(err1**2 + err2**2)
        return scorr, scorr_err
    else:
        #if error spectrum is not given
        mag1 = filter_mag(filt1, zerof1, syn_spec, wrange=wrange1)
        mag2 = filter_mag(filt2, zerof2, syn_spec, wrange=wrange2)
        scorr = mag2 - mag1
        return scorr

#function: evaluate K correction on spec from 0 to z
def Kcorr(z, filt, zerof, syn_spec, syn_err=None, wrange=None):
    #Kcorr = 2.5log(flux(z)/flux1(z=0)) such that mag(z) = mag(z=0) - Kcorr
    if syn_err is not None:
        #if error spectrum is given
        mag1, err = filter_mag(filt, zerof, syn_spec, syn_err,
                               z=z, wrange=wrange)
        mag2, err = filter_mag(filt, zerof, syn_spec, syn_err,
                               z=0, wrange=wrange)
        kcorr = -2.5*np.log10(1.+z) + mag2 - mag1
        #kcorr = mag2 - mag1
        kcorr_err = err
        return kcorr, kcorr_err
    else:
        #if error spectrum is not given
        mag1 = filter_mag(filt, zerof, syn_spec, z=z, wrange=wrange)
        mag2 = filter_mag(filt, zerof, syn_spec, z=0, wrange=wrange)
        kcorr = -2.5*np.log10(1.+z) + mag2 - mag1
        #kcorr = mag2 - mag1
        return kcorr

#################################################################
# Spectral Line Fitting                                         #
#################################################################

#function: bin spectrum to resolution
def bin_spec(wave, spec, filt='boxcar', R=400, filtwin=180.0):
    if filt == 'boxcar':
        spec_filt = np.zeros(len(spec))
        #bin resolved intervals
        for i in range(len(spec)):
            dw = wave[i]/R
            mask = np.logical_and(wave >= wave[i]-dw/2., wave < wave[i]+dw/2)
            spec_filt[i] = np.mean(spec[mask])
    elif filt == 'savgol':
        from scipy.signal import savgol_filter

        #filtering window [A] R~wave/filtwin
        dwave = np.median(wave[1:] - wave[:-1])
        window = int(np.ceil(filtwin/dwave) // 2 * 2 + 1)
        #2nd order savitsky-golay filter
        spec_filt = savgol_filter(spec, window, 2)
    return spec_filt

#function: Doppler velocity from line shift
def Dopp_v(zdop, zdop_err=None):
    #speed of light
    c = 299.792458 #1e3 km/s
    #relativistic Doppler equation
    beta = ((zdop+1)**2 - 1)/((zdop+1)**2 + 1)
    if zdop_err is None:
        #return Doppler velocity
        return beta*c
    else:
        #calculate error
        beta_err = zdop_err*2*(zdop+1)*(1 - beta)/((zdop+1)**2 + 1)
        #return Doppler velocity
        return beta*c, beta_err*c
#function: Doppler z from velocity
def Dopp_z(vdop, vdop_err=None):
    #speed of light
    c = 299.792458 #1e3 km/s
    #relativistic Doppler equation
    z = np.sqrt((c+vdop)/(c-vdop)) - 1
    if vdop_err is None:
        #return redshift
        return z
    else:
        #calculate error
        zerr = vdop_err*c/(c-vdop)**2/z
        #return redshift
        return z, zerr

#function: Skewed Gaussian line profile
def lpNorm(x, A, mu, sig, skew, b):
    from scipy.stats import skewnorm, beta

    return A*skewnorm.pdf((x-mu)/sig, skew) + b
    
#function: fit the center of a line feature
def fitLine(center, width, spec, skew=None, spec_err=None, plot=False):
    from scipy.optimize import curve_fit, minimize
    
    spec_wave = spec.waveset.value
    b_est = max(spec(spec_wave, flux_unit='flam').value*1.e14)
    if skew is None:
        skew = 0.5
    mask = np.logical_and(spec_wave < center+skew*width,
                          spec_wave > center-(1.-skew)*width)
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

def poly_gaus(x, *params):
    from scipy.stats import norm

    #count number of gaussians [A, mu, sig]
    ngaus = int(len(params)/3)
    #construction number of gaussians [A, mu, sig]
    y = np.zeros(len(x))
    for i in range(ngaus):
        y += params[i*3]*norm(params[i*3+1], params[i*3+2]).pdf(x)
    return y

def poly_caii3(x, m, b, *params):
    from scipy.stats import norm
    #Ca II triplet separation
    x0 = 8498.02
    x1 = 8542.09
    x2 = 8662.14
    #EW of lines propto line strengths propto gi*fik
    #this is opticaly thin limit (Silverman et al. 2015)
    #values below are log(gi*fik) from NIST absorption oscillator strength
    log_gifik0, gi0 = -1.318, 2*1.5+1
    log_gifik1, gi1 = -0.36, 2*2.5+1
    log_gifik2, gi2 = -0.622, 2*1.5+1
    fik0 = np.power(10, log_gifik0)/gi0
    fik1 = np.power(10, log_gifik1)/gi1
    fik2 = np.power(10, log_gifik2)/gi2
    #note sigma is from temp/velocity, not strength

    #count number of Ca II triplets [A, mu, sig]
    ngaus = int(len(params)/3)
    #construction number of gaussians [A, mu, sig]
    y = np.zeros(len(x))
    for i in range(ngaus):
        #line continuum
        z0 = (params[i*3+1]-x0)/x0
        c0 = m*(params[i*3+1])+b
        c1 = m*(z0*x1+x1)+b
        c2 = m*(z0*x2+x2)+b
        #line absoption ratio
        y0 = params[i*3]*norm(0,params[i*3+2]).pdf(0)
        r0 = (1 + y0/c0) #I/Ic
        r1 = r0*np.exp(-fik1+fik0)
        r2 = r0*np.exp(-fik2+fik0)
        A1 = -(1-r1)*c1/norm(0,params[i*3+2]).pdf(0)
        A2 = -(1-r2)*c2/norm(0,params[i*3+2]).pdf(0)
        #Sum Ca II triplet
        y += params[i*3]*norm(params[i*3+1], params[i*3+2]).pdf(x)
        y += A1*norm(z0*x1+x1, params[i*3+2]).pdf(x)
        y += A2*norm(z0*x2+x2, params[i*3+2]).pdf(x)
    return y

#function: concave down quadratic
def quadpeak(x,a,b,c):
    return -1 * abs(a) * x**2 + b * x + c 

#function: fit HVF and PVF of Si II line (based on Silverman et al. 2015)
def fit_SiII(center, lim1, lim2, spec, spec_err=None,
             HVF=True, C=True, plot=False):
    from scipy.optimize import curve_fit, minimize
    
    #decompose spectrum
    spec_wave = spec.waveset.value
    spec_flux = spec(spec_wave, flux_unit='flam').value

    win=20
    #measure left limit
    limmask = np.logical_and(spec_wave > lim1-win, spec_wave < lim1+win)
    p, pcov = curve_fit(quadpeak, spec_wave[limmask], spec_flux[limmask])
    #p = np.polyfit(spec_wave[limmask], spec_flux[limmask], 2)
    lim1 = p[1]/(2*abs(p[0]))
    cont1 = (p[1]**2 + 4*abs(p[0])*p[2])/(4*abs(p[0]))
    #measure right limit
    limmask = np.logical_and(spec_wave > lim2-win, spec_wave < lim2+win)
    p, pcov = curve_fit(quadpeak, spec_wave[limmask], spec_flux[limmask])
    #p = np.polyfit(spec_wave[limmask], spec_flux[limmask], 2)
    lim2 = p[1]/(2*abs(p[0]))
    cont2 = (p[1]**2 + 4*abs(p[0])*p[2])/(4*abs(p[0]))
    #represent continuum
    slope = (cont2 - cont1)/(lim2-lim1)

    #cut out feature
    mask = np.logical_and(spec_wave > lim1, spec_wave < lim2)
    feat_wave = spec_wave[mask]
    #subtract out continuum
    cont_flux = slope*(feat_wave - lim1) + cont1
    feat_flux = spec_flux[mask] - cont_flux

    #PVF estimate
    est = [-0.20, center[0], center[0]*1e4/3e5]
    #HVF + C estimate
    if HVF:
        #add to fit list
        est += [-0.10, center[1], center[1]*1e4/3e5]
    if C:
        #add to fit list
        est += [-0.02, center[2], center[2]*1e4/3e5]

    if plot:
        import matplotlib.pyplot as plt
        
        #plot feature and continuum definition
        plt.plot(spec_wave, spec_flux, 'k')
        plt.plot(feat_wave, cont_flux, 'r')
        plt.title("Feature and continuum definition")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux [flam]")
        plt.tight_layout()
        plt.show()

    #fit poly gaussian
    popt, pcov = curve_fit(poly_gaus, feat_wave, feat_flux, p0=est)
    perr = np.sqrt(np.diag(pcov))
    #order gaussians by wavelength
    if HVF:
        ngaus = len(popt)/3
        centers = popt[1::3]
        #locate HVF
        hvf_i = np.argmin(centers)
        hvf_popt = popt[hvf_i*3:hvf_i*3+3]
        hvf_perr = perr[hvf_i*3:hvf_i*3+3]
        popt = np.delete(popt, np.s_[hvf_i*3:hvf_i*3+3])
        perr = np.delete(perr, np.s_[hvf_i*3:hvf_i*3+3])
        centers = np.delete(centers, hvf_i)
        #locate PVF
        pvf_i = np.argmin(centers)
        pvf_popt = popt[pvf_i*3:pvf_i*3+3]
        pvf_perr = perr[pvf_i*3:pvf_i*3+3]
        popt = np.delete(popt, np.s_[pvf_i*3:pvf_i*3+3])
        perr = np.delete(perr, np.s_[pvf_i*3:pvf_i*3+3])
        #concatenate everything together
        popt = np.concatenate([pvf_popt, hvf_popt, popt])
        perr = np.concatenate([pvf_perr, hvf_perr, perr])
    
    #plot best fit
    if plot:
        from scipy.stats import norm

        #plot HVF and PVF fitting
        plt.plot(feat_wave, feat_flux, 'k')
        #do each gaussian separately
        ngaus = int(len(popt)/3)
        colors = ['b', 'c', 'm']
        labels = ['Si II PVF', 'Si II HVF', 'C II']
        for i in range(ngaus):
            gaus_flux = popt[i*3]*norm(popt[i*3+1], popt[i*3+2]).pdf(feat_wave)
            plt.plot(feat_wave, gaus_flux, color=colors[i], label=labels[i])
        #do the gaussians together
        poly_flux = poly_gaus(feat_wave, *popt)
        plt.plot(feat_wave, poly_flux, color='y', label='Combined')
        #continuum level
        plt.plot(feat_wave, np.zeros(len(feat_wave)), 'r')
        plt.title("Fit of HVF and PVF")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux [flam]")
        plt.legend()
        plt.tight_layout()
        plt.show()

    #return line centers
    return popt, perr

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

#function: measure equivalent width of Na I D lines
def fitVoigt(line0, spec, spec_err=None, r=15, z=0, plot=False, params=True):
    #r is annulus size around lines for fitting
    from scipy.optimize import curve_fit
    
    #load spectrum data
    spec_wave = spec.waveset.value
    spec_flux = spec(spec_wave, flux_unit='flam').value

    #line locations
    line = line0*(1.+z)

    filt=False
    #estimate spectrum error
    if isinstance(spec_err, list):
        #load spectrum data
        filt_flux = bin_spec(spec_wave, spec_flux, filt='savgol', filtwin=200)
        nmask = np.logical_and(spec_wave > spec_err[0]*(1.+z),
                               spec_wave < spec_err[1]*(1.+z))
        spec_err = np.std(spec_flux[nmask] - filt_flux[nmask])
        filt=True

    if line0 == 6564.614:
        print "Halpha"
        r1 = r
        r2 = r
    else:
        r1 = r
        r2 = r

    #crop a fitting window around the Na I D
    mask = np.logical_and(spec_wave > line-r1, spec_wave < line+r2)
    #estimate background level
    b_est = 0.5*(spec_flux[mask][0] + spec_flux[mask][-1])

    #fitting doublet profile
    est = [z, 0.4, 0.8, 100.0, 0, b_est]
    print est
    LineVoigt = lambda nu, z, aD1, aL1, A1, a, b: Voigt(nu, aD1, aL1, line0*(1.+z), A1) + lin(nu, a, b)
    if spec_err is None:
        popt, pcov = curve_fit(LineVoigt, spec_wave[mask], spec_flux[mask],
                               p0=est, maxfev=1000000)
        popt, pcov = curve_fit(LineVoigt, spec_wave[mask], spec_flux[mask],
                               p0=popt, maxfev=1000000)
    elif isinstance(spec_err, float):
        spec_flux_err = np.ones(len(spec_wave[mask]))*spec_err
        popt, pcov = curve_fit(LineVoigt, spec_wave[mask], spec_flux[mask],
                               p0=est, sigma=spec_flux_err,
                               absolute_sigma=True, maxfev=1000000)
        popt, pcov = curve_fit(LineVoigt, spec_wave[mask], spec_flux[mask],
                               p0=popt, sigma=spec_flux_err,
                               absolute_sigma=True, maxfev=1000000)
    else: #error spectrum given
        spec_flux_err = spec_err(spec_wave, flux_unit='flam').value
        popt, pcov = curve_fit(LineVoigt, spec_wave[mask], spec_flux[mask],
                               p0=est, sigma=spec_flux_err[mask],
                               absolute_sigma=True, maxfev=1000000)
        popt, pcov = curve_fit(LineVoigt, spec_wave[mask], spec_flux[mask],
                               p0=popt, sigma=spec_flux_err[mask],
                               absolute_sigma=True, maxfev=1000000)
    perr = np.sqrt(np.diag(pcov))
    print popt
    print perr
    wave = np.linspace(line-r1, line+r2, 1000)
    V = LineVoigt(wave, *popt)
    
    z, zerr = popt[0], perr[0]
    print "Doppler shift:", z, zerr

    #plotting the fit
    if plot:
        import matplotlib.pyplot as plt
        
        plt.plot(spec_wave, spec_flux, c='g')
        #plt.plot(spec_wave[mask], DVoigt(spec_wave[mask], *popt), c='r')
        plt.plot(wave, V, c='b', linestyle=':')
        if filt:
            plt.plot(spec_wave[nmask], filt_flux[nmask],
                     c='orange', linestyle=':')
        plt.title("Voigt Line Fit")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux [flam]")
        plt.show()

    if params:
        return popt, perr

def DVoigt(nu, aD1, aL1, nu1, A1, aD2, aL2, nu2, A2, a, b):
    # The Voigt line shape in terms of its physical parameters
    # aD, aL are half widths at half max for Doppler and Lorentz(not FWHM)
    # A - scaling factor
    #first Voigt
    V1 = Voigt(nu, aD1, aL1, nu1, A1)
    #second Voigt
    V2 = Voigt(nu, aD2, aL2, nu2, A2)
    #linear background
    backg = lin(nu, a, b)
    V = V1+V2+backg
    return V

#function: measure equivalent width of Na I D lines
def fitNaID(spec, spec_err=None, r=15, z=0, plot=False, params=True):
    #r is annulus size around lines for fitting
    from scipy.optimize import curve_fit
    
    #load spectrum data
    spec_wave = spec.waveset.value
    spec_flux = spec(spec_wave, flux_unit='flam').value

    filt=False
    #estimate spectrum error
    if isinstance(spec_err, list):
        #load spectrum data
        filt_flux = bin_spec(spec_wave, spec_flux, filt='savgol', filtwin=200)
        nmask = np.logical_and(spec_wave > spec_err[0]*(1.+z),
                               spec_wave < spec_err[1]*(1.+z))
        spec_err = np.std(spec_flux[nmask] - filt_flux[nmask])
        filt=True

    #Na I D locations
    line20 = 5889.950
    line2 = line20*(1.+z)
    line10 = 5895.924
    line1 = line10*(1.+z)

    #crop a fitting window around the Na I D
    mask = np.logical_and(spec_wave > line2-r, spec_wave < line1+r)
    #mask = np.logical_and(spec_wave > 6864, spec_wave <6882)
    #estimate background level
    b_est = 0.5*(spec_flux[mask][0] + spec_flux[mask][-1])

    #fitting doublet profile
    est = [z, 1, 0.1, -0.2*b_est,
           1, 0.1, -0.2*b_est,
           0.0, b_est]
    print est
    NaVoigt = lambda nu, z, aD1, aL1, A1, aD2, aL2, A2, a, b: DVoigt(nu, aD1, aL1, line20*(1.+z), A1, aD2, aL2, line10*(1.+z), A2, a, b)
    if spec_err is None:
        popt, pcov = curve_fit(NaVoigt, spec_wave[mask], spec_flux[mask],
                               p0=est, maxfev=1000000)
        popt, pcov = curve_fit(NaVoigt, spec_wave[mask], spec_flux[mask],
                               p0=popt, maxfev=1000000)
    elif isinstance(spec_err, float):
        spec_flux_err = np.ones(len(spec_wave[mask]))*spec_err
        popt, pcov = curve_fit(NaVoigt, spec_wave[mask], spec_flux[mask],
                               p0=est, sigma=spec_flux_err,
                               absolute_sigma=True, maxfev=1000000)
        popt, pcov = curve_fit(NaVoigt, spec_wave[mask], spec_flux[mask],
                               p0=popt, sigma=spec_flux_err,
                               absolute_sigma=True, maxfev=1000000)
    else:
        spec_flux_err = spec_err(spec_wave, flux_unit='flam').value
        popt, pcov = curve_fit(NaVoigt, spec_wave[mask], spec_flux[mask],
                               p0=est, sigma=spec_flux_err[mask],
                               absolute_sigma=True, maxfev=1000000)
        popt, pcov = curve_fit(NaVoigt, spec_wave[mask], spec_flux[mask],
                               p0=popt, sigma=spec_flux_err[mask],
                               absolute_sigma=True, maxfev=1000000)
    perr = np.sqrt(np.diag(pcov))
    print popt
    print perr
    wave = np.linspace(line2-r, line1+r, 1000)
    V = NaVoigt(wave, *popt)
    
    z, zerr = popt[0], perr[0]
    print "Doppler shift:", z, zerr

    #calculate equivalent widths
    wave2 = np.copy(wave)
    V2 = Voigt(wave2, popt[1], popt[2], line20*(1.+z), popt[3])
    backg2 = popt[8] + popt[7] * wave2
    EW2 = np.trapz(-V2/backg2, wave2)
    
    wave1 = np.copy(wave)
    V1 = Voigt(wave1, popt[4], popt[5], line10*(1.+z), popt[6])
    backg1 = popt[8] + popt[7] * wave1
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
        plt.plot(wave, V, c='k', linestyle=':')
        if filt:
            plt.plot(spec_wave[nmask], filt_flux[nmask],
                     c='orange', linestyle=':')
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
    spec_flux = spec(spec_wave, flux_unit='flam').value

    #Na I D locations
    line2 = 5889.950*(1.+z)
    line1 = 5895.924*(1.+z)

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
    n = 1000
    A = np.linspace(-0.02, -1., n)
    
    SNR1 = np.zeros(n)
    wave1= np.linspace(line1-r, line1+r, 1000)
    backg1 = lin(wave1, *plin)
    for i in range(n):
        V1 = Voigt(wave1, popt[1], popt[2], line1, A[i])
        EW1 = np.trapz(-V1/backg1, wave1)
        #signal to noise ratio
        EWmask1 = np.logical_and(wave1>=line1-0.5*EW1, wave1<=line1+0.5*EW1)
        SNR1[i] = np.mean(V1[EWmask1]/noise)
    i1 = np.argmin(np.absolute(SNR1 + sn))
    A1 = A[i1]
    SNR1 = SNR1[i1]
    V1 = Voigt(wave1, popt[1], popt[2], line1, A1)
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
    err2 = EBV2*np.log(10)*0.15 #scatter of D2 relation
    EBV1 = np.power(10, 2.47*EW1-1.76)
    err1 = EBV1*np.log(10)*0.17 #scatter of D1 relation
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
def EWdirect(spec,lim1,lim2,r=20,spec_err=None, fit=True,plot=False,ret=False):
    #lim1 and lim2 are limits surrounding the line
    #r is annulus size for fitting
    
    #load spectrum data
    spec_wave = spec.waveset.value
    spec_flux = spec(spec_wave, flux_unit='flam').value
    
    #crop a window around the line
    lmask = np.logical_and(spec_wave > lim1, spec_wave < lim2)

    if fit == True:
        #crop an annulus around the line
        allmask = np.logical_and(spec_wave > lim1-30, spec_wave < lim2+30)
        amask = np.logical_and(allmask, np.logical_not(lmask))
        #fit continuum emission using 2nd order polynomial
        cont_poly = np.polyfit(spec_wave[amask], spec_flux[amask], 3)
        cont_flux = np.polyval(cont_poly, spec_wave[lmask])
        #continuum representation
        rep_w = spec_wave[allmask]
        rep_f = np.polyval(cont_poly, spec_wave[allmask])
        #measure equivalent width
        ratio = (cont_flux - spec_flux[lmask])/cont_flux
    elif fit == False:
        #measure continuum emission
        amask1 = np.logical_and(spec_wave>lim1-30, spec_wave<lim1)
        cont_flux1 = np.mean(spec_flux[amask1])
        amask2 = np.logical_and(spec_wave>lim2, spec_wave<lim2+30)
        cont_flux2 = np.mean(spec_flux[amask2])
        slope = (cont_flux2 - cont_flux1)/(lim2-lim1)
        #continuum representation
        rep_w = spec_wave[lmask]
        rep_f = slope*(rep_w - lim1) + cont_flux1
        #measure equivalent width
        ratio = (rep_f - spec_flux[lmask])/rep_f

    #integrate equivalent width
    EW = np.trapz(ratio, spec_wave[lmask])
    print "EW:", EW

    #plotting the fit
    if plot:
        import matplotlib.pyplot as plt

        plt.plot(spec_wave[lmask], ratio)
        plt.show()
        
        plt.plot(spec_wave, spec_flux, c='g')
        plt.plot(rep_w, rep_f, c='r')
        plt.title("Direct Integration")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Flux [flam]")
        plt.show()

    if ret:
        return EW
