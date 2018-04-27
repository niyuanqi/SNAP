#################################################################
# Name:     Cosmology.py                                        #
# Author:   Yuan Qi Ni                                          #
# Version:  July, 18, 2016                                      #
# Function: Program contains various routines for calculating   #
#           cosmological quantities. Also contains constants.   #
#################################################################

#essential modules
import numpy as np

#global constants
#Planck 2015 results. XIII. Cosmological parameters from CMB
H0 = 73.24 #km/s/Mpc Hubble constant (1.74)
H = H0*1e3/1e6 #m/s/pc Hubble Constant
c = 2.9979e8 #m/s
Wm = 0.27 #matter density parameter
Wl = 0.73 #vacuum density parameter
#Reiss  Cosmological parameters from SN1a

#calibration information
#band zero fluxes [Jy] from Bessel 1998 (outdated)
#flux_0 = [1790, 4063, 3636, 3064, 2416]
#band zero fluxes [Jy] in UBVri with mix from AB system
flux_0 = [1790, 4063, 3636, 3631, 3631]

#telescope information
#band isophotal wavelengths [Angstrom] from Bessel 2005
wave_0 = [3663, 4361, 5448, 6407, 7980]
#band widths [Angstrom] from Bessel 2005
widths = [650, 890, 840, 1580, 1540]
#band numbering convention
bands = {'U':0, 'B':1, 'V':2, 'R':3, 'I':4, 'i':4}

#function: cosmological infinitesimal comoving distance [pc]
def comov(z):
    return (c/H)/np.sqrt(Wm*np.power(1.0+z,3)+Wl)

#function: integrate comoving distance to some redshift [pc]
def intDc(z):

    import scipy.integrate as integrate
    
    return integrate.quad(comov,0,z)[0]*(1.0+z)

#function: integrate angular diameter distance to some redshift
def intDa(z):
    return intDc(z)/(1.0+z)

#function: integrate luminosity distance to some redshift
def intDl(z):
    return intDc(z)*(1.0+z)

#function: projected dist (pc) between angle separated objects at redshift z
def sepDistProj(sepRad, z):
    return intDa(z)*sepRad

#function: calculate dereddened magnitudes
def deredMag(appMag, EBV, Coef):
    return appMag - EBV*Coef
#function: calculate dereddened fluxes
def deredFlux(appFlux, EBV, Coef):
    return appFlux*np.power(10,EBV*Coef/2.512)
    
#function: calculate absolute magnitude at redshift z
def absMag(appMag, z, appMag_err=None, z_err=None, Kcorr=None):
    dl = intDl(z)
    if Kcorr is None:
        #estimate absolute magnitude using luminosity distance and naive K
        Mabs = appMag - 5.024*np.log10(dl/10.0) - 2.512*np.log10(1+z)
    else:
        #compute absolute magnitude using luminosity distance and K correction
        Mabs = appMag - 5.024*np.log10(dl/10.0) - Kcorr
    if appMag_err is not None:
        #calculate corresponding errors
        dl_err = (intDl(z+z_err)-intDl(z-z_err))/2.0
        if Kcorr is None:
            Mabs_err = np.sqrt(np.square(appMag_err) + np.square((5.024/(10*np.log(10)))*(dl_err/dl)) + np.square((2.512/np.log(10))*(z_err/(1.0+z))))
        else:
            Mabs_err = np.sqrt(np.square(appMag_err) + np.square((5.024/(10*np.log(10.0)))*(dl_err/dl)))
        return Mabs, Mabs_err
    else:
        #don't calculate errors
        return Mabs

#function: calculate absolute magnitude at redshift z
def absFlux(appFlux, z, appFlux_err=None, z_err=None, Kcorr=None):
    dl = intDl(z)
    if Kcorr is None:
        #estimate absolute flux using luminosity distance and naive K
        Fabs = appFlux*(1+z)*(dl/10.0)**2
    else:
        #compute absolute magnitude using luminosity distance and K correction
        Fabs = appFlux*np.power(10,Kcorr/2.512)*(dl/10.0)**2
    if appFlux_err is not None:
        #calculate corresponding errors
        dl_err = (intDl(z+z_err)-intDl(z-z_err))/(2.0*np.sqrt(3))
        if Kcorr is None:
            Fabs_err = np.absolute(Fabs)*np.sqrt(np.square(appFlux_err/appFlux) + np.square(2*(dl_err/dl)) + np.square(z_err/(1.0+z)))
        else:
            Fabs_err = np.absolute(Fabs)*np.sqrt(np.square(appFlux_err/appFlux) + np.square(2*(dl_err/dl)))
        return Fabs, Fabs_err
    else:
        #don't calculate errors
        return Fabs

#function: calculate flux [Jy] from mag
def Mag_toFlux(band, mag, mag_err=None):
    flux = flux_0[bands[band]]*np.power(10,mag/-2.512)
    if mag_err is not None:
        #calculate errors
        flux_err = flux*np.log(10)*mag_err/-2.512
        return flux, flux_err
    else:
        #don't calculate errors
        return flux
#function: calculate flux [Jy] from mag
def Flux_toMag(band, flux, flux_err=None):
    mod = flux/flux_0[bands[band]]
    mag =  -2.512 * np.array([np.log10(m) for m in mod])
    if flux_err is not None:
        #calculate errors
        mag_err = -2.512*flux_err/(np.log(10)*flux)
        return mag, mag_err
    else:
        #don't calculate errors
        return mag

#function: calculate absolute time at redshift z
def absTime(appTime, z):
    #apply a simple time dilation
    return appTime/(1+z)
