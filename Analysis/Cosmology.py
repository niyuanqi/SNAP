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

#Reiss 2016 Cosmological parameters from SN1a
H0 = 73.24 #km/s/Mpc Hubble constant (1.74)
H = H0*1e3/1e6 #m/s/pc Hubble Constant
c = 2.99792485e8 #m/s
Wm = 0.27 #matter density parameter
Wl = 0.73 #vacuum density parameter
Tcmb=2.725 #CMB current temperature

#calibration information
#band zero fluxes [Jy] from Bessell 1998 (outdated)
#flux_0 = [1790, 4063, 3636, 3064, 2416]
#band zero fluxes [Jy] in UBVri with mix from AB system
flux_0 = [1790, 4063, 3636, 3631, 3631]

#telescope information
#band isophotal wavelengths [Angstrom] from Bessell 2005
wave_0 = [3663, 4361, 5448, 6407, 7980]
#band widths [Angstrom] from Bessel 2005
widths = [650, 890, 840, 1580, 1540]
#band numbering convention
bands = {'U':0, 'B':1, 'V':2, 'R':3, 'I':4, 'i':4}

#function: cosmological infinitesimal comoving distance [pc]
def comov(z):
    return (c/H)/np.sqrt(Wm*np.power(1.0+z,3)+Wl)

#function: integrate comoving distance to some redshift [pc]
def intDc(z, H1=H0):
    #import scipy.integrate as integrate
    #return integrate.quad(comov,0,z)[0]
    from astropy.cosmology import FlatLambdaCDM

    cosmo = FlatLambdaCDM(H0=H1, Om0=Wm, Tcmb0=2.725)
    return cosmo.comoving_distance(z).value*1e6

#function: integrate angular diameter distance to some redshift
def intDa(z, H1=H0):
    return intDc(z, H1=H1)/(1.0+z)

#function: integrate luminosity distance to some redshift
def intDl(z, H1=H0):
    return intDc(z, H1=H1)*(1.0+z)

#function: integrate distance modulus to some redshift
def intDM(z, H1=H0):
    return 5.*np.log10(intDl(z, H1=H1)/10)

#function: invert DM to equivalent hubble flow redshift
def DM_toz(DM, DMerr, z0=0, ra=None, dec=None):
    from scipy.optimize import fsolve
    diff_func = lambda z: MCDM(z, 0, ra=ra, dec=dec)[0] - DM
    zDM = fsolve(diff_func, z0)[0]
    return zDM

#function: bootstrap DM error
def MCDM(z, zerr, n=1000, ra=None, dec=None, Herr=None, verbose=False):
    if ra is not None:
        if dec is None:
            print "COORDINATE ERROR"
        #Mould et al. (2000) local velocity field correction
        z, zerr = M00_corrz(z, zerr, ra, dec, verbose)
    if Herr is None:
        Hs = np.ones(n)*H0
    else:
        Hs = np.random.normal(H0, Herr, n)
    dzs = np.random.normal(0,zerr,n)
    DMs = []
    for i in range(n):
        #get distance modulus
        DM = intDM(z+dzs[i], H1=Hs[i])
        DMs.append(DM)
    DM = intDM(z)
    DMerr = np.std(DMs)
    return DM, DMerr

#function: Mould et al. (2000) infall velocity component
def r_infall(Vo, Va, theta):
    gma = 2.
    roa = np.sqrt(Vo**2 + Va**2 - 2*Vo*Va*np.cos(theta))
    return np.cos(theta)+np.power(roa/Va, 1-gma)*(Vo-Va*np.cos(theta))/roa
#function: Transform heliocentric to LG coordinates
def Hel_toLG(v, l, b):
    #(Karachentsev & Makarov 1996)
    vap = 316
    lap, bap = 93*np.pi/180., -4*np.pi/180.
    dvLG = vap*(np.sin(b)*np.sin(bap) + np.cos(b)*np.cos(bap)*np.cos(l-lap))
    return v + dvLG

#function: Mould et al. (2000) local velocity field correction
def M00_corrz(z, zerr, ra, dec, verbose=False):
    from SNAP.Astrometry import Hel_toGal, sepAngle
    #observed redshift z
    v, verr = c*1e-3*z, c*1e-3*zerr
    
    l, b = Hel_toGal(ra, dec)
    l, b = l*np.pi/180., b*np.pi/180.
    #correction for local group velocity 
    vLG = Hel_toLG(v, l, b)
    
    #Virgo infall
    Virgo_pos = (187.0791667, 12.0666667)
    #Virgo_l, Virgo_b = Hel_toGal(*Virgo_pos)
    Virgo_fid = 200 #km/s
    #Virgo_v = 1035 #km/s
    #Virgo_vLG = Hel_toLG(Virgo_v, Virgo_l, Virgo_b)
    Virgo_vLG = 957 #km/s
    theta = sepAngle((ra, dec), Virgo_pos)*np.pi/180.
    Virgo_vin = Virgo_fid*r_infall(vLG, Virgo_vLG, theta)

    #GA infall
    GA_pos = (200.0, -44.0)
    #GA_l, GA_b = Hel_toGal(*GA_pos)
    GA_fid = 400 #km/s
    #GA_v = 4600 #km/s
    #GA_vLG = Hel_toLG(GA_v, GA_l, GA_b)
    GA_vLG = 4380 #km/s
    theta = sepAngle((ra, dec), GA_pos)*np.pi/180.
    GA_vin = GA_fid*r_infall(vLG, GA_vLG, theta)

    #Shapley infall
    Shapley_pos = (202.5, -31.0)
    #Shapley_l, Shapley_b = Hel_toGal(*Shapley_pos)
    Shapley_fid = 85 #km/s
    #Shapley_v = 13800 #km/s
    #Shapley_vLG = Hel_toLG(Shapley_v, Shapley_l, Shapley_b)
    Shapley_vLG = 13600 #km/s
    theta = sepAngle((ra, dec), Shapley_pos)*np.pi/180.
    Shapley_vin = Shapley_fid*r_infall(vLG, Shapley_vLG, theta)

    #corrected to cosmic velocity
    vin = Virgo_vin + GA_vin + Shapley_vin
    vC = vLG + vin
    vCerr = np.sqrt(verr**2 + (0.06*(vLG-v))**2 + (0.07*Virgo_vin)**2
                    +(0.07*GA_vin)**2 + (0.07*Shapley_vin)**2)
    if verbose:
        print "vhelio", v
        print "vLG", vLG
        print "Virgo", vLG+Virgo_vin
        print "Virgo + GA", vLG+Virgo_vin+GA_vin
        print "Virgo + GA + Shapley", vC, vCerr
    
    #convert to cosmic redshift
    zC, zCerr = vC*1e3/c, vCerr*1e3/c
    return zC, zCerr

#function: convert luminosity distance (Mpc) to DM
def Dl_toDM(Dl, Dlerr):
    DM = 5.*np.log10(Dl*1e6) - 5.
    DMerr = 5.*Dlerr/(Dl*np.log(10.))
    return DM, DMerr

#function: convert DM to luminosity distance
def DM_toDl(DM, DMerr):
    Dl = np.power(10., 0.2*(DM+5))/1.e6
    Dlerr = Dl*np.log(10.)*DMerr
    return Dl, Dlerr #Mpc

#function: projected dist (pc) between angle separated objects at redshift z
def sepDistProj(sepRad, z):
    return intDa(z)*sepRad

#function: calculate dereddened magnitudes
def deredMag(appMag, EBV, Coef):
    return appMag - EBV*Coef
#function: calculate dereddened fluxes
def deredFlux(appFlux, EBV, Coef):
    return appFlux*np.power(10,EBV*Coef/2.5)
    
#function: calculate absolute magnitude at redshift z
def absMag(appMag, z, appMag_err=None, z_err=None, DM=None, DMerr=None, Kcorr=None):
    if DM is None:
        dl = intDl(z)
        DM = 5.*np.log10(dl/10.0)
    
    if Kcorr is None:
        #estimate absolute magnitude using luminosity distance and naive K
        Mabs = appMag - DM + 2.5*np.log10(1+z)
    else:
        #compute absolute magnitude using luminosity distance and K correction
        Mabs = appMag - DM - Kcorr
    if appMag_err is not None:
        #calculate corresponding errors
        if DM is None:
            dl_err = (intDl(z+z_err)-intDl(z-z_err))/2.0
            DMerr = (5./np.log(10))*(dl_err/dl)
        if Kcorr is None:
            Mabs_err = np.sqrt(np.square(appMag_err) + np.square(DMerr) + np.square((2.5/np.log(10))*(z_err/(1.0+z))))
        else:
            Mabs_err = np.sqrt(np.square(appMag_err) + np.square(DMerr))
        return Mabs, Mabs_err
    else:
        #don't calculate errors
        return Mabs

#function: calculate absolute magnitude at redshift z
def absFlux(appFlux, z, appFlux_err=None, z_err=None, DM=None, DMerr=None, Kcorr=None):
    if DM is None:
        dl = intDl(z)
    else:
        dl = np.power(10, 1.+DM/5.0)
    
    if Kcorr is None:
        #estimate absolute flux using luminosity distance and naive K
        Fabs = appFlux*(1+z)*(dl/10.0)**2
    else:
        #compute absolute magnitude using luminosity distance and K correction
        Fabs = appFlux*np.power(10,Kcorr/2.5)*(dl/10.0)**2
    if appFlux_err is not None:
        if DM is None:
            #calculate corresponding errors
            dl_err = (intDl(z+z_err)-intDl(z-z_err))/(2.0*np.sqrt(3))
        else:
            dl_err = (np.log(10)/5.0)*DMerr*dl
        if Kcorr is None:
            appFlux[appFlux == 0] = 1e-10
            Fabs_err = np.absolute(Fabs)*np.sqrt(np.square(appFlux_err/appFlux) + np.square(2*(dl_err/dl)) + np.square(z_err/(1.0+z)))
        else:
            Fabs_err = np.absolute(Fabs)*np.sqrt(np.square(appFlux_err/appFlux) + np.square(2*(dl_err/dl)))
        return Fabs, Fabs_err
    else:
        #don't calculate errors
        return Fabs

#function: calculate flux [Jy] from mag
def Mag_toFlux(band, mag, mag_err=None):
    if isinstance(band, str):
        flux = flux_0[bands[band]]*np.power(10,mag/-2.5)
    else:
        flux = band*np.power(10,mag/-2.5)
    if mag_err is not None:
        #calculate errors
        flux_err = flux*np.log(10)*mag_err/2.5
        return flux, flux_err
    else:
        #don't calculate errors
        return flux
#function: calculate mag from flux [Jy]
def Flux_toMag(band, flux, flux_err=None):
    if isinstance(band, str):
        mod = flux/flux_0[bands[band]]
    else:
        mod = flux/band
    if hasattr(mod, '__iter__'):
        mag =  -2.5 * np.array([np.log10(m) for m in mod])
    else:
        mag =  -2.5 * np.log10(mod)
    if flux_err is not None:
        #calculate errors
        mag_err = 2.5*flux_err/(np.log(10)*flux)
        return mag, mag_err
    else:
        #don't calculate errors
        return mag

#function: transform magnitude to AB system
def Mag_toAB(band, mag):
    flux = Mag_toFlux(band, mag) #Jy
    mag_AB = Flux_toMag(3631, flux)
    return mag_AB
#function: transform AB magnitude to other system
def AB_toMag(band, mag):
    flux = Mag_toFlux(3631, mag) #Jy
    mag_AB = Flux_toMag(band, flux)
    return mag_AB

#function: calculate absolute time at redshift z
def absTime(appTime, z):
    #apply a simple time dilation
    return appTime/(1+z)

#function: subtract a constant magnitude object from mags
def Mag_subMag(mag1, err1, mag2, err2):
    flux1 = np.power(10,mag1/-2.5)
    flux2 = np.power(10,mag2/-2.5)
    relflux_sub = flux1-flux2
    mag_sub = -2.5*np.log10(relflux_sub)
    err_sub = np.sqrt(np.square(err1*flux1/relflux_sub) + np.square(err2*flux2/relflux_sub))
    return mag_sub, err_sub
#function: add a constant magnitude object to limiting mags
def Mag_addMag(mag1, mag2):
    flux1 = np.power(10,mag1/-2.5)
    flux2 = np.power(10,mag2/-2.5)
    relflux_sub = flux1+flux2
    mag_add = -2.5*np.log10(relflux_sub)
    return mag_add
#function: add a source of error to limiting mag
def Lim_addMag(mag1, mag2, err2, SNR=3):
    flux1 = np.power(10,mag1/-2.5)
    flux2 = np.power(10,mag2/-2.5)
    flux2err = flux2*np.log(10)*err2/2.5
    relflux_sub = np.sqrt(np.square(flux1/SNR)+np.square(flux2err))*SNR
    mag_add = -2.5*np.log10(relflux_sub)
    return mag_add
#function: subtract a constant magnitude object from fluxes
def Flux_subMag(flux1, err1, mag2, err2, band='i'):
    flux2 = flux_0[bands[band]]*np.power(10,mag2/-2.5)
    flux_sub = flux1 - flux2*1e6
    err_sub = np.sqrt(np.square(err1) + np.square(err2*flux2*np.log10(10)/(-2.5)))
    return flux_sub, err_sub

#function: subtract a constant flux object from mags
def Mag_subFlux(mag1, flux2, band='i'):
    flux1 = Mag_toFlux(band, mag1)
    relflux_sub = flux1-flux2*1e-6
    mag_sub = Flux_toMag(band, relflux_sub)
    return mag_sub
#function: subtract a constant magnitude object from fluxes
def Flux_subFlux(flux1, flux2):
    return flux1-flux2
