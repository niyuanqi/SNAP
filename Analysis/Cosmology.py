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

#function: calculate absolute magnitude at redshift z
def absMag(appMag, z, appMag_err=None, z_err=None, Kcorr=None):
    if appMag_err is None:
        dl = intDl(z)
        if Kcorr is None:
            #estimate absolute magnitude using luminosity distance and naive K
            Mabs = appMag - 5.024*np.log10(dl/10.0) - 2.512*np.log10(1+z)
        else:
            #compute absolute magnitude using luminosity distance and K correction
            Mabs = appMag - 5.024*np.log10(dl/10.0) - Kcorr
        return Mabs
    else:
        dl = intDl(z)
        dl_err = (intDl(z+z_err)-intDl(z-z_err))/2.0
        if Kcorr is None:
            #estimate absolute magnitude using luminosity distance and naive K
            Mabs = appMag - 5.024*np.log10(dl/10.0) - 2.512*np.log10(1.0+z)
            Mabs_err = np.sqrt(np.square(appMag_err) + np.square((5.024/(10*np.log(10)))*(dl_err/dl)) + np.square((2.512/np.log(10))*(z_err/(1.0+z))))
        else:
            #compute absolute magnitude using luminosity distance and K correction
            Mabs = appMag - 5.024*np.log10(dl/10.0) - Kcorr
            Mabs_err = np.sqrt(np.square(appMag_err) + np.square((5.024/(10*np.log(10.0)))*(dl_err/dl)))
        return Mabs, Mabs_err

#function: calculate absolute time at redshift z
def absTime(appTime, z):
    return appTime/(1+z)



