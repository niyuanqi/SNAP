#################################################################
# Name:     Cosmology.py                                        #
# Author:   Yuan Qi Ni                                          #
# Version:  July, 18, 2016                                      #
# Function: Program contains various routines for calculating   #
#           cosmological quantities. Also contains constants.   #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

#global constants
H0 = 69.6 #km/s/Mpc Hubble constant
H = H0*1e3/1e6 #m/s/pc Hubble Constant
c = 2.9979e8 #m/s
Wm = 0.286 #matter density parameter
Wl = 0.714 #vacuum density parameter

#function: cosmological infinitesimal comoving distance [pc]
def comov(z):
    return (c/H)/np.sqrt(Wm*np.power(1+z,3)+Wl)

#function: integrate comoving distance to some redshift [pc]
def intDc(z):
    return integrate.quad(comov,0,z)[0]*(1+z)

#function: integrate angular diameter distance to some redshift
def intDa(z):
    return intDc(z)/(1+z)

#function: integrate luminosity distance to some redshift
def intDl(z):
    return intDc(z)*(1+z)

#function: projected dist (pc) between angle separated objects at redshift z
def sepDistProj(sepRad, z):
    return intDa(z)*sepRad

#function: calculate absolute magnitude at redshift z
def absMag(appMag, z, appMag_err=None, z_err=None, Kcorr=None):
    if appMag_err == None:
        dl = intDl(z)
        if Kcorr == None:
            #estimate absolute magnitude using luminosity distance and naive K
            Mabs = appMag - 5.024*np.log10(dl/10.0) + 2.512*np.log10(1+z)
        else:
            #compute absolute magnitude using luminosity distance and K correction
            Mabs = appMag - 5.024*np.log10(dl/10.0) + Kcorr
        return Mabs
    else:
        dl = intDl(z)
        dl_err = (intDl(z+z_err)-intDl(z-z_err))/2.0
        if Kcorr == None:
            #estimate absolute magnitude using luminosity distance and naive K
            Mabs = appMag - 5.024*np.log10(dl/10.0) + 2.512*np.log10(1+z)
            Mabs_err = np.sqrt(np.square(appMag_err) + np.square((5.024/(10*np.log(10)))*(dl_err/dl)) + np.square((2.512/np.log(10))*(z_err/(1+z))))
        else:
            #compute absolute magnitude using luminosity distance and K correction
            Mabs = appMag - 5.024*np.log10(dl/10.0) + Kcorr
            Mabs_err = np.sqrt(np.square(appMag_err) + np.square((5.024/(10*np.log(10)))*(dl_err/dl)))
        return Mabs, Mabs_err

#function: calculate absolute time at redshift z
def absTime(appTime, z):
    return appTime/(1+z)



