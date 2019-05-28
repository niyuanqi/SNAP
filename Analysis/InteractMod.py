#################################################################
# Name:     InteractMod.py                                      #
# Author:   Yuan Qi Ni                                          #
# Version:  May 28, 2018                                        #
# Function: Program contains ejecta interaction models.         #
#################################################################

#essential modules
import numpy as np

#################################################################
# Kasen Companion-Ejecta Interaction (CEI) model.               #
#################################################################

#function: Kasen model of shock interaction with companion
def Kasen2010(t_day,a13,m_c=1,e_51=1,kappa=1.0):
    """This calculates the luminosity, Liso, and Teff for the Kasen2010 analytic models.
    This incorporates the parameterization of viewing angle from Olling 2015
    
    :param t_day: time (days since explosion in rest frame)
    :param a13: semi-major axis of binary separation (10^13 cm)
    :param theta: viewing angle (degrees) minimum at 180.
    :param m_c: ejecta mass in units of M_chandra. default = 1
    :param e_51: explosion energy in units of 10^51 ergs. default=1
    :param kappa: opacity. default = 0.2 cm^2/g
    :return: luminosity (erg/s) (isotropic, angular), Teff (K)
    """

    #offset t_day to account for time it takes for interaction to begin
    L_u = 1.69 # constant related to ejecta density profile.
    vt = 6.0 * 10**8 * L_u * np.sqrt(e_51/m_c) # transition velocity
    v9 = vt / 10**9
    
    ti = (1.0e4 * a13 / v9) / 86400.0
    t_day = t_day - ti

    #check validity of kasen
    if t_day > 0 and e_51/m_c > 0:
        # Equations for Luminosity and Teff
        Lc_iso = 10**43 * a13 * m_c * v9**(7./4.) * kappa**(-3./4.) * t_day**(-1./2.) # (erg/s)
        Teff = 2.5 * 10**4 * a13**(1./4.) * kappa**(-35./36) * t_day**(-37./72.)
    else:
        Lc_iso = 0
        Teff = 1000
    return Lc_iso,Teff #erg/s

#function: Kasen fitting functyion
def KasenFit(t_day,a13,kappa,wave,m_c,e_51,z,DM,t0):
    #essential imports
    from SEDAnalysis import BBflux
    
    #shift time to rest frame
    t_rest = (t_day)/(1+z) - t0
    #calculate Kasen luminosity in rest frame
    Lk, Tk = Kasen2010(t_rest,a13,m_c,e_51,kappa)
    #shift luminosity to observer frame flux in band
    Fk = BBflux(Lk,Tk,wave,z,DM)
    #return predicted flux in band
    return Fk

#function: Kasen isotropic correction for viewing angle
def Kasen_isocorr(theta):
    #param theta: viewing angle (degrees) minimum at 180.
    return 0.982 * np.exp(-((theta % 180.0)/99.7)**2) + 0.018

#function: Error function for multi-band early light curve leastsq fitting
def kasenMultiErr(p, t, L, L_err, z, DM, m_c, e_51):
    from Cosmology import wave_0, bands
    #Kasen component p0=epoch, p1=a13, p2=theta
    B_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['B']],
                                m_c, e_51, z, DM, p[0])
                       for ti in t[0]])*Kasen_isocorr(p[2])
    V_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['V']],
                                m_c, e_51, z, DM, p[0])
                       for ti in t[1]])*Kasen_isocorr(p[2])
    I_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['i']],
                                m_c, e_51, z, DM, p[0])
                       for ti in t[2]])*Kasen_isocorr(p[2])
    #Power law component
    B_pred = np.array(B_pred) + earlyFit(t[0], p[0], p[3], p[6]) 
    V_pred = np.array(V_pred) + earlyFit(t[1], p[0], p[4], p[7]) 
    I_pred = np.array(I_pred) + earlyFit(t[2], p[0], p[5], p[8]) 
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#function: rule out Kasen model to sig at angle theta
def ruleout(F, Ferr, Fk, Fkerr, theta, sig, lims):
    #angle corrected Kasen luminosity
    Fk_theta = Fk*Kasen_isocorr(theta)
    Fk_theta_err = Fkerr*Kasen_isocorr(theta)
    #total error
    Err = np.sqrt(np.square(Ferr)+np.square(Fk_theta_err))
    #which is more constraining? datapoint or limit?
    level = F + sig*Err
    level[level < lims] = lims[level < lims]
    #check if any points rule out angle with conf
    if any(Fk_theta > level):
        return True
    else:
        return False

#function: rule out Kasen model to sig at angle theta (using both distributions)
def sym_ruleout(F, Ferr, Fk, Fkerr, JN, theta, sig):
    #noise in data number
    N = np.sqrt(np.square(Ferr/JN) - np.absolute(F)/JN)
    #angle corrected Kasen luminosity
    Fk_theta = Fk*Kasen_isocorr(theta)
    Fk_theta_err = Fkerr*Kasen_isocorr(theta)
    #total error
    FN_err = np.sqrt(np.square(Fk_theta_err/JN)+np.square(N)+Fk_theta)*JN
    #check if any points rule out angle with conf
    if any(Fk_theta - sig*FN_err > F + sig*Ferr):
        return True
    else:
        return False

#################################################################
# Piro CSM-Ejecta Interaction (CSM) model.                      #
#################################################################

def CSMpeak(Eej, Mej, Mext, Rext):
    #Eej is ejecta kinetic energy in 10^51 ergs 
    #Mej is ejecta mass in solar masses
    #Mext is mass of extended material in 0.01 solar masses
    #Rext is radius of extended material in 10^13 cm

    #opacity
    k_opt = 0.1 #0.1g/cm^2
    k034 = k_opt/0.34
    
    #peak scalings from Nakar and Piro 2014
    #with fit to numerical work on Woosley et al. 1994
    # and Bersten et al. 2012
    t_peak = 0.9 * (k034**0.5)*(Eej**-0.25)*(Mej**0.17)*(Mext**0.57) #days
    L_peak = 2.0e43 * (k034**-1)*Eej*Rext*(Mej**-0.7)*(Mext**-0.3) #erg/s
    Teff_peak = 3.0e4 * (k034**-0.25)*(t_peak**-0.5)*(Rext**0.25) #K
 
    return t_peak, L_peak, Teff_peak

def CSMmod(t_day, Eej, Mej, Mext, Rext):
    #t_day is time in days since explosion
    #Eej is ejecta kinetic energy in 10^51 ergs 
    #Mej is ejecta mass in solar masses
    #Mext is mass of extended material in 0.01 solar masses
    #Rext is radius of extended material in 10^13 cm

    #time in seconds
    ts = t_day*86400.0

    #sB constant
    sb_const = 5.6704e-5 #erg/cm^2/s/K^4
    #opacity
    k_opt = 0.1 #0.1g/cm^2
    k034 = k_opt/0.34

    #velocity imparted on extended material
    vext = 2.0e9 * (Eej**0.5)*(Mej**-0.35)*(Mext**-0.15) #cm/s
    #expansion timescale
    te = Rext*1.0e13/vext #s
    #Energy imparted on extended material
    Eext = 4.0e49 * Eej*(Mej**-0.7)*(Mext**0.7) #erg

    #peak scaling quantities
    t_peak, L_peak, Teff_peak = CSMpeak(Eej, Mej, Mext, Rext)
    t_peak = t_peak*86400.0 #s

    #luminosity evolution of CSM interaction Piro 2015
    Lcsm = (te*Eext/t_peak**2)*np.exp(-ts*(ts+2*te)/(2.0*t_peak**2)) #erg/s
    #Radius evolution of CSM
    Rcsm = Rext*1.0e13 + vext*t #cm
    #Temperature evolution of CSM given blackbody
    Tcsm = np.power(Lcsm/(4*sb_const*np.pi*Rcsm**2), 0.25) #K
    
    return Lcsm, Tcsm

def CSMFit(t_day, wave, z, DM, Eej, Mej, Mext, Rext, t0):
    from LCFitting import BBflux

    #shift time to rest frame
    t_rest = (t_day)/(1+z) - t0
    #calculate CSM luminosity in rest frame
    Lcsm, Tcsm = CSMmod(t_day, Eej, Mej, Mext, Rext)
    #shift luminosity to observer frame flux in band
    Fcsm = BBflux(Lcsm,Tcsm,wave,z,DM)
    #return predicted flux in band
    return Fcsm
