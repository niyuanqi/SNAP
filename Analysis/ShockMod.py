#################################################################
# Name:     ShockMod.py                                         #
# Author:   Yuan Qi Ni                                          #
# Version:  May 28, 2018                                        #
# Function: Program contains ejecta shock interaction models.   #
#################################################################

#essential modules
import numpy as np

#################################################################
# Piro, Chang, Weinberg 2010 White Dwarf Shock Breakout model.  #
#################################################################

#function: Piro, Chang, Weinberg 2010 model of shock breakout cooling
def BreakoutFlash(R8,T5,m_c=1):
    """This calculates the luminosity, L, and Teff for the Piro, Chang, and Weinberg 2010 analytic prompt breakout flash model.
    
    :param R8: radius of envelope surface [10^8 cm].
    :param T5: effective temperature of envelope [10^5 K]
    :param m_c: ejecta mass in units of M_chandra. default = 1
    :return: energy (erg), Teff (K)
    """

    #Errata Bloom et al. 2012
    EL = 1./(7.**(4./3.))
    ET = 1./(7.**(1./3.))

    #constants
    G = 6.674e-8 #cm^3/g/s^2
    Msun = 1.988e33 #g
    mu_e = 2.
    #surface gravity
    M = m_c*1.4*Msun #g
    g = G*M/(R8*1.e8)**2 #cm/s^2
    
    #scalings
    t = t_day * 86400. #s
    t4 = t/(1.e4) #10^4 s
    R85 = R8 / (3.) #3x10^8 cm
    g9 = g/1.e9 #x10^9 cm/s^2
    
    #shock runaway characteristic scales
    V9 = 0.6 #x10^9 cm/s
    rho6 = 2.0*(g9**0.11) #x10^6 g/cm^3
    #critical density
    rho_c = 9.e3 #g/cm^3

    #polytropic constants
    K1 = 6.1e13 * (g9**(-1./3.))*(T5**(4./3.)) #shallow depth (relativistic)
    #scalings
    K1 = K1/(6.e13)

    #Energy and temperature model
    Efl = EL* (4.e40) * ((g9/K1)**-0.16) * v9 * (rho6**0.22) * (R85**2.) #erg
    Tfl = ET* (2.e8) * ((g9/K1)**0.14) * (v9**0.32) * (rho6**0.068) #K

    return Efl, Tfl #erg, K

#function: Piro, Chang, Weinberg 2010 model of shock breakout cooling
def ShockCoolingMod(t_day,R8,T5,m_c=1,late=True):
    """This calculates the luminosity, L, and Teff for the Piro, Chang, and Weinberg 2010 analytic shock cooling model.
    
    :param t_day: time (days since explosion in rest frame)
    :param R8: radius of envelope surface [10^8 cm].
    :param T5: effective temperature of envelope [10^5 K]
    :param m_c: ejecta mass in units of M_chandra. default = 1
    :return: luminosity (erg/s), Teff (K)
    """

    #Errata Bloom et al. 2012
    EL = 1./(7.**(4./3.))
    ET = 1./(7.**(1./3.))
    
    #constants
    G = 6.674e-8 #cm^3/g/s^2
    Msun = 1.988e33 #g
    mu_e = 2.
    #surface gravity
    M = m_c*1.4*Msun #g
    g = G*M/(R8*1.e8)**2 #cm/s^2
    
    #scalings
    t = t_day * 86400. #s
    t4 = t/(1.e4) #10^4 s
    R85 = R8*1.e8 / (10**8.5) #3x10^8 cm
    g9 = g/1.e9 #x10^9 cm/s^2
    
    #shock runaway characteristic scales
    v9 = 0.6 #x10^9 cm/s
    rho6 = 2.0*(g9**0.11) #x10^6 g/cm^3
    #rho6 = 2.0
    #critical density
    rho_c = 9.e3 #g/cm^3

    #polytropic constants
    K1 = 6.1e13 * (g9**(-1./3.))*(T5**(4./3.)) #shallow depth (relativistic)
    K2 = 9.91e12 / (mu_e**(5./3.)) #deep depth (non-relativistic)
    #scalings
    K1 = K1/(10**13.8)
    K2 = K2/(1.e13)
    
    #time to reach non-relativistic material
    t_c = (rho_c / (2.*((v9*g9/K1)**0.66)*(rho6**0.12)/(R85**1.3)))**(1./1.3)

    #late, t_c =False, 100000000
    #Luminosity and temperature model
    if t < 0:
        Lsh = 0
        Tsh = 100
    elif t < t_c and not late:
        #shallow diffusion depth (relativistic)
        Lsh = EL* (3.e41) * ((g9/K1)**-0.5) * (
            v9**1.8) * (rho6**0.42) * R85 * (t**-0.34) #erg/s 
        Tsh = ET* (1.e6) * ((g9/K1)**-0.065) * (
            v9**0.019) * (rho6**0.0035) * (R85**0.13) * (t**-0.46) #K
    else:
        #deep diffusion depth (non-relativistic)
        Lsh = EL* (2.e40) * ((g9/K2)**-0.41) * (
            v9**1.9) * (rho6**0.36) * (R85**0.83) * (t4**-0.16) #erg/s
        Tsh = ET* (2.e4) * ((g9/K2)**-0.058) * (
            v9**0.030) * (rho6**0.0058) * (R85**0.11) * (t4**-0.44) #K
        
    return Lsh, Tsh #erg/s, K

def ShockCoolingFit(t_day, wave, z, DM, m_c, R8, T5, t0, late=True):
    from SEDAnalysis import BBflux
    #shift time to rest frame
    t_rest = t_day/(1+z) - t0
    #calculate shock cooling luminosity in rest frame
    Lsh, Tsh = ShockCoolingMod(t_rest, R8, T5, m_c, late=late)
    #shift luminosity to observer frame flux in band
    Fsh = BBflux(Lsh,Tsh,wave,z,DM)
    #return predicted flux in band
    return Fsh

#function: Error function for multi-band early light curve leastsq fitting
def ShockCoolingMultiErr(p, t, L, L_err, z, DM, Mej):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #Shock cooling component p0=epoch (in rest frame), p1=R8
    #Fix t5=0.1 temperature, because assume late time
    B_pred = np.array([ShockCoolingFit(ti, wave_0[bands['B']], z, DM, Mej,
                                       p[1], 0.1, p[0], late=True)
                       for ti in t[0]])
    V_pred = np.array([ShockCoolingFit(ti, wave_0[bands['V']], z, DM, Mej,
                                       p[1], 0.1, p[0], late=True)
                       for ti in t[1]])
    I_pred = np.array([ShockCoolingFit(ti, wave_0[bands['i']], z, DM, Mej,
                                       p[1], 0.1, p[0], late=True)
                       for ti in t[2]])
    #Power law component, p2=epoch (in rest frame)
    B_pred = np.array(B_pred) + earlyFit(t[0], p[2]*(1.+z), p[3], p[6]) 
    V_pred = np.array(V_pred) + earlyFit(t[1], p[2]*(1.+z), p[4], p[7])
    I_pred = np.array(I_pred) + earlyFit(t[2], p[2]*(1.+z), p[5], p[8]) 
    #Error
    B_err = (B_pred - L[0])/L_err[0]/1000.
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#function: Error function for multi-band early light curve leastsq fitting
def ShockCoolingViErr(p, t, L, L_err, z, DM, Mej):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #Shock cooling component p0=epoch (in rest frame), p1=R8
    #Fix t5=0.1 temperature, because assume late time
    V_pred = np.array([ShockCoolingFit(ti, wave_0[bands['V']], z, DM, Mej,
                              p[1], 0.1, p[0]) for ti in t[0]])
    I_pred = np.array([ShockCoolingFit(ti, wave_0[bands['i']], z, DM, Mej,
                              p[1], 0.1, p[0]) for ti in t[1]])
    #Power law component, p2=epoch (in rest frame) 
    V_pred = np.array(V_pred) + earlyFit(t[0], p[2]*(1.+z), p[3], p[5])
    I_pred = np.array(I_pred) + earlyFit(t[1], p[2]*(1.+z), p[4], p[6]) 
    #Error
    V_err = (V_pred - L[0])/L_err[0]
    I_err = (I_pred - L[1])/L_err[1]
    return np.concatenate([V_err, I_err],axis=0)

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
def KasenFit(t_day,a13,kappa,wave,z,m_c,e_51,DM,t0):
    #essential imports
    from SEDAnalysis import BBflux
    
    #shift time to rest frame
    t_rest = t_day/(1+z) - t0
    
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
    from LCFitting import earlyFit
    #Kasen component p0=epoch (in rest frame), p1=a13, p2=theta
    B_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['B']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[0]])*Kasen_isocorr(p[2])
    V_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['V']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[1]])*Kasen_isocorr(p[2])
    I_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['i']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[2]])*Kasen_isocorr(p[2])
    #Power law component, p3=epoch
    B_pred = np.array(B_pred) + earlyFit(t[0], p[3]*(1.+z), p[4], p[7]) 
    V_pred = np.array(V_pred) + earlyFit(t[1], p[3]*(1.+z), p[5], p[8]) 
    I_pred = np.array(I_pred) + earlyFit(t[2], p[3]*(1.+z), p[6], p[9]) 
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#function: Error function for multi-band early light curve leastsq fitting
def kasenFixedMultiErr(p, t, L, L_err, z, DM, m_c, e_51, angle):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #Kasen component p0=epoch (in rest frame), p1=a13, p2=theta
    B_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['B']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[0]])*Kasen_isocorr(angle)
    V_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['V']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[1]])*Kasen_isocorr(angle)
    I_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['i']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[2]])*Kasen_isocorr(angle)
    #Power law component, p3=epoch
    B_pred = np.array(B_pred) + earlyFit(t[0], p[2]*(1.+z), p[3], p[6]) 
    V_pred = np.array(V_pred) + earlyFit(t[1], p[2]*(1.+z), p[4], p[7]) 
    I_pred = np.array(I_pred) + earlyFit(t[2], p[2]*(1.+z), p[5], p[8]) 
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#function: Error function for multi-band early light curve leastsq fitting
def kasent0MultiErr(p, t, L, L_err, z, DM, m_c, e_51, t0):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #Kasen component p0=epoch (in rest frame), p1=a13, p2=theta
    B_pred = np.array([KasenFit(ti, p[0], 1.0, wave_0[bands['B']], z,
                                m_c, e_51, DM, t0)
                       for ti in t[0]])*Kasen_isocorr(p[1])
    V_pred = np.array([KasenFit(ti, p[0], 1.0, wave_0[bands['V']], z,
                                m_c, e_51, DM, t0)
                       for ti in t[1]])*Kasen_isocorr(p[1])
    I_pred = np.array([KasenFit(ti, p[0], 1.0, wave_0[bands['i']], z,
                                m_c, e_51, DM, t0)
                       for ti in t[2]])*Kasen_isocorr(p[1])
    #Power law component, p3=epoch
    B_pred = np.array(B_pred) + earlyFit(t[0], p[2]*(1.+z), p[3], p[6]) 
    V_pred = np.array(V_pred) + earlyFit(t[1], p[2]*(1.+z), p[4], p[7]) 
    I_pred = np.array(I_pred) + earlyFit(t[2], p[2]*(1.+z), p[5], p[8]) 
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#function: Error function for multi-band early light curve leastsq fitting
def kasenPowMultiErr(p, t, L, L_err, z, DM, m_c, e_51, sep, angle):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #Kasen component p0=epoch (in rest frame), p1=a13, p2=theta
    B_pred = np.array([KasenFit(ti, sep, 1.0, wave_0[bands['B']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[0]])*Kasen_isocorr(angle)
    V_pred = np.array([KasenFit(ti, sep, 1.0, wave_0[bands['V']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[1]])*Kasen_isocorr(angle)
    I_pred = np.array([KasenFit(ti, sep, 1.0, wave_0[bands['i']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[2]])*Kasen_isocorr(angle)
    #Power law component, p3=epoch
    B_pred = np.array(B_pred) + earlyFit(t[0], p[1]*(1.+z), p[2], p[5]) 
    V_pred = np.array(V_pred) + earlyFit(t[1], p[1]*(1.+z), p[3], p[6]) 
    I_pred = np.array(I_pred) + earlyFit(t[2], p[1]*(1.+z), p[4], p[7]) 
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#function: Error function for multi-band early light curve leastsq fitting
def kasenViErr(p, t, L, L_err, z, DM, m_c, e_51):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #Kasen component p0=epoch (in rest frame), p1=a13, p2=theta
    V_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['V']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[0]])*Kasen_isocorr(p[2])
    I_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['i']], z,
                                m_c, e_51, DM, p[0])
                       for ti in t[1]])*Kasen_isocorr(p[2])
    #Power law component
    V_pred = np.array(V_pred) + earlyFit(t[0], p[3]*(1.+z), p[4], p[6]) 
    I_pred = np.array(I_pred) + earlyFit(t[1], p[3]*(1.+z), p[5], p[7]) 
    #Error
    V_err = (V_pred - L[0])/L_err[0]
    I_err = (I_pred - L[1])/L_err[1]
    return np.concatenate([V_err, I_err],axis=0)

#function: rule out Kasen model to sig at angle theta
def ruleout(F, Ferr, Fk, Fkerr, theta, sig, lims=None):
    #angle corrected Kasen luminosity
    Fk_theta = Fk*Kasen_isocorr(theta)
    Fk_theta_err = Fkerr*Kasen_isocorr(theta)
    #total error
    Err = np.sqrt(np.square(Ferr)+np.square(Fk_theta_err))
    #which is more constraining? datapoint or limit?
    level = F + sig*Err
    if lims is not None:
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
    #k_opt = 0.1 #0.1g/cm^2
    k_opt = 0.2
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

    if t_day > 0 and Eej/Mej > 0:
        #time in seconds
        ts = t_day*86400.0
        
        #sB constant
        sb_const = 5.6704e-5 #erg/cm^2/s/K^4
        #opacity
        #k_opt = 0.1 #0.1g/cm^2
        k_opt = 0.2
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
        Rcsm = Rext*1.0e13 + vext*ts #cm
        #Temperature evolution of CSM given blackbody
        Tcsm = max(np.power(Lcsm/(4*sb_const*np.pi*Rcsm**2), 0.25), 100) #K
    else:
        Lcsm = 0
        Tcsm = 1000
    return Lcsm, Tcsm

def CSMFit(t_day, wave, z, DM, Mej, Eej, Mext, Rext, t0):
    from SEDAnalysis import BBflux
    #shift time to rest frame
    t_rest = t_day/(1+z) - t0
    #calculate CSM luminosity in rest frame
    Lcsm, Tcsm = CSMmod(t_rest, Eej, Mej, Mext, Rext)
    #shift luminosity to observer frame flux in band
    Fcsm = BBflux(Lcsm,Tcsm,wave,z,DM)
    #return predicted flux in band
    return Fcsm

#function: Error function for multi-band early light curve leastsq fitting
def CSMMultiErr(p, t, L, L_err, z, DM, Mej, Eej):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #CSM component p0=epoch (in rest frame), p1=Mext, p2=Rext
    B_pred = np.array([CSMFit(ti, wave_0[bands['B']], z, DM, Mej, Eej,
                              p[1], p[2], p[0]) for ti in t[0]])
    V_pred = np.array([CSMFit(ti, wave_0[bands['V']], z, DM, Mej, Eej,
                              p[1], p[2], p[0]) for ti in t[1]])
    I_pred = np.array([CSMFit(ti, wave_0[bands['i']], z, DM, Mej, Eej,
                              p[1], p[2], p[0]) for ti in t[2]])
    #Power law component
    B_pred = np.array(B_pred) + earlyFit(t[0], p[3]*(1.+z), p[4], p[7]) 
    V_pred = np.array(V_pred) + earlyFit(t[1], p[3]*(1.+z), p[5], p[8])
    I_pred = np.array(I_pred) + earlyFit(t[2], p[3]*(1.+z), p[6], p[9]) 
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#function: Error function for multi-band early light curve leastsq fitting
def CSMt0MultiErr(p, t, L, L_err, z, DM, Mej, Eej, t0):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #CSM component p0=epoch (in rest frame), p1=Mext, p2=Rext
    B_pred = np.array([CSMFit(ti, wave_0[bands['B']], z, DM, Mej, Eej,
                              p[0], p[1], t0) for ti in t[0]])
    V_pred = np.array([CSMFit(ti, wave_0[bands['V']], z, DM, Mej, Eej,
                              p[0], p[1], t0) for ti in t[1]])
    I_pred = np.array([CSMFit(ti, wave_0[bands['i']], z, DM, Mej, Eej,
                              p[0], p[1], t0) for ti in t[2]])
    #Power law component
    B_pred = np.array(B_pred) + earlyFit(t[0], p[2]*(1.+z), p[3], p[6]) 
    V_pred = np.array(V_pred) + earlyFit(t[1], p[2]*(1.+z), p[4], p[7])
    I_pred = np.array(I_pred) + earlyFit(t[2], p[2]*(1.+z), p[5], p[8]) 
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#function: Error function for multi-band early light curve leastsq fitting
def CSMViErr(p, t, L, L_err, z, DM, Mej, Eej):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #CSM component p0=epoch (in rest frame), p1=Mext, p2=Rext
    V_pred = np.array([CSMFit(ti, wave_0[bands['V']], z, DM, Mej, Eej,
                              p[1], p[2], p[0]) for ti in t[0]])
    I_pred = np.array([CSMFit(ti, wave_0[bands['i']], z, DM, Mej, Eej,
                              p[1], p[2], p[0]) for ti in t[1]])
    #Power law component 
    V_pred = np.array(V_pred) + earlyFit(t[0], p[3]*(1.+z), p[4], p[6])
    I_pred = np.array(I_pred) + earlyFit(t[1], p[3]*(1.+z), p[5], p[7]) 
    #Error
    V_err = (V_pred - L[0])/L_err[0]
    I_err = (I_pred - L[1])/L_err[1]
    return np.concatenate([V_err, I_err],axis=0)

#################################################################
# Composite models.                                             #
#################################################################

#function: Error function for multi-band early light curve leastsq fitting
def CompMultiErr(p, t, L, L_err, z, DM, m_c, e_51, tbase, Lbase):
    from Cosmology import wave_0, bands
    #Kasen component p0=epoch (in rest frame), p1=a13, p2=theta
    B_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['B']],
                                m_c, e_51, z, DM, p[0])
                       for ti in t[0]])*Kasen_isocorr(p[2])
    V_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['V']],
                                m_c, e_51, z, DM, p[0])
                       for ti in t[1]])*Kasen_isocorr(p[2])
    I_pred = np.array([KasenFit(ti, p[1], 1.0, wave_0[bands['i']],
                                m_c, e_51, z, DM, p[0])
                       for ti in t[2]])*Kasen_isocorr(p[2])
    #Base component
    B_base = np.interp(t[0], tbase[0]-p[3], Lbase[0])
    V_base = np.interp(t[1], tbase[1]-p[3], Lbase[1])
    I_base = np.interp(t[2], tbase[2]-p[3], Lbase[2])
    #Add components
    B_pred = B_pred + B_base
    V_pred = V_pred + V_base
    I_pred = I_pred + I_base
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)
