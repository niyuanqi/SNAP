#################################################################
# Name:     LCFitting.py                                        #
# Author:   Yuan Qi Ni                                          #
# Version:  Oct, 19, 2016                                       #
# Function: Program contains various models to analyse and fit  #
#           single-band light curves.                           #
#################################################################

#essential modules
import numpy as np

#################################################################
# Light Curve Analysis and Fitting Functions                    #
#################################################################

#function: compute monte carlo polynomial fit and parameters
def LCpolyFit(t, M, M_err=None, order=6, N=None, plot=False):

    import warnings
    import matplotlib.pyplot as plt

    #crop out peak section of light curve
    #mask = np.logical_and(t<20, t>-10)
    #t = t[mask]
    #M = M[mask]
    #M_err = M_err[mask]
    
    #ignore least square coefficient matrix rank deficiency
    warnings.simplefilter('ignore', np.RankWarning)
    if M_err is None:
        #fit light curve
        popt = np.polyfit(t, M, order)
        Mt = np.polyval(popt, t)
            
        #generate continuous light curve
        ta = np.linspace(min(t),max(t),1000)
        Ma = np.polyval(popt, ta)
        if plot:
            #plot fit
            plt.scatter(t, M, c='r')
            plt.plot(ta, Ma)
            plt.show()
        #get point of maximum light
        t_max = ta[np.argmin(Ma)]
        M_max = min(Ma)
        #get deltaM15
        dM15 = np.polyval(popt, t_max+15.0) - min(Ma)
        #add to params list
        params = [t_max, M_max, dM15]
        #return parameters
        return popt, params
        
    else:
        #generate set of N Light curves using M and M_err
        LCtrials = np.zeros((len(t),N))
        #for each time, draw N monte carlo values
        for i in range(len(t)):
            #draw from gaussian centered at M, with sigma M_err
            LCtrials[i] = np.random.normal(M[i],np.absolute(M_err[i]),N)
        LCtrials = LCtrials.T
    
        #initializations
        fits = np.zeros((N,order+1))
        t_maxes, M_maxes, dM15s = np.zeros(N), np.zeros(N), np.zeros(N)
        #For each light curve, extract polynomial fit
        for j, LC in enumerate(LCtrials):
            #fit light curve
            popt = np.polyfit(t, LC, order, w=1.0/M_err)
            Mt = np.polyval(popt, t)
            x2dof = np.sum(np.square((Mt-LC)/M_err))/(len(LC)-len(popt))
            #check fit
            if x2dof > 10:
                print j, x2dof
                print ','.join([str(term) for term in popt])
        
            #generate continuous light curve
            ta = np.linspace(min(t),max(t),1000)
            Ma = np.polyval(popt, ta)
            #get point of maximum light
            t_max = ta[np.argmin(Ma)]
            M_max = min(Ma)
            #get deltaM15
            dM15 = np.polyval(popt, t_max+15.0) - min(Ma)
            #append values to lists
            fits[j] = popt
            t_maxes[j] = t_max
            M_maxes[j] = M_max
            dM15s[j] = dM15
                           
        #average parameters among monte carlo datasets
        fit_err = np.std(fits,0)/np.sqrt(N)
        fit = np.mean(fits,0)
        t_max_err = np.std(t_maxes)/np.sqrt(N)
        t_max = np.mean(t_maxes)
        M_max_err = np.std(M_maxes)/np.sqrt(N)
        M_max = np.mean(M_maxes)
        dM15_err = np.std(dM15s)/np.sqrt(N)
        dM15 = np.mean(dM15s)
        if plot:
            #generate analytic curve
            ta = np.linspace(min(t),max(t),1000)
            Ma = np.polyval(fit, ta)
            #plot fit
            plt.errorbar(t, M, yerr=M_err, fmt='r+')
            plt.plot(ta, Ma)
            plt.plot([t_max,t_max],[M_max,M_max+dM15],c='g')
            plt.show()
        #add to params list
        params = [t_max, M_max, dM15]
        params_err = [t_max_err, M_max_err, dM15_err]
        #return parameters
        return fit, fit_err, params, params_err

#function: 10 parameter Supernova 1a fit function
def SN1aLC(t, g0, t0, sigma0, g1, t1, sigma1, gamma, f0, tau, theta):
    gaus0 = g0*np.exp(-np.square((t-t0)/sigma0)/2)
    gaus1 = g1*np.exp(-np.square((t-t1)/sigma1)/2)
    lin = gamma*(t-t0)
    factor = 1.0 - np.exp((tau-t)/theta)
    return (f0 + lin + gaus0 + gaus1)/factor

#function: compute monte carlo 10 parameter SN1a fit and parameters
def LCSN1aFit(t, M, M_err=None, p0=None, N=30, plot=False):

    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    
    #expected parameters
    t0 = t[np.argmin(M)]
    t1 = t[np.argmin(M)] + 25.0
    if p0 is None:
        p0 = [-1.0, t0, 10.0, -0.5, t1, 8.0, 0.05, 1.0, -25.0, 5.0]
    #check for given errors
    if M_err is None:
        #fit light curve
        popt, pcov = curve_fit(SN1aLC, t, M, p0=p0)
        Mt = SN1aLC(t, *popt)
            
        #generate continuous light curve
        ta = np.linspace(min(t),max(t),1000)
        Ma = SN1aLC(ta, *popt)
        if plot:
            #plot fit
            plt.scatter(t, M, c='r')
            plt.plot(ta, Ma)
            plt.show()
        #get point of maximum light
        t_max = ta[np.argmin(Ma)]
        M_max = min(Ma)
        #get deltaM15
        dM15 = SN1aLC(t_max+15.0, *popt) - min(Ma)
        #add to params list
        params = [t_max, M_max, dM15]
        #return parameters
        return popt, params
        
    else:
        #generate set of N Light curves using M and M_err
        LCtrials = np.zeros((len(t),N))
        #for each time, draw N monte carlo values
        for i in range(len(t)):
            #draw from gaussian centered at M, with sigma M_err
            LCtrials[i] = np.random.normal(M[i],np.absolute(M_err[i]),N)
        LCtrials = LCtrials.T
    
        #initializations
        fits = np.zeros((N,10))
        t_maxes, M_maxes, dM15s = np.zeros(N), np.zeros(N), np.zeros(N)
        #For each light curve, extract polynomial fit
        for j, LC in enumerate(LCtrials):
            #fit light curve
            popt, pcov = curve_fit(SN1aLC, t, LC, sigma=M_err, p0=p0, maxfev=100000)
            Mt = SN1aLC(t, *popt)
            x2dof = np.sum(np.square((Mt-LC)/M_err))/(len(LC)-len(popt))
            #check fit
            if x2dof > 10:
                print j, x2dof
                print ','.join([str(term) for term in popt])
        
            #generate continuous light curve
            ta = np.linspace(min(t),max(t),1000)
            Ma = SN1aLC(ta, *popt)
            #get point of maximum light
            t_max = ta[np.argmin(Ma)]
            M_max = min(Ma)
            #get deltaM15
            dM15 = SN1aLC(t_max+15.0, *popt) - min(Ma)
            #append values to lists
            fits[j] = popt
            t_maxes[j] = t_max
            M_maxes[j] = M_max
            dM15s[j] = dM15
                           
        #average parameters among monte carlo datasets
        fit_err = np.std(fits,0)/np.sqrt(N)
        fit = np.mean(fits,0)
        t_max_err = np.std(t_maxes)/np.sqrt(N)
        t_max = np.mean(t_maxes)
        M_max_err = np.std(M_maxes)/np.sqrt(N)
        M_max = np.mean(M_maxes)
        dM15_err = np.std(dM15s)/np.sqrt(N)
        dM15 = np.mean(dM15s)
        if plot:
            #generate analytic curve
            ta = np.linspace(min(t),max(t),1000)
            Ma = SN1aLC(ta, *fit)
            #plot fit
            plt.errorbar(t, M, yerr=M_err, fmt='r+')
            plt.plot(ta, Ma)
            plt.show()
        #add to params list
        params = [t_max, M_max, dM15]
        params_err = [t_max_err, M_max_err, dM15_err]
        #return parameters
        return fit, fit_err, params, params_err

#function: computes chi squared error between template and LC
def SN1aX2fit(t, M, M_err, template, stretch, t0, M0):
    #typical t0, M0 = 260, -19
    #test data
    MT = SN1aLC(t*stretch - t0, *template) + M0
    #return chi2
    return np.sum(np.square((M-MT)/M_err))

#function: SN1a template (10 parameter function) fit and parameters
def LCtemplateFit(t, M, M_err=None, template=None, plot=False):
    #check if template is given, if none just fit 10 parameter function
    if template is None:
        return LCSN1aFit(t, M, M_err, plot=plot)
    #if no errors given, make all errors one
    if M_err is None:
        M_err = np.ones(len(t))
    
    #perform chi2 fitting of templates parameterized by stretch factor
    stretch = np.linspace(0.5,2.0,100)
    x2list = np.zeros(len(stretch))
    for i, s in enumerate(stretch):
        #synthesize template at stretch
        temp = SN1aLC(t/s, *template)
        x2list[i] = np.sum(np.square((M-temp)/M_err))

    #get least chi2 stretch factor
    s_opt = stretch[np.argmin(x2list)]

    if plot:
        plt.errorbar(t/s_opt, M, yerr=M_err, fmt='r+')
        plt.scatter(t, SN1aLC(t, *template), c='b')
        plt.show()
    return s_opt

#function: Fit Arnett Ni56 mass to bolometric light curve
def ArnettFit(t, M_N, MejE):
    #Inputs
    #################################
    #M_N = Mass of Nickel (Solar Mass)
    #MejE = (Mej^3/Ek)^(1/4)
    #Mej = Mass of ejecta (Solar mass)
    #Ek = Kinetic energy of ejecta (*10^51 ergs)
    #t = time from epoch in days

    #Outputs
    #################################
    #array including time since explosion (days) and luminosity (erg/s)

    from scipy.integrate import simps
    
    #Constants
    M_sun=2.e33
    c=3.e10
    #parameters to be fitted
    M_Ni=M_N*M_sun
    M_ejE_K = MejE*((M_sun)**3/(1.e51))**(0.25)
    #time axis (sec)
    #dt=(np.arange(103*4)/4.+0.25)*86400.
    #dt = np.arange(0.25,103.25,0.25)*86400.
    dt = t*86400.

    beta=13.8 #constant of integration (Arnett 1982)
    k_opt=0.1 #g/cm^2 optical opacity (this corresponds to electron scattering)

    tau_Ni=8.8*86400. #decay time of Ni56 in sec
    tau_Co=9.822e6 #decay time of Co56 in sec

    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co

    #tau_m is the timescale of the light-curve
    tau_m=((k_opt/(beta*c))**0.5)*((10./3.)**(0.25))*M_ejE_K

    #integrate up the A(z) factor where z goes from 0 to x
    int_A=np.zeros(len(dt)) 
    int_B=np.zeros(len(dt)) 
    L_ph=np.zeros(len(dt))

    x=dt/tau_m
    y=tau_m/(2.*tau_Ni)
    s=tau_m*(tau_Co-tau_Ni)/(2.*tau_Co*tau_Ni)

    for i in range(len(dt)):
	z=np.arange(100)*x[i]/100.
	Az=2.*z*np.exp(-2.*z*y+np.square(z))
	Bz=2.*z*np.exp(-2.*z*y+2.*z*s+np.square(z))
	int_A[i]=simps(Az,z)
	int_B[i]=simps(Bz,z)
	L_ph[i]=(M_Ni*np.exp(-1.*np.square(x[i])))*((e_Ni-e_Co)*int_A[i]+e_Co*int_B[i])

    #return results
    return L_ph

#function: Arnett max light error for leastsq fitting
def ArnettMaxErr(p, tmax, tmax_err, Lmax, Lmax_err):
    t_arnett = np.arange(0.25,103.25,0.1)
    L_arnett = ArnettFit(t_arnett, abs(p[0]), abs(p[1]))
    La_max = max(L_arnett)
    ta_max = t_arnett[np.argmax(L_arnett)]
    LX2 = (La_max-Lmax)/Lmax_err
    tX2 = (ta_max-tmax)/tmax_err
    err = np.sqrt(np.square(LX2) + np.square(tX2))
    #helping to choose step size in leastsq
    #print p, np.sqrt(np.sum(np.square(err)))
    return err

#function: break Arnett degeneracy
def ArnettMejE(MejE, MejEerr, vej, vejerr):
    #Constants
    M_sun=2.e33
    #convert to cgs
    MejEK = MejE*((M_sun)**3/(1.e51))**(0.25)
    MejEKerr = MejEerr*((M_sun)**3/(1.e51))**(0.25)
    #calculate ejecta mass, kinetric energy
    Mej = (3.0/10.0)**0.5*MejEK**2*vej
    Kej = (3.0/10.0)*Mej*vej**2
    #calculate errors
    Mejerr = Mej*((2*MejEKerr/MejEK)**2 + (vejerr/vej)**2)**0.5
    Kejerr = Kej*((2*vejerr/vej)**2 + (Mejerr/Mej)**2)**0.5
    #return M[Msun], K[ergs]
    return Mej/M_sun, Mejerr/M_sun, Kej/10**51, Kejerr/10**51

#function: Fit power law to early light curve
def earlyFit(t, t0, C, a):
    return np.concatenate((np.zeros(len(t[t<t0])),C*np.power(t[t>=t0]-t0,a)),axis=0)

#function: Error function for multi-band early light curve leastsq fitting
def earlyMultiErr(p, t, L, L_err):
    B_err = (earlyFit(t[0], p[0], p[1], p[4]) - L[0])/L_err[0]
    V_err = (earlyFit(t[1], p[0], p[2], p[5]) - L[1])/L_err[1]
    I_err = (earlyFit(t[2], p[0], p[3], p[6]) - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#function: normalized planck distribution (wavelength)
def planck(x, T):
    #constants
    wave = x*1e-8 #angstrom to cm
    sb = 5.67051e-5 #erg/s/cm2/K4
    h = 6.6260755e-27 #erg*s
    c = 2.99792458e10 #cm/s
    k = 1.380658e-16 #erg/k
    freq = c/wave #Hz
    #black body luminosity
    integ = sb*np.power(T,4) #ergs/s/cm2
    #planck distribution
    p_rad = (2*h*freq**3/c**2)/(np.exp(h*freq/(k*T))-1.0) #erg/s/cm2/rad2/Hz power per area per solid angle per frequency
    #integrated planck over solid angle
    p_int = np.pi*p_rad #erg/s/cm2/Hz #power per area per frequency
    #normalized planck distribution
    return (p_int/integ) #1/Hz, luminosity density

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
    v9 = 6.0e8 * 1.69 * np.sqrt(e_51/m_c) / 1.0e9
    ti = (1.0e4 * a13 / v9) / 86400.0
    t_day = t_day - ti

    #check validity of kasen
    if t_day > 0:
        # Set-up other parameters from inputs
        L_u = 1.69 # constant related to ejecta density profile.
        vt = 6.0 * 10**8 * L_u * np.sqrt(e_51/m_c) # transition velocity
        v9 = vt / 10**9

        # Equations for Luminosity and Teff
        Lc_iso = 10**43 * a13 * m_c * v9**(7./4.) * kappa**(-3./4.) * t_day**(-1./2.) # (erg/s)
        Teff = 2.5 * 10**4 * a13**(1./4.) * kappa**(-35./36) * t_day**(-37./72.)
        
        #Lc_iso[np.isnan(Lc_iso)] = 0
        #Teff[np.isnan(Teff)] = 1
    else:
        Lc_iso = 0
        Teff = 1000
    return Lc_iso,Teff #erg/s

#function: Observable Kasen model from rest frame
def KasenShift(Lc_angle,Teff,wave,z):
    #give wave in observer frame

    from Cosmology import intDl
    
    #luminosity distance [pc -> cm]
    dl = intDl(z)*3.086*10**18
    Area = 4.0*np.pi*np.square(dl) #cm^2
    #kasen model in observer band
    Lc_angle_wave = planck(wave/(1.0+z),Teff)*Lc_angle/Area
    Lc_angle_wave = np.nan_to_num(Lc_angle_wave)
    #ergs/s/Hz/cm^2, luminosity density in observer frame
    #return in uJy, 10**29 uJy = 1 ergs/s/Hz/cm^2
    return Lc_angle_wave*10**29

#function: Kasen fitting functyion
def KasenFit(t_day,a13,kappa,wave,m_c,e_51,z,t0):
    #shift time to rest frame
    t_rest = (t_day)/(1+z) - t0
    #calculate Kasen luminosity in rest frame
    Lk, Tk = Kasen2010(t_rest,a13,m_c,e_51,kappa)
    #shift luminosity to observer frame flux in band
    Fk = KasenShift(Lk,Tk,wave,z)
    #return predicted flux in band
    return Fk

#function: Kasen isotropic correction for viewing angle
def Kasen_isocorr(theta):
    #param theta: viewing angle (degrees) minimum at 180.
    return 0.982 * np.exp(-(theta/99.7)**2) + 0.018

#function: rule out Kasen model to sig at angle theta
def ruleout(F, Ferr, Fk, Fkerr, theta, sig):
    #angle corrected Kasen luminosity
    Fk_theta = Fk*Kasen_isocorr(theta)
    Fk_theta_err = Fkerr*Kasen_isocorr(theta)
    #total error
    Err = np.sqrt(np.square(Ferr)+np.square(Fk_theta_err))
    #check if any points rule out angle with conf
    if any(Fk_theta > F + sig*Err):
        return True
    else:
        return False

#function: Monte Carlo Error Analysis (independent gaussian errors)
def MCerr(func, ins, params, errs, nums):
    #func : function taking in parameters
    #ins : list of inputs to function
    #params : list of parameters to put into function
    #err : list of error associated with parameters
    #nums: list of number of trials to compute for each parameter
    #np.random.seed(0)
    n = len(params)
    val_errs = np.zeros(n)
    val_means = np.zeros(n)
    #for each parameter
    for i in range(n):
        #perturb parameter N times by STD
        trials = np.random.normal(params[i], errs[i], nums[i])

        vals = np.zeros(nums[i])
        #for each perturbation
        for j in range(nums[i]):
            #calculate value using perturbed perameter
            trial_params = params
            trial_params[i] = trials[j]
            #perform processes in parallel
            vals[j] = func(*(ins+trial_params))
        #error associated with perturbation of parameter
        val_errs[i] = vals.std()
        val_means[i] = vals.mean()
    print ins[1], ins[0], val_errs
    #total summed error associated with all perturbation
    val_err = np.sqrt(np.square(val_errs).sum())
    #function evaluated at central values
    val = val_means.mean()
    #return value and error
    return val, val_err

#function: least chi2 fitting method
def fit_leastchi2(p0, datax, datay, yerr, function, errfunc=False):

    from scipy.optimize import leastsq
    
    if not errfunc:
        #define error function for leastsq
        errfunc = lambda p, x, y, yerr: (function(x,p) - y)/yerr
    else:
        errfunc = function
        
    # Fit
    pfit, perr = leastsq(errfunc, p0, args=(datax, datay, yerr), full_output=0, maxfev=100000)
    return pfit

#function: bootstrap fitting method (Pedro Duarte)
def fit_bootstrap(p0, datax, datay, yerr, function, errfunc=False, perturb=True, n=3000, nproc=4):

    from multiprocessing import Pool
    pool = Pool(nproc)

    popt = fit_leastchi2(p0, datax, datay[0], yerr, function, errfunc)

    if perturb:
        # n random data sets are generated and fitted
        randomDelta = np.random.normal(0., yerr, (n, len(datay)))
        randomdataY = datay + randomDelta
    else:
        randomdataY = datay
    #perform processes asynchronously
    procs = [pool.apply_async(fit_leastchi2, [p0, datax, randY, yerr, function, errfunc]) for randY in randomdataY]
    ps = np.array([proc.get(timeout=10) for proc in procs])
    pool.terminate()
    
    #mean fit parameters
    #mean_pfit = np.mean(ps,0)

    # You can choose the confidence interval that you want for your
    # parameter estimates: 
    Nsigma = 1. # 1sigma gets approximately the same as methods above
                # 1sigma corresponds to 68.3% confidence interval
                # 2sigma corresponds to 95.44% confidence interval
    err_pfit = Nsigma * np.std(ps,0)
    mean_pfit = np.mean(ps,0)

    #pfit_bootstrap = mean_pfit
    pfit_bootstrap = popt
    perr_bootstrap = err_pfit
    return pfit_bootstrap, perr_bootstrap 
