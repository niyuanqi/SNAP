#################################################################
# Name:     Ni56Mod.py                                          #
# Author:   Yuan Qi Ni                                          #
# Version:  Apr. 12, 2018                                       #
# Function: Program contains Ni56 powered light curve models.   #
#################################################################

#essential modules
import numpy as np

#################################################################
# Arnett 1982 based Ni56 model (fixed errata).                  #
#################################################################

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
    #array of luminosity (erg/s)

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
    dt = np.array([dt]) if (isinstance(dt, np.float64) or isinstance(dt, float)) else dt
    n = len(dt)

    beta=13.8 #constant of integration (Arnett 1982)
    k_opt=0.1 #g/cm^2 optical opacity (this corresponds to electron scattering)
    #k_opt = 0.08 #as in Arnett et al. or Li et al. 2019
    tau_Ni=8.8*86400. #decay time of Ni56 in sec
    tau_Co=9.822e6 #decay time of Co56 in sec

    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co

    #tau_m is the timescale of the light-curve
    #tau_m=((k_opt/(beta*c))**0.5)*((10./3.)**(0.25))*M_ejE_K
    tau_m=((k_opt/(beta*c))**0.5)*((6./5.)**(0.25))*M_ejE_K
    #print "tau_m=", tau_m

    #integrate up the A(z) factor where z goes from 0 to x
    int_A=np.zeros(n) 
    int_B=np.zeros(n) 
    L_ph=np.zeros(n)

    x=dt/tau_m
    y=tau_m/(2.*tau_Ni)
    s=tau_m*(tau_Co-tau_Ni)/(2.*tau_Co*tau_Ni)

    for i in range(n):
	z=np.arange(100)*x[i]/100.
	Az=2.*z*np.exp(-2.*z*y+np.square(z))
	Bz=2.*z*np.exp(-2.*z*y+2.*z*s+np.square(z))
	int_A[i]=simps(Az,z)
	int_B[i]=simps(Bz,z)
	L_ph[i]=(M_Ni*np.exp(-1.*np.square(x[i])))*((e_Ni-e_Co)*int_A[i]+e_Co*int_B[i])

    #return results
    return L_ph

#function: break Arnett degeneracy
def ArnettMejE(MejE, MejEerr, vej, vejerr):
    #Constants
    M_sun=2.e33
    #convert to cgs
    MejEK = MejE*((M_sun)**3/(1.e51))**(0.25)
    MejEKerr = MejEerr*((M_sun)**3/(1.e51))**(0.25)
    #calculate ejecta mass, kinetic energy
    Mej = (3.0/10.0)**0.5*MejEK**2*vej
    Kej = (3.0/10.0)*Mej*vej**2
    #calculate errors
    Mejerr = Mej*((2*MejEKerr/MejEK)**2 + (vejerr/vej)**2)**0.5
    Kejerr = Kej*((2*vejerr/vej)**2 + (Mejerr/Mej)**2)**0.5
    #return M[Msun], K[ergs]
    return Mej/M_sun, Mejerr/M_sun, Kej/10**51, Kejerr/10**51

#function: easily get nickel mass from arnett model using max
def ArnettNi56(p, tmax, Lmax):
    return Lmax/ArnettFit(tmax, 1.0, np.absolute(p))

#function: get Ni26 mass from arnett model using max and errors
def ArnettNi56MC(p, tmax, Lmax, Lmax_err, n=100):
    #bootstrap sample intercept
    y = np.random.normal(Lmax, Lmax_err, n)
    Nis = np.absolute(ArnettNi56(p, tmax, y))
    return np.mean(Nis), np.std(Nis)

#function: Arnett error function for determining MejEk parameter
def ArnettMaxErr1(p, tmax, tmax_err):

    from scipy.optimize import fmin

    #get maximum of 
    ta_max = fmin(lambda t: -1*ArnettFit(t, 1.0, np.absolute(p)), 50, (), 0.01)[0]
    tX2 = np.absolute(ta_max-tmax)/tmax_err
    return tX2

#function: fit 2 parameter function to a 2D intercept using Monte Carlo
def ArnettIntercept(tmax, Lmax, tmax_err, Lmax_err, p0=1.2, n=100, nproc=4):

    from scipy.optimize import fmin
    from multiprocessing import Pool

    #to pickle properly, you need dill.
    from multi import apply_async

    pool = Pool(nproc)
    
    #bootstrap sample in x-direction to determine MejEk
    xs = np.random.normal(tmax, tmax_err, n)
    #For each point, solve function
    procs = []
    for i, x in enumerate(xs):
        print str(i+1)+'/'+str(n)
        errfunc = lambda p: ArnettMaxErr1(p, x, tmax_err)
        procs.append(apply_async(pool, fmin, [errfunc, p0, (), 0.001, 0.01]))
    #retrieve processes
    popt = [proc.get()[0] for proc in procs]
    pool.terminate()
    #interpret results
    ME = np.mean(popt, axis=0)
    MEerr = np.std(popt, axis=0)
    #use vertical error to get Ni56
    Ni56, Ni56err = ArnettNi56MC(ME, tmax, Lmax, Lmax_err, n=n)
    #return parameters
    return Ni56, ME, Ni56err, MEerr

#################################################################
# Piro and Nakar 2014 based Ni56 model.                         #
#################################################################

#function: Piro and Nakar model
def PN13Fit(t, t_diff, L_diff, Mej, Ek, beta, x_2, plot=False):
    #Inputs
    #################################
    #t_diff = time for diffusion wave to reach core of ejecta (day)
    #L_diff = luminosity at t_diff (should be only L_direct), (ergs/s)
    #MEej = optional, ejecta mass-energy parameter

    #Fit parameters:
    #Mej = Mass of ejecta (in 1.4 solar masses)
    #Ek = Kinetic energy of ejecta (*10^51 ergs)
    #x_2 = t/t_peak at which half of Ni56 is reached by diffusion wave
    #beta = slope of x56 distribution
    
    #t = time from epoch in days

    #Outputs
    #################################
    #array of luminosity (erg/s)

    from scipy.integrate import simps
    from scipy.special import erfc

    #Constants
    M_sun=2.e33

    #t_diff = diff_time(Mej, Ek)
    
    #time axis, in days
    t = np.array([t]) if (isinstance(t, np.float64) or
                            isinstance(t, float)) else t
    n = len(t)

    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    #k_opt=0.8
    
    #diffusion wave depth in solar masses
    dM = t**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    #diffusion wave depth at peak in solar masses
    dM_diff = t_diff**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    
    tau_Ni=8.8 #decay time of Ni56 in day
    tau_Co=9.822e6/86400. #decay time of Co56 in day
    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    #specific heating rate from Ni56 decay in erg/s/g
    eps = e_Ni*np.exp(-t/tau_Ni) +e_Co*(np.exp(-t/tau_Co) - 
                                        np.exp(-t/tau_Ni))
    #specific heating rate at peak
    eps_diff = e_Ni*np.exp(-t_diff/tau_Ni) +e_Co*(np.exp(-t_diff/tau_Co) - 
                                                  np.exp(-t_diff/tau_Ni))

    #time normalized by rise time
    x = t/t_diff

    #normalize Ni56 distribution
    x_range=np.linspace(0.0005,1,2000, endpoint=True)
    intg_x56 = x_range**0.76/(1.+np.exp(-beta*(x_range-x_2)))
    #Ni56 mass
    M_ni = L_diff/eps_diff/M_sun
    #normalization factor
    norm = M_ni/(1.76*dM_diff)/simps(intg_x56, x_range)
    #Ni56 mass fraction at depth
    X56_x = norm/(1+np.exp(-beta*(x-x_2)))
    #total Ni56 distribution
    X56_range = norm/(1+np.exp(-beta*(x_range-x_2)))
    
    #approx. local heating from Ni56 in erg/s
    L56 = X56_x*dM*M_sun*eps
    
    #integrate lambda function
    intg_lambd = simps(X56_range*x_range**0.76, x_range)
    lambd = (eps/eps_diff)/(intg_lambd/(X56_x*x**1.76))
    #integrate from 0 to t at each epoch
    L_direct = np.zeros(n)
    L_tail = np.zeros(n)
    L_ph = np.zeros(n)
    for i in range(n):
        #before diffusion time is reached
        if x[i] < 1: #there is Ni56 deeper than diffusion depth
            #luminosity due to Ni56 shallower than diffusion depth.
            mask_direct = x_range<=x[i]
            x_direct = x_range[mask_direct]
            X56_direct = X56_range[mask_direct]
            intg_direct = (X56_direct/X56_x[i])*(x_direct/x[i])**1.76/x_direct
            L_direct[i] = L_diff * lambd[i] * simps(intg_direct, x_direct)
            
            #luminosity due to Ni56 deeper than diffusion depth.
            mask_tail = x_range>=x[i]
            x_tail = x_range[mask_tail]
            X56_tail = X56_range[mask_tail]
            intg_tail =  (X56_tail/X56_x[i])*(x_tail/x[i])**1.76/x_tail
            diff_corr = erfc(x_tail/x[i]/np.sqrt(2.))/erfc(1./np.sqrt(2.))
            L_tail[i] = L_diff * lambd[i] * simps(intg_tail*diff_corr, x_tail)
        else: #diffusion depth has exposed all nickel
            #L_direct[i] = 0.35*M_sun*eps[i]
            
            #luminosity due to Ni56 shallower than diffusion depth.
            x_direct = x_range
            X56_direct = X56_range
            intg_direct = (X56_direct/X56_x[i])*(x_direct/x[i])**1.76/x_direct
            L_direct[i] = L_diff * lambd[i] * simps(intg_direct, x_direct)
            
            #no Ni56 is deeper than diffusion depth
            L_tail[i] = 0
        
        #total luminosity
        L_ph[i] = L_direct[i] + L_tail[i]
        
    if plot:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6,6))
        ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=2)
        ax2 = plt.subplot2grid((4, 1), (2, 0))
        ax3 = plt.subplot2grid((4, 1), (3, 0))
        ax = [ax1, ax2, ax3]
        
        #plot luminosity over time
        ax[0].plot(x, L_tail/L_diff, 'k', linestyle='--', label="$L_{tail}$")
        ax[0].plot(x, L_direct/L_diff, 'k', linestyle=':', label="$L_{direct}$")
        ax[0].plot(x, L_ph/L_diff, 'k', label="$L = L_d + L_t$")
        ax[0].set_ylabel("$L/L_{peak}$")
        ax[0].set_ylim([0.0011,2])
        ax[0].set_xlim([-0.05,1])
        ax[0].set_yscale('log')
        ax[0].axes.get_xaxis().set_ticklabels([])

        ax[1].plot(x, L56/L_ph, 'k')
        ax[1].set_ylabel("$L_{56}/L$")
        ax[1].set_ylim([0.0011, 3])
        ax[1].set_xlim([-0.05,1])
        ax[1].set_yscale('log')
        ax[1].axes.get_xaxis().set_ticklabels([])

        X56_diff = norm/(1.+np.exp(-beta*(1.-x_2)))
        print "Ni56 mass:", M_ni
        ax[2].plot(x, X56_x/X56_diff, 'k')
        ax[2].set_yscale('log')
        ax[2].set_ylim([0.001,2.0])
        ax[2].set_ylabel("$X_{56}/X_{56, peak}$")
        
        ax[-1].set_xlabel("$x=t/t_{peak}$")
        ax[-1].set_xlim([-0.05,1])
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.0)
        plt.show()
    
    #return results
    return L_ph

#function: plot Ni56 model
def plotNi56mod(tB, tfit, LB, LBerr, t_diff, L_diff, Mni, Mej, Ek, beta, x_2, etc):
    from scipy.integrate import simps
    from scipy.special import erfc
    
    t = np.arange(t_diff/1000,tfit,0.01)
    L = PN13Fit(t, t_diff, L_diff, Mej, Ek, beta, x_2)

    #Constants
    M_sun=2.e33
    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    #k_opt=0.8
    sb_const = 5.6704e-5 #erg/cm^2/s/K^4
    
    #diffusion wave depth in solar masses
    dM = t**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    #diffusion wave depth at peak in solar masses
    dM_diff = t_diff**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)

    #photospheric radius for calculating temperature [cm]
    rph = 3.0e14 * (t**0.78)*(k_opt**0.11)*(Ek**0.39)/(Mej*1.40)**0.28
    
    tau_Ni=8.8 #decay time of Ni56 in day
    tau_Co=9.822e6/86400. #decay time of Co56 in day
    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    #specific heating rate from Ni56 decay in erg/s/g
    eps = e_Ni*np.exp(-t/tau_Ni) +e_Co*(np.exp(-t/tau_Co) - 
                                        np.exp(-t/tau_Ni))
    #specific heating rate at peak
    eps_diff = e_Ni*np.exp(-t_diff/tau_Ni) +e_Co*(np.exp(-t_diff/tau_Co) - 
                                                  np.exp(-t_diff/tau_Ni))

    #time normalized by rise time
    x = t/t_diff

    #normalize Ni56 distribution
    x_range=np.linspace(0.0005,1,2000, endpoint=True)
    intg_x56 = x_range**0.76/(1.+np.exp(-beta*(x_range-x_2)))
    #Ni56 mass
    M_ni = L_diff/eps_diff/M_sun
    #normalization factor
    norm = M_ni/(1.76*dM_diff)/simps(intg_x56, x_range)
    #Ni56 mass fraction at depth
    X56_x = norm/(1+np.exp(-beta*(x-x_2)))
    #total Ni56 distribution
    X56_range = norm/(1+np.exp(-beta*(x_range-x_2)))

    #quantity of nickel above diffusion depth
    M_ni = np.zeros(len(t))
    for i in range(len(t)):
        mask_ni = x_range<=x[i]
        x_ni = x_range[mask_ni]
        X56_ni = X56_range[mask_ni]
        intg_ni = X56_ni*x_ni**0.76
        M_ni[i] = 1.76*dM[i]*simps(intg_ni, x_ni)
    
    #approx. local heating from Ni56 in erg/s
    M56 = X56_x*dM
    L56 = M56*M_sun*eps
    #color temperature (approximate)
    Tc = np.power(L/(4.*np.pi*rph**2*sb_const), 0.25)

    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    plt.rcParams.update({'font.size': 11})
    plt.rcParams.update({'axes.linewidth': 1.})

    plt.figure(figsize=(6,6))
    ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    ax2 = plt.subplot2grid((3, 1), (2, 0))
    ax = [ax1, ax2]
    
    #plot luminosity over time
    #ax[0].plot(t, L56, 'k', linestyle=':')
    ax[0].plot(t, L, 'k', label=r'L$_{\rm PN14}$')
    ax[0].errorbar(tB, LB, LBerr, fmt='ko', ms=4, mfc='k')
    mask = tB < 5.5
    #mask = tB < 4.2
    ax[0].errorbar(tB[mask], LB[mask], LBerr[mask], ms=4,
                   fmt='ko', mfc='white')
    #ax[0].plot(t, L56, 'k', linestyle="--", label=r'L$_{56}$')
    ax[0].plot(etc[0], etc[1], 'k', linestyle=":", label=r'L$_{\rm Arnett}$')
    ax[0].set_ylabel("Luminosity [erg s$^{-1}$]", fontsize=14)
    ax[0].set_ylim([1.5e41,1.2e43])
    #ax[0].set_ylim([1.5e39,1e43])
    ax[0].set_xlim([-0.1,t[-1]])
    ax[0].set_yscale('log')
    ax[0].axes.get_xaxis().set_ticklabels([])
    ax[0].legend(framealpha=0, fontsize=12, loc='upper left')
    """
    ax[1].plot(t, Tc, 'k')
    arrs = t[::100][1:]
    arrows = Tc[::100][1:]
    for i in range(len(arrs)):
        ax[1].arrow(arrs[i], arrows[i], 0, 1000,
                    head_width=0.1, head_length=500, fc='k', ec='k')
    ax[1].set_ylabel("$T_c$", fontsize=14)
    ax[1].set_ylim([3e3+100, 1.5e4])
    ax[1].set_yscale('log')
    ax[1].set_xlim([-0.1,t[-1]])
    ax[1].axes.get_xaxis().set_ticklabels([])
    ax[1].yaxis.set_minor_formatter(NullFormatter())
    """
    ax[1].plot(t, dM, 'k', label="M$_{diff}$")
    ax[1].plot(t, M56, 'k--', label="M$_{56}$")
    ax[1].set_ylabel("Mass [M$_{\odot}$]", fontsize=14)
    ax[1].set_ylim([0.0011, 0.9])
    ax[1].set_xlim([-0.1,t[-1]])
    ax[1].set_yscale('log')
    ax[1].legend(framealpha=0, fontsize=12)
    
    ax[-1].set_xlabel("Days since explosion", fontsize=14)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.show()

#function: get Ni56 distributions
def Ni56dist(t, t_diff, L_diff, Mej, Ek, beta, x_2):
    #Inputs
    #################################
    #t_diff = time for diffusion wave to reach core of ejecta (day)
    #L_diff = luminosity at t_diff (should be only L_direct), (ergs/s)
    #MEej = optional, ejecta mass-energy parameter

    #Fit parameters:
    #Mej = Mass of ejecta (in 1.4 solar masses)
    #Ek = Kinetic energy of ejecta (*10^51 ergs)
    #x_2 = t/t_peak at which half of Ni56 is reached by diffusion wave
    #beta = slope of x56 distribution
    
    #t = time from epoch in days

    #Outputs
    #################################
    #array of Ni56 fraction

    from scipy.integrate import simps
    from scipy.special import erfc

    #Constants
    M_sun=2.e33

    #t_diff = diff_time(Mej, Ek)
    
    #time axis, in days
    t = np.array([t]) if (isinstance(t, np.float64) or
                            isinstance(t, float)) else t
    n = len(t)

    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    #k_opt=0.8
    
    #diffusion wave depth in solar masses
    dM = t**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    #diffusion wave depth at peak in solar masses
    dM_diff = t_diff**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    
    tau_Ni=8.8 #decay time of Ni56 in day
    tau_Co=9.822e6/86400. #decay time of Co56 in day
    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    #specific heating rate from Ni56 decay in erg/s/g
    eps = e_Ni*np.exp(-t/tau_Ni) +e_Co*(np.exp(-t/tau_Co) - 
                                        np.exp(-t/tau_Ni))
    #specific heating rate at peak
    eps_diff = e_Ni*np.exp(-t_diff/tau_Ni) +e_Co*(np.exp(-t_diff/tau_Co) - 
                                                  np.exp(-t_diff/tau_Ni))

    #time normalized by rise time
    x = t/t_diff

    #normalize Ni56 distribution
    x_range=np.linspace(0.0005,1,2000, endpoint=True)
    intg_x56 = x_range**0.76/(1.+np.exp(-beta*(x_range-x_2)))
    #Ni56 mass
    M_ni = L_diff/eps_diff/M_sun
    #normalization factor
    norm = M_ni/(1.76*dM_diff)/simps(intg_x56, x_range)
    #Ni56 mass fraction at depth
    X56_x = norm/(1+np.exp(-beta*(x-x_2)))
    #Ni56 distribution
    X56 = norm/(1+np.exp(-beta*(x-x_2)))
    return X56

#function: predict observations in some band
def predNi56mod(t, wave, z, DM, taus, t_diff, L_diff, Mej, Ek, beta, x_2):
    from scipy.integrate import simps
    from SEDAnalysis import BBflux

    #t = np.arange(t_diff/1000,twin,0.01)
    tr = t/(1.+z)
    L = PN13Fit(tr, t_diff, L_diff, Mej, Ek, beta, x_2)

    #Constants
    M_sun=2.e33
    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    #k_opt=0.8
    sb_const = 5.6704e-5 #erg/cm^2/s/K^4

    #diffusion wave depth in solar masses
    dM = tr**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    #diffusion wave depth at peak in solar masses
    dM_diff = t_diff**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)

    #photospheric radius for calculating temperature [cm]
    rph = 3.0e14 * (tr**0.78)*(k_opt**0.11)*(Ek**0.39)/(Mej*1.40)**0.28

    tau_Ni=8.8 #decay time of Ni56 in day
    tau_Co=9.822e6/86400. #decay time of Co56 in day
    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    #specific heating rate from Ni56 decay in erg/s/g
    eps = e_Ni*np.exp(-tr/tau_Ni) +e_Co*(np.exp(-tr/tau_Co) - 
                                         np.exp(-tr/tau_Ni))
    #specific heating rate at peak
    eps_diff = e_Ni*np.exp(-t_diff/tau_Ni) +e_Co*(np.exp(-t_diff/tau_Co) - 
                                                  np.exp(-t_diff/tau_Ni))

    #time normalized by rise time
    x = tr/t_diff

    #normalize Ni56 distribution
    x_range=np.linspace(0.0005,1,2000, endpoint=True)
    intg_x56 = x_range**0.76/(1.+np.exp(-beta*(x_range-x_2)))
    #Ni56 mass
    M_ni = L_diff/eps_diff/M_sun
    #normalization factor
    norm = M_ni/(1.76*dM_diff)/simps(intg_x56, x_range)
    #Ni56 mass fraction at depth
    X56_x = norm/(1+np.exp(-beta*(x-x_2)))
    #total Ni56 distribution
    X56_range = norm/(1+np.exp(-beta*(x_range-x_2)))

    #quantity of nickel above diffusion depth
    M_ni = np.zeros(len(t))
    for i in range(len(t)):
        mask_ni = x_range<=x[i]
        x_ni = x_range[mask_ni]
        X56_ni = X56_range[mask_ni]
        intg_ni = X56_ni*x_ni**0.76
        M_ni[i] = 1.76*dM[i]*simps(intg_ni, x_ni)
    
    #approx. local heating from Ni56 in erg/s
    M56 = X56_x*dM
    L56 = M56*M_sun*eps
    #color temperature (approximate)
    Tc = np.power(L*taus/(4.*np.pi*rph**2*sb_const), 0.25)
    #mask = t>0.5
    #print t[mask][::5]
    #print Tc[mask][::5]

    return BBflux(L, Tc, wave, z, DM) #uJy

#function: error function for multi-band fitting of shallow Ni model
def Ni56Err(ts, Ls, L_errs, waves, z, DM, taus, t_diff, L_diff, Mej, Ek, beta, x_2):
    n_bands = len(waves)
    #Compute error in each band
    chi2err = 0
    for i in range(n_bands):
        #mask out positive time
        tmask =  ts[i] > 0
        #predict Ni model luminosity
        L_pred_neg = np.zeros(len(ts[i][np.logical_not(tmask)]))
        L_pred_pos = predNi56mod(ts[i][tmask], waves[i], z, DM, taus[i], 
                                      t_diff, L_diff, Mej, Ek, beta, x_2)
        L_pred = np.concatenate([L_pred_neg, L_pred_pos])
        chi2err += np.sum(np.square((Ls[i] - L_pred*1e-6)/L_errs[i]))
    return chi2err

#function: plot multiband shallow Ni model
def Ni56Plot(ts, Ms, M_errs, bs, z, DM, t_pred, taus, t_diff, L_diff, 
                  Mej, Ek, beta, x_2, rest=True):
    #essential import
    import matplotlib.pyplot as plt
    from SNAP.Analysis.Cosmology import Flux_toMag, bands, wave_0
    
    nbands = len(bs)
    #create Figure
    fig, ax = plt.subplots(nbands, sharex=True)
    
    #plot observations and model fits
    for i in range(nbands):
        L_pred = predNi56mod(t_pred, wave_0[bands[bs[i]]], z, DM, taus[i], 
                             t_diff, L_diff, Mej, Ek, beta, x_2)
        M_pred = Flux_toMag(bs[i], L_pred*1e-6)
        if rest:
            #plot model in SN restframe
            ax[i].errorbar(ts[i]/(1.+z), Ms[i] - DM, M_errs[i], fmt='k+')
            ax[i].plot(t_pred/(1.+z), M_pred - DM, color='r')
            ax[i].set_ylim([max(Ms[i]-DM)+1.0, min(Ms[i]-DM)])
            if i == np.ceil(nbands/2):
                ax[i].set_ylabel("Absolute Magnitude")
        else:
            #plot model in observer frame
            ax[i].errorbar(ts[i], Ms[i], M_errs[i], fmt='k+')
            ax[i].plot(t_pred, M_pred, color='r')
            ax[i].set_ylim([max(Ms[i])+0.2, min(Ms[i])])
            if i == np.ceil(nbands/2):
                ax[i].set_ylabel("Apparent Magnitude")
        ax[i].set_xlim([-1, max(t_pred)])
        ax[i].set_xlabel("Days since explosion.")
            
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.show()

#################################################################
# Shallow layer Ni56 model.                                     #
#################################################################

#function: Piro and Nakar model
def ShallowNiFit(t, t_diff, L_diff, Mej, Ek, beta, x_2, x_s, a_s, plot=False):
    #Inputs
    #################################
    #t_diff = time for diffusion wave to reach core of ejecta (day)
    #L_diff = luminosity at t_diff (should be only L_direct), (ergs/s)
    #MEej = optional, ejecta mass-energy parameter

    #Fit parameters:
    #Mej = Mass of ejecta (in 1.4 solar masses)
    #Ek = Kinetic energy of ejecta (*10^51 ergs)
    #x_2 = t/t_peak at which half of Ni56 is reached by diffusion wave
    #beta = slope of x56 distribution
    
    #t = time from epoch in days

    #Outputs
    #################################
    #array of luminosity (erg/s)

    from scipy.integrate import simps
    from scipy.special import erfc

    #Constants
    M_sun=2.e33

    #t_diff = diff_time(Mej, Ek)
    
    #time axis, in days
    t = np.array([t]) if (isinstance(t, np.float64) or
                            isinstance(t, float)) else t
    n = len(t)

    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    #k_opt=0.8
    
    #diffusion wave depth in solar masses
    dM = t**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    #diffusion wave depth at peak in solar masses
    dM_diff = t_diff**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    
    tau_Ni=8.8 #decay time of Ni56 in day
    tau_Co=9.822e6/86400. #decay time of Co56 in day
    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    #specific heating rate from Ni56 decay in erg/s/g
    eps = e_Ni*np.exp(-t/tau_Ni) +e_Co*(np.exp(-t/tau_Co) - 
                                        np.exp(-t/tau_Ni))
    #specific heating rate at peak
    eps_diff = e_Ni*np.exp(-t_diff/tau_Ni) +e_Co*(np.exp(-t_diff/tau_Co) - 
                                                  np.exp(-t_diff/tau_Ni))

    #time normalized by rise time
    x = t/t_diff

    #normalize Ni56 distribution
    x_range=np.linspace(0.0005,1,2000, endpoint=True)
    intg_x56 = x_range**0.76/(1.+np.exp(-beta*(x_range-x_2)))
    #Ni56 mass
    M_ni = L_diff/eps_diff/M_sun
    #normalization factor
    norm = M_ni/(1.76*dM_diff)/simps(intg_x56, x_range)
    #Ni56 mass fraction at depth
    X56_x = norm/(1+np.exp(-beta*(x-x_2)))
    #X56_x[x < x_s] = a_s*X56_x[x < x_s]
    X56_x[x < x_s] = a_s
    #total Ni56 distribution
    X56_range = norm/(1+np.exp(-beta*(x_range-x_2)))
    #X56_range[x_range < x_s] = a_s*X56_range[x_range < x_s]
    X56_range[x_range < x_s] = a_s
    #import matplotlib.pyplot as plt
    #plt.plot(x, X56_x)
    #plt.show()
    
    #approx. local heating from Ni56 in erg/s
    L56 = X56_x*dM*M_sun*eps
    
    #integrate lambda function
    intg_lambd = simps(X56_range*x_range**0.76, x_range)
    lambd = (eps/eps_diff)/(intg_lambd/(X56_x*x**1.76))
    #integrate from 0 to t at each epoch
    L_direct = np.zeros(n)
    L_tail = np.zeros(n)
    L_ph = np.zeros(n)
    for i in range(n):
        #before diffusion time is reached
        if x[i] < 1: #there is Ni56 deeper than diffusion depth
            #luminosity due to Ni56 shallower than diffusion depth.
            mask_direct = x_range<=x[i]
            x_direct = x_range[mask_direct]
            X56_direct = X56_range[mask_direct]
            intg_direct = (X56_direct/X56_x[i])*(x_direct/x[i])**1.76/x_direct
            L_direct[i] = L_diff * lambd[i] * simps(intg_direct, x_direct)
            
            #luminosity due to Ni56 deeper than diffusion depth.
            mask_tail = x_range>=x[i]
            x_tail = x_range[mask_tail]
            X56_tail = X56_range[mask_tail]
            intg_tail =  (X56_tail/X56_x[i])*(x_tail/x[i])**1.76/x_tail
            diff_corr = erfc(x_tail/x[i]/np.sqrt(2.))/erfc(1./np.sqrt(2.))
            L_tail[i] = L_diff * lambd[i] * simps(intg_tail*diff_corr, x_tail)
        else: #diffusion depth has exposed all nickel
            #L_direct[i] = 0.35*M_sun*eps[i]
            
            #luminosity due to Ni56 shallower than diffusion depth.
            x_direct = x_range
            X56_direct = X56_range
            intg_direct = (X56_direct/X56_x[i])*(x_direct/x[i])**1.76/x_direct
            L_direct[i] = L_diff * lambd[i] * simps(intg_direct, x_direct)
            
            #no Ni56 is deeper than diffusion depth
            L_tail[i] = 0
        
        #total luminosity
        L_ph[i] = L_direct[i] + L_tail[i]
        
    if plot:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6,6))
        ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=2)
        ax2 = plt.subplot2grid((4, 1), (2, 0))
        ax3 = plt.subplot2grid((4, 1), (3, 0))
        ax = [ax1, ax2, ax3]
        
        #plot luminosity over time
        ax[0].plot(x, L_tail/L_diff, 'k', linestyle='--', label="$L_{tail}$")
        ax[0].plot(x, L_direct/L_diff, 'k', linestyle=':', label="$L_{direct}$")
        ax[0].plot(x, L_ph/L_diff, 'k', label="$L = L_d + L_t$")
        ax[0].set_ylabel("$L/L_{peak}$")
        ax[0].set_ylim([0.0011,2])
        ax[0].set_xlim([-0.05,1])
        ax[0].set_yscale('log')
        ax[0].axes.get_xaxis().set_ticklabels([])

        ax[1].plot(x, L56/L_ph, 'k')
        ax[1].set_ylabel("$L_{56}/L$")
        ax[1].set_ylim([0.0011, 3])
        ax[1].set_xlim([-0.05,1])
        ax[1].set_yscale('log')
        ax[1].axes.get_xaxis().set_ticklabels([])

        X56_diff = norm/(1.+np.exp(-beta*(1.-x_2)))
        print "Ni56 mass:", M_ni
        ax[2].plot(x, X56_x/X56_diff, 'k')
        ax[2].set_yscale('log')
        ax[2].set_ylim([0.001,2.0])
        ax[2].set_ylabel("$X_{56}/X_{56, peak}$")
        
        ax[-1].set_xlabel("$x=t/t_{peak}$")
        ax[-1].set_xlim([-0.05,1])
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.0)
        plt.show()
    
    #return results
    return L_ph

#function: predict observations in some band
def predShallowNimod(t, wave, z, DM, taus, t_diff, L_diff, Mej, Ek, beta, x_2, x_s, a_s, prnt=False):
    from scipy.integrate import simps
    from SEDAnalysis import BBflux
    
    #t = np.arange(t_diff/1000,twin,0.01)
    tr = t/(1.+z)
    L = ShallowNiFit(tr, t_diff, L_diff, Mej, Ek, beta, x_2, x_s, a_s)
    
    #Constants
    M_sun=2.e33
    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    #k_opt=0.8
    sb_const = 5.6704e-5 #erg/cm^2/s/K^4
    
    #diffusion wave depth in solar masses
    dM = tr**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    #diffusion wave depth at peak in solar masses
    dM_diff = t_diff**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    if prnt:
        print "Time when clump exposed", tr[tr<x_s*t_diff][-1]
        print "Mass in clump", dM[tr<x_s*t_diff][-1]
        tspec = 0.9
        print "Diffusion mass depth when t=0.9", tspec**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    
    #photospheric radius for calculating temperature [cm]
    rph = 3.0e14 * (tr**0.78)*(k_opt**0.11)*(Ek**0.39)/(Mej*1.40)**0.28
    if prnt:
        print "tau when clump exposed", taus[tr<x_s*t_diff][-1]
        print "Rphot when clump exposed", rph[tr<x_s*t_diff][-1]
        print "Rphot when t=0.9", 3.0e14*(tspec**0.78)*(k_opt**0.11)*(Ek**0.39)/(Mej*1.40)**0.28
    
    tau_Ni=8.8 #decay time of Ni56 in day
    tau_Co=9.822e6/86400. #decay time of Co56 in day
    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    #specific heating rate from Ni56 decay in erg/s/g
    eps = e_Ni*np.exp(-tr/tau_Ni) +e_Co*(np.exp(-tr/tau_Co) - 
                                         np.exp(-tr/tau_Ni))
    #specific heating rate at peak
    eps_diff = e_Ni*np.exp(-t_diff/tau_Ni) +e_Co*(np.exp(-t_diff/tau_Co) - 
                                                  np.exp(-t_diff/tau_Ni))
    
    #time normalized by rise time
    x = tr/t_diff
    
    #normalize Ni56 distribution
    x_range=np.linspace(0.0005,1,2000, endpoint=True)
    intg_x56 = x_range**0.76/(1.+np.exp(-beta*(x_range-x_2)))
    #Ni56 mass
    M_ni = L_diff/eps_diff/M_sun
    #normalization factor
    norm = M_ni/(1.76*dM_diff)/simps(intg_x56, x_range)
    #Ni56 mass fraction at depth
    X56_x = norm/(1+np.exp(-beta*(x-x_2)))
    #X56_x[x < x_s] = a_s*X56_x[x < x_s]
    X56_x[x < x_s] = a_s
    #total Ni56 distribution
    X56_range = norm/(1+np.exp(-beta*(x_range-x_2)))
    #X56_range[x_range < x_s] = a_s*X56_range[x_range < x_s]
    X56_range[x_range < x_s] = a_s
    
    #quantity of nickel above diffusion depth
    M_ni = np.zeros(len(tr))
    for i in range(len(tr)):
        mask_ni = x_range<=x[i]
        x_ni = x_range[mask_ni]
        X56_ni = X56_range[mask_ni]
        intg_ni = X56_ni*x_ni**0.76
        M_ni[i] = 1.76*dM[i]*simps(intg_ni, x_ni)
    
    #approx. local heating from Ni56 in erg/s
    M56 = X56_x*dM
    if prnt:
        print "Ni56 in clump", M56[x<x_s][-1]
    L56 = M56*M_sun*eps
    #color temperature (approximate)
    Tc = np.power(L*taus/(4.*np.pi*rph**2*sb_const), 0.25)
    
    return BBflux(L, Tc, wave, z, DM) #uJy

#function: error function for multi-band fitting of shallow Ni model
def ShallowNiErr(ts, Ls, L_errs, waves, z, DM, taus, t_diff, L_diff, Mej, Ek, beta, x_2, x_s, a_s):
    n_bands = len(waves)
    #Compute error in each band
    chi2err = 0
    for i in range(n_bands):
        #mask out positive time
        tmask =  ts[i] > 0
        #predict Ni model luminosity
        L_pred_neg = np.zeros(len(ts[i][np.logical_not(tmask)]))
        L_pred_pos = predShallowNimod(ts[i][tmask], waves[i], z, DM, taus[i], 
                                      t_diff, L_diff, Mej, Ek, beta, x_2, x_s, a_s)
        L_pred = np.concatenate([L_pred_neg, L_pred_pos])
        chi2err += np.sum(np.square((Ls[i] - L_pred*1e-6)/L_errs[i]))
    return chi2err

#function: plot multiband shallow Ni model
def ShallowNiPlot(ts, Ms, M_errs, bs, z, DM, t_pred, taus, t_diff, L_diff, 
                  Mej, Ek, beta, x_2, x_s, a_s, rest=True, flux=False, prnt=False):
    #essential import
    import matplotlib.pyplot as plt
    from SNAP.Analysis.Cosmology import Flux_toMag, Mag_toFlux, bands, wave_0
    
    nbands = len(bs)
    #create Figure
    fig, ax = plt.subplots(nbands, sharex=True)
    
    #plot observations and model fits
    for i in range(nbands):
        L_pn14 = predNi56mod(t_pred, wave_0[bands[bs[i]]], z, DM, taus[i], 
                             t_diff, L_diff, Mej, Ek, beta, x_2)
        M_pn14 = Flux_toMag(bs[i], L_pn14*1e-6)
        L_pred = predShallowNimod(t_pred, wave_0[bands[bs[i]]], z, DM, taus[i], 
                                  t_diff, L_diff, Mej, Ek, beta, x_2, x_s, a_s, prnt=prnt)
        M_pred = Flux_toMag(bs[i], L_pred*1e-6)

        L_pred2 = predShallowNimod(t_pred, wave_0[bands[bs[i]]], z, DM, taus[i], 
                                   t_diff, L_diff, Mej, Ek, beta, x_2, x_s, a_s/2)
        M_pred2 = Flux_toMag(bs[i], L_pred2*1e-6)
        L_pred3 = predShallowNimod(t_pred, wave_0[bands[bs[i]]], z, DM, taus[i], 
                                   t_diff, L_diff, Mej, Ek, beta, x_2, x_s, a_s*2)
        M_pred3 = Flux_toMag(bs[i], L_pred3*1e-6)
        if not flux:
            #plot magnitudes
            if rest:
                #plot model in SN restframe
                ax[i].errorbar(ts[i]/(1.+z), Ms[i] - DM, M_errs[i], fmt='k+')
                ax[i].plot(t_pred/(1.+z), M_pn14 - DM, color='b', linestyle='--')
                ax[i].plot(t_pred/(1.+z), M_pred - DM, color='r')
                ax[i].plot(t_pred/(1.+z), M_pred2 - DM, color='g')
                ax[i].plot(t_pred/(1.+z), M_pred3 - DM, color='m')
                ax[i].set_ylim([max(Ms[i]-DM)+1.0, min(Ms[i]-DM)])
                if i == np.ceil(nbands/2):
                    ax[i].set_ylabel("Absolute Magnitude")
            else:
                #plot model in observer frame
                ax[i].errorbar(ts[i], Ms[i], M_errs[i], fmt='k+')
                ax[i].plot(t_pred, M_pn14, color='b', linestyle='--')
                ax[i].plot(t_pred, M_pred, color='r')
                ax[i].set_ylim([max(Ms[i])+0.2, min(Ms[i])])
                if i == np.ceil(nbands/2): 
                    ax[i].set_ylabel("Apparent Magnitude")
        else:
            if not rest: 
                Ls, L_errs = Mag_toFlux(bs[i], Ms[i], M_errs[i])
                #plot janskies in observer frame
                ax[i].errorbar(ts[i], Ls*1e6, L_errs*1e6, fmt='k+')
                ax[i].plot(t_pred, L_pn14, color='b', linestyle='--')
                ax[i].plot(t_pred, L_pred, color='r')
                ax[i].plot(t_pred, L_pred2, color='g')
                ax[i].plot(t_pred, L_pred3, color='m')
                #ax[i].set_ylim([max(Ms[i]-DM)+1.0, min(Ms[i]-DM)])
                ax[i].set_yscale('log')
                if i == np.ceil(nbands/2):
                    ax[i].set_ylabel("Flux [uJy]")
            else:
                if isinstance(rest, bool):
                    DM2 = 0
                else:
                    DM2 = 5*np.log10(rest/3.086e19)
                Ls, L_errs = Mag_toFlux(bs[i], Ms[i] - DM + DM2, M_errs[i])
                L_pn14 = Mag_toFlux(bs[i], M_pn14 - DM + DM2)
                L_pred = Mag_toFlux(bs[i], M_pred - DM + DM2)
                L_pred2 = Mag_toFlux(bs[i], M_pred2 - DM + DM2)
                L_pred3 = Mag_toFlux(bs[i], M_pred3 - DM + DM2)
                
                #get suppression (errors) by Janskies
                masko = np.logical_and(ts[i]>0.4, ts[i]<0.7)
                maski = np.logical_and(t_pred>ts[i][masko][0], t_pred<ts[i][masko][-1])
                print bs[i]
                print np.mean(L_pred[maski]) - np.mean(Ls[masko])
                print np.mean(L_pred2[maski]) - np.mean(Ls[masko])
                print np.mean(L_pred3[maski]) - np.mean(Ls[masko])
                
                #plot janskies in rest frame (or certain radius)
                ax[i].errorbar(ts[i], Ls, L_errs, fmt='k+')
                ax[i].plot(t_pred/(1.+z), L_pn14, color='b', linestyle='--')
                ax[i].plot(t_pred/(1.+z), L_pred, color='r')
                ax[i].plot(t_pred/(1.+z), L_pred2, color='g')
                ax[i].plot(t_pred/(1.+z), L_pred3, color='m')
                #ax[i].set_ylim([max(Ms[i]-DM)+1.0, min(Ms[i]-DM)])
                ax[i].set_yscale('log')
                if i == np.ceil(nbands/2):
                    ax[i].set_ylabel("Flux [Jy]")
                    
        ax[i].set_xlim([-1, max(t_pred)])
        ax[i].set_xlabel("Days since explosion.")
        
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.show()

#function: plot multiband shallow Ni model
def ShallowNiColorPlot(ts, Ms, M_errs, bs, z, DM, t_pred, taus, t_diff, L_diff, 
                       Mej, Ek, beta, x_2, x_s, a_s):
    #essential import
    import matplotlib.pyplot as plt
    from SNAP.Analysis.Cosmology import Flux_toMag, bands, wave_0
    from SNAP.Analysis.LCRoutines import LCcolors

    tc, C, C_err = LCcolors(ts, Ms, M_errs)
    
    nbands = len(bs)
    #create Figure
    fig, ax = plt.subplots(nbands-1, sharex=True)
    
    M_14s = []
    M_preds = []
    #plot observations and model fits
    for i in range(nbands):
        L_pn14 = predNi56mod(t_pred, wave_0[bands[bs[i]]], z, DM, taus[i], 
                             t_diff, L_diff, Mej, Ek, beta, x_2)
        M_pn14 = Flux_toMag(bs[i], L_pn14*1e-6)
        L_pred = predShallowNimod(t_pred, wave_0[bands[bs[i]]], z, DM, taus[i], 
                                  t_diff, L_diff, Mej, Ek, beta, x_2, x_s, a_s)
        M_pred = Flux_toMag(bs[i], L_pred*1e-6)
        M_14s.append(M_pn14)
        M_preds.append(M_pred)
    
    C_14s = [M_14s[0]-M_14s[1], M_14s[1]-M_14s[2]]
    C_preds = [M_preds[0]-M_preds[1], M_preds[1]-M_preds[2]]

    ax[0].errorbar(tc[0]/(1.+z), C[0], C_err[0], fmt='k+')
    ax[0].plot(t_pred/(1.+z), C_14s[0], color='b')
    ax[0].plot(t_pred/(1.+z), C_preds[0], color='r')
    ax[0].set_ylabel("B-V")
    ax[1].errorbar(tc[1]/(1.+z), C[1], C_err[1], fmt='k+')
    ax[1].plot(t_pred/(1.+z), C_14s[1], color='b')
    ax[1].plot(t_pred/(1.+z), C_preds[1], color='r')
    ax[1].set_ylabel("V-I")
    ax[1].set_xlim([-1, max(t_pred)])
    ax[1].set_xlabel("Days since explosion")
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.show()

#################################################################
# Ni56 pure-shell model.                                        #
#################################################################

#function: Piro and Nakar model
def NiSFit(t, Mej, Ek, t_s, a_s, plot=False):
    #Inputs
    #################################
    #Parameters:
    #Mej = Mass of ejecta (in 1.4 solar masses)
    #Ek = Kinetic energy of ejecta (*10^51 ergs)
    
    #t = time from epoch in days

    #Outputs
    #################################
    #array of luminosity (erg/s)

    from scipy.integrate import simps
    from scipy.special import erfc

    #Constants
    M_sun=2.e33
    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    tau_Ni=8.8 #decay time of Ni56 in day
    tau_Co=9.822e6/86400. #decay time of Co56 in day
    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    
    #time axis, in days
    t = np.array([t]) if (isinstance(t, np.float64) or
                            isinstance(t, float)) else t
    n = len(t)

    #diffusion wave depth in solar masses
    dM = t**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    #total mass and Ni56 mass in Ni56 shell
    dM_s = t_s**1.76*(2.0e-2*Ek**0.44)/(k_opt**0.88*Mej**0.32)
    M_ni = dM_s*a_s
    #specific heating rate from Ni56 decay in erg/s/g
    eps = e_Ni*np.exp(-t/tau_Ni) +e_Co*(np.exp(-t/tau_Co) - 
                                        np.exp(-t/tau_Ni))
    
    #integration range
    t_range=np.linspace(0.0005,t_s,2000, endpoint=True)
    #Ni56 mass fraction at depth
    X56_t = np.zeros(len(t))
    X56_t[t < t_s] = a_s
    #total Ni56 distribution
    X56_range = np.zeros(len(t_range))
    X56_range[t_range < t_s] = a_s
    #approx. local heating from Ni56 in erg/s
    L56 = X56_t*dM*M_sun*eps

    #integrate from 0 to t_s at each epoch
    L_direct = np.zeros(n)
    L_tail = np.zeros(n)
    L_ph = np.zeros(n)
    for i in range(n):
        #before t_s is reached
        if t[i] < t_s: #there is Ni56 deeper than diffusion depth
            #luminosity due to Ni56 shallower than diffusion depth.
            mask_direct = t_range<=t[i]
            t_direct = t_range[mask_direct]
            X56_direct = X56_range[mask_direct]
            intg_direct = (X56_direct/X56_t[i])*(t_direct/t[i])**1.76/t_direct
            L_direct[i] = 1.76 * L56[i] * simps(intg_direct, t_direct)
            
            #luminosity due to Ni56 deeper than diffusion depth.
            mask_tail = t_range>=t[i]
            t_tail = t_range[mask_tail]
            X56_tail = X56_range[mask_tail]
            intg_tail =  (X56_tail/X56_t[i])*(t_tail/t[i])**1.76/t_tail
            diff_corr = erfc(t_tail/t[i]/np.sqrt(2.))/erfc(1./np.sqrt(2.))
            L_tail[i] = 1.76 * L56[i] * simps(intg_tail*diff_corr, t_tail)
        else: #diffusion depth has exposed all nickel
            #luminosity due to Ni56 shallower than diffusion depth.
            t_direct = t_range
            X56_direct = X56_range
            intg_direct = (X56_direct)*(t_direct/t[i])**1.76/t_direct
            L_direct[i] = 1.76 * dM[i] * M_sun * eps[i] * simps(intg_direct, t_direct)
            
            #no Ni56 is deeper than diffusion depth
            L_tail[i] = 0
        
        #total luminosity
        L_ph[i] = L_direct[i] + L_tail[i]
        
    if plot:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6,6))
        ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=2)
        ax2 = plt.subplot2grid((4, 1), (2, 0))
        ax3 = plt.subplot2grid((4, 1), (3, 0))
        ax = [ax1, ax2, ax3]
        
        #plot luminosity over time
        ax[0].plot(t, L_tail, 'k', linestyle='--', label="$L_{tail}$")
        ax[0].plot(t, L_direct, 'k', linestyle=':', label="$L_{direct}$")
        ax[0].plot(t, L_ph, 'k', label="$L = L_d + L_t$")
        ax[0].set_ylabel("$L$ [ergs/s]")
        #ax[0].set_ylim([0.0011,2])
        ax[0].set_ylim([1e39,1e41])
        ax[0].set_xlim([-0.5,10])
        ax[0].set_yscale('log')
        ax[0].axes.get_xaxis().set_ticklabels([])

        ax[1].plot(t, L56/L_ph, 'k')
        ax[1].set_ylabel("$L_{56}/L$")
        ax[1].set_ylim([0.0011, 3])
        ax[1].set_xlim([-0.5,10])
        ax[1].set_yscale('log')
        ax[1].axes.get_xaxis().set_ticklabels([])

        print "Clump mass:", dM_s
        print "Ni56 mass:", M_ni
        ax[2].plot(t, X56_t, 'k')
        ax[2].set_yscale('log')
        ax[2].set_ylim([0.001,2.0])
        ax[2].set_ylabel("$X_{56}$")
        
        ax[-1].set_xlabel("$x=t$ [days]")
        ax[-1].set_xlim([-0.5,10])
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.0)
        plt.show()
    
    #return results
    return L_ph

#function: predict observations in some band
def predNiSmod(t, wave, z, DM, taus, Mej, Ek, t_s, a_s, t0, prnt=False):
    #Inputs
    #################################
    #Parameters:
    #Mej = Mass of ejecta (in 1.4 solar masses)
    #Ek = Kinetic energy of ejecta (*10^51 ergs)
    
    #t = time from epoch in days

    #Outputs
    #################################
    #array of flux (uJy)
    
    from scipy.integrate import simps
    from SEDAnalysis import BBflux
    
    tr = t/(1.+z) - t0
    L = np.zeros(len(tr))
    L[tr>0] = NiSFit(tr[tr>0], Mej, Ek, t_s, a_s)
    
    #Constants
    M_sun=2.e33
    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    sb_const = 5.6704e-5 #erg/cm^2/s/K^4
    
    #photospheric radius for calculating temperature [cm]
    rph = 3.0e14 * (tr[tr>0]**0.78)*(k_opt**0.11)*(Ek**0.39)/(Mej*1.40)**0.28
    if prnt:
        print "Rphot when clump exposed", rph[tr < t_s/(1.+z)][-1]
    #color temperature (approximate)

    Tc = np.ones(len(tr))
    Tc[tr>0] = np.power(L[tr>0]*taus[tr>0]/(4.*np.pi*rph**2*sb_const), 0.25)
    if prnt:
        print "Tc when clump exposed", Tc[tr < t_s/(1.+z)][-1]
    
    return BBflux(L, Tc, wave, z, DM) #uJy

#function: error function for multi-band fitting of shallow Ni model
def NiSMultiErr(p, t, L, L_err, z, DM, taus, Mej, Ek):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #Ni shell component p0=epoch (in rest frame), p1=t_s, p2=a_s
    B_pred = predNiSmod(t[0], wave_0[bands['B']], z, DM, taus[0], Mej, Ek, p[1], p[2], p[0])
    V_pred = predNiSmod(t[1], wave_0[bands['V']], z, DM, taus[1], Mej, Ek, p[1], p[2], p[0])
    I_pred = predNiSmod(t[2], wave_0[bands['i']], z, DM, taus[2], Mej, Ek, p[1], p[2], p[0])
    #Power law component
    B_pred = B_pred + earlyFit(t[0], p[3]*(1.+z), p[4], p[7]) 
    V_pred = V_pred + earlyFit(t[1], p[3]*(1.+z), p[5], p[8])
    I_pred = I_pred + earlyFit(t[2], p[3]*(1.+z), p[6], p[9])
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    err = np.concatenate([B_err, V_err, I_err],axis=0)
    return np.array(err, dtype=float)

#function: error function for multi-band fitting of shallow Ni model
def NiSTMultiErr(p, t, L, L_err, z, DM, taus, Mej, Ek):
    from Cosmology import wave_0, bands
    from LCFitting import earlyFit
    #Ni shell component p0=epoch (in rest frame), p1=t_s, p2=a_s, p3=tau_s
    B_pred = predNiSmod(t[0], wave_0[bands['B']], z, DM, p[3]*taus[0], Mej, Ek, p[1], p[2], p[0])
    V_pred = predNiSmod(t[1], wave_0[bands['V']], z, DM, p[3]*taus[1], Mej, Ek, p[1], p[2], p[0])
    I_pred = predNiSmod(t[2], wave_0[bands['i']], z, DM, p[3]*taus[2], Mej, Ek, p[1], p[2], p[0])
    #Power law component
    B_pred = B_pred + earlyFit(t[0], p[4]*(1.+z), p[5], p[8]) 
    V_pred = V_pred + earlyFit(t[1], p[4]*(1.+z), p[6], p[9])
    I_pred = I_pred + earlyFit(t[2], p[4]*(1.+z), p[7], p[10])
    #Error
    B_err = (B_pred - L[0])/L_err[0]
    V_err = (V_pred - L[1])/L_err[1]
    I_err = (I_pred - L[2])/L_err[2]
    err = np.concatenate([B_err, V_err, I_err],axis=0)
    return np.array(err, dtype=float)

