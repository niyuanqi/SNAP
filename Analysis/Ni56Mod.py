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

    tau_Ni=8.8*86400. #decay time of Ni56 in sec
    tau_Co=9.822e6 #decay time of Co56 in sec

    e_Ni=3.90e10 #erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co

    #tau_m is the timescale of the light-curve
    #tau_m=((k_opt/(beta*c))**0.5)*((10./3.)**(0.25))*M_ejE_K
    tau_m=((k_opt/(beta*c))**0.5)*((6./5.)**(0.25))*M_ejE_K

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
    #calculate ejecta mass, kinetric energy
    Mej = (3.0/10.0)**0.5*MejEK**2*vej
    Kej = (3.0/10.0)*Mej*vej**2
    #calculate errors
    Mejerr = Mej*((2*MejEKerr/MejEK)**2 + (vejerr/vej)**2)**0.5
    Kejerr = Kej*((2*vejerr/vej)**2 + (Mejerr/Mej)**2)**0.5
    #return M[Msun], K[ergs]
    return Mej/M_sun, Mejerr/M_sun, Kej/10**51, Kejerr/10**51


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
    x_range=np.linspace(0.01,1,1000, endpoint=True)
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
def plotNi56mod(tB, LB, LBerr, t_diff, L_diff, Mni, Mej, Ek, beta, x_2, etc):
    from scipy.integrate import simps
    from scipy.special import erfc

    t = np.arange(t_diff/100,8,0.01)
    L = PN13Fit(t, t_diff, L_diff, Mej, Ek, beta, x_2)

    #Constants
    M_sun=2.e33
    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    
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
    x_range=np.linspace(0.01,1,1000, endpoint=True)
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

    import matplotlib.pyplot as plt

    plt.figure(figsize=(6,6))
    ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    ax2 = plt.subplot2grid((3, 1), (2, 0))
    ax = [ax1, ax2]
    
    #plot luminosity over time
    #ax[0].plot(t, L56, 'k', linestyle=':')
    ax[0].plot(t, L, 'k')
    ax[0].errorbar(tB, LB, LBerr, fmt='k+')
    ax[0].plot(etc[0], etc[1], 'k', linestyle=":")
    ax[0].set_ylabel("$Luminosity$ [ergs/s]")
    ax[0].set_ylim([1.5e41,1e43])
    ax[0].set_xlim([-0.1,t[-1]])
    ax[0].set_yscale('log')
    ax[0].axes.get_xaxis().set_ticklabels([])
    
    ax[1].plot(t, dM/(Mej*1.4), 'k')
    ax[1].plot(t, M_ni/(Mni), 'k--')
    ax[1].set_ylabel("$M_{diff}$ / $M_{total}$")
    ax[1].set_ylim([0.0011, 1])
    ax[1].set_xlim([-0.1,t[-1]])
    ax[1].set_yscale('log')
    
    ax[-1].set_xlabel("Days since explosion")
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.show()


"""
#function: diffusion time through some ejecta mass
def diff_time(Mej, Ek):
    #Mej in chandrasekhar masses
    #Ek in 10^51ergs

    k_opt=1.0 #x0.1g/cm^2 (this corresponds to electron scattering)
    t_diff = np.power(Mej*1.4*(k_opt**0.88*Mej**0.32)/(2.0e-2*Ek**0.44),
                      1./1.76)
    return t_diff #in days
"""
