#################################################################
# Name:     LCFitting.py                                        #
# Author:   Yuan Qi Ni                                          #
# Version:  Oct, 19, 2016                                       #
# Function: Program contains various routines to fit multi-band #
#           supernova light curves.                             #
#################################################################

#essential modules
import numpy as np

#################################################################
# Basic Light Curve Models                                      #
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
            dM15 = np.polyval(popt, t_max+15.0) - M_max
            #append values to lists
            fits[j] = popt
            t_maxes[j] = t_max
            M_maxes[j] = M_max
            dM15s[j] = dM15
                           
        #average parameters among monte carlo datasets
        fit_err = np.std(fits,0)
        fit = np.mean(fits,0)
        t_max_err = np.std(t_maxes)
        M_max_err = np.std(M_maxes)
        dM15_err = np.std(dM15s)
        #generate analytic curve
        ta = np.linspace(min(t),max(t),5000)
        Ma = np.polyval(fit, ta)
        #M_max = np.mean(M_maxes)
        t_max = ta[np.argmin(Ma)]
        M_max = min(Ma)
        #dM15 = np.mean(dM15s)
        dM15 = np.polyval(fit, t_max+15.0) - M_max
        if plot:
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

#function: linear
def linfunc(t, a, b):
    return a*t + b

#function: sBV fit from Burns 2014
def sBVfit(t, s0, s1, tmax, tau, c):
    y = 0.5*(s0 + s1)*(t-tmax)+0.5*tau*(s0-s1)*np.log(np.cosh((t-tmax)/tau))+c
    return y
    
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
        fit_err = np.std(fits,0)
        fit = np.mean(fits,0)
        t_max_err = np.std(t_maxes)
        t_max = np.mean(t_maxes)
        M_max_err = np.std(M_maxes)
        M_max = np.mean(M_maxes)
        dM15_err = np.std(dM15s)
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

#function: Fit power law to early light curve
def earlyFit(t, t0, C, a):
    return np.concatenate((np.zeros(len(t[t<t0])),C*np.power(t[t>=t0]-t0,a)),axis=0)

#function: Error function for multi-band early light curve leastsq fitting
def earlyMultiErr(p, t, L, L_err):
    B_err = (earlyFit(t[0], p[0], p[1], p[4]) - L[0])/L_err[0]
    V_err = (earlyFit(t[1], p[0], p[2], p[5]) - L[1])/L_err[1]
    I_err = (earlyFit(t[2], p[0], p[3], p[6]) - L[2])/L_err[2]
    return np.concatenate([B_err, V_err, I_err],axis=0)

#################################################################
# Optimization and Bootstrap routines                           #
#################################################################

#function: Monte Carlo Error Analysis (independent gaussian errors)
def MCerr(func, ins, params, errs, nums, conf, nproc=1):
    #func : function taking in parameters
    #ins : list of inputs to function
    #params : list of parameters to put into function
    #err : list of error associated with parameters
    #nums: list of number of trials to compute for each parameter
    #np.random.seed(0)

    from scipy.stats import norm

    #val = func(*(ins+params))
    n = len(params)
    val_errs = np.zeros(n)
    val_means = np.zeros(n)
    #val_means = np.zeros(n)
    #for each parameter
    for i in range(n):
        #print "computing parameter "+str(i+1)+"/"+str(n)
        #perturb parameter N times by STD
        trials = np.random.normal(params[i], errs[i], nums[i])
        #confidence interval
        conf_int = norm.interval(conf, loc=params[i], scale=errs[i])
        trials = trials[np.logical_and(trials>conf_int[0], trials<conf_int[1])]

        if nproc > 1:
            from multiprocessing import Pool
            pool = Pool(nproc)
            procs = []
            #vals = np.zeros(nums[i])
            #for each perturbation
            for j in range(len(trials)):
                #calculate value using perturbed perameter
                trial_params = np.copy(params)
                trial_params[i] = trials[j]
                #perform processes in parallel
                #vals[j] = func(*(ins+trial_params))
                procs.append(pool.apply_async(func, ins+list(trial_params)))
            vals = np.array([proc.get(timeout=10) for proc in procs])
            pool.terminate()
        else:
            vals = np.zeros(len(trials))
            #for each perturbation
            for j in range(len(trials)):
                #calculate value using perturbed perameter
                trial_params = np.copy(params)
                trial_params[i] = trials[j]
                #perform process
                vals[j] = func(*(ins+list(trial_params)))
        
        #error associated with perturbation of parameter
        val_errs[i] = vals.std()
        val_means[i] = vals.mean()
        #val_means[i] = vals.mean()
    #total summed error associated with all perturbation
    val_err = np.sqrt(np.square(val_errs).sum())
    val = val_means.mean()
    #return value and error
    return val, val_err

#function: least chi2 fitting method
def fit_leastchi2(p0, datax, datay, yerr, function, bounds, errfunc=False):

    from scipy.optimize import least_squares
    
    if not errfunc:
        #define error function for leastsq
        errfunc = lambda p, x, y, yerr: (function(x,p) - y)/yerr
    else:
        errfunc = function
       
    # Fit
    res = least_squares(errfunc, p0, args=(datax, datay, yerr),
                        bounds=bounds, max_nfev=900, verbose=2)
    #pfit, ier = leastsq(errfunc, p0, args=(datax, datay, yerr), full_output=0, maxfev=100000)
    return res.x

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

    

    
