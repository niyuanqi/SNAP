#################################################################
# Name:     LCRoutines.py                                       #
# Author:   Yuan Qi Ni                                          #
# Version:  Oct, 19, 2016                                       #
# Function: Program contains various routines for manipulating  #
#           light curve files and deriving fit parameters.      #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import simps
import warnings

#################################################################
# String Formatting Functions                                   #
#################################################################

#function: places right-aligned string in a string of space-padded chars
def padstr(string, l):
    ######################################################
    # Input                                              #
    # -------------------------------------------------- #
    # string: str to be padded                           #
    #      l: length of padded str                       #
    # -------------------------------------------------- #
    # Output                                             #
    # -------------------------------------------------- #
    # string: right-aligned space padded str of length l #
    ######################################################
    pad = l - len(string)
    for i in range(pad):
        string = " "+string
    return string

#function: split value(error) format into value
def splitval(string):
    ######################################
    # Input                              #
    # ---------------------------------- #
    # string: value(error) in str format #
    # ---------------------------------- #
    # Output                             #
    # ---------------------------------- #
    #       : value in float format      #
    ######################################
    try :
        return float(string.split('(')[0])
    except ValueError:
        return string
    
#function: split value(error) format into error
def spliterr(string):
    ######################################
    # Input                              #
    # ---------------------------------- #
    # string: value(error) in str format #
    # ---------------------------------- #
    # Output                             #
    # ---------------------------------- #
    #       : error in float format      #
    ######################################
    try:
        return float(string.split('(')[1].split(')')[0])/100.0
    except IndexError:
        return string
    
#function: format value and error as value(error)
def valerr(val, err):
    #########################################
    # Input                                 #
    # ------------------------------------- #
    #    val: value in float format         #
    #    err: error in float format         #
    # ------------------------------------- #
    # Output                                #
    # ------------------------------------- #
    #       : value(error) in string format #
    #########################################
    return ('%.3f'%val)+'('+str(int(err*1000))+')'

#################################################################
# Light Curve Processing Functions                              #
#################################################################

#function: split data formatting in light curve mag(err)
def LCsplit(valerrs):
    #####################################################################
    # Input                                                             #
    # ----------------------------------------------------------------- #
    # valerrs: list of light curves (eg. in different bands [B, V, I])  #
    #          where each is an array of value(error) in string format. #
    # ----------------------------------------------------------------- #
    # Output                                                            #
    # ----------------------------------------------------------------- #
    #    mags: array of light curves where each is an array of floats   #
    #          representing the magnitude value.                        #
    #    errs: array of light curves where each is an array of floats   #
    #          representing the magnitude error.                        #
    #####################################################################
    mags, errs = [[]]*len(valerrs), [[]]*len(valerrs)
    for i, valerr in enumerate(valerrs):
        mags[i] = np.array([splitval(string) for string in valerr])
        errs[i] = np.array([spliterr(string) for string in valerr])
    return mags, errs

#function: filter bad data from light curve
def LCpurify(ts, mags, errs, strs=None, lims=None, flags=['_'], aflag='sn'):
    #######################################################################
    # Input                                                               #
    # ------------------------------------------------------------------- #
    #   mags: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitudes in float.              #
    #                                                                     #
    #   errs: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitude errors in float.        #
    #                                                                     #
    #     ts: list of time arrays (eg. [tB, tV, tI]) where each is an     #
    #         array of time (in float) corresponding to the light curve.  #
    #                                                                     #
    #   strs; list of string arrays (eg. [sB, sV, sI]) where each is an   #
    #         array of str comments accompanying the light curve.         #
    #                                                                     #
    #   lims; list of limiting magnitude arrays (eg. [Blim, Vlim, Ilim])  #
    #         where each is an array of limiting magnitudes in float      #
    #         evaluated at each sample of its corresponding light curve.  #
    #         If given, routine would purge all corresponding light curve #
    #         values to magnitudes measured below detection threshold.    #
    #                                                                     #
    #  flags; list of str flags defining what values in strs would cause  #
    #         all corresponding light curve values to be purged.          #
    #                                                                     #
    #  aflag; antiflag defining what value in strs would protect light    #
    #         curve values corresponding to a magnitude below the         #
    #         detection limit from purging.                               #
    #         "" indicates don't purge at all based on lims.              #
    # ------------------------------------------------------------------- #
    # Output                                                              #
    # ------------------------------------------------------------------- #
    #     ts: list of light curve time float arrays where bad elements    #
    #         defined by flags and lims were purged.                      #
    #                                                                     #
    #   mags: list of float magnitude light curves where bad elements     #
    #         defined by flags and lims were purged.                      #
    #                                                                     #
    #   errs: list of float error light curves where bad elements         #
    #         defined by flags and lims were purged.                      #
    #                                                                     #
    #   lims; list of float limiting magnitude arrays (if given) where    #
    #         bad elements defined by flags and lims were purged.         #
    #######################################################################
    #for each band
    for i in range(len(ts)):
        #remove elements with flag value for each flag
        for flag in flags:
            index = np.logical_and(mags[i]!=flag,errs[i]!=flag)
            mags[i] = mags[i][index]
            errs[i] = errs[i][index]
            ts[i] = ts[i][index]
            if strs is not None:
                strs[i] = strs[i][index]
            if lims is not None:
                lims[i] = lims[i][index]
        #check for string comments
        if strs is not None:
            #str flag values
            sflags = ['BAD_IMAGE', 'SATURATED', 'FALSE_DET']
            #remove elements with bad sflag value
            for sflag in sflags:
                #apply filter based on strs
                index = strs[i]!=sflag
                mags[i] = mags[i][index]
                errs[i] = errs[i][index]
                ts[i] = ts[i][index]
                strs[i] = strs[i][index]
                if lims is not None:
                    lims[i] = lims[i][index]
        #check for limiting magnitudes
        if lims is not None:
            #remove elements with bad lim value or beyond detection limit
            index = mags[i].astype(float) < lims[i]
            if aflag is not None:
                #antiflag prevents deletion of points below det lim
                matches = np.array([strs[i][j][:len(aflag)]==aflag for j in range(len(strs[i]))])
                index = np.logical_or(index, matches)
            mags[i] = mags[i][index]
            errs[i] = errs[i][index]
            ts[i] = ts[i][index]
            lims[i] = lims[i][index]
            if strs is not None:
                strs[i] = strs[i][index]
    if lims is not None:
        return ts, mags, errs, lims
    else:
        return ts, mags, errs

#function: crop time segment from dataset
def LCcrop(t, t1, t2, M, M_err=None, Mlim=None):
    ##############################################
    # Input                                      #
    # ------------------------------------------ #
    #      M: array of magnitudes in light curve #
    #  M_err; array of magnitude errors          #
    #   Mlim; array of limiting magnitudes       #
    #      t: array of time                      #
    #     t1: start of time segment crop         #
    #     t2: end of time segment crop           #
    # ------------------------------------------ #
    # Output                                     #
    # ------------------------------------------ #
    #      t: cropped time array                 #
    #      M: cropped magnitude array            #
    #  M_err; cropped error array                #
    #   Mlim; cropped limiting magnitude array   #
    ##############################################
    index = np.logical_and(t<t2,t>t1)
    if M_err is not None:
        if Mlim is not None:
            return t[index], M[index], M_err[index], Mlim[index]
        else:
            return t[index], M[index], M_err[index]
    else:
        if Mlim is not None:
            return t[index], M[index], Mlim[index]
        else:
            return t[index], M[index]

#function: return first difference (color) between a set of light curves
def LCcolors(ts, mags, errs):
    #######################################################################
    # Input                                                               #
    # ------------------------------------------------------------------- #
    #   mags: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitudes in float.              #
    #                                                                     #
    #   errs: list of light curves (eg. in different bands [B, V, I])     #
    #         where each is an array of magnitude errors in float.        #
    #                                                                     #
    #     ts: list of time arrays (eg. [tB, tV, tI]) where each is an     #
    #         array of time (in float) corresponding to the light curve.  #
    # ------------------------------------------------------------------- #
    # Output                                                              #
    # ------------------------------------------------------------------- #
    #  tdiff: list of color curve time float arrays.                      #
    #                                                                     #
    #  diffs: list of float first difference color curves.                #
    #                                                                     #
    #  derrs: list of float errors color curves.                          #
    #######################################################################
    tdiff, diffs, derrs = [], [], []
    #for each adjacent pair of light curve in mags
    for i in range(len(ts)-1):
        #create common axis (ordered union)
        tdiff.append(np.sort(np.concatenate((ts[i],ts[i+1]))))
        #interpolate light curves on matching axis
        interp1 = np.interp(tdiff[i],ts[i],mags[i])
        interp1_err2 = np.interp(tdiff[i],ts[i],np.square(errs[i]))
        interp2 = np.interp(tdiff[i],ts[i+1],mags[i+1])
        interp2_err2 = np.interp(tdiff[i],ts[i+1],np.square(errs[i+1]))
        #take first difference with light curves
        diffs.append(interp1 - interp2)
        derrs.append(np.sqrt(interp1_err2 + interp2_err2))
    #return first difference light curves
    return tdiff, diffs, derrs

#function: load light curve from text file
def LCload(filenames, tcol, magcols, errcols=None, limcol=None, scol=None, flags=['_'], aflag='sn', mode='single'):
    #######################################################################
    # Input                                                               #
    # ------------------------------------------------------------------- #
    #      mode: 'single' => read from single txt containing all bands.   #
    #             'multi' => read from many txt each containing one band. #
    #                                                                     #
    # filenames: string file name (single), or list of file names (multi) #
    #                                                                     #
    #      tcol: int location of time column.                             #
    #                                                                     #
    #   magcols: int location of magnitude column (multi) or location of  #
    #            magnitude columns (single).                              #
    #                                                                     #
    #   errcols; int location of error column (multi) or location of      #
    #            error columns (single).                                  #
    #                                                                     #
    #    limcol; int location of limiting magnitude column (multi) or     #
    #            location of limiting magnitude columns (single).         #
    #            If given, rows of light curve will be purged when        #
    #            measured magnitude is below detection threshold.         #
    #                                                                     #
    #      scol; int location of comments column.                         #
    #                                                                     #
    #     flags; list of str flags defining what values in scol would     #
    #            cause a row in the file to be purged.                    #
    #                                                                     #
    #     aflag; antiflag defining what value in strs would protect rows  #
    #            corresponding to a measured magnitude below the          #
    #            detection limit from purging.                            #
    #            "" indicates don't purge at all based on lims.           #
    # ------------------------------------------------------------------- #
    # Output                                                              #
    # ------------------------------------------------------------------- #
    #     ts: list of light curve time float arrays where bad elements    #
    #         defined by flags and lims were purged.                      #
    #                                                                     #
    #   mags: list of float magnitude light curves where bad elements     #
    #         defined by flags and lims were purged.                      #
    #                                                                     #
    #   errs: list of float error light curves where bad elements         #
    #         defined by flags and lims were purged.                      #
    #                                                                     #
    #   lims; list of float limiting magnitude arrays (if given) where    #
    #         bad elements defined by flags and lims were purged.         #
    #######################################################################
    #check load mode, multifile or singlefile
    if mode == 'single':
        #load time column
        t = np.loadtxt(filenames,usecols=(tcol,),comments=';',unpack=True)
        ts = [t]*len(magcols)
        #load mag columns
        mags = list(np.loadtxt(filenames,dtype=str,usecols=magcols,comments=';',unpack=True))
        #retrieve mag errors
        if errcols == 'valerr': #error columns in formatting
            #extract errors from magnitude
            mags, errs = LCsplit(mags)
        elif errcols is not None: #error columns in file
            #load err columns
            errs = list(np.loadtxt(filenames,dtype=str,usecols=errcols,comments=';',unpack=True))
        else: #no errors given
            errs = [np.array(['1.0']*len(t)) for t in ts]

        #check if comment strings are given
        if scol is not None:
            #load comment column
            s = np.loadtxt(filenames,dtype=str,usecols=(scol,),comments=';',unpack=True)
            strs = [s]*len(magcols)
        else:
            strs = [np.array(['_']*len(t)) for t in ts]
        #check if limiting magnitudes are given
        if limcol is not None:
            #load Mlims column
            l = np.loadtxt(filenames,usecols=(limcol,),comments=';',unpack=True)
            lims = [l]*len(magcols)
            #filter out bad data with all information
            ts, mags, errs, lims = LCpurify(ts, mags, errs, strs=strs, lims=lims, flags=flags, aflag=aflag)
        else:
            lims = [np.array([-1.0]*len(t)) for t in ts]
            #filter out bad data with all information
            ts, mags, errs = LCpurify(ts, mags, errs, strs=strs, flags=flags, aflag=aflag)      
            
    elif mode == 'multi':
        #load from each file
        ts, mags = [], []
        for filename in filenames:
            #load time column
            ts.append(np.loadtxt(filename,usecols=(tcol,),comments=';',unpack=True))
            #load mag column
            mags.append(np.loadtxt(filename,dtype=str,usecols=(magcols,),comments=';',unpack=True))
        #retrieve mag errors
        if errcols == 'valerr': #error columns in formatting
            #extract errors from magnitude
            mags, errs = LCsplit(mags)
        elif errcols is not None: #error columns in files
            #extract errors from files
            errs = []
            for filename in filenames:
                #load err columns
                errs.append(np.loadtxt(filename,dtype=str,usecols=(errcols,),comments=';',unpack=True))
        else: #no errors given
            errs = [np.array(['1.0']*len(t)) for t in ts]

        #check if comment strings are given
        if scol is not None:
            #extract comments from files
            strs = []
            for filename in filenames:
                #load comment column
                s = np.loadtxt(filename,dtype=str,usecols=(scol,),comments=';',unpack=True)
                strs.append(s)
        else:
            strs = [np.array(['_']*len(t)) for t in ts]
        #check if limiting magnitudes are given
        if limcol is not None:
            #extract mlims from files
            lims = []
            for filename in filenames:
                #load comment column
                l = np.loadtxt(filename,usecols=(limcol,),comments=';',unpack=True)
                lims.append(l)
            #filter out bad data with all information
            ts, mags, errs, lims = LCpurify(ts, mags, errs, strs=strs, lims=lims, flags=flags, aflag=aflag)
        else:
            lims = [np.array([1.0]*len(t)) for t in ts]
            #filter out bad data with all information
            ts, mags, errs = LCpurify(ts, mags, errs, strs=strs, flags=flags, aflag=aflag)
            
    #convert mags to float
    mags = [mag.astype(float) for mag in mags]
    if errcols is not None:
        if limcol is not None:
            #convert errors to float
            errs = [err.astype(float) for err in errs]
            return ts, mags, errs, lims
        else:
            #convert errors to float
            errs = [err.astype(float) for err in errs]
            return ts, mags, errs
    else:
        if limcol is not None:
            return ts, mags, lims
        else:
            return ts, mags

#################################################################
# Light Curve Analysis and Fitting Functions                    #
#################################################################

#function: compute monte carlo polynomial fit and parameters
def LCpolyFit(t, M, M_err=None, order=6, N=None, plot=False):
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
def ArnettFit(M_N, MejE):
    #Inputs
    #################################
    #M_N = Mass of Nickel (Solar Mass)
    #MejE = (Mej^3/Ek)^(1/4)
    #Mej = Mass of ejecta (Solar mass)
    #Ek = Kinetic energy of ejecta (*10^51 ergs)

    #Outputs
    #################################
    #array including time since explosion (days) and luminosity (erg/s)

    #Constants
    M_sun=2.e33
    c=3.e10
    #parameters to be fitted
    M_Ni=M_N*M_sun
    M_ejE_K = MejE*((M_sun)**3/(1.e51))**(0.25)
    #time axis (sec)
    dt=(np.arange(103*4)/4.+0.25)*86400.

    beta=13.8 #constant of integration (Arnett 1982)
    k_opt=0.07 #g/cm^2 optical opacity (this corresponds to electron scattering)

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
    return dt/86400., L_ph

#function: Fit power law to early light curve
def earlyFit(t, t0, C, a):
    return np.concatenate((np.zeros(len(t[t<t0])),C*np.power(t[t>=t0]-t0,a)),axis=0)

#function: normalized planck distribution (wavelength)
def planck(x, T): 
    sb_const = 5.6704e-5 #erg/s/cm2/K4
    integ = sb_const*np.power(T,4) #ergs/s/cm2
    p_rad = (1.1910430e-5)*np.power(x,-5)/(np.exp(1.439/(x*T))-1.0) #erg/s/cm2/rad2/cm power per area per solid angle per wavelength
    p_int = np.pi*p_rad #erg/s/cm2/cm
    return (1.0/integ)*(p_int) #1/cm

#function: Kasen isotropic correction for viewing angle
def isoAngle(theta):
    return 0.982*np.exp(-np.square(theta/99.7))+0.018

#function: Fit Kasen companion shock model to early bolometric light curve
def KasenFit(dt, wave, Mc, a13, v9, theta):
    #Inputs
    #################################
    #dt [days]
    #wave [angstroms]
    #Mc [chandrasekhar mass]
    #a13 [10^13 cm]
    #v9 [10^9 cm/s]
    #theta [degree]
    
    #Outputs
    #################################
    #array of luminosity density at band (erg/s/cm)

    #Constants
    M_sun=2.e33
    c=3.e10

    k_opt=0.2 #g/cm^2 optical opacity (this corresponds to electron scattering)

    #Luminosity evolution [ergs/s]
    L = 1e43*a13*Mc**(1.0/4)*v9**(7.0/4)*k_opt**(-3.0/4)*np.power(dt,-1.0/2)
    
    #Temperature evolution [K]
    T = 2.5e4*a13**(1.0/4)*k_opt**(-35.0/36)*np.power(dt,-37.0/72)
    print T
    #Angle parameter, isotropic correction []
    f = 0.982*np.exp(-np.square(theta/99.7))+0.018
    
    #Observed Luminosity [ergs/s]
    L_obs = L*f

    #Luminosity density at wavelength [ergs/s/cm]
    L_s = L_obs*planck(wave*1.0e-8,T)
    #luminosity density at wavelength [erg/s/m]
    L_s = L_s*1.0e2
    
    #give band luminosity evolution
    return L_s
