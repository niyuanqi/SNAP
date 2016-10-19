#################################################################
# Name:     LCRoutines.py                                       #
# Author:   Yuan Qi Ni                                          #
# Version:  Oct, 19, 2016                                       #
# Function: Program contains various routines for manipulating  #
#           multi-band light curve files and arrays.            #
#################################################################

#essential modules
import numpy as np

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
