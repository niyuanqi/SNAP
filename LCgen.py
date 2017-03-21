#################################################################
# Name:     LCgen.py                                            #
# Author:   Yuan Qi Ni                                          #
# Date:     March, 21, 2017                                     #
# Function: Program uses MagCalc routine and Source Extractor   #
#           generated light curves to generate light curve file #
#           with both magnitudes for cross comparison.          #
#           Compatible currently with KMTNet nomenclature.      #
#################################################################

#essential modules
import numpy as np
import os
from glob import glob
import math

#essential files
from Analysis.LCRoutines import *
from MagCalc import*
from Catalog import*
from Photometry import*

#object position
RA = 140.92247
DEC = -21.969278
#object name
name = 'KSP-OT-1'
#file prefix
prefix = 'N2784-7.Q1.'
#catalog to use
catname = 'N2784-7.Q1.DPRS.cat'
#year observed
year = 2015
#current time
t_now = "170321_0303"
#user running this code
user = "Chris Ni"
#noise level
SNRnoise = 3.0
#number of reference stars used in each band
nrefs = [1,1,1]
#SExtractor time series data files at source
Bfile = "N2784-7.Q1.B.092341D394-215809D4.150301_0131-150626_0842.0046-0375.0.0066.var.lc.txt"
Vfile = "N2784-7.Q1.V.092341D394-215809D4.150301_0131-150626_0842.0046-0375.0.0066.var.lc.txt"
Ifile = "N2784-7.Q1.I.092341D394-215809D4.150301_0131-150626_0842.0046-0375.0.0066.var.lc.txt"
files = [Bfile, Vfile, Ifile]
#number of rows to skip
skip = 3
#output light curve filenames
outBname = "N2784-7.Q1.B.092341D394-215809D4.150301_0131-150626_0842.0046-0375.0.0066.var.lc.CN_170321.txt"
outVname = "N2784-7.Q1.V.092341D394-215809D4.150301_0131-150626_0842.0046-0375.0.0066.var.lc.CN_170321.txt"
outIname = "N2784-7.Q1.I.092341D394-215809D4.150301_0131-150626_0842.0046-0375.0.0066.var.lc.CN_170321.txt"

#observation filters
bands = ['B','V','I']
bindex = {'B':0, 'V':1, 'I':2}
#observatory positions
observatories = {'A':[210.9383,-31.2712,1143.0], 'S':[339.8104,-32.3789,1762.0], 'C':[70.8040,-30.1672,2167.0]}

#function which fills a row with column entries
def rowGen(to,fo,RAs,DECs,Ms,Ms_err,RAo,DECo,Io,SNo,Mo,Mo_err,Mlim,so):
    sto = padstr("%.5f"%to,10)
    sfo = padstr(fo,18)
    sRAs = padstr("%.1f"%RAs,10)
    sDECs = padstr("%.1f"%DECs,10)
    sMs = padstr("%.3f"%Ms,10)
    sMs_err = padstr("%.3f"%Ms_err,10)
    sRAo = padstr("%.1f"%RAo,10)
    sDECo = padstr("%.1f"%DECo,10)
    sIo = padstr(str(Io)[:9],10)
    sSNo = padstr(str(SNo)[:5],10)
    sMo = padstr("%.3f"%Mo,10)
    sMo_err = padstr("%.3f"%Mo_err,10)
    sMlim = padstr("%.3f"%Mlim,10)
    ss = "   "+so
    out = '\n  '+sto+sf+sRAs+sDECs+sMs+sMs_err+sRAo+sDECo+sIo+sSNo+sMo+sMo_err+sMlim+ss
    return out
#fills first row with column headers
def headGen():
    sto = padstr("OBSDAY"+str(year),10)
    sfo = padstr("STRTXT",18)
    sRAs = padstr("RA_SE(\")",10)
    sDECs = padstr("DEC_SE(\")",10)
    sMs = padstr("MAG_SE",10)
    sMs_err = padstr("MAGERR_SE",10)
    sRAo = padstr("RA_MC(\")",10)
    sDECo = padstr("DEC_MC(\")",10)
    sIo = padstr("Flux_MC",10)
    sSNo = padstr("SNR_MC",10)
    sMo = padstr("MAG_MC",10)
    sMo_err = padstr("MAGERR_MC",10)
    sMlim = padstr("LIM_MC",10)
    ss = "   "+"NOTE"
    out = "\n; "+sto+sf+sRAs+sDECs+sMs+sMs_err+sRAo+sDECo+sIo+sSNo+sMo+sMo_err+sMlim+ss
    return out

#read SExtractor files
#time series
tB, RAB, DECB, MB, MBerr = np.loadtxt(Bfile, usecols=(0,2,3,4,5), skiprows=skip, unpack=True)
fB = np.loadtxt(Bfile, dtype=str, usecols=(1,), skiprows=skip)
tV, RAV, DECV, MV, MVerr = np.loadtxt(Vfile, usecols=(0,2,3,4,5), skiprows=skip, unpack=True)
fV = np.loadtxt(Vfile, dtype=str, usecols=(1,), skiprows=skip)
tI, RAI, DECI, MI, MIerr = np.loadtxt(Ifile, usecols=(0,2,3,4,5), skiprows=skip, unpack=True)
fI = np.loadtxt(Ifile, dtype=str, usecols=(1,), skiprows=skip)
fs = [fB, fV, fI]
ts = [tB, tV, tI]
RAs = [RAB, RAV, RAI]
DECs = [DECB, DECV, DECI]
Ms = [MB, MV, MI]
Ms_err = [MBerr, MVerr, MIerr]

#generate output files
os.system('touch '+outBname)
os.system('touch '+outVname)
os.system('touch '+outIname)
outB = open(outBname, 'a')
outV = open(outVname, 'a')
outI = open(outIname, 'a')
outs = [outB, outV, outI]
#write headers for BVI light curve output files
for i in range(len(bands)):
    outs[i].write("; SOURCE_RA_DEC\t"+str(RA)+"\t"+str(DEC))
    outs[i].write("\n; NUMBER_OF_REFERENCES\t"+str(nrefs[bindex[bands[i]]]))
    outs[i].write("\n; "+str(user)+"\t"+str(t_now))
    outs[i].write(headGen())

#search for fits files with which to construct light curve
files = sorted(glob(prefix+'*.fits'))

#generate light curve
for i in range(len(files)):
    filename = files[i]
    print "Computing file "+str(i+1)+"/"+str(len(files))+": "+filename
    #decipher information from KMTNet filename convention
    fo = filename[10:25]
    band = fo[0]
    #check if SrcExt processed this image
    if band in bands: #valid band selection
        b = bindex[band]
        if fo in fs[bindex[band]]:
            #source extractor did process image
            j = float(np.argwhere(fB==fo))
            ts, RAs, DECs, Ms, Mserr = ts[b][j], RAs[b][j], DECs[b][j], Ms[b][j], Ms_err[b][j]
            Stest = True
        else:
            #source extractor did not process image
            ts, RAs, DECs, Ms, Ms_err = -99.99999, -99.9, -99.9, -99.999, -99.999
            Stest = False        
            
    Mtest = True
    so = "_"
    try: #try to load image
        image, to, wcs = loadFits(filename, verbosity=0)
    except FitsError:
        #image critically failed to load
        Mtest = False
        so = "FITS_ERROR"
        t0 = 0

    if Mtest:
        #get moon ephemeris
        obs = fo[-1]
        loc = observatories[obs]
        time = day_isot(to,year)
        RAmoon, DECmoon = moonEQC(time,loc)
        ALTmoon, AZmoon = moonLC(time,loc)
        #check if moon bright
        if ALTmoon > 15.0:
            so = "MOON_BRIGHT"
        elif ALTmoon > 0.0 and sepAngle((RA,DEC),(RAmoon,DECmoon)) < 90.0:
            so = "MOON_BRIGHT"

    if Mtest:
        try:
            RAo, DECo, Io, SNo, Mo, Mo_err, Mlim = magnitude(image, wcs, 'dprs', catname, (RA,DEC), radius=1000.0, name=name, band=band, fwhm=5.0, limsnr=SNRnoise, satmag=14.0, verbosity=0)
            #check if MagCalc returns nonsense
            if any([math.isnan(Mo),math.isinf(Mo),math.isnan(Mo_err),math.isinf(Mo_err),math.isnan(Io),math.isinf(Io),math.isnan(SNo),math.isinf(SNo)]):
                Io, SNo = -99.99999, -99.99
                Mo, Mo_err = -99.999, -99.999
                if any([math.isnan(Mlim),math.isinf(Mlim)]):
                    Mlim = -99.999
                    RAo, DECo = -99.9, -99.9
                    Mtest = False
            
            if any([math.isnan(Mlim),math.isinf(Mlim)]):
                Mlim = -99.999
                if any([math.isnan(Mo),math.isinf(Mo),math.isnan(Mo_err),math.isinf(Moerr),math.isnan(Io),math.isinf(Io),math.isnan(SNo),math.isinf(SNo)]):
                    Io, SNo = -99.99999, -99.99
                    Mo, Mo_err = -99.999, -99.999
                    RAo, DECo = -99.9, -99.9
                    Mtest = False
                else:
                    RAo = RAo - RA
                    DECo = DECo - DEC
            else:
                RAo = RAo - RA
                DECo = DECo - DEC
                
        except PSFError: #if image PSF cant be extracted
            RAo, DECo, Io, SNo, Mo, Mo_err, Mlim  = -99.9, -99.9, -99.99999, -99.99, -99.999, -99.999, -99.999
            so = "PSF_ERROR"
            Mtest = False
        except: #General catastrophic failure
            RAo, DECo, Io, SNo, Mo, Mo_err, Mlim  = -99.9, -99.9, -99.99999, -99.99, -99.999, -99.999, -99.999
            Mtest = False
    else:
        RAo, DECo, Io, SNo, Mo, Mo_err, Mlim  = -99.9, -99.9, -99.99999, -99.99, -99.999, -99.999, -99.999

    #check for total failure
    if (not Stest) and (not Mtest):
        so = so + "_BAD_IMAGE"
    if Mtest and Mlim < 0:
        so = "INCONV"

    #format output
    out = rowGen(to,fo,RAs,DECs,Ms,Ms_err,RAo,DECo,Io,SNo,Mo,Mo_err,Mlim,so)
    print out+'\n'

    if band in bands:
        outs[b].write(out)
for out in outs:
    out.close()
