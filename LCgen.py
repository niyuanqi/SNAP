#################################################################
# Name:     LCgen.py                                            #
# Author:   Yuan Qi Ni                                          #
# Date:     March 21, 2017                                      #
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
RA = 14.263303
DEC = -37.039900
#object name
name = 'N300-1.Q0.SN'
#file prefix
prefix = 'N300-1.Q0.'
#catalog to use
catname = 'N300_1_Q0_SN.csv'
#type of catalog
cattype = 'phot'
#year observed
year = 2015
#current time
t_now = "170503_1500"
#user running this code
user = "Chris Ni"
#noise level
SNRnoise = 3.0
#saturation level
satlvl = 14.0
#number of reference stars used in each band
nrefs = [1,1,1]
#SExtractor time series data files at source
Bfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lc.txt"
Vfile = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lc.txt"
Ifile = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lc.txt"
files = [Bfile, Vfile, Ifile]
#photometric radius
radphot = 1000.0
#output light curve filenames
outBname = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lc.CN_170503.txt"
outVname = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lc.CN_170503.txt"
outIname = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lc.CN_170503.txt"

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
    out = '\n  '+sto+sfo+sRAs+sDECs+sMs+sMs_err+sRAo+sDECo+sIo+sSNo+sMo+sMo_err+sMlim+ss
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
    out = "\n; "+sto+sfo+sRAs+sDECs+sMs+sMs_err+sRAo+sDECo+sIo+sSNo+sMo+sMo_err+sMlim+ss
    return out

#read SExtractor files
#time series
tB, RAB, DECB, MB, MBerr = np.loadtxt(Bfile, usecols=(0,2,3,4,5), comments=';', unpack=True)
fB = np.loadtxt(Bfile, dtype=str, usecols=(1,), comments=';')
tV, RAV, DECV, MV, MVerr = np.loadtxt(Vfile, usecols=(0,2,3,4,5), comments=';', unpack=True)
fV = np.loadtxt(Vfile, dtype=str, usecols=(1,), comments=';')
tI, RAI, DECI, MI, MIerr = np.loadtxt(Ifile, usecols=(0,2,3,4,5), comments=';', unpack=True)
fI = np.loadtxt(Ifile, dtype=str, usecols=(1,), comments=';')
fSs = [fB, fV, fI]
tSs = [tB, tV, tI]
RASs = [RAB, RAV, RAI]
DECSs = [DECB, DECV, DECI]
MSs = [MB, MV, MI]
MSs_err = [MBerr, MVerr, MIerr]

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
    fo = '.'.join(filename.split('.')[2:5])
    band = fo[0]
    #check if SrcExt processed this image
    if band in bands: #valid band selection
        b = bindex[band]
        if fo in fSs[b]:
            #source extractor did process image
            j = int(np.argwhere(fSs[b]==fo))
            ts, RAs, DECs, Ms, Ms_err = tSs[b][j], RASs[b][j], DECSs[b][j], MSs[b][j], MSs_err[b][j]
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
        print "Critical error loading image!"

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
            RAo, DECo, Io, SNo, Mo, Mo_err, Mlim = magnitude(image, wcs, cattype, catname, (RA,DEC), radius=radphot, name=name, band=band, fwhm=5.0, limsnr=SNRnoise, satmag=satlvl, verbosity=0)
            #check if MagCalc returns nonsense
            if any([math.isnan(Mo),math.isinf(Mo),math.isnan(Mo_err),math.isinf(Moerr)]):
                Mo, Mo_err = -99.999, -99.999
            if any([math.isnan(Io),math.isinf(Io),math.isnan(SNo),math.isinf(SNo)]):
                Io, SNo = -99.999, -99.999
                if any([math.isnan(Mlim),math.isinf(Mlim)]):
                    Mlim = -99.999
                    RAo, DECo = -99.9, -99.9
                    Mtest = False
            
            if any([math.isnan(Mlim),math.isinf(Mlim)]):
                Mlim = -99.999
                if any([math.isnan(Io),math.isinf(Io),math.isnan(SNo),math.isinf(SNo)]):
                    Io, SNo = -99.999, -99.999
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
            print "PSF can't be extracted!"
        except: #General catastrophic failure
            RAo, DECo, Io, SNo, Mo, Mo_err, Mlim  = -99.9, -99.9, -99.99999, -99.99, -99.999, -99.999, -99.999
            Mtest = False
            print "Unknown catastrophic failure!"
    else:
        RAo, DECo, Io, SNo, Mo, Mo_err, Mlim  = -99.9, -99.9, -99.99999, -99.99, -99.999, -99.999, -99.999

    #check for total failure
    if (not Stest) and (not Mtest):
        so = so + "_BAD_IMAGE"
    elif Mtest:
        if any([math.isnan(RAo),math.isinf(RAo),math.isnan(DECo),math.isinf(DECo)]):
            RAo, DECo = 0.0, 0.0
        if Mlim < 0:
            so = "INCONV"

    #format output
    out = rowGen(to,fo,RAs,DECs,Ms,Ms_err,RAo,DECo,Io,SNo,Mo,Mo_err,Mlim,so)
    print out+'\n'

    if band in bands:
        outs[b].write(out)
for out in outs:
    out.close()
