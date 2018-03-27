#################################################################
# Name:     LCbgen.py                                           #
# Author:   Yuan Qi Ni                                          #
# Date:     May 23, 2017                                        #
# Function: Program uses MagCalc routine and BinIm routine      #
#           to add binned data to a light curve file.           #
#################################################################

#essential modules
import numpy as np
import os
from glob import glob
import math
import subprocess

#essential files
from SNAP.Analysis.LCRoutines import *
from SNAP.MagCalc import*
from SNAP.Catalog import*
from SNAP.Photometry import*

#N300-1.Q0.SN time series data files
Bfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcAp.CN_170727.txt"
Vfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcAp.CN_170727.txt"
Ifile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcAp.CN_170727.txt"
files = [Bfile, Vfile, Ifile]

#get N300-1.Q0.SN light curve
t, M, M_err, F, SN, Mlim, f = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, fcols=1, scols=9, flags=['-99.99999'], mode='multi')
print f
#time intervals
t_ints = [262.8, 263.0, 264.5, 264.65, 265.0, 266.7, 266.975, 267.06,
          267.5, 268.64, 268.8, 270.167, 270.28, 270.4, 270.647, 270.8,
          272.4, 272.8, 274.22, 274.313, 274.4, 274.6, 274.8]

#lowest limiting magnitude (worst image quality)
lim_lim = 19.0

#file to contain bin files
bindir = '../N300-1.Q0.SN.bin/'

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
t_now = "180327_1100"
#user running this code
user = "Chris Ni"
#noise level
SNRnoise = 1.0
#saturation level
satlvl = 14.0
#number of reference stars used in each band
nrefs = [6,3,2]
#photometric radius
radphot = 1000.0
#observation filters
bands = ['B','V','I']
bindex = {'B':0, 'V':1, 'I':2}

#N300-1.Q0.SN time series data files
outBname = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S1.txt"
outVname = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S1.txt"
outIname = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S1.txt"

#function which fills a row with column entries
def rowGen(to,fo,RAo,DECo,Io,SNo,Mo,Mo_err,Mlim,so):
    sto = padstr("%.5f"%to,10)
    sfo = padstr(fo,27)
    sRAo = padstr("%.1f"%RAo,10)
    sDECo = padstr("%.1f"%DECo,10)
    sIo = padstr(str(Io)[:9],10)
    sSNo = padstr(str(SNo)[:5],10)
    sMo = padstr("%.3f"%Mo,10)
    sMo_err = padstr("%.3f"%Mo_err,10)
    sMlim = padstr("%.3f"%Mlim,10)
    ss = "   "+so
    out = '\n  '+sto+sfo+sRAo+sDECo+sIo+sSNo+sMo+sMo_err+sMlim+ss
    return out
#fills first row with column headers
def headGen():
    sto = padstr("OBSDAY"+str(year),10)
    sfo = padstr("STRTXT",27)
    sRAo = padstr("RA_MC(\")",10)
    sDECo = padstr("DEC_MC(\")",10)
    sIo = padstr("Flux(uJy)",10)
    sSNo = padstr("SNR",10)
    sMo = padstr("MAG_MC",10)
    sMo_err = padstr("MAGERR_MC",10)
    sMlim = padstr("LIM_MC",10)
    ss = "   "+"NOTE"
    out = "\n; "+sto+sfo+sRAo+sDECo+sIo+sSNo+sMo+sMo_err+sMlim+ss
    return out

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

#for each band

bin_ts = []
for i in range(len(bands)):
    bin_ts.append([])
    #cycle through the intervals
    for j in range(len(t_ints)-1):
        #interval boundaries
        t1 = t_ints[j]
        t2 = t_ints[j+1]
        mask = np.logical_and(t[i]>t1, t[i]<t2)
        #filter out bad images
        mask = np.logical_and(mask, Mlim[i]>lim_lim)
        #bin all that remains
        t_bin = t[i][mask].mean()
        bin_ts[i].append(t[i][mask])
        bin_names = f[i][mask]
        print "Binning the following files:"
        print i
        print bin_names
        bin_files = [glob(prefix+name+'*.fits')[0] for name in bin_names]
        out_base = bindir+prefix+bands[i]+'.'+bin_names[0][2:-2]+'-'+bin_names[-1][2:-2]+".coadd."
        out_name = out_base+'fits'
        wt_name = out_base+'weight.fits'
        xml_name = out_base+'xml'
        #swarp files between t1 and t2
        #subprocess.call(['swarp','-COMBINE_TYPE','SUM','-IMAGEOUT_NAME',
        #                 out_name,'-WEIGHTOUT_NAME',wt_name,'-XML_NAME',
        #                 xml_name]+bin_files)

        filename = out_name
        to = t_bin
        fo = bands[i]+bin_names[0][2:-2]+'-'+bin_names[-1][2:-2]
        #compute magnitude at 
        Mtest = True
        so = "_"
        try: #try to load image
            image, to, wcs = loadFits(filename, year=2015, getwcs=True, verbosity=0)
            to = t_bin
        except FitsError:
            #image critically failed to load
            Mtest = False
            so = "FITS_ERROR"
            to = 0
            print "Critical error loading image!"

        if Mtest:
            try:
                RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo = magnitude(image, image, wcs, cattype, catname, (RA,DEC), radius=radphot, aperture=0, name=name, band=bands[i], fwhm=5.0, limsnr=SNRnoise, satmag=satlvl, verbosity=0)
                
                #check if MagCalc returns nonsense
                if any([math.isnan(Mo),math.isinf(Mo),math.isnan(Mo_err),math.isinf(Mo_err)]):
                    Mo, Mo_err = -99.999, -99.999
                    
                if any([math.isnan(Io),math.isinf(Io),math.isnan(SNo),math.isinf(SNo)]):
                    Io, SNo = -99.99999, -99.99
                    if any([math.isnan(Mlimo),math.isinf(Mlimo)]):
                        Mlimo = -99.999
                        RAo, DECo = -99.9, -99.9
                        Mtest = False
            
                if any([math.isnan(Mlimo),math.isinf(Mlimo)]):
                    Mlim = -99.999
                    if any([math.isnan(Io),math.isinf(Io),math.isnan(SNo),math.isinf(SNo)]):
                        Io, SNo = -99.99999, -99.99
                        RAo, DECo = -99.9, -99.9
                        Mtest = False
                    else:
                        RAo = RAo - RA
                        DECo = DECo - DEC
                else:
                    RAo = RAo - RA
                    DECo = DECo - DEC
            
            except PSFError: #if image PSF cant be extracted
                RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo  = -99.9, -99.9, -99.99999, -99.99, -99.999, -99.999, -99.999
                so = "PSF_ERROR"
                Mtest = False
                print "PSF can't be extracted!"
            except: #General catastrophic failure
                RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo  = -99.9, -99.9, -99.99999, -99.99, -99.999, -99.999, -99.999
                Mtest = False
                print "Unknown catastrophic failure!"
        else:
            RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo  = -99.9, -99.9, -99.99999, -99.99, -99.999, -99.999, -99.999

        #check for total failure
        if not Mtest:
            so = so + "_BAD_IMAGE"
        else:
            if any([math.isnan(RAo),math.isinf(RAo),math.isnan(DECo),math.isinf(DECo)]):
                RAo, DECo = 0.0, 0.0
            if Mlimo < 0:
                so = "INCONV"

        #format output
        out = rowGen(to,fo,RAo,DECo,Io,SNo,Mo,Mo_err,Mlimo,so)
        print out+'\n'

        outs[i].write(out)
for out in outs:
    out.close()

print bin_ts
