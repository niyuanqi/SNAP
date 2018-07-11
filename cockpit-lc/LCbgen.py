#################################################################
# Name:     LCbgen.py                                           #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
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

#essential data
from ObjData import*

#generate names using suffix
outBname = name+'.B.'+suffix
outVname = name+'.V.'+suffix
outIname = name+'.I.'+suffix
#time series data files
files = [outBname, outVname, outIname]

#get light curve
t, M, M_err, F, SN, Mlim, f = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, fcols=1, scols=9, flags=['-99.99999'], mode='multi')
t1 = 85
t2 = 94.5
for i in range(len(M)):
    t[i], M[i], M_err[i], Mlim[i], f[i] = LCcrop(t[i], t1, t2, M[i], M_err[i], Mlim[i], f[i])

#time intervals
t_ints = [85.4, 85.6, 86.0, 87.2, 87.6]

#lowest limiting magnitude (worst image quality)
lim_lim = 19.0

#file to contain bin files
bindir = '../bin/'

#observation filters
bands = ['B','V','I']
bindex = {'B':0, 'V':1, 'I':2}

#generate names using suffix
outBname = name+'.B.'+binsuffix
outVname = name+'.V.'+binsuffix
outIname = name+'.I.'+binsuffix

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

#generate output files if they don't already exist
if os.path.exists(outBname) and os.path.exists(outVname) and os.path.exists(outIname):
    print "Continuing "+outBname
    print "Continuing "+outVname
    print "Continuing "+outIname
    #load list of already processed files
    fB_done = np.loadtxt(outBname, dtype=str, comments=';', usecols=[1])
    fV_done = np.loadtxt(outVname, dtype=str, comments=';', usecols=[1])
    fI_done = np.loadtxt(outIname, dtype=str, comments=';', usecols=[1])
    f_done = [fB_done, fV_done, fI_done]
    outs = [outBname, outVname, outIname]
else:
    print "Starting "+outBname
    print "Starting "+outVname
    print "Starting "+outIname
    os.system('touch '+outBname)
    os.system('touch '+outVname)
    os.system('touch '+outIname)
    outB = open(outBname, 'a')
    outV = open(outVname, 'a')
    outI = open(outIname, 'a')
    outs = [outB, outV, outI]
    #already processed files is empty
    f_done = [[],[],[]]
    #write headers for BVI light curve output files
    for i in range(len(bands)):
        outs[i].write("; SOURCE_RA_DEC\t"+str(ra)+"\t"+str(dec))
        outs[i].write("\n; NUMBER_OF_REFERENCES\t"+str(nrefs[bindex[bands[i]]]))
        outs[i].write("\n; "+str(user)+"\t"+str(t_now))
        outs[i].write(headGen())
        outs[i].close()
    outs = [outBname, outVname, outIname]

#for each band
for i in range(len(bands)):
    #cycle through the intervals
    bin_ts = []
    for j in range(len(t_ints)-1):
        #interval boundaries
        t1 = t_ints[j]
        t2 = t_ints[j+1]
        mask = np.logical_and(t[i]>t1, t[i]<t2)
        #filter out bad images
        mask = np.logical_and(mask, Mlim[i]>lim_lim)
        #bin all that remains
        t_bin = t[i][mask].mean()
        bin_names = f[i][mask]
        print "Binning the following files:"
        bin_files = [glob('../raw/'+prefix+name+'*.fits')[0] for name in bin_names]
        print bin_files
        out_base = bindir+prefix+bands[i]+'.'+bin_names[0][2:-2]+'-'+bin_names[-1][2:-2]+".coadd."
        out_name = out_base+'fits'
        wt_name = out_base+'weight.fits'
        xml_name = out_base+'xml'

        #file information
        filename = out_name
        to = t_bin
        fo = bands[i]+'.'+bin_names[0][2:-2]+'-'+bin_names[-1][2:-2]

        if fo in f_done[bindex[band]]:
            print "Already processed "+fo
        else:
            print "Processing "+fo
            #swarp files between t1 and t2
            subprocess.call(['swarp','-COMBINE_TYPE','SUM','-IMAGEOUT_NAME',
                             out_name,'-WEIGHTOUT_NAME',wt_name,'-XML_NAME',
                             xml_name]+bin_files)
            #compute magnitude at 
            Mtest = True
            so = "_"
            try: #try to load image
                image, to, wcs = loadFits(filename, year=year, getwcs=True, verbosity=0)
                to = t_bin
            except FitsError:
                #image critically failed to load
                Mtest = False
                so = "FITS_ERROR"
                to = 0
                print "Critical error loading image!"

            if Mtest:
                try:
                    RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo = magnitude(image, image, wcs, cattype, catname, (ra,dec), radius=size, aperture=0, psf=1, name=name, band=bands[i], fwhm=5.0, limsnr=SNRnoise, satmag=satlvl, refmag=rellvl, fitsky=True, satpix=1000000000.0, verbosity=0)
                
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
                            RAo = RAo - ra
                            DECo = DECo - dec
                    else:
                        RAo = RAo - ra
                        DECo = DECo - dec
            
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

            if band in bands:
                outfile = open(outs[bindex[band]], 'a')
                outfile.write(out)
                outfile.close()
