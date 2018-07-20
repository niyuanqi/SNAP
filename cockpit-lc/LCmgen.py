#################################################################
# Name:     LCgen.py                                            #
# Author:   Yuan Qi Ni                                          #
# Date:     July 10, 2017                                       #
# Function: Program uses MagCalc routine to generate light      #
#           curve file of magnitudes and limiting magnitudes.   #
#################################################################

#essential modules
import numpy as np
import os
from glob import glob
import math

#essential files from SNAP
from SNAP.Analysis.LCRoutines import *
from SNAP.MagCalc import*
from SNAP.Catalog import*
from SNAP.Photometry import*
from SNAP.Astrometry import*
#essential data
from ObjData import *

#Single object? Generalize to multiple object. (Compatibility)
if not hasattr(name, '__iter__'):
    #Make all listable objects into list
    name = [name]
    ra = [ra]
    dec = [dec]
    psftype = [psftype]
    fitsky = [fitsky]
#number of sources to perform photometry on
Nobj = len(name)

#observation filters
bands = ['B','V','I']
bindex = {'B':0, 'V':1, 'I':2}
#observatory positions
observatories = {'A':[210.9383,-31.2712,1143.0], 'S':[339.8104,-32.3789,1762.0], 'C':[70.8040,-30.1672,2167.0]}

#function which fills a row with column entries
def rowGen(to,fo,RAo,DECo,Io,SNo,Mo,Mo_err,Mlimo,so):
    sto = padstr("%.5f"%to,10)
    sfo = padstr(fo,27)
    sRAo = padstr("%.7f"%RAo,13)
    sDECo = padstr("%.7f"%DECo,13)
    sIo = padstr(str(Io)[:9],10)
    sSNo = padstr(str(SNo)[:7],10)
    sMo = padstr("%.3f"%Mo,10)
    sMo_err = padstr("%.3f"%Mo_err,10)
    sMlimo = padstr("%.3f"%Mlimo,10)
    ss = "   "+so
    out = '\n  '+sto+sfo+sRAo+sDECo+sIo+sSNo+sMo+sMo_err+sMlimo+ss
    return out
#fills first row with column headers
def headGen():
    sto = padstr("OBSDAY"+str(year),10)
    sfo = padstr("STRTXT",27)
    sRAo = padstr("RA_MC(\")",13)
    sDECo = padstr("DEC_MC(\")",13)
    sIo = padstr("Flux(uJy)",10)
    sSNo = padstr("SNR",10)
    sMo = padstr("MAG_MC",10)
    sMo_err = padstr("MAGERR_MC",10)
    sMlimo = padstr("LIM_MC",10)
    ss = "   "+"NOTE"
    out = "\n; "+sto+sfo+sRAo+sDECo+sIo+sSNo+sMo+sMo_err+sMlimo+ss
    return out

#generate names using suffix
outBname, outVname, outIname = [], [], []
for i in range(Nobj):
    outBname.append(name[i]+'.B.'+suffix)
    outVname.append(name[i]+'.V.'+suffix)
    outIname.append(name[i]+'.I.'+suffix)

outs = []
f_done = []
for i in range(Nobj):
    #generate output files if they don't already exist
    if os.path.exists(outBname[i]) and os.path.exists(outVname[i]) and os.path.exists(outIname[i]):
        print "Continuing "+outBname[i]
        print "Continuing "+outVname[i]
        print "Continuing "+outIname[i]
        #load list of already processed files
        fB_done = np.loadtxt(outBname[i], dtype=str, comments=';', usecols=[1])
        fV_done = np.loadtxt(outVname[i], dtype=str, comments=';', usecols=[1])
        fI_done = np.loadtxt(outIname[i], dtype=str, comments=';', usecols=[1])
        f_done.append([fB_done, fV_done, fI_done])
        print "Already done list:"
        print f_done
        outs.append([outBname[i], outVname[i], outIname[i]])
    else:
        print "Starting "+outBname[i]
        print "Starting "+outVname[i]
        print "Starting "+outIname[i]
        os.system('touch '+outBname[i])
        os.system('touch '+outVname[i])
        os.system('touch '+outIname[i])
        outB = open(outBname[i], 'a')
        outV = open(outVname[i], 'a')
        outI = open(outIname[i], 'a')
        outfiles = [outB, outV, outI]
        #already processed files is empty
        f_done.append([[],[],[]])
        #write headers for BVI light curve output files
        for j in range(len(bands)):
            outfiles[j].write("; SOURCE_NAME : "+name[i])
            outfiles[j].write("\n; SOURCE_RA_DEC\t"+str(ra)+"\t"+str(dec))
            outfiles[j].write("\n; NUMBER_OF_REFERENCES\t"+str(nrefs[bindex[bands[j]]]))
            outfiles[j].write("\n; "+str(user)+"\t"+str(t_now))
            outfiles[j].write(headGen())
            outfiles[j].close()
        outs.append([outBname[i], outVname[i], outIname[i]])

#search for fits files with which to construct light curve
files = sorted(glob('../crop/'+prefix+'*.fits'))

#generate light curve
for i in range(len(files)):
    filename = files[i].split('/')[-1]
    print "Computing file "+str(i+1)+"/"+str(len(files))+": "+filename
    #decipher information from KMTNet filename convention
    fo = '.'.join(filename.split('.')[2:5])
    band = fo[0]

    if fo in f_done[0][bindex[band]]:
        print "Already processed "+fo
    else:
        print "Processing "+fo
        #compute magnitude
        Mtest = True
        so = "_"
        try: #try to load image
            image, to, wcs = loadFits("../crop/"+filename, year=year, getwcs=True, verbosity=0)
        except FitsError:
            #image critically failed to load
            Mtest = False
            so = "FITS_ERROR"
            to = 0
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
            elif ALTmoon > 0.0 and sepAngle((ra[0],dec[0]),(RAmoon,DECmoon)) < 90.0:
                so = "MOON_BRIGHT"

        if Mtest:
            try:
                
                RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo = magnitude(image, image, wcs, cattype, catname, (ra,dec), radius=size, psf=psftype, name=name, band=band, fwhm=5.0, limsnr=SNRnoise, satmag=satlvl, refmag=rellvl, fitsky=fitsky, satpix=satpix, verbosity=0)
                    
                """
                #Replace with this block for more stable measurements
                #when psftype[0]=2
                print "Let's try photometry with fixed centroid. psf=1"
                psfmod = list(psftype)
                psfmod[0] = 1
                RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo = magnitude(image, image, wcs, cattype, catname, (ra,dec), radius=size, psf=psfmod, name=name, band=band, fwhm=5.0, limsnr=SNRnoise, satmag=satlvl, refmag=rellvl, fitsky=fitsky, satpix=satpix, verbosity=0)
                if SNo[0]>SNRnoise:
                    print "Source is bright, get a fix on centroid. psf=2"
                    RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo = magnitude(image, image, wcs, cattype, catname, (ra,dec), radius=size, psf=psftype, name=name, band=band, fwhm=5.0, limsnr=SNRnoise, satmag=satlvl, refmag=rellvl, fitsky=fitsky, satpix=satpix, verbosity=0)
                """
            
                #check if MagCalc returns nonsense
                for j in range(Nobj):
                    if any([math.isnan(Mo[j]),math.isinf(Mo[j]),math.isnan(Mo_err[j]),math.isinf(Mo_err[j])]):
                        Mo[j], Mo_err[j] = -99.999, -99.999
            
                    if any([math.isnan(Io[j]),math.isinf(Io[j]),math.isnan(SNo[j]),math.isinf(SNo[j])]):
                        Io[j], SNo[j] = -99.99999, -99.99
                        if any([math.isnan(Mlimo),math.isinf(Mlimo)]):
                            Mlimo = -99.999
                            RAo[j], DECo[j] = -99.9999999, -99.9999999
                            Mtest = False
            
                    if any([math.isnan(Mlimo),math.isinf(Mlimo)]):
                        Mlimo = -99.999
                        if any([math.isnan(Io[j]),math.isinf(Io[j]),math.isnan(SNo[j]),math.isinf(SNo[j])]):
                            Io[j], SNo[j] = -99.99999, -99.99
                            RAo[j], DECo[j] = -99.9999999, -99.9999999
                            Mtest = False
                        else:
                            RAo[j] = RAo[j]
                            DECo[j] = DECo[j]
                    else:
                        RAo[j] = RAo[j]
                        DECo[j] = DECo[j]

                    #Maybe magnitude error is insane
                    if Mo_err[j] > 1000:
                        Mo_err[j] = -99.999

            except PSFError: #if image PSF cant be extracted
                RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo  = [-99.9999999]*Nobj, [-99.9999999]*Nobj, [-99.99999]*Nobj, [-99.99]*Nobj, [-99.999]*Nobj, [-99.999]*Nobj, -99.999
                so = "PSF_ERROR"
                Mtest = False
                print "PSF can't be extracted!"
            except: #General catastrophic failure
                RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo  = [-99.9999999]*Nobj, [-99.9999999]*Nobj, [-99.99999]*Nobj, [-99.99]*Nobj, [-99.999]*Nobj, [-99.999]*Nobj, -99.999
                Mtest = False
                print "Unknown catastrophic failure!"

        else:
            RAo, DECo, Io, SNo, Mo, Mo_err, Mlimo  = [-99.9999999]*Nobj, [-99.9999999]*Nobj, [-99.99999]*Nobj, [-99.99]*Nobj, [-99.999]*Nobj, [-99.999]*Nobj, -99.999

        #check for total failure
        if not Mtest:
            so = so + "_BAD_IMAGE"
        else:
            for j in range(Nobj):
                if any([math.isnan(RAo[j]),math.isinf(RAo[j]),math.isnan(DECo[j]),math.isinf(DECo[j])]):
                    RAo[j], DECo[j] = -99.9999999, -99.9999999
            if Mlimo < 0:
                so = "INCONV"
        
        print ""
        for j in range(Nobj):
            #format output
            out = rowGen(to,fo,RAo[j],DECo[j],Io[j],SNo[j],Mo[j],Mo_err[j],Mlimo,so)
            print name[j]+" : "+out

            if band in bands:
                outfile = open(outs[j][bindex[band]], 'a')
                outfile.write(out)
                outfile.close()
        print ""
