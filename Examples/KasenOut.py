#################################################################
# Name:     Kasen.py                                            #
# Author:   Yuan Qi Ni                                          #
# Date:     May, 31, 2017                                       #
# Function: Program computes early flux light curve calculates  #
#           a kasen result for some given parameters.           #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.stats import norm
from multiprocessing import Pool
import sys

#essential files
from SNAP.Analysis.LCRoutines import*
from SNAP.Analysis.LCFitting import*
from SNAP.Analysis.Cosmology import*

plot = False #plot polynomial fits to light curves

#Epoch
t0 = -18.74
t0err = 0.64 #in rest frame
#t0 = t0 + t0err
#redshift of N300-1.Q0.SN
z = 0.057
zerr = 0.003
#Extinction coefficient (galactic) in each band S & F (2011)
EBVgal = 0.107
Coefs = np.array([3.641, 2.682, 1.516])
band = ['B','V','i']
#max time
Tmax = 281.982
#explosion parameters
m_c = 0.900 #1.26/1.4 Mchandra
e_51 = 0.92 #x10^51 ergs
m_c_err = 0.086
e_51_err = 0.19
#bands
band = ['B','V','i']
Band = ['B','V','I']

print "loading binned data"
#N300-1.Q0.SN binned time series data files
Bfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcbin.CN_170804.txt"
Vfile = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lcbin.CN_170804.txt"
Ifile = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lcbin.CN_170804.txt"
binfiles = [Bfile, Vfile, Ifile] 
print "Loading binned early light curve."
#get N300-1.Q0.SN binned light curve
t, M, M_err, F, SN, Mlim = LCload(binfiles, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, scols=9, flags=['-99.99999'], mode='multi')
#get noise in flux
F_err = [F[i]/SN[i] for i in range(3)]

#deredden flux and get templates
print "Correcting for galactic reddening"
tlim = []
Flim = []
for i in range(len(F)):
    #correct fluxes for galactic reddening
    F[i] = deredFlux(F[i], EBVgal, Coefs[i])
    F_err[i] = deredFlux(F_err[i], EBVgal, Coefs[i])
    #shift to zero centered time
    t[i] = t[i] - Tmax
    #get explosion centered time in observer frame
    t[i] = t[i] - t0

    #get limiting magnitudes
    Mlim[i] = deredMag(Mlim[i], EBVgal, Coefs[i])
    Flim.append(Mag_toFlux(band[i], Mlim[i])*10**6)
    #tlim.append(t[i][SN[i]<2.0])
    #Flim[i] = Flim[i][SN[i]<2.0]

    #get useful times
    t[i], F[i], F_err[i], Flim[i] = LCcrop(t[i], -10,6, F[i], F_err[i], Flim[i])

print "Computing viewing angles at given separation distance"
#list of sample models
a13s = 0.15
a13s = float(sys.argv[-2])
outfilename = sys.argv[-1]
confs = [68.27, 99.54, 99.73]
print [norm.ppf(conf/100.0) for conf in confs]
print a13s
#list of viewing angles
thetas = np.linspace(0,180,100)

#function: test a given a13
def test_a13(a13, sig):
    print a13
    #boolean mask for whether angle is ruled out to given confidence
    mask = np.array([True]*len(thetas))
    #for each band
    for i in range(len(t)):
        #print band[i]
        Fk = np.zeros(len(t[i]))
        Fk_err = np.zeros(len(t[i]))
        for r in range(len(t[i])):
            #get theoretical light curve at each time
            Fk[r], Fk_err[r] = MCerr(KasenFit, [t[i][r], a13, 1.0,
                                                wave_0[bands[band[i]]]],
                                     [m_c, e_51, z, 0],
                                     [m_c_err, e_51_err, zerr, t0err],
                                     [1000000,1000000,1000000,10000000],
                                     nproc=45)
            #print Fk[r], Fk_err[r]
        #if i == 0:
            #print a13, max(Fk), Fk_err[np.argmax(Fk)]
            #tt = np.linspace(0,40,1000)
            #Ft = [KasenFit(ti, a13, 1.0, wave_0[bands[band[i]]], m_c, e_51, z, 0) for ti in tt]
            #print a13, tt[np.argmax(Ft)]
        #for each angle
        for k, theta in enumerate(thetas):
            #check if any points rule out angle with conf
            if ruleout(F[i], F_err[i], Fk, Fk_err, theta, sig):
                mask[k] = False
            #else:
                #print "Consistent!", band[i], a13, norm.cdf(sig), theta
                
    #At this confidence level, we rule out some percent of angles
    outangles = 180.0*float(len(thetas)-len(thetas[mask]))/len(thetas)
    return outangles

style = ['k:', 'k--', 'k-']
outangles = []
#for each confidence interval
for n, conf in enumerate(confs):
    #sigma needed to establish confidence below LC
    sig = norm.ppf(conf/100.0)
    #calculate test of a13
    outangles = test_a13(a13s, sig)
    print outangles, conf
    #output
    outfile = open(outfilename, 'a')
    outfile.write("\t".join([conf, outangles])+"\n")
    outfile.close()
