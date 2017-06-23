#################################################################
# Name:     Kasen.py                                            #
# Author:   Yuan Qi Ni                                          #
# Date:     May, 31, 2017                                       #
# Function: Program computes early flux light curve and fits    #
#           for explosion epoch.                                #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from snpy import *

#essential files
from SNAP.Analysis.LCRoutines import*
from SNAP.Analysis.LCFitting import*
from SNAP.Analysis.Cosmology import*

plot = False #plot polynomial fits to light curves

#Epoch
t0 = -16.71
t0err = 0.52
#redshift of N300-1.Q0.SN
z = 0.056
zerr = 0.005
#Extinction coefficient (galactic) in each band S & F (2011)
EBVgal = 0.107
Coefs = np.array([3.641, 2.682, 1.516])
band = ['B','V','i']
#max time
Tmax = 281.3
#explosion parameters
m_c = 1.17/1.4 #Mchandra
e_51 = 0.85 #x10^51 ergs
#bands
band = ['B','V','i']

print "loading binned data"
#N300-1.Q0.SN binned time series data files
Bfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcbin.CN_170505.txt"
Vfile = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lcbin.CN_170505.txt"
Ifile = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lcbin.CN_170505.txt"
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
    #get epoch centered dilated time
    t[i] = absTime(t[i], z) - t0

    #get limiting magnitudes
    Mlim[i] = deredMag(Mlim[i], EBVgal, Coefs[i])
    Flim.append(Mag_toFlux(band[i], Mlim[i])*10**6)
    tlim.append(t[i][SN[i]<2.0])
    Flim[i] = Flim[i][SN[i]<2.0]

    #get reliable detections
    t[i] = t[i][SN[i]>=2.0]
    F[i] = F[i][SN[i]>=2.0]
    F_err[i] = F_err[i][SN[i]>=2.0]
    
print "plotting early data"
#plot
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("L/Lmax")
for i in range(len(t)):
    #fit early light curve
    ax[i].errorbar(t[i],F[i],F_err[i],fmt='g+')
    ax[i].scatter(tlim[i], Flim[i], marker='v', c='r')
plt.show()

print "Trial thetas:", np.linspace(0,180,10)
tk = np.arange(0.001, 10, 0.001)
###################################################################

print "Computing 1Ms RG Kasen model"
M_comp = 1.0/1.4 #Mchandra
a13 = 2.0 #10^13cm
#plot comparison
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("Flux [uJy]")
for i in range(len(t)):
    #plot I band fluxes
    thetas = np.linspace(0,180,10)
    for theta in thetas:
        #rest time angle corrected Kasen luminosity
        Lk, Lk_theta, Tk = Kasen2010(tk, a13, theta, m_c, e_51)
        #rest time observer frame angle corrected Kasen flux
        Fk, Fk_err = KasenFit(tk,Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
        ax[i].plot(tk, Fk)
    #plot red corrected observer frame rest time observed fluxes
    ax[i].errorbar(t[i],F[i],F_err[i],fmt='g+')
    ax[i].scatter(tlim[i], Flim[i], marker='v', c='r')
plt.show()

print "Computing 6Ms MS Kasen model"
M_comp = 6.0/1.4 #Mchandra
a13 = 0.2 #10^13cm
#plot comparison
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("Flux [uJy]")
for i in range(len(t)):
    #fit early light curve
    #plot I band fluxes
    thetas = np.linspace(0,180,10)
    for theta in thetas:
        #angle corrected Kasen luminosity
        Lk, Lk_theta, Tk = Kasen2010(tk, a13, theta, m_c, e_51)
        #observer frame angle corrected Kasen flux
        Fk, Fk_err = KasenFit(tk,Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
        ax[i].plot(tk, Fk)
    #plot red corrected observer frame observed fluxes
    ax[i].errorbar(t[i],F[i],F_err[i],fmt='g+')
    ax[i].scatter(tlim[i], Flim[i], marker='v', c='r')
plt.show()

print "Computing 2Ms MS Kasen model"
M_comp = 2.0/1.4 #Mchandra
a13 = 0.05 #10^13cm
#plot comparison
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("Flux [uJy]")
for i in range(len(t)):
    #fit early light curve
    #plot I band fluxes
    thetas = np.linspace(0,180,10)
    for theta in thetas:
        #angle corrected Kasen luminosity
        Lk, Lk_theta, Tk = Kasen2010(tk, a13, theta, m_c, e_51)
        #observer frame angle corrected Kasen flux
        Fk, Fk_err = KasenFit(tk,Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
        ax[i].plot(tk, Fk)
    #plot red corrected observer frame observed fluxes
    ax[i].errorbar(t[i],F[i],F_err[i],fmt='g+')
    ax[i].scatter(tlim[i], Flim[i], marker='v', c='r')
plt.show()
