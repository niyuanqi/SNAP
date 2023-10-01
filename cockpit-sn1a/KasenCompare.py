#################################################################
# Name:     KasenCompare.py                                     #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2017                                       #
# Function: Program compares early light curve to Kasen model.  #
#           Update ObjData.py first.                            #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

#essential files
from SNAP.Analysis.LCRoutines import*
from SNAP.Analysis.LCFitting import*
from SNAP.Analysis.Cosmology import*
from ObjData import *

#Ejecta mass in chandrasekhar masses
m_c = m_ej/1.4
m_c_err = m_ej_err/1.4
plot = False #plot polynomial fits to light curves

print "Loading binned early light curve."
#get binned light curve
t, M, M_err, F, SN, Mlim = LCload(binfiles, tcol=0, magcols=7, errcols=8, fluxcols=5, SNcols=6, limcols=9, SNthres=-10.0, scols=10, flags=['-99.99999'], mode='multi')
#get noise in flux
F_err = [F[i]/SN[i] for i in range(3)]

#crop window in data
for i in range(len(M)):
    t[i], M[i], M_err[i], F[i], SN[i], Mlim[i] = LCcrop(t[i], t1_early, t2_early, M[i], M_err[i], F[i], SN[i], Mlim[i])

#deredden flux and shift light curves
print "Correcting for galactic reddening"
Flim = []
for i in range(len(F)):
    #correct fluxes for galactic reddening
    F[i] = deredFlux(F[i], EBVgal, Coefs[i])
    F_err[i] = deredFlux(F_err[i], EBVgal, Coefs[i])
    #shift to zero centered time
    t[i] = t[i] - Tmax
    #get explosion centered time in observer frame
    t[i] = t[i] - t0obs

    #get limiting magnitudes
    Mlim[i] = deredMag(Mlim[i], EBVgal, Coefs[i])
    Flim.append(Mag_toFlux(band[i], Mlim[i])*10**6)

###################################################################

print "Trial thetas:", np.linspace(0,180,10)
tk = np.arange(0.001, 7, 0.001)
Fks = [[],[],[]]

print "Computing 1Ms RG Kasen model"
M_comp = 1.0/1.4 #Mchandra
a13 = 2.0 #10^13cm
#plot comparison
f, ax = plt.subplots(3, sharex=True)
ax[0].set_title("Comparison to 1RG model")
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("Flux [uJy]")
for i in range(len(t)):
    #plot I band fluxes
    thetas = np.linspace(0,180,10)
    Lk = np.array([KasenFit(ti, a13, 1.0, wave_0[bands[band[i]]], m_c, e_51, z, 0) for ti in tk])
    for theta in thetas:
        #rest time angle corrected Kasen luminosity
        Fk = Lk*Kasen_isocorr(theta)
        ax[i].plot(tk, Fk)
        if theta == thetas[0]:
            Fks[i].append(Fk)
    #plot red corrected observer frame rest time observed fluxes
    ax[i].errorbar(t[i],F[i],F_err[i],fmt='g+')
    ax[i].scatter(t[i], Flim[i], marker='v', c='r')
plt.show()

print "Computing 6Ms MS Kasen model"
M_comp = 6.0/1.4 #Mchandra
a13 = 0.2 #10^13cm
#plot comparison
f, ax = plt.subplots(3, sharex=True)
ax[0].set_title("Comparison to 6MS model")
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("Flux [uJy]")
for i in range(len(t)):
    #fit early light curve
    #plot I band fluxes
    thetas = np.linspace(0,180,10)
    Lk = np.array([KasenFit(ti, a13, 1.0, wave_0[bands[band[i]]], m_c, e_51, z, 0) for ti in tk])
    for theta in thetas:
        #angle corrected Kasen luminosity
        Fk = Lk*Kasen_isocorr(theta)
        ax[i].plot(tk, Fk)
        if theta == thetas[0]:
            Fks[i].append(Fk)
    #plot red corrected observer frame observed fluxes
    ax[i].errorbar(t[i],F[i],F_err[i],fmt='g+')
    ax[i].scatter(t[i], Flim[i], marker='v', c='r')
plt.show()

print "Computing 2Ms MS Kasen model"
M_comp = 2.0/1.4 #Mchandra
a13 = 0.05 #10^13cm
#plot comparison
f, ax = plt.subplots(3, sharex=True)
ax[0].set_title("Comparison to 2MS model")
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("Flux [uJy]")
for i in range(len(t)):
    #fit early light curve
    #plot I band fluxes
    thetas = np.linspace(0,180,10)
    Lk = np.array([KasenFit(ti, a13, 1.0, wave_0[bands[band[i]]], m_c, e_51, z, 0) for ti in tk])
    for theta in thetas:
        #angle corrected Kasen luminosity
        Fk = Lk*Kasen_isocorr(theta)
        ax[i].plot(tk, Fk)
        if theta == thetas[0]:
            Fks[i].append(Fk)
    #plot red corrected observer frame observed fluxes
    ax[i].errorbar(t[i],F[i],F_err[i],fmt='g+')
    ax[i].scatter(t[i], Flim[i], marker='v', c='r')
plt.show()

print "Plotting zero theta models"
#Make plot of zero theta models over light curve
f, ax = plt.subplots(len(t), sharex=True) 
ax[-1].set_xlabel("Days since first light", fontsize=16)
#for each band
for i in range(len(t)):
    #plot light curve
    Mlim = Flux_toMag(Band[i],Flim[i]*1e-6)
    
    detmask = F[i]>Flim[i]
    Fp, Fp_err = F[i][detmask], F_err[i][detmask]
    tp = t[i][detmask]
    M, M_err = Flux_toMag(Band[i],Fp*1e-6,Fp_err*1e-6)
    
    nonmask = F[i]<Flim[i]
    Mlimn = Mlim[nonmask]
    tn = t[i][nonmask]
    
    ax[i].errorbar(tp,M, yerr=2*M_err,fmt='k+')
    ax[i].scatter(tn, Mlimn, marker='v', c='k')
    ax[i].set_ylabel(Band[i], fontsize=16, fontstyle='italic', fontweight='bold')
    ax[i].set_ylim(25.0,19.0)
    ax[i].yaxis.set_ticks([24,22,20])
    ax[i].tick_params(labelsize=14)
    ax[i].set_xlim(t1_early - Tmax - t0obs, t2_early - Tmax - t0obs)
    #ax[i].scatter(0, 19.9, marker='v', c='k')
    #ax[i].plot([0,0], [18.5, 20], 'k')
    #ax[i].plot([t0obserr,t0obserr], [18.5, 20], 'k', ls=':')
    #ax[i].plot([-t0obserr,-t0obserr], [18.5, 20], 'k', ls=':')
    #plot each zero theta model
    for Fk in Fks[i]:
        ts, Fs = tk[Fk>0], Fk[Fk>0]
        Ms = Flux_toMag(Band[i],Fs*1e-6)
        ax[i].plot(ts, Ms)
f.subplots_adjust(hspace=0)
plt.show()

