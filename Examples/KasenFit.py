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
from scipy.stats import norm
from snpy import *

#essential files
from SNAP.Analysis.LCRoutines import*
from SNAP.Analysis.LCFitting import*
from SNAP.Analysis.Cosmology import*

plot = False #plot polynomial fits to light curves

#Epoch
t0 = -17.66
t0err = 0.89
#redshift of N300-1.Q0.SN
z = 0.057
zerr = 0.005
#Extinction coefficient (galactic) in each band S & F (2011)
EBVgal = 0.107
Coefs = np.array([3.641, 2.682, 1.516])
band = ['B','V','i']
#max time
Tmax = 281.982
#explosion parameters
m_c = 1.24/1.4 #Mchandra
e_51 = 0.90 #x10^51 ergs
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
    #get epoch centered dilated time
    t[i] = absTime(t[i], z) - t0

    #get limiting magnitudes
    Mlim[i] = deredMag(Mlim[i], EBVgal, Coefs[i])
    Flim.append(Mag_toFlux(band[i], Mlim[i])*10**6)
    #tlim.append(t[i][SN[i]<2.0])
    #Flim[i] = Flim[i][SN[i]<2.0]

    #get reliable detections
    #t[i] = t[i][SN[i]>=2.0]
    #F[i] = F[i][SN[i]>=2.0]
    #F_err[i] = F_err[i][SN[i]>=2.0]
    
print "plotting early data"
#plot
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("L/Lmax")
for i in range(len(t)):
    #fit early light curve
    ax[i].errorbar(t[i],F[i],F_err[i],fmt='g+')
    ax[i].scatter(t[i], Flim[i], marker='v', c='r')
    ax[i].plot(t[i],[0]*len(t[i]),'k')
plt.show()

print "Trial thetas:", np.linspace(0,180,10)
tk = np.arange(0.001, 10, 0.001)
Fks = [[],[],[]]
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
        Lk, Tk = Kasen2010(tk, a13, m_c, e_51)
        Lk_theta = Lk*Kasen_isocorr(theta)
        #rest time observer frame angle corrected Kasen flux
        Fk, Fk_err = KasenFit(Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
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
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("Flux [uJy]")
for i in range(len(t)):
    #fit early light curve
    #plot I band fluxes
    thetas = np.linspace(0,180,10)
    for theta in thetas:
        #angle corrected Kasen luminosity
        Lk, Tk = Kasen2010(tk, a13, m_c, e_51)
        Lk_theta = Lk*Kasen_isocorr(theta)
        #observer frame angle corrected Kasen flux
        Fk, Fk_err = KasenFit(Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
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
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("Flux [uJy]")
for i in range(len(t)):
    #fit early light curve
    #plot I band fluxes
    thetas = np.linspace(0,180,10)
    for theta in thetas:
        #angle corrected Kasen luminosity
        Lk, Tk = Kasen2010(tk, a13, m_c, e_51)
        Lk_theta = Lk*Kasen_isocorr(theta)
        #observer frame angle corrected Kasen flux
        Fk, Fk_err = KasenFit(Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
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
ax[-1].set_xlabel("Time since epoch (Days)", fontsize=14)
#for each band
for i in range(len(t)):
    #plot light curve
    Fp, Fp_err = F[i][F[i]>0], F_err[i][F[i]>0]
    tp = t[i][F[i]>0]
    M, M_err = Flux_toMag(Band[i],Fp*1e-6,Fp_err[i]*1e-6)
    Mlim = Flux_toMag(Band[i],Flim[i]*1e-6)
    ax[i].errorbar(tp,M,M_err,fmt='g+')
    ax[i].scatter(t[i], Mlim, marker='v', c='r')
    ax[i].set_ylabel(Band[i], fontsize=14, fontstyle='italic', fontweight='bold')
    ax[i].set_ylim(25,18)
    ax[i].yaxis.set_ticks([24,22,20])
    #plot each zero theta model
    for Fk in Fks[i]:
        ts, Fs = tk[Fk>0], Fk[Fk>0]
        Ms = Flux_toMag(Band[i],Fs*1e-6)
        ax[i].plot(ts, Ms)
f.subplots_adjust(hspace=0)
plt.show()


print "Computing viewing angles at each confidence interval"
#list of confidence intervals to traverse
confs = np.linspace(65,99.99,200)
#list of sample models
a13s = [2.0, 0.2, 0.05] #1RG, 6MS, 2MS
#list of viewing angles
thetas = np.linspace(0,180,100)

Lk_terr = []
for i in range(len(t)):
    tk = t[i][t[i]>0]
    #perturb time axis by error in epoch
    dts = np.random.normal(0.0, t0err, 100)
    Lk_dif = np.zeros([len(dts), len(tk)])
    for r, dt in enumerate(dts):
        Lk_dif[r] = Kasen2010(tk+dt, a13, m_c, e_51)[0]
    Lk_terr.append(Lk_dif.std(axis=0))
    
#for each sample model
for a13 in a13s:
    #array to hold percent of viewing angles ruled out at each conf
    outangles = np.zeros(len(confs))
    #for each confidence interval
    for j, conf in enumerate(confs):
        #sigma needed to establish confidence below LC
        sig = norm.ppf(conf/100.0)
        print conf
        #boolean mask for whether angle is ruled out to given confidence
        mask = np.array([True]*len(thetas))
        #for each band
        for i in range(len(t)):
            tk = t[i][t[i]>0]
            Fobs = F[i][t[i]>0]
            Ferr = F_err[i][t[i]>0]
            #get theoretical light curve
            Lk, Tk = Kasen2010(tk, a13, m_c, e_51)
            #for each angle
            for k, theta in enumerate(thetas):
                #angle corrected Kasen luminosity
                Lk_theta = Lk*Kasen_isocorr(theta)
                #observer frame angle corrected Kasen flux
                Fk, Fk_zerr = KasenFit(Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
                #epoch error
                Lk_theta_terr = Lk_terr[i]*Kasen_isocorr(theta)
                Fk_terr = KasenFit(Lk_theta_terr,Tk,wave_0[bands[band[i]]],z,zerr)[0]
                #total error
                Fk_err = np.square(Fk_terr)+np.square(Fk_zerr)
                Fk_err = np.array([np.sqrt(f) for f in Fk_err])
                #total error
                Err = np.square(Ferr)+np.square(Fk_err)
                Err = np.array([np.sqrt(E) for E in Err])
                #check if any points rule out angle with conf
                if any(Fk > Fobs + sig*Err):
                    mask[k] = False
                else:
                    print "Consistent!", conf, theta
        #At this confidence level, we rule out some percent of angles
        outangles[j] = 100.0*float(len(thetas)-len(thetas[mask]))/len(thetas)
    #At this a13, we can plot ruled out angles vs confidence
    print confs, outangles
    plt.plot(confs, outangles)
    plt.xlim(0,100)
    plt.ylim(0,101)
plt.show()


print "Computing viewing angles at each separation distance"
#list of sample models
a13s = np.arange(0.05,5.0,0.05) #1RG, 6MS, 2MS
confs = [68.0, 95.0, 99.7]
#list of viewing angles
thetas = np.linspace(0,180,100)

Lk_terr = []
for i in range(len(t)):
    tk = t[i][t[i]>0]
    #perturb time axis by error in epoch
    dts = np.random.normal(0.0, t0err, 100)
    Lk_dif = np.zeros([len(dts), len(tk)])
    for r, dt in enumerate(dts):
        Lk_dif[r] = Kasen2010(tk+dt, a13, m_c, e_51)[0]
    Lk_terr.append(Lk_dif.std(axis=0))
    
#for each confidence interval
for conf in confs:
    #array to hold percent of viewing angles ruled out at each conf
    outangles = np.zeros(len(a13s))
    #sigma needed to establish confidence below LC
    sig = norm.ppf(conf/100.0)
    #for each sample model
    for j, a13 in enumerate(a13s):
        print a13
        #boolean mask for whether angle is ruled out to given confidence
        mask = np.array([True]*len(thetas))
        #for each band
        for i in range(len(t)):
            tk = t[i][t[i]>0]
            Fobs = F[i][t[i]>0]
            Ferr = F_err[i][t[i]>0]
            #get theoretical light curve
            Lk, Tk = Kasen2010(tk, a13, m_c, e_51)
            #for each angle
            for k, theta in enumerate(thetas):
                #angle corrected Kasen luminosity
                Lk_theta = Lk*Kasen_isocorr(theta)
                #observer frame angle corrected Kasen flux
                Fk, Fk_zerr = KasenFit(Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
                #epoch error
                Lk_theta_terr = Lk_terr[i]*Kasen_isocorr(theta)
                Fk_terr = KasenFit(Lk_theta_terr,Tk,wave_0[bands[band[i]]],z,zerr)[0]
                #total error
                Fk_err = np.square(Fk_terr)+np.square(Fk_zerr)
                Fk_err = np.array([np.sqrt(f) for f in Fk_err])
                #total error
                Err = np.square(Ferr)+np.square(Fk_err)
                Err = np.array([np.sqrt(E) for E in Err])
                #check if any points rule out angle with conf
                if any(Fk > Fobs + sig*Err):
                    mask[k] = False
                else:
                    print "Consistent!", a13, conf, theta
        #At this confidence level, we rule out some percent of angles
        outangles[j] = 100.0*float(len(thetas)-len(thetas[mask]))/len(thetas)
    #At this conf, we can plot ruled out angles vs a13
    print conf, outangles
    plt.plot(a13s, outangles)
    plt.xlim(0,5)
    plt.ylim(0,101)
plt.show()

"""
#Load Wang distribution
Masses, Probs = np.loadtxt("Wang2010.txt", unpack=True)
print Masses, Probs.sum()
#For each mass, compute a13
Radii = []
for mass in Masses:
    if mass > 1.0:
        Radii.append(6.963*10**10*np.power(mass,0.57))
    else:
        Radii.append(6.963*10**10*np.power(mass,0.8))
Radii = np.array(Radii)
a13s = 5.0*Radii/10**13
print a13s

Probtheta = []
for j in range(len(Masses)):
    a13 = a13s[j]
    print "Computing model:",Masses[j], a13
    thetas = np.linspace(0,180,100)
    mask = np.array([True]*len(thetas))
    for i in range(len(t)):
        #fit early light curve
        for k, theta in enumerate(thetas):
            #angle corrected Kasen luminosity
            tk = t[i]
            Lk, Lk_theta, Tk = Kasen2010(tk, a13, theta, m_c, e_51)
            #observer frame angle corrected Kasen flux
            Fk, Fk_err = KasenFit(tk,Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
            #light curve constraint
            Fk[0] = 0
            if any(Fk>F[i]+F_err[i]):
                mask[k] = False
                
            tk = t[i]
            Lk, Lk_theta, Tk = Kasen2010(tk, a13, theta, m_c, e_51)
            #observer frame angle corrected Kasen flux
            Fk, Fk_err = KasenFit(tk,Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
            #light curve constraint
            Fk[0] = 0
            if any(Fk>Flim[i]):
                mask[k] = False
    Probtheta.append(float(len(thetas[mask]))/len(thetas))
Probtheta = np.array(Probtheta)
Probfinal = Probtheta*Probs
print Probtheta
print Probfinal.sum()


index = np.arange(len(Masses))
bar_width = 1.0

f, ax = plt.subplots()
ax.bar(index, Probfinal, bar_width, alpha=0.5, edgecolor='k', color='b')
ax.xaxis.set_ticks(index[:10]*2)
ax.set_xticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.6, 1.8, 2.0])
ax.set_yticks([0.0, 0.1, 0.2])
ax.set_xlabel("$M_{2} (M_{\odot})$", fontsize=14)
ax.set_ylabel("Observing Probability", fontsize=14)
ax.set_xlim(0.0,21)
ax.set_ylim(0.0,0.25)
plt.tight_layout()
plt.show()
"""
