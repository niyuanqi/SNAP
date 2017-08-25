#################################################################
# Name:     Epoch.py                                            #
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

#max magnitudes
maxes = np.array([-18.8394,-19.0216,-18.2893])
#redshift of N300-1.Q0.SN
z = 0.057
zerr = 0.003
#Extinction coefficient (galactic) in each band S & F (2011)
EBVgal = 0.107
Coefs = np.array([3.641, 2.682, 1.516])
band = ['B','V','i']

#N300-1.Q0.SN binned time series data files
binBfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcbin.CN_170804.txt"
binVfile = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lcbin.CN_170804.txt"
binIfile = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lcbin.CN_170804.txt"
binfiles = [binBfile, binVfile, binIfile] 
print "Loading binned early light curve."
#get N300-1.Q0.SN binned light curve
t, M, M_err, F, SN, Mlim = LCload(binfiles, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, scols=9, flags=['-99.99999'], mode='multi')
#get noise in flux
F_err = [F[i]/SN[i] for i in range(3)]

print "Loading reliable total light curve"
#get absolute fluxes
s = get_sn("N300-1.Q0.SN.txt")
#don't plot fit
s.replot = 0
#fit model
s.choose_model("EBV_model2", stype="st")
print "Performing SNpy fit and conversion of LC to rest frame"
s.fit(band)
Flim = Mlim
for i in range(len(F)):
    #correct fluxes for galactic reddening
    F[i] = deredFlux(F[i], EBVgal, Coefs[i])
    F_err[i] = deredFlux(F_err[i], EBVgal, Coefs[i])
    #need to correct for K correction
    kcorr, mask = s.model.kcorr(band[i],t[i]-s.Tmax)
    kcorr_val = kcorr[mask]
    kcorr_inval = np.zeros(len(kcorr[np.invert(mask)]))
    kcorr = np.concatenate([kcorr_inval, kcorr_val], axis=0)
    F[i], F_err[i] = absFlux(F[i]*10**-6, z, appFlux_err=F_err[i]*10**-6, z_err=zerr, Kcorr=kcorr)
    F[i], F_err[i] = F[i]*10**6, F_err[i]*10**6
    Flim[i] = Mag_toFlux(band[i], absMag(deredMag(Mlim[i], EBVgal, Coefs[i]), z, Kcorr=kcorr))*10**6 #in uJy
    #peak abs fluxes
    maxes[i] = Mag_toFlux(band[i],maxes[i])*10**6 #in uJy
    #get Max centered dilated time
    t[i] = absTime(t[i]-s.Tmax, z)

print "Converting early LC to rest frame normalized luminosity"
#scale luminosity as fraction of max lum
L_err = [l/maxes[i] for i, l in enumerate(F_err)]
L = [l/maxes[i] for i, l in enumerate(F)]
Llim = [l/maxes[i] for i, l in enumerate(Flim)]

print "plotting early data"
#plot
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("L/Lmax")
for i in range(len(t)):
    #fit early light curve
    #plot I band fluxes
    ax[i].errorbar(t[i],L[i],L_err[i],fmt='g+')
    #ax[i].scatter(t[i],Llim[i],c='r',marker='v')
    #ax[i].set_xlim(-20.0,-5.0)
    #ax[i].set_ylim(-0.1,1)
plt.show()


print "Fitting power law to section of early light curve"
#for each band crop to right section
f = 0.52
t = [time[L[i]<f] for i, time in enumerate(t)]
L_err = [lerr[L[i]<f] for i, lerr in enumerate(L_err)]
L = [l[l<f] for l in L]

#crop out early light section for reliable power indices
#for i in range(len(t)):
#    mask = np.logical_or(t[i]<-17.0,t[i]>-15.0)
#    t[i] = t[i][mask]
#    L_err[i] = L_err[i][mask]
#    L[i] = L[i][mask]
    
t1 = -20
t2 = -8
for i in range(len(t)):
    t[i], L[i], L_err[i] = LCcrop(t[i], t1, t2, L[i], L_err[i])

#plot
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("Days from peak", fontsize = 14)
ax[1].set_ylabel("Normalized Flux", fontsize = 14)
tT = np.linspace(t[0][0]-10, t[0][-1]+10, 1000)
#fit for early light curve using leastsq
p0 = [-17.8, 0.003,0.003,0.004, 2.2,2.3,2.3]

#for each band, perturb flux by flux errors
n = 100000
randomdataY = [L]
for j in range(n):
    L_pert = []
    for i in range(len(t)):
        L_pert.append(L[i] + np.random.normal(0., L_err[i], len(L[i])))
    randomdataY.append(L_pert)
x, err_x = fit_bootstrap(p0, t, randomdataY, L_err, earlyMultiErr, errfunc=True, perturb=False, n=3000, nproc=4)

#x, cov_x, infodict, mesg, ier = leastsq(earlyMultiErr, p0, args=(t, L, L_err), full_output=1)
#err_x = np.sqrt(np.diag(cov_x))

t0, t0_err = x[0], err_x[0]
C = [x[1],x[2],x[3]]
C_err = [err_x[1],err_x[2],err_x[3]]
a = [x[4],x[5],x[6]]
a_err = [err_x[4],err_x[5],err_x[6]]
LT = [earlyFit(tT,t0,C[0],a[0]),earlyFit(tT,t0,C[1],a[1]),earlyFit(tT,t0,C[2],a[2])]
x2dof = np.sqrt(np.square(earlyMultiErr(x, t, L, L_err)).sum())/(len(L[0])+len(L[1])+len(L[2])-len(x))

print "Epoch:", t0, t0_err
print "Coefficient", C, C_err
print "Power:", a, a_err
print "Fit Chi2", x2dof

#plot fit
for i in range(len(t)):
    #plot fluxes
    ax[i].errorbar(t[i],L[i],L_err[i],fmt='k+')
    ax[i].plot(tT, LT[i], 'k-')
    ax[i].set_xlim(-19.0,-7.0)
    ax[i].set_ylim(-0.1,0.8)
ax[0].text(-18.5,0.6,'B', fontsize = 14, fontstyle='italic', fontweight='bold')
ax[1].text(-18.5,0.6,'V', fontsize = 14, fontstyle='italic', fontweight='bold')
ax[2].text(-18.5,0.6,'I', fontsize = 14, fontstyle='italic', fontweight='bold')
f.subplots_adjust(hspace=0)
plt.tight_layout()
plt.show()

#plot residuals
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("Days from peak", fontsize = 14)
ax[1].set_ylabel("Normalized Flux", fontsize = 14)
for i in range(len(t)):
    #plot residuals
    ax[i].errorbar(t[i],L[i]-earlyFit(t[i],t0,C[i],a[i]),L_err[i],fmt='g+')
    #ax[i].scatter(t[i],Llim[i]-earlyFit(t[i],*popts[i]),c='r',marker='v')
    ax[i].plot(tT,[0]*len(tT),label='power fit')
    ax[i].set_xlim(-19.0,-7.0)
    ax[i].set_ylim(-0.09,0.09)
ax[0].text(-18.5,0.05,'B', fontsize = 14, fontstyle='italic', fontweight='bold')
ax[1].text(-18.5,0.05,'V', fontsize = 14, fontstyle='italic', fontweight='bold')
ax[2].text(-18.5,0.05,'I', fontsize = 14, fontstyle='italic', fontweight='bold')
f.subplots_adjust(hspace=0)
plt.tight_layout()
plt.show()
