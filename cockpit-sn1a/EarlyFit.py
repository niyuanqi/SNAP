#################################################################
# Name:     EarlyFit.py                                         #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2017                                       #
# Function: Program computes early flux light curve and fits    #
#           for explosion epoch. Update ObjData.py first.       #
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
from ObjData import *

plot = False #plot polynomial fits to light curves

print "Loading binned early light curve."
#get N300-1.Q0.SN binned light curve
t, M, M_err, F, SN, Mlim = LCload(binfiles, tcol=0, magcols=7, errcols=8, fluxcols=5, SNcols=6, limcols=9, SNthres=-10.0, scols=10, flags=['-99.99999'], mode='multi')

#crop window in data
for i in range(len(M)):
    t[i], M[i], M_err[i], F[i], SN[i], Mlim[i] = LCcrop(t[i], t1_early, t2_early, M[i], M_err[i], F[i], SN[i], Mlim[i])

#get noise in flux
F_err = [F[i]/SN[i] for i in range(len(t))]

#Where theres no SN, estimate background level.
bg = np.zeros(len(M))
bg_err = np.zeros(len(M))
for i in range(len(M)):
    mask = t[i] < 418.0
    w = 1/np.square(F_err[i][mask])
    bg[i] = np.sum(F[i][mask]*w)/np.sum(w)
    bg_err[i] = np.sqrt(1/np.sum(w))

    #Warning, you can't do that because data number =/= jansky!
    #F[i] = F[i] - bg[i]
    #F_err[i] = np.sqrt(np.square(F_err[i])+np.square(bg_err[i]))

print "Loading reliable total light curve"
#get absolute fluxes
s = get_sn(sn_file)
#don't plot fit
s.replot = 0
#fit model
s.choose_model("EBV_model2", stype="st")
print "Performing SNpy fit and conversion of LC to rest frame"
s.fit(band)
print "Tmax ", s.Tmax, s.e_Tmax
Flim = Mlim
for i in range(len(F)):
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
    #early light curve
    ax[i].errorbar(t[i],L[i],L_err[i],fmt='g+')
    #ax[i].scatter(t[i],Llim[i],c='r',marker='v')
plt.show()

print "Fitting power law to section of early light curve"
#for each band crop to right section
f = 0.4 #fit up to 40% of maximum flux (Olling 2015)
t = [time[L[i]<f] for i, time in enumerate(t)]
L_err = [lerr[L[i]<f] for i, lerr in enumerate(L_err)]
L = [l[l<f] for l in L]

#fit for early light curve using leastsq
p0 = [-16.7120900068, 0.0057208154898343525, 0.0074170018284233702, 0.019945558712668135, 2.0016550276367129, 1.9739100509713501, 1.5029759166486247]

#for each band, perturb flux by flux errors
n = 100000 #number of perturbations
randomdataY = [L]
for j in range(n):
    print str(j+1)+"/"+str(n)
    L_pert = []
    for i in range(len(t)):
        L_pert.append(L[i] + np.random.normal(0., L_err[i], len(L[i])))
    randomdataY.append(L_pert)
print "Fitting bootstrap using 4 processes"
x, err_x = fit_bootstrap(p0, t, randomdataY, L_err, earlyMultiErr, errfunc=True, perturb=False, n=3000, nproc=4)
print "Fitting bootstrap using better initial parameters"
x, err_x = fit_bootstrap(x, t, randomdataY, L_err, earlyMultiErr, errfunc=True, perturb=False, n=3000, nproc=4)

#interpret results
t0, t0_err = x[0], err_x[0]
C = [x[1],x[2],x[3]]
C_err = [err_x[1],err_x[2],err_x[3]]
a = [x[4],x[5],x[6]]
a_err = [err_x[4],err_x[5],err_x[6]]
x2dof = np.sqrt(np.square(earlyMultiErr(x, t, L, L_err)).sum())/(len(L[0])+len(L[1])+len(L[2])-len(x))
#output data
print ""
print "Epoch of first light in rest frame:", t0, t0_err
print "Epoch of first light in obs frame:", t0*(1.0+z), np.absolute(t0*(1.0+z)*np.sqrt(np.square(zerr/(1.0+z))+np.square(t0_err/t0)))
print "Coefficient", C, C_err
print "Power:", a, a_err
print "Fit Chi2/dof", x2dof
print ""

print "Plotting early fit"
#plot
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("Days from peak", fontsize = 14)
ax[1].set_ylabel("Normalized Flux", fontsize = 14)
#best fit curves
tT = np.linspace(t[0][0]-10, t[0][-1]+10, 1000)
LT = [earlyFit(tT,t0,C[0],a[0]),earlyFit(tT,t0,C[1],a[1]),earlyFit(tT,t0,C[2],a[2])]

#plot fit
for i in range(len(t)):
    #plot fluxes
    ax[i].errorbar(t[i],L[i],2*L_err[i],fmt='k+')
    ax[i].plot(tT, LT[i], 'k-')
    ax[i].set_xlim(t1_early-s.Tmax, t2_early-s.Tmax)
    ax[i].set_ylim(-0.1,0.7)
    ax[i].plot([t0,t0], [-10,30], 'k:', linewidth=1)
ax[0].text(-18,0.5,'B', fontsize = 14, fontstyle='italic', fontweight='bold')
ax[1].text(-18,0.5,'V', fontsize = 14, fontstyle='italic', fontweight='bold')
ax[2].text(-18,0.5,'I', fontsize = 14, fontstyle='italic', fontweight='bold')
f.subplots_adjust(hspace=0)
plt.tight_layout()
plt.show()

#plot residuals
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("Days from peak", fontsize = 14)
ax[1].set_ylabel("Normalized Flux", fontsize = 14)
for i in range(len(t)):
    #plot residuals
    ax[i].errorbar(t[i],L[i]-earlyFit(t[i],t0,C[i],a[i]),2*L_err[i],fmt='g+')
    #ax[i].scatter(t[i],Llim[i]-earlyFit(t[i],*popts[i]),c='r',marker='v')
    ax[i].plot(tT,[0]*len(tT),label='power fit')
    ax[i].set_xlim(t1_early-s.Tmax, t2_early-s.Tmax)
    ax[i].set_ylim(-0.095,0.095)
ax[0].text(-18,0.05,'B', fontsize = 14, fontstyle='italic', fontweight='bold')
ax[1].text(-18,0.05,'V', fontsize = 14, fontstyle='italic', fontweight='bold')
ax[2].text(-18,0.05,'I', fontsize = 14, fontstyle='italic', fontweight='bold')
f.subplots_adjust(hspace=0)
plt.tight_layout()
plt.show()
