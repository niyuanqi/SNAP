#################################################################
# Name:     KasenF.py                                           #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 27, 2017                                       #
# Function: Program computes kasen models at various separation #
#           distances and viewing angles, comparing them to     #
#           early light curve. Uses flux errors to estimate     #
#           confidence interval of analysis.                    #
#           Update ObjData.py with the right files first.       #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from multiprocessing import Pool

#essential files
from SNAP.Analysis.LCRoutines import*
from SNAP.Analysis.LCFitting import*
from SNAP.Analysis.Cosmology import*
from ObjData import *

#ejecta mass in chandrasekhar masses
m_c = m_ni/1.4 #*1.4 Mchandra
m_c_err = m_ni_err/1.4

print "loading binned data"
#load binned early LC
t, M, M_err, F, SN, Mlim = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, scols=9, flags=['-99.99999'], mode='multi')
#get noise in flux
F_err = [F[i]/SN[i] for i in range(len(t))]
#crop to relevant window
for i in range(len(t)):
    t[i], F[i], F_err[i], Mlim[i] = LCcrop(t[i], t1_early, t2_early, F[i], F_err[i], Mlim[i])
#limiting magnitudes at different SNR
limconf = []
#load each set of limiting magnitudes
for i in range(len(limSNs)):
    tc, Mc, Mc_err, Fc, SNc, Mlimc = LCload(conffiles.T[i], tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, scols=9, flags=['-99.99999'], mode='multi')
    limconf.append(Mlim)

#deredden flux and shift light curve
print "Correcting for galactic reddening"
flimconf = np.copy(limconf)
for i in range(len(t)):
    #correct fluxes for galactic reddening
    F[i] = deredFlux(F[i], EBVgal, Coefs[i])
    F_err[i] = deredFlux(F_err[i], EBVgal, Coefs[i])
    #shift to zero centered time
    t[i] = t[i] - Tmax
    #get explosion centered time in observer frame
    t[i] = t[i] - t0obs

    #get limiting magnitudes into flux
    for n in range(len(limconf)):
        limconf[n][i] = deredMag(limconf[n][i], EBVgal, Coefs[i])
        flimconf[n][i] = Mag_toFlux(band[i], limconf[n][i])*10**6

###################################################################

print "Computing viewing angles at each separation distance"
#list of sample a13 models
a13s = np.concatenate((np.arange(0.001,0.05,0.001), np.arange(0.05,0.2,0.01), np.arange(0.2, 2.0, 0.01), np.arange(2.0,10.1,0.1)))
#list of viewing angles
thetas = np.linspace(0,180,100)
#SNR of 1, 2, 3, 4, 5
confs = [norm.cdf(sn) for sn in limSNs]
#[0.84134474606854293, 0.97724986805182079, 0.9986501019683699, 0.99996832875816688, 0.99999971334842808]
print "Confidence levels:",confs
print "Sigma levels:",limSNs
print "Trial a13s:", a13s
print "Trial thetas:",thetas

#for each confidence interval
for n, conf in enumerate(confs):
    #array to hold percent of viewing angles ruled out for each a13 at this conf
    outangles = np.zeros(len(a13s))
    #sigma needed to establish confidence below LC
    sig = norm.ppf(conf/100.0)
    #for each a13 model
    for j, a13 in enumerate(a13s):
        print "Testing model at a13:",a13
        #boolean mask for whether angle is ruled out to given confidence
        mask = np.array([True]*len(thetas))
        #for each band
        for i in range(len(t)):
            #generate theoretical light curve
            Lk = np.array([KasenFit(ti, a13, 1.0, wave_0[bands[band[i]]], m_c, e_51, z, 0) for ti in t[i]])

            #compare to observed for each viewing angle
            for k, theta in enumerate(thetas):
                #angle corrected Kasen luminosity
                Fk = Lk*Kasen_isocorr(theta)
                #check if any points rule out angle with conf
                level = F[i] + sig*F_err[i]
                #if not above limit at this conf, replace with limit
                level[level < flimconf[n][i]] = flimconf[n][i][level < flimconf[n][i]]
                if any(Fk > level):
                    #ruled out
                    mask[k] = False
                else:
                    #not ruled out
                    print "Consistent!", a13, conf, theta

        """
        #Diagnostic
        if n == 0 and a13 == 0.05 and False:
            print "plotting section"
            #plot section
            f, ax = plt.subplots(len(t), sharex=True) 
            for i in range(len(t)):
                Lk = np.array([KasenFit(ti, a13, 1.0, wave_0[bands[band[i]]], m_c, e_51, z, 0) for ti in t[i]])
                Fk = Lk*Kasen_isocorr(0)
                ax[i].errorbar(t[i], F[i], yerr=sig*F_err[i], fmt="k+")
                ax[i].scatter(t[i], Fk)
                ax[i].scatter(t[i], flimconf[n][i], color='r', marker='v')
            plt.subplots_adjust(hspace=None)
            plt.show() 
        """
                
        #At this confidence level, we rule out some angles at each a13
        outangles[j] = 180.0*float(len(thetas)-len(thetas[mask]))/len(thetas)
    #At this conf, we can plot ruled out angles vs a13
    print "Confidence level:",conf
    print "Angles ruled out at each a13:",outangles
    plt.plot(a13s, outangles, style[n])
plt.xlim(0,10.0)
plt.ylim(0,185)
#plot positions of 1RG, 6MS, 2MS
plt.plot([0.05,0.05], [175,185], 'k:', linewidth=1)
plt.plot([0.2,0.2], [175,185], 'k:', linewidth=1)
plt.plot([2.0,2.0], [175,185], 'k:', linewidth=1)
plt.ylabel("Lower Limit of Acceptable \nViewing Angles (Degrees)", fontsize=16)
plt.xlabel("Separation Distance ($10^{13}$ cm)", fontsize=16)
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.show()
