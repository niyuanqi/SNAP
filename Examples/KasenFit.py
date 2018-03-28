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
from multiprocessing import Pool
import matplotlib.pyplot as plt

#essential files
from SNAP.Analysis.LCRoutines import*
from SNAP.Analysis.LCFitting import*
from SNAP.Analysis.Cosmology import*

plot = False #plot polynomial fits to light curves

#Epoch
t0 = -19.09 #in observer frame
t0err = 0.67 #in rest frame
#t0 = t0 + t0err
#redshift of N300-1.Q0.SN
z = 0.063
zerr = 0.005
#Extinction coefficient (galactic) in each band S & F (2011)
EBVgal = 0.0107
Coefs = np.array([3.641, 2.682, 1.516])
band = ['B','V','i']
#max time
Tmax = 282.01
#explosion parameters
m_c = 0.97 #1.36/1.4 Mchandra
e_51 = 0.99 #x10^51 ergs
m_c_err = 0.12
e_51_err = 0.20
#bands
band = ['B','V','i']
Band = ['B','V','I']

print "loading binned data"
#N300-1.Q0.SN binned time series data files
Bfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S1.txt"
Vfile = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S1.txt"
Ifile = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S1.txt"
binfiles = [Bfile, Vfile, Ifile] 
#get N300-1.Q0.SN binned light curve
t, M, M_err, F, SN, Mlim1 = LCload(binfiles, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, scols=9, flags=['-99.99999'], mode='multi')
#N300-1.Q0.SN binned time series data files
Bfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S2.txt"
Vfile = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S2.txt"
Ifile = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S2.txt"
binfiles = [Bfile, Vfile, Ifile] 
#get N300-1.Q0.SN binned light curve
t, M, M_err, F, SN, Mlim2 = LCload(binfiles, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, scols=9, flags=['-99.99999'], mode='multi')
#get noise in flux
F_err = [F[i]/SN[i] for i in range(3)]
#N300-1.Q0.SN binned time series data files
Bfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S3.txt"
Vfile = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S3.txt"
Ifile = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lcbin.CN_180327.S3.txt"
binfiles = [Bfile, Vfile, Ifile] 
print "Loading binned early light curve."
#get N300-1.Q0.SN binned light curve
t, M, M_err, F, SN, Mlim3 = LCload(binfiles, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, scols=9, flags=['-99.99999'], mode='multi')
#get noise in flux
F_err = [F[i]/SN[i] for i in range(3)]
#limiting magnitudes at different SNR
limconf = [Mlim1, Mlim2, Mlim3]

#deredden flux and get templates
print "Correcting for galactic reddening"
flimconf = limconf
for i in range(len(F)):
    #correct fluxes for AB calibration
    if i == 2:
        F[i] = F[i]*1.42392
        F_err[i] = F_err[i]*1.42392
    #correct fluxes for galactic reddening
    F[i] = deredFlux(F[i], EBVgal, Coefs[i])
    F_err[i] = deredFlux(F_err[i], EBVgal, Coefs[i])
    #shift to zero centered time
    t[i] = t[i] - Tmax
    #get explosion centered time in observer frame
    t[i] = t[i] - t0

    #get limiting magnitudes
    for n in range(len(limconf)):
        limconf[n][i] = deredMag(limconf[n][i], EBVgal, Coefs[i])
        flimconf[n][i] = Mag_toFlux(band[i], limconf[n][i])*10**6
        #tlim.append(t[i][SN[i]<2.0])
        #Flim[i] = Flim[i][SN[i]<2.0]

        #get useful times
        t[i], F[i], F_err[i], flimconf[n][i] = LCcrop(t[i], -10, 10, F[i], F_err[i], flimconf[n][i])
    
"""   
print "plotting early data"
#plot
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("t [days]")
ax[1].set_ylabel("L/Lmax")
for i in range(len(t)):
    #fit early light curve
    ax[i].errorbar(t[i],F[i],F_err[i],fmt='g+')
    ax[i].scatter(t[i], Flim[i], marker='v', c='r')
    ax[i].plot(t[i],[0]*len(t[i]),'k')
plt.show()
"""

#print "Trial thetas:", np.linspace(0,180,10)
tk = np.arange(0.001, 10, 0.001)
Fks = [[],[],[]]
###################################################################
"""
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
    Fp, Fp_err = F[i][F[i]>0], F_err[i][F[i]>0]
    tp = t[i][F[i]>0]
    M, M_err = Flux_toMag(Band[i],Fp*1e-6,Fp_err[i]*1e-6)
    Mlim = Flux_toMag(Band[i],Flim[i]*1e-6)
    ax[i].errorbar(tp,M,M_err,fmt='g+')
    ax[i].scatter(t[i], Mlim, marker='v', c='r')
    ax[i].set_ylabel(Band[i], fontsize=16, fontstyle='italic', fontweight='bold')
    ax[i].set_ylim(25,18)
    ax[i].yaxis.set_ticks([24,22,20])
    ax[i].tick_params(labelsize=14)
    #plot each zero theta model
    for Fk in Fks[i]:
        ts, Fs = tk[Fk>0], Fk[Fk>0]
        Ms = Flux_toMag(Band[i],Fs*1e-6)
        ax[i].plot(ts, Ms)
f.subplots_adjust(hspace=0)
plt.show()
"""

"""
f, ax = plt.subplots(1, 3, sharey=True)
ax[0].set_ylabel("Percentage of ruled out angles", fontsize=16)
ax[1].set_xlabel("Confidence Level (%)", fontsize=16)
epochs = [t0-t0err, t0, t0+t0err]
confidences = [np.arange(99.50,99.99,0.01), np.arange(99.50,99.99,0.01), np.arange(97.0,99.99,0.05)]
lims = [(99.51,100.0),(99.51,100.0),(97.1,100.0)]
for n in range(3):
    t0 = epochs[n]
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
    
    print "Computing viewing angles at each confidence interval"
    #list of confidence intervals to traverse
    confs = confidences[n]
    #list of sample models
    a13s = [2.0, 0.2, 0.05] #1RG, 6MS, 2MS
    #list of viewing angles
    thetas = np.linspace(0,180,100)

    #Lk_terr = []
    #for i in range(len(t)):
    #    tk = t[i][t[i]>0]
    #    #perturb time axis by error in epoch
    #    dts = np.random.normal(0.0, t0err, 100)
    #    Lk_dif = np.zeros([len(dts), len(tk)])
    #    for r, dt in enumerate(dts):
    #        Lk_dif[r] = Kasen2010(tk+dt, a13, m_c, e_51)[0]
    #    Lk_terr.append(Lk_dif.std(axis=0))
    
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
                    #Lk_theta_terr = Lk_terr[i]*Kasen_isocorr(theta)
                    #Fk_terr = KasenFit(Lk_theta_terr,Tk,wave_0[bands[band[i]]],z,zerr)[0]
                    #total error
                    #Fk_err = np.square(Fk_terr)+np.square(Fk_zerr)
                    #Fk_err = np.array([np.sqrt(f) for f in Fk_err])
                    Fk_err = Fk_zerr
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
        ax[n].plot(confs, outangles)
        ax[n].set_xlim(lims[n])
        ax[n].set_ylim(0,101)
plt.show()
"""

"""
print "Computing viewing angles at each separation distance"
#list of sample models
#a13s = np.concatenate([np.arange(0.001,0.2,0.001), np.arange(5.7,5.9,0.001)])#1RG, 6MS, 2MS
#a13s = np.arange(0.001,0.5,0.001) #1RG, 6MS, 2MS
a13s = np.concatenate((np.arange(0.001,0.05,0.001), np.arange(0.05,0.2,0.01), np.arange(0.2, 2.0, 0.01), np.arange(2.0,10.1,0.1)))
#confs = [68.27, 95.45, 99.73]
confs = [99, 99.5, 99.73]
#list of viewing angles
thetas = np.linspace(0,180,100)

#print "Calculating t err"
#Lk_terr = []
#perturb time axis by error in epoch
#dts = np.random.normal(0, t0err, 10000000)
#perturb z by error
#dzs = np.random.normal(0, zerr, 100000)
#account for epoch error
#for i in range(len(t)):
#    tk = t[i][t[i]>0]
    #dts = np.random.normal(0, t0err, 10)
#    Lk_dif = np.zeros([len(dts), len(tk)])
#    for r, dt in enumerate(dts):
#        Lk_dif[r] = Kasen2010(absTime(tk+dt,z), a13, m_c, e_51)[0]
#    Lk_terr.append(Lk_dif.std(axis=0))


#Ts = []
#Ls = []
#Lf = []
#Es = []
style = ['k:', 'k-.', 'k-']
#for each confidence interval
for n, conf in enumerate(confs):
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
            #get theoretical light curve
            Lk = np.array([KasenFit(ti, a13, 1.0, wave_0[bands[band[i]]], m_c, e_51, z, 0) for ti in t[i]])
            #Lk, Tk = Kasen2010(absTime(tk,z), a13, m_c, e_51)

            #if n == 1 and i == 2:
            #    Ts.append(max(Tk))
            #    Ls.append(max(Lk))
            #    Lf.append(max(KasenFit(Lk,Tk,wave_0[bands[band[i]]],z,zerr)[0]))
            #    Es.append(KasenFit(Lk_terr[i],Tk,wave_0[bands[band[i]]],z,zerr)[0][np.argmax(KasenFit(Lk,Tk,wave_0[bands[band[i]]],z,zerr)[0])])
            
            #for each angle
            for k, theta in enumerate(thetas):
                #angle corrected Kasen luminosity
                Lk_theta = Lk*Kasen_isocorr(theta)
                #observer frame angle corrected Kasen flux
                #Fk, Fk_zerr = KasenFit(Lk_theta,Tk,wave_0[bands[band[i]]],z,zerr)
                Fk = Lk*Kasen_isocorr(theta)
                #epoch error
                #Lk_theta_terr = Lk_terr[i]*Kasen_isocorr(theta)
                #Fk_terr = KasenFit(Lk_theta_terr,Tk,wave_0[bands[band[i]]],z,zerr)[0]
                #total error
                #Fk_err = np.square(Fk_terr)+np.square(Fk_zerr)
                #Fk_err = np.array([np.sqrt(f) for f in Fk_err])
                #Fk_err = Fk_zerr
                #total error
                #Err = np.square(Ferr)+np.square(Fk_err)
                #Err = np.array([np.sqrt(E) for E in Err])
                #check if any points rule out angle with conf
                if any(Fk > F[i] + sig*F_err[i]):
                    mask[k] = False
                else:
                    print "Consistent!", a13, conf, theta
        #At this confidence level, we rule out some percent of angles
        outangles[j] = 180.0*float(len(thetas)-len(thetas[mask]))/len(thetas)
    #At this conf, we can plot ruled out angles vs a13
    print conf, outangles
    plt.plot(a13s, outangles, style[n])
plt.xlim(0,10.0)
plt.ylim(0,185)
plt.plot([0.05,0.05], [0,5], 'k', linewidth=1)
plt.plot([0.2,0.2], [0,5], 'k', linewidth=1)
plt.plot([2.0,2.0], [0,5], 'k', linewidth=1)
#plt.plot([0.05,0.05], [175,185], 'k', linewidth=1)
#plt.plot([0.2,0.2], [175,185], 'k', linewidth=1)
#plt.plot([2.0,2.0], [175,185], 'k', linewidth=1)
plt.ylabel("Unacceptable viewing angles (deg)", fontsize=16)
plt.xlabel("Separation Distance ($10^{13}$ cm)", fontsize=16)
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.show()
"""


print "Computing viewing angles at each separation distance"
#list of sample models
#a13s = np.arange(6.01,10.01,0.1) #1RG, 6MS, 2MS
a13s = np.concatenate((np.arange(0.001,0.05,0.005), np.arange(0.05,0.2,0.02), np.arange(0.2, 2.0, 0.05), np.arange(2.0,11.0,1.0)))
#SNR of 1, 2 3 respectively
confs = [84.1345, 97.7250, 99.8650]
print [norm.ppf(conf/100.0) for conf in confs]
print a13s
#list of viewing angles
thetas = np.linspace(0,180,100)

#function: test a given a13
def gen_a13(a13):
    print a13
    #for each band
    Fks = []
    Fk_errs = []
    for i in range(len(t)):
        #print band[i]
        Fk = np.zeros(len(t[i]))
        Fk_err = np.zeros(len(t[i]))
        for r in range(len(t[i])):
            #get theoretical light curve at each time
            #1000000,1000000,1000000,1000000
            #200000,200000,200000,200000
            Fk[r], Fk_err[r] = MCerr(KasenFit, [t[i][r], a13, 1.0,
                                                wave_0[bands[band[i]]]],
                                     [m_c, e_51, z, 0],
                                     [m_c_err, e_51_err, zerr, t0err],
                                     [10000,10000,10000,10000], confs[-1])
            #Assumptions here:
            #Independent parameters is a good assumption (MCerr uses this)
            #No covariance simulation needed.
            #for max 3sig confidence, only ~1000 trials needed for MC.
            
            #print Fk[r], Fk_err[r]
        #if i == 0:
            #print a13, max(Fk), Fk_err[np.argmax(Fk)]
            #tt = np.linspace(0,40,1000)
            #Ft = [KasenFit(ti, a13, 1.0, wave_0[bands[band[i]]], m_c, e_51, z, 0) for ti in tt]
            #print a13, tt[np.argmax(Ft)]
        #for each angle
        Fks.append(Fk)
        Fk_errs.append(Fk_err)
    print "done", a13
    return Fks, Fk_errs

def test_a13(gen, gen_err, sig, flim):
    #boolean mask for whether angle is ruled out to given confidence
    print sig
    mask = np.array([True]*len(thetas))
    
    for i in range(len(t)):
        for k, theta in enumerate(thetas):
            #check if any points rule out angle with conf
            if ruleout(F[i], F_err[i], gen[i], gen_err[i], theta, sig, flim[i]):
                mask[k] = False
            #else:
                #print "Consistent!", band[i], a13, norm.cdf(sig), theta

    #At this confidence level, we rule out some percent of angles
    outangles = 180.0*float(len(thetas)-len(thetas[mask]))/len(thetas)
    
    outfile = open("log.txt", 'a')
    outfile.write(str(norm.cdf(sig))+"\t"+str(a13)+"\t"+str(outangles)+"\n")
    outfile.close()
    
    return outangles

nproc = 32

print "Generating Synthetic Light Curves"
#generate synthetic light curves
pool = Pool(nproc)
procs = []
#for each sample model
for j, a13 in enumerate(a13s):
    procs.append(pool.apply_async(gen_a13, [a13]))
#array to hold percent of viewing angles ruled out at each conf
genlcs = []
print "Processes", len(procs)
for n, proc in enumerate(procs):
    print "getting proc",n
    genlcs.append(proc.get())
    print "got proc",n
pool.terminate()
print "Generated Light Curves"

print "Checking Against Observations"
style = ['k:', 'k--', 'k-']
outangles = []
#for each confidence interval
for n, conf in enumerate(confs):
    #sigma needed to establish confidence below LC
    sig = norm.ppf(conf/100.0)
    
    pool = Pool(nproc)
    procs = []
    #for each sample model
    for j, a13 in enumerate(a13s):
        Fks = genlcs[j][0]
        Fk_errs = genlcs[j][1]
        procs.append(pool.apply_async(test_a13, [Fks, Fk_errs, sig, flimconf[n]]))
        if n = 2 and a13 > 0.555 and a13 < 0.639:
            #plot section
            f, ax = plt.subplots(len(t), sharex=True) 
            for i in range(len(t)):
                #ax[i].errorbar(t1[i], M1[i], yerr=M1_err[i], fmt='r+', label='SrcExt')
                ax[i].errorbar(t[i], M[i], yerr=sig*M_err[i], fmt="k+")
                ax[i].errorbar(t[i], Fks[i], yerr=sig*Fk_errs[i], fmt="g+")
                ax[i].scatter(t[i], flimconf[n][i], color='r', marker='v')
            plt.show()
        
    #array to hold percent of viewing angles ruled out at each conf
    outangles.append([proc.get() for proc in procs])
    pool.terminate()
print "Checked Observations"
    
    #At this conf, we can plot ruled out angles vs a13
    #print conf, outangles
outangles = np.array(outangles)
out = np.concatenate(([a13s], outangles), axis=0)
    
#for n, conf in enumerate(confs):
    #plt.plot(a13s, outangles[n], style[n])
    #print "DONE!"
    #print "confidence", conf
    #print outangles[n]

np.savetxt("kasen.txt", out.T)

    
#plt.xlim(0,1.0)
#plt.ylim(0,185)
#plt.plot([0.05,0.05], [0,5], 'k', linewidth=1)
#plt.plot([0.2,0.2], [0,5], 'k', linewidth=1)
#plt.plot([2.0,2.0], [0,5], 'k', linewidth=1)
#plt.plot([0.05,0.05], [175,185], 'k', linewidth=1)
#plt.plot([0.2,0.2], [175,185], 'k', linewidth=1)
#plt.plot([2.0,2.0], [175,185], 'k', linewidth=1)
#plt.ylabel("Unacceptable viewing angles (deg)", fontsize=16)
#plt.xlabel("Separation Distance ($10^{13}$ cm)", fontsize=16)
#plt.tick_params(labelsize=14)
#plt.tight_layout()
#plt.savefig("/home/chrisni/trials/kasen.pdf")
#plt.show()


"""
print Ls
print Ts
print Lf
plt.plot(a13s, Ls)
plt.show()
plt.plot(a13s, Lf)
plt.show()
plt.plot(a13s, Es)
plt.show()
plt.plot(a13s, 2900000.0/np.array(Ts))
plt.plot(a13s, [445]*len(a13s))
plt.plot(a13s, [551]*len(a13s))
plt.plot(a13s, [806]*len(a13s))
plt.show()
"""

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
