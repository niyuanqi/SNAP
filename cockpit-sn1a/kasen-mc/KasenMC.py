#################################################################
# Name:     KasenMC.py                                          #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 27, 2017                                       #
# Function: Program computes kasen models at various separation #
#           distances and viewing angles, comparing them to     #
#           early light curve. Uses Monte Carlo to estimate     #
#           uncertainties caused by errors in estimated         #
#           parameters ejecta mass and kinetic energy, redshift #
#           and epoch of first light. Uses those along with     #
#           flux errors to estimate confidence intervals.       #
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
    limconf.append(Mlimc)

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
#list of sample models
a13s = np.concatenate((np.arange(0.001,0.05,0.005), np.arange(0.05,0.2,0.02), np.arange(0.2, 2.0, 0.05), np.arange(2.0,11.0,1.0)))
#list of viewing angles
thetas = np.linspace(0,180,100)
#SNR of 1, 2, 3, 4, 5
confs = [norm.cdf(sn) for sn in limSNs]
print "Confidence levels:",confs
print "Sigma levels:",limSNs
print "Trial a13s:", a13s
print "Trial thetas:",thetas

#Note in log for reference
outfile = open(logfile, 'a')
outfile.write(" ; Confs "+str(confs))
outfile.write(" ; Sigs "+str(limSNs))
outfile.close()

#function: test a given a13
def gen_a13(a13, conf):
    print "Generating model at a13:",a13
    #for each band
    Fks = []
    Fk_errs = []
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
                                     [nmc,nmc,nmc,nmc], conf)
            #Assumptions here:
            #Independent parameters is a good assumption (MCerr uses this)
            #No covariance simulation needed.
        Fks.append(Fk)
        Fk_errs.append(Fk_err)
    print "done model at a13:", a13
    return Fks, Fk_errs

def test_a13(gen, gen_err, sig, flim):
    #boolean mask for whether angle is ruled out to given confidence
    print "Testing model at sig:", sig
    mask = np.array([True]*len(thetas))
    
    for i in range(len(t)):
        for k, theta in enumerate(thetas):
            #check if any points rule out angle with conf
            if ruleout(F[i], F_err[i], gen[i], gen_err[i], theta, sig, flim[i]):
                mask[k] = False
    
    #At this confidence level, we rule out some percent of angles
    outangles = 180.0*float(len(thetas)-len(thetas[mask]))/len(thetas)
    #record some stuff into log for some intermediary reference
    outfile = open(logfile, 'a')
    outfile.write(str(norm.cdf(sig))+"\t"+str(a13)+"\t"+str(outangles)+"\n")
    outfile.close()
    #Angles ruled out at confidence sig 
    return outangles

nproc = 32

print "Generating Synthetic Light Curves"
#generate synthetic light curves
pool = Pool(nproc)
procs = []
#for each model
for j, a13 in enumerate(a13s):
    #start process to generate model, taking highest confidence interval
    procs.append(pool.apply_async(gen_a13, [a13, confs[-1]]))
#array to hold generated light curves
genlcs = []
#get processes
print "Processes", len(procs)
for n, proc in enumerate(procs):
    print "getting proc",n
    genlcs.append(proc.get())
    print "got proc",n
pool.terminate()
print "Generated Light Curves"

print "Checking Against Observations"
outangles = []
#for each confidence interval
for n, conf in enumerate(confs):
    #sigma needed to establish confidence below LC
    sig = norm.ppf(conf)
    
    pool = Pool(nproc)
    procs = []
    #for each generated model
    for j, a13 in enumerate(a13s):
        Fks = genlcs[j][0]
        Fk_errs = genlcs[j][1]
        #start process to test model against data
        procs.append(pool.apply_async(test_a13, [Fks, Fk_errs, sig, flimconf[n]]))
        """
        #Diagnostic tool
        if n == 2 and a13 > 0.195 and a13 < 0.205:
            print "plotting section"
            #plot section
            f, ax = plt.subplots(len(t), sharex=True) 
            for i in range(len(t)):
                c = Kasen_isocorr(80)
                ax[i].errorbar(t[i], F[i], yerr=sig*F_err[i], fmt="k+")
                ax[i].errorbar(t[i], Fks[i]*c, yerr=sig*Fk_errs[i]*c, fmt="g+")
                ax[i].scatter(t[i], flimconf[n][i], color='r', marker='v')
            plt.subplots_adjust(hspace=None)
            plt.show()
        """
    #percent of viewing angles ruled out for each a13 model at this conf
    outangles.append([proc.get() for proc in procs])
    pool.terminate()
print "Checked Against Observations."
print "DONE!"    

print "Saving..."
#Output ruled out angles vs a13
outangles = np.array(outangles)
out = np.concatenate(([a13s], outangles), axis=0)
np.savetxt(kasfile, out.T)
