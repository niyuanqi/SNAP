#################################################################
# Name:     ArnettFit.py                                        #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Program constructs UVOIR light curve using SNooPy   #
#           SN1a template. Fits using Arnett model.             #
#           Update ObjData.py first.                            #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.integrate as integrate
from scipy.optimize import leastsq
from snpy import *

#essential files
from SNAP.Analysis.LCRoutines import*
from SNAP.Analysis.Cosmology import*
from SNAP.Analysis.LCFitting import*
from ObjData import *

plot = True #plot polynomial fits to light curves
print "Loading SN File"
s = get_sn(sn_file)
#don't plot fit
s.replot = 0
#fit SN
s.choose_model("EBV_model2", stype='st')
s.fit(band)

print "Bolometric integration range (angstrom):", 3663.0/(1+s.z), 8750/(1+s.z)
#calculate Bolometric Luminosity
n=500 #use some monte carlo to get bolometric luminosity error
sBs = []
dzs = np.random.normal(0,zerr,n)
for i in range(n):
    s.z = z+dzs[i]
    print i, s.z
    #luminosity distance
    dl = intDl(s.z)
    #distance modulus
    DM = 5.024*np.log10(dl/10)
    
    tB, sB, fB, lB = s.bolometric(['B','V','i'], method='SED',
                                  EBVhost=0, Rv=0, SED='H3',
                                  lam1=3663.0/(1+s.z),
                                  lam2=8750/(1+s.z),
                                  DM=DM,
                                  use_stretch=True)
    print len(sB)
    sBs.append(sB)
sBs = np.array(sBs)
sB = sBs.mean(axis=0)
sBerr = sBs.std(axis=0)

tB, sB, sBerr = LCcrop(tB, -20, 50, sB, M_err=sBerr)
#window in which to perform Arnett fit
tph1 = -10 #from -10 days to 15 days, this is the "photospheric phase"
tph2 = 15
tfit, sfit, sfiterr = LCcrop(tB, tph1, tph2, sB, M_err=sBerr)

#perform fit to find peak of function
fit, fit_err, params, params_err = LCpolyFit(tfit, -sfit, sfiterr, order=4, N=1000, plot=plot)
print params, params_err
tmax = -t0+params[0]
tmax_err = np.sqrt(np.square(t0err)+np.square(params_err[0]))

#Use monte carlo to find optimal parameters of MNi and MejEk
n=5000
p0 = 1.2 #initial guess for ejecta MejEk parameter
print ""
print "Performing Monte Carlo to find optimal MNi, and MejEk"
print "This takes a while, go for a walk."
MNi, MejEk, MNi_err, MejEk_err  = ArnettIntercept(tmax, -params[1], tmax_err, params_err[1], p0, n=n, nproc=4)
print "DONE!"

vej, vej_err = 11.0, 1.0
#break degeneracy in MejEk
Mej, Mej_err, Eej, Eej_err = ArnettMejE(MejEk,MejEk_err,vej*10**8,vej_err*10**8)
print ""
print "Nickel mass:", MNi, MNi_err, "Msun"
print "Ejecta mass-energy parameter:", MejEk, MejEk_err
print "Assumed ejecta velocity:", vej, vej_err, "km/s"
print "Ejecta mass:", Mej, Mej_err, "Msun"
print "Ejecta Ek:", Eej, Eej_err, "x10^51ergs"

t_arnett = np.arange(0.25,103.25,0.25)
L_arnett = ArnettFit(t_arnett, MNi, MejEk)
#plt.plot(t_arnett, L_arnett)
#plt.errorbar(tmax, Lmax, xerr=tmax_err, yerr=Lmax_err, fmt='r+')
#plt.show()
logL_arnett = np.log10(L_arnett)
print "Maximum luminosity (ergs/s):",max(L_arnett)

print "plotting peak fit"
#plot Arnett fit
#gotta get error of fit somehow
#plt.title('Arnett Fitted to UVOIR Bolometric Light Curve')
plt.errorbar(tB, np.log10(sB), yerr=sBerr/np.log(10)/sB, marker='o', c='g')
plt.plot(t_arnett+t0, logL_arnett, c='r')
#plt.errorbar(tB, sB, yerr=sBerr, marker='o', c='g')
#plt.plot(t_arnett, L_arnett, c='r')
plt.ylabel('$log_{10}(L)$ [erg/s]', fontsize = 14)
plt.xlabel('Time [rest-frame days]', fontsize = 14)
plt.xlim(-20,40)
plt.ylim(40.5,43.5)
#plt.legend()
plt.show()
