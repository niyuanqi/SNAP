#################################################################
# Name:     ArnettFit.py                                        #
# Author:   Yuan Qi Ni                                          #
# Date:     May, 31, 2016                                      #
# Function: Program corrects light curve for galactic reddening,#
#           fits light curve with 10 parameter function, and    #
#           constructs UVOIR light curve using SN1a template.   #
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

sn_file = "N300-1.Q0.SN.txt"
plot = False #plot polynomial fits to light curves

print "Loading SN File"
s = get_sn(sn_file)
#s.z = z
print s.z
z = s.z
zerr = 0.003
n = 50 #number of monte carlo trials
print 3663.0/(1+s.z), 8750/(1+s.z)

#Bands in light curve
band = ['B','V','i']
#Extinction coefficient (galactic) in each band S & F (2011)
EBVgal = 0.107
Coefs = np.array([3.641, 2.682, 1.516])
#Epoch of explosion
t0 = -17.73
t0_err = 0.79

s = get_sn("N300-1.Q0.SN.txt")
#don't plot fit
s.replot = 0
#fit SN
s.choose_model("EBV_model2", stype='st')
s.fit(['B','V','i'])

#calculate Bolometric Luminosity
sBs = []
dzs = np.random.normal(0,zerr,n)
for i in range(n):
    s.z = z+dzs[i]
    print i, s.z
    #luminosity distance
    dl = intDl(s.z)
    #distance modulus
    DM = 5.024*np.log10(dl/10)
    
    tB, sB, fB, lB = s.bolometric(['B','V','i'], method='SED', refband='i',
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

tB, sB, sBerr = LCcrop(tB, -15, 30, sB, M_err=sBerr)
#window in which to perform Arnett fit
t1 = -10
t2 = 15
tfit, sfit, sfiterr = LCcrop(tB, t1, t2, sB, M_err=sBerr)
#plt.errorbar(tB, np.log10(sB), yerr=sBerr/np.log(10)/sB, marker='o', c='g')
#plt.show()
#perform fit to find peak of function

fit, fit_err, params, params_err = LCpolyFit(tfit, -sfit, sfiterr, order=4, N=1000, plot=True)
print params, params_err

#draw some points around max by monte carlo
N=10000
tmax = np.random.normal(-t0, t0_err, N)
Lmax = np.random.normal(-params[1], params_err[1], N)
tmax_err = t0_err
Lmax_err = params_err[1]

x2cent = (np.square((-params[1]-Lmax)/Lmax_err) + np.square((-t0-tmax)/tmax_err)).sum()
x2edge = (np.square((-params[1]+Lmax_err-Lmax)/Lmax_err) + np.square((-t0+t0_err-tmax)/tmax_err)).sum()
x2diff = x2edge - x2cent

print np.sqrt(x2cent)

print "Fitting Arnett Function to Peak of LC"
p1 = np.linspace(0.275,0.285,10)
p2 = np.linspace(1.15,1.25,10)
p1, p2 = np.meshgrid(p1, p2)
p1 = p1.flatten()
p2 = p2.flatten()
p = np.array([p1,p2]).T

chi2 = np.zeros(len(p))
for i in range(len(p)):
    chi2[i] = np.square(ArnettMaxErr(p[i], tmax, tmax_err, Lmax, Lmax_err)).sum()
    print i, chi2[i]
chiarg = (chi2-x2cent)<x2diff
chimin = np.argmin(chi2)
popt = p[chimin]
perr = p[chiarg].std(axis=0)
print popt, perr

#p0 = [0.27, 1.2]
#popt, cov_x, infodict, mesg, ier = leastsq(ArnettMaxErr, p0, args=(tmax, tmax_err, Lmax, Lmax_err), full_output=1, ftol=1.49012e-10, xtol=1.49012e-10, factor=0.01)
#perr = np.sqrt(np.diag(cov_x))

MNi, MNi_err = popt[0], perr[0]
vej, vej_err = 11.0, 1.0
Mej, Mej_err, Eej, Eej_err = ArnettMejE(popt[1],perr[1],vej*10**8,vej_err*10**8)
print popt, perr
print "Assumed ejecta velocity:", vej, vej_err, "km/s"
print "Nickel mass:", MNi, MNi_err, "Msun"
print "Ejecta mass:", Mej, Mej_err, "Msun"
print "Ejecta Ek:", Eej, Eej_err, "x10^51ergs"

t_arnett = np.arange(0.25,103.25,0.25)
L_arnett = ArnettFit(t_arnett, *popt)
plt.plot(t_arnett, L_arnett)
plt.errorbar(tmax, Lmax, xerr=tmax_err, yerr=Lmax_err, fmt='r+')
plt.show()
logL_arnett = np.log10(L_arnett)
print max(logL_arnett)

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
plt.ylim(41,43)
#plt.legend()
plt.show()
