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
from scipy.optimize import curve_fit
from snpy import *

#essential files
from SNAP.Analysis.LCRoutines import*
from SNAP.Analysis.Cosmology import*
from SNAP.Analysis.LCFitting import*

sn_file = "N300-1.Q0.SN.txt"
plot = False #plot polynomial fits to light curves

#redshift of N300-1.Q0.SN
z = 0.065
z_err = 0.005

print "Loading SN File"
s = get_sn(sn_file)
#s.z = z
print s.z
print 3663.0/(1+s.z), 8750/(1+s.z)

#Bands in light curve
band = ['B','V','i']
#Extinction coefficient (galactic) in each band S & F (2011)
EBVgal = 0.107
Coefs = np.array([3.641, 2.682, 1.516])

s = get_sn("N300-1.Q0.SN.txt")
#don't plot fit
s.replot = 0
#fit SN
s.choose_model("EBV_model2", stype='st')
s.fit(['B','V','i'])

#luminosity distance
dl = intDl(s.z)
dlerr1 = intDl(s.z+z_err)
dlerr2 = intDl(s.z-z_err)

#distance modulus
DM = 5.024*np.log10(dl/10)
DMerr1 = 5.024*np.log10(dlerr1/10)
DMerr2 = 5.024*np.log10(dlerr2/10)

#Calculate Bolometric Luminosity
tB, sB, fB, lB = s.bolometric(['B','V','i'], method='SED', refband='i',
                                  EBVhost=0, Rv=0, SED='H3',
                                  lam1=3663.0/(1+s.z),
                                  lam2=8750/(1+s.z),
                                  DM=DM,
                                  use_stretch=True)
#calculate Bolometric Luminosity error
t, sBerr1, f, l = s.bolometric(['B','V','i'], method='SED', refband='i',
                                  EBVhost=0, Rv=0, SED='H3',
                                  lam1=3663.0/(1+s.z+z_err),
                                  lam2=8750/(1+s.z-z_err),
                                  DM=DMerr1,
                                  use_stretch=True)
t, sBerr2, f, l = s.bolometric(['B','V','i'], method='SED', refband='i',
                                  EBVhost=0, Rv=0, SED='H3',
                                  lam1=3663.0/(1+s.z-z_err),
                                  lam2=8750/(1+s.z-z_err),
                                  DM=DMerr2,
                                  use_stretch=True)
sBerr = np.absolute((sBerr1-sBerr2)/2.0)
tB, sB, sBerr = LCcrop(tB, -15, 30, sB, M_err=sBerr)
#window in which to perform Arnett fit
t1 = -10
t2 = 15
tfit, sfit, sfiterr = LCcrop(tB, t1, t2, sB, M_err=sBerr)
#plt.errorbar(tB, np.log10(sB), yerr=sBerr/np.log(10)/sB, marker='o', c='g')
#plt.show()

print "Fitting Arnett Function to Peak of LC"
p0 = [-16.0, 0.2, 1.0]
popt, pcov = curve_fit(ArnettFit, tfit, sfit, sigma=sfiterr, p0=p0)
perr = np.sqrt(np.diag(pcov))
print popt, perr

t_arnett = np.arange(0.25,103.25,0.25)+popt[0]
L_arnett = ArnettFit(t_arnett, *popt)
t_trial, L_trial = LCcrop(t_arnett, t1, t2, L_arnett)
plt.plot(t_trial, L_trial)
plt.plot(tfit, sfit)
plt.show()
logL_arnett = np.log10(L_arnett) 

print "plotting peak fit"
#plot Arnett fit
#gotta get error of fit somehow
#plt.title('Arnett Fitted to UVOIR Bolometric Light Curve')
plt.errorbar(tB, np.log10(sB), yerr=sBerr/np.log(10)/sB, marker='o', c='g')
plt.plot(t_arnett, logL_arnett, c='r')
#plt.errorbar(tB, sB, yerr=sBerr, marker='o', c='g')
#plt.plot(t_arnett, L_arnett, c='r')
plt.ylabel('$log_{10}(L)$ [erg/s]', fontsize = 14)
plt.xlabel('Time to maximum [rest-frame days]', fontsize = 14)
plt.xlim(-20,40)
plt.ylim(41,43)
#plt.legend()
plt.show()
