#################################################################
# Name:     LCplot.py                                           #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Program generates light curve plot using input      #
#           light curve data. Light curve data produced by      #
#           LCgen.py routine. Update ObjData.py first.          #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt

#essential files
from SNAP.Analysis.LCRoutines import*
from ObjData import *

#get light curve for SN>SNthres=2.0
t, M, M_err, F, SN, Mlim = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=2.0, scols=9, flags=['-99.99999'], mode='multi')
print "Maximum measured magnitudes:",[min(m) for m in M]
print "Mean limiting magnitude:",[l.mean() for l in Mlim]
print "Std of limiting magnitude:",[np.std(lims) for lims in Mlim]
print "Average cadence:",[np.mean(time[1:] - time[:-1]) for time in t]
print ""

#crop window in data
for i in range(len(M)):
    t[i], M[i], M_err[i], F[i], SN[i], Mlim[i] = LCcrop(t[i], t1, t2, M[i], M_err[i], F[i], SN[i], Mlim[i])
F_err = [F[i]/SN[i] for i in range(len(t))]

#get color curve
tc, C, C_err = LCcolors(t, M, M_err)

print "Plotting data"
#plot
f, ax = plt.subplots(len(t)+len(t)-1, sharex=True) 
ax[-1].set_xlabel("Days into 2018")

for i in range(len(t)):
    ax[i].errorbar(t[i], M[i], yerr=M_err[i], fmt="k", label='MagCalc')
    ax[i].scatter(t[i], M[i])
    ax[i].scatter(t[i],Mlim[i],color='r',marker='v', label='det. limit')
    ax[i].set_ylim(24.5,11.5)
    ax[i].set_xlim(t1, t2)
    #ax[i].legend()
    if i == 0:
        ax[i].set_ylabel("B mag")
    elif i == 1:
        ax[i].set_ylabel("V mag")
    elif i == 2:
        ax[i].set_ylabel("I mag")

for i in range(len(t)-1):
    ax[i+len(t)].errorbar(tc[i], C[i], C_err[i], fmt='k', label='1st Diff')
    ax[i+len(t)].set_xlim(t1, t2)
    ax[i+len(t)].set_ylim(-2.6, 2.6)
    if i == 0:
        ax[i+len(t)].set_ylabel("B-V mag")
    elif i == 1:
        ax[i+len(t)].set_ylabel("V-I mag")

f.subplots_adjust(hspace=0)
plt.show()

print "plotting early data"
#get light curve
t, M, M_err, F, SN, Mlim = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-10.0, scols=9, flags=['-99.99999'], mode='multi')

for i in range(len(M)):
    t[i], M[i], M_err[i], F[i], SN[i], Mlim[i] = LCcrop(t[i], t1_early, t2_early, M[i], M_err[i], F[i], SN[i], Mlim[i])
F_err = [F[i]/SN[i] for i in range(len(t))]
#plot early light curve
f, ax = plt.subplots(3, sharex=True)
ax[-1].set_xlabel("t rest [days]")
ax[1].set_ylabel("Flux [uJy]")
for i in range(len(t)):
    ax[i].errorbar(t[i], F[i], yerr=F_err[i], fmt='g+')
    ax[i].set_xlim(t1_early,t2_early)
plt.show()
