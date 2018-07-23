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
t, M, M_err, F, SN, Mlim, ra, dec = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=2.0, racols=2, deccols=3, scols=9, flags=['-99.99999'], mode='multi')
print "Maximum measured magnitudes:",[min(m) for m in M]
print "Mean limiting magnitude:",[l.mean() for l in Mlim]
print "Std of limiting magnitude:",[np.std(lims) for lims in Mlim]
print "Average cadence:",[np.mean(time[1:] - time[:-1]) for time in t]
print ""

#crop window in data
for i in range(len(M)):
    ttemp, ra[i], dec[i] = LCcrop(t[i], t1, t2, ra[i], dec[i])
    t[i], M[i], M_err[i], F[i], SN[i], Mlim[i] = LCcrop(t[i], t1, t2, M[i], M_err[i], F[i], SN[i], Mlim[i])
F_err = [F[i]/SN[i] for i in range(len(t))]
for i in range(len(M)):
    limtemp, ra[i], dec[i] = LCcrop(Mlim[i], limcuts[i], 100, ra[i], dec[i])
    Mlim[i], t[i], M[i], M_err[i], F[i], SN[i] = LCcrop(Mlim[i], limcuts[i], 100, t[i], M[i], M_err[i], F[i], SN[i])

#Get average position of source
ra = np.concatenate(ra)
dec = np.concatenate(dec)
ra_ave = '{:3.9f}'.format(ra.mean())
dec_ave = '{:3.9f}'.format(dec.mean())
ra_err = '{:3.9f}'.format(ra.std())
dec_err = '{:3.9f}'.format(dec.std())
print "Mean position of source: {: >14} {: >14}".format(ra_ave, dec_ave)
print "Errors in mean position: {: >14} {: >14}".format(ra_err, dec_err)
print ""

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
    ax[i].set_ylim(24.0,18.5)
    #ax[i].set_xlim(t1, t2)
    #ax[i].legend()
    if i == 0:
        ax[i].set_ylabel("B mag")
    elif i == 1:
        ax[i].set_ylabel("V mag")
    elif i == 2:
        ax[i].set_ylabel("I mag")

for i in range(len(t)-1):
    ax[i+len(t)].errorbar(tc[i], C[i], C_err[i], fmt='k', label='1st Diff')
    #ax[i+len(t)].set_xlim(t1, t2)
    ax[i+len(t)].set_ylim(-1.3, 2.6)
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
#In case you want to test binning intervals
#t_ints = [279.0, 279.21, 279.35, 279.6, 279.8, 280.0, 281.22, 281.4, 281.8, 283.24, 283.5, 284.0, 284.5, 285.42, 285.73, 286.24, 286.76, 287.5, 288.35, 288.85, 291.5, 292.3, 292.8, 293.2, 293.7, 294.5, 295.0, 295.7, 296.0, 296.53, 296.7, 297.52, 297.8, 298.15, 298.3, 298.54, 298.8, 299.12, 299.3, 299.8, 300.1, 300.4, 302.0, 302.51, 302.8, 303.11, 304.0, 305.45, 305.7, 306.0, 306.5, 307.0, 307.5, 308.0, 308.53, 308.7, 308.9, 310.0, 310.5, 310.7, 310.9, 311.9, 312.07, 312.5]
for i in range(len(t)):
    #ax[i].scatter(t_ints, [10]*len(t_ints), color='b', marker='.')
    ax[i].errorbar(t[i], F[i], yerr=F_err[i], fmt='g+')
    mask = Mlim[i] < limcuts[i]
    ax[i].errorbar(t[i][mask], F[i][mask], yerr=F_err[i][mask], fmt='r+')
    ax[i].set_xlim(t1_early,t2_early)
plt.show()
