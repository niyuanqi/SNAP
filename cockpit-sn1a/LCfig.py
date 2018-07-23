#################################################################
# Name:     LCfig.py                                            #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Program generates light curve plot using input      #
#           light curve data. Light curve data produced by      #
#           LCgen.py routine. Update ObjData.py first.          #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
from scipy import interp
from scipy.optimize import curve_fit

#essential files
from SNAP.Analysis.LCRoutines import*
from ObjData import *

#get light curve for SN>SNthres=2.0
t, M, M_err, F, SN, Mlim = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=2.0, scols=9, flags=['-99.99999', "INCONV"], mode='multi')

for i in range(len(t)):
    t[i], M[i], M_err[i], F[i], Mlim[i] = LCcrop(t[i], t1, t2, M[i], M_err[i], F[i], Mlim[i])
    t[i] = t[i] - Tmax

for i in range(len(t)):
    m = np.argmin(M[i])
    print M[i][m], M_err[i][m]

#get KSP-OT-1 color curve
tc, C, C_err = LCcolors(t, M, M_err)

#plot
f, ax = plt.subplots(len(t)+len(t)-1, sharex=True)
ax[-1].set_xlabel("Days from peak", fontsize = 14)

tw1=t1 - Tmax + 15
tw2=t2 - Tmax 

for i in range(len(t)):
    ax[i].errorbar(t[i], M[i], yerr=M_err[i], fmt="k+", label='MagCalc')
    ax[i].plot(t[i], M[i], c='k')
    #ax[i].plot([t0,t0], [0,30], 'k:', linewidth=1)
    #ax[i].plot([t0+t0err,t0+t0err], [0,30], 'k:')
    #ax[i].plot([t0-t0err,t0-t0err], [0,30], 'k:')
    ax[i].plot([0,0], [0,30], 'k--', linewidth=1)
    #ax[i].plot([0+Tmaxerr,0+Tmaxerr], [0,30], 'k:')
    #ax[i].plot([0-Tmaxerr,0-Tmaxerr], [0,30], 'k:')
    #ax[i].scatter(t[i],Mlim[i],color='r',marker='v', label='det. limit')
    ax[i].set_ylim(24.0, 18.5)
    ax[i].yaxis.set_ticks([19,21,23])
    ax[i].tick_params(labelsize=12)
    ax[i].set_xlim(tw1,tw2)
    #ax[i].legend()
    if i == 0:
        ax[i].set_ylabel("B", fontsize = 14, fontstyle='italic', fontweight='bold')
        ax[i].text(tw1+0.5, 20, "(a)", fontsize = 14)
    elif i == 1:
        ax[i].set_ylabel("V", fontsize = 14, fontstyle='italic', fontweight='bold')
        ax[i].text(tw1+0.5, 20, "(b)", fontsize = 14)
    elif i == 2:
        ax[i].set_ylabel("I", fontsize = 14, fontstyle='italic', fontweight='bold')
        ax[i].text(tw1+0.5, 20, "(c)", fontsize = 14)

for i in range(len(t)-1):
    ax[i+len(t)].errorbar(tc[i], C[i], C_err[i], fmt='k', label='1st Diff')
    ax[i+len(t)].plot(tc[i], C[i], c='k')
    #ax[i+len(t)].plot([t0,t0], [-5,5], 'k:', linewidth=1)
    #ax[i+len(t)].plot([t0+t0err,t0+t0err], [-5,5], 'k:')
    #ax[i+len(t)].plot([t0-t0err,t0-t0err], [-5,5], 'k:')
    ax[i+len(t)].plot([0,0], [-5,5], 'k--', linewidth=1)
    #ax[i+len(t)].plot([0+Tmaxerr,0+Tmaxerr], [-5,5], 'k:')
    #ax[i+len(t)].plot([0-Tmaxerr,0-Tmaxerr], [-5,5], 'k:')
    ax[i+len(t)].set_ylim(-1.6, 2.6)
    ax[i+len(t)].yaxis.set_ticks([-1.0,0.0,1.0])
    ax[i+len(t)].tick_params(labelsize=12)
    ax[i+len(t)].set_xlim(tw1,tw2)
    if i == 0:
        ax[i+len(t)].set_ylabel("B$-$V", fontsize = 14, fontstyle='italic', fontweight='bold')
        ax[i+len(t)].text(tw1+0.5, 1.5, "(d)", fontsize = 14)
    elif i == 1:
        ax[i+len(t)].set_ylabel("V$-$I", fontsize = 14, fontstyle='italic', fontweight='bold')
        ax[i+len(t)].text(tw1+0.5, 1.5, "(e)", fontsize = 14)

f.subplots_adjust(hspace=0)
plt.show()
