#################################################################
# Name:     LCplot.py                                           #
# Author:   Yuan Qi Ni                                          #
# Date:     May, 29, 2016                                       #
# Function: Program generates light curve plot using input      #
#           light curve data. Light curve data produced by      #
#           LCgen.py routine.                                   #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt

#essential files
from SNAP.Analysis.LCRoutines import*

#time series data files
Bfile = "N300-1.Q0.B.005703D193-370223D6.150625-160111.var.lc.CN_170505.txt"
Vfile = "N300-1.Q0.V.005703D193-370223D6.150625-160111.var.lc.CN_170505.txt"
Ifile = "N300-1.Q0.I.005703D193-370223D6.150625-160111.var.lc.CN_170505.txt"
files = [Bfile, Vfile, Ifile]

#get N300-1.Q0.SN light curve
t, M, M_err, F, SN, Mlim = LCload(files, tcol=0, magcols=10, errcols=11, fluxcols=8, SNcols=9, limcols=12, SNthres=2.0, scols=13, flags=['-99.99999'], mode='multi')
print [min(m) for m in M]
print [l.mean() for l in Mlim]

#'-99.99999'
#'-99.999'

#get N300-1.Q0.SN color curve
tc, C, C_err = LCcolors(t, M, M_err)

#plot
f, ax = plt.subplots(len(t)+len(t)-1, sharex=True) 
ax[-1].set_xlabel("Days into 2016")

for i in range(len(t)):
    #ax[i].errorbar(t1[i], M1[i], yerr=M1_err[i], fmt='r+', label='SrcExt')
    ax[i].errorbar(t[i], M[i], yerr=M_err[i], fmt="k", label='MagCalc')
    #ax[i].plot(t[i], M[i])
    #ax[i].scatter(t[i],Mlim[i],color='r',marker='v', label='det. limit')
    #ax[i].errorbar(tn[i],Mn[i],fmt='kx', label='no moon')
    #ax[i].errorbar(tm[i],Mm[i],fmt='bx', label='moon bright')
    ax[i].set_ylim(23.5,17.5)
#    ax[i].set_xlim(260,370)
    #ax[i].set_xlim(170, 370)
    #ax[i].legend()
    if i == 0:
        ax[i].set_ylabel("B mag")
    elif i == 1:
        ax[i].set_ylabel("V mag")
    elif i == 2:
        ax[i].set_ylabel("I mag")

for i in range(len(t)-1):
    ax[i+len(t)].errorbar(tc[i], C[i], C_err[i], fmt='k', label='1st Diff')
    #ax[i+len(t)].set_xlim(260,370)
    if i == 0:
        ax[i+len(t)].set_ylabel("B-V mag")
    elif i == 1:
        ax[i+len(t)].set_ylabel("V-I mag")

f.subplots_adjust(hspace=0)
#plt.savefig("N300-1.Q2.Nova.005509D412-374216D5.150625-151103.var.moon.CN_final.pdf")
plt.show()
