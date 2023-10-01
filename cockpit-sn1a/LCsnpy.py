#################################################################
# Name:     LCsnpy.py                                           #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Program generates light curve in snpy format using  #
#           input light curve data. Light curve data produced   #
#           by LCgen.py routine in associated format.           #
#           Update light curves and ObjData.py first.           #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt

#essential files  
from SNAP.Analysis.LCRoutines import*
from ObjData import *

#make output file
out = open(sn_file, 'w')
#write header
out.write(name+" "+str(z)+" "+str(RA)+" "+str(DEC)+"\n")

#get light curve for SN>SNthres=3.0
t, M, M_err, F, SN, Mlim = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=3.0, scols=9, flags=['-99.99999'], mode='multi')

#crop window in data
for i in range(len(M)):
    t[i], M[i], M_err[i], F[i], SN[i], Mlim[i] = LCcrop(t[i], t1, t2, M[i], M_err[i], F[i], SN[i], Mlim[i])
for i in range(len(M)):
    Mlim[i], t[i], M[i], M_err[i], F[i], SN[i] = LCcrop(Mlim[i], limcuts[i], 100, t[i], M[i], M_err[i], F[i], SN[i])

#for each band, write light curve
for i in range(len(t)):
    #write header for band data
    out.write("filter "+band[i]+"\n")
    print "Writing "+band[i]
    #for each line in light curve, write to file
    for j in range(len(t[i])):
        out.write(padstl(str(t[i][j]),11)+"\t"+padstl(str(M[i][j]),8)+"\t"+padstl(str(M_err[i][j]),7)+"\n")
print "DONE"



