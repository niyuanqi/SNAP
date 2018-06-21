#################################################################
# Name:     DiffFits.py                                         #
# Author:   Yuan Qi Ni                                          #
# Date:     July 14, 2017                                       #
# Function: Program uses DiffIm routine to subtract images.     #
#           Update /raw files and ObjData.py before running.    #
#################################################################

#essential modules
import numpy as np
import os
from glob import glob
from multiprocessing import Pool

#essential imports
from SNAP.DiffIm import make_diff_image
from ContextManager import cd
from ObjData import *

#number of processors to use
nproc = 7

#reference files
band = ['B','V','I']
refs = ['../ref/'+Brefname, '../ref/'+Vrefname, '../ref/'+Irefname] 

#current working directory
wd = os.getcwd()
#make directory for diff images
with cd(wd+"/../"):
    if not os.path.isdir("diff"): os.mkdir('diff')
    if not os.path.isdir("conv"): os.mkdir('conv')

#initialize multiprocess
pool = Pool(nproc)
queue = []
#for each band
for i in range(len(band)):
    #get all band files
    files = sorted(glob('../raw/'+prefix+band[i]+'*.fits'))
    for n, filename in enumerate(files):  
        #output filename
        diffname = '.'.join(filename.split('.')[:-1])+".diff.fits"
        diffname = '../diff/'+'/'.join(diffname.split('/')[2:])
        convname = '.'.join(filename.split('.')[:-1])+".conv.fits"
        convname = '../conv/'+'/'.join(convname.split('/')[2:])

        #subtract if not already subtracted
        if os.path.exists(diffname) and os.path.exists(convname):
            print "Already subtracted "+filename
        else:
            print "subtracting "+filename
            #subtract reference image
            args = [filename, refs[i], diffname, convname, "DITemp"+str(n)]
            queue.append(pool.apply_async(make_diff_image, args))

#retrieve results
results = [proc.get() for proc in queue]


