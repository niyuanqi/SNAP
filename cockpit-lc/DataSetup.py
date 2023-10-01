#################################################################
# Name:     DataSetup.py                                        #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Set up files for MagCalc data analysis.             #
#           Update ObjData.py beforehand.                       #
#           Synchronizes .fz files from remote.                 #
#           Makes .reg file for easy ds9.                       #
#################################################################

#essential modules
import os
from glob import glob
import subprocess

#essential imports
from ContextManager import cd
from MakeReg import makeReg
from ObjData import *

#current working directory
wd = os.getcwd()

#make directories
with cd(wd+"/../"):
    if not os.path.isdir("raw"): os.mkdir('raw')
    if not os.path.isdir("ref"): os.mkdir('ref')

#make reg file
makeReg(name+".reg", ra, dec)

#synchronize reference files from remote
os.system("rsync -tv "+rawfiles+"REF_Images/*.fits ../ref/")

#synchronize raw files from remote
os.system("rsync -tv "+rawfiles+"*.fz ../raw/")

#unpack files
with cd(wd+"/../raw/"):
    filenames = sorted(glob(prefix+'*.fits.fz'))
    #if not already unpacked, unpack
    for filename in filenames:
        outname = filename[:-3]
        if os.path.exists(outname):
            print "Already unpacked "+outname
        else:
            subprocess.call(['funpack',filename])
