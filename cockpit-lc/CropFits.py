#################################################################
# Name:     CropFits.py                                         #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr., 26, 2018                                      #
# Function: Program crops raw fits files into given size.       #
#           Update /raw files and ObjData.py before running.    #
#################################################################

#essential modules
import os
from glob import glob

#essential imports
from SNAP.CropIm import make_crop_image
from ObjData import *
from ContextManager import cd

#current working directory
wd = os.getcwd()

#make directory for cropped images
with cd(wd+"/../"):
    if not os.path.isdir("crop"): os.mkdir('crop')

#raw image files
filenames = sorted(glob('../raw/'+prefix+'*.fits'))

#crop files write path
outpath = '../crop/'

#crop each file if not already cropped
for filename in filenames:
    #input file
    imname = filename.split("/")[-1]
    #output file
    outname = outpath+imname[:-4]+'crop.fits'

    if os.path.exists(outname):
        print outname+" already exists"
    else:
        make_crop_image(filename, outname, size/2)
