#################################################################
# Name:     ObjData.py                                          #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Important configuration data.                       #
#################################################################

##########################################
# Setup section:
##########################################
t_now = "180426_1100"
user = "Chris Ni"

#year of first detection
year = 2018

#You need to fill out the following before running DataSetup.py
#name of object
name = "KSP-N3923-2-2018ku"
#files used
prefix = "N3923-2.Q1."
rawfiles = "chrisni@sn1987a.astro.utoronto.ca:/data/ksp/data/PROCESSED/" + prefix.split("-")[0] + "/" + prefix.split(".")[0] + "/" + prefix.split(".")[1] + "/"

##########################################
# Image processing section:
##########################################

#You need to fill out the following before cropping
#position of source
ra = 177.757506 #deg
dec = -28.744022 #deg
#radius (pixels) of reference stars (also size of cropped image)
size = 2000.0

"""
#Uncomment this section if you want to perform image subtraction
#reference files
Brefname = 'N2292-1.Q1.B.161023_1631-161025_0722.XCXA.064318N2550.00005.00005.FM48.BS0512.coadd.REF.fits'
Vrefname = 'N2292-1.Q1.V.161024_1637-161025_1638.XCXA.064318N2550.00005.00005.FM42.BS0512.coadd.REF.fits'
Irefname = 'N2292-1.Q1.I.161023_1745-161025_0835.XCXA.064317N2550.00005.00005.FM37.BS0512.coadd.REF.fits'
#reference image fwhm (measure using MagCalc and flag diagnosis=True)
ref_fwhms = [3.325, 3.215, 2.266]
"""

##########################################
# Photometry section:
##########################################

#Fill out the following before making light curves
#write some notes here on what each of your suffixes  means
#e.g.: CN_180426 chris' saturated star photometry with satpix=40000, psf=2
suffix = "lc.CN_180426.txt"
binsuffix = "lcbin.CN_180426.txt"

#catalog to use
cattype = 'aavso'
catname = prefix+'.AAVSO.cat'
#psf to fit
psftype = 2 #default. Also change the sequence in LC*gen.py
fitsky = 1
#signal to noise of detection limits
SNRnoise = 3.0
#limits for reliable reference star magnitudes
satlvl = 15.0
rellvl = 16.0
#saturation pixel count
satpix = 40000.0
#number of reference stars used in each band
nrefs = [11,13,19]
