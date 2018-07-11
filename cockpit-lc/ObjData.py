#################################################################
# Name:     ObjData.py                                          #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Important configuration data.                       #
#################################################################


#Always fill the following
t_now = "180426_1100"
user = "Chris Ni"
suffix = "lc.CN_180426.txt"
binsuffix = "lcbin.CN_180426.txt"

#Fill out the following before setup
#year of observation
year = 2018
#name of object
name = "KSP-N3923-2-2018ku"
#files used
rawfiles = "chrisni@sn1987a.astro.utoronto.ca:/data/ksp/data/PROCESSED/N3923/N3923-2/Q1/"
prefix = "N3923-2.Q1."

#Fill out the following before cropping
#position of source
ra = 177.757506
dec = -28.744022
#size of cropped image
size = 2000.0

#Fill out the following before making light curves
#psf to fit
psftype = 2
fitsky = 1
#catalog to use
cattype = 'aavso'
catname = prefix+'.AAVSO.cat'
#signal to noise of detection limits
SNRnoise = 2.0
#limits for reliable reference star magnitudes
satlvl = 14.0
rellvl = 16.0
#saturation pixel count
satpix = 40000.0
#number of reference stars used in each band
nrefs = [11,13,19]

#output light curve filenames
#180421 normal photometry
#180423 saturated star photometry with satpix=25000
#180425 saturated star photometry with satpix=40000, relax centroid
#180426 same method. Update with new data.

#output binned light curve filenames
#180423 saturated star photometry with satpix=25000
#180425 saturated star photometry with satpix=40000, relax centroid for SNR>2
#180426 don't bin saturated images. Use aperture photometry for very early.

"""
#Uncomment this section if you want to perform adv. subtraction
#reference files (keep masks in /ref/)
Brefname = 'N2292-1.Q1.B.161023_1631-161025_0722.XCXA.064318N2550.00005.00005.FM48.BS0512.coadd.REF.fits'
Vrefname = 'N2292-1.Q1.V.161024_1637-161025_1638.XCXA.064318N2550.00005.00005.FM42.BS0512.coadd.REF.fits'
Irefname = 'N2292-1.Q1.I.161023_1745-161025_0835.XCXA.064317N2550.00005.00005.FM37.BS0512.coadd.REF.fits'
#reference images [saturation levels, lower limits]
reflims = [[40000, -36.1005803567],
           [40000, -52.6813824473],
           [40000, -106.147897866]]
"""
