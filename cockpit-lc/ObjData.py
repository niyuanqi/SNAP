#################################################################
# Name:     ObjData.py                                          #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Important configuration data.                       #
#################################################################


#Always fill the following
t_now = "180426_1100"
user = "Chris Ni"

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
#catalog to use
cattype = 'aavso'
catname = 'N3923-2.Q1.AAVSO.cat'
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
outBname = "KSP-N3923-2-2018ku.B.lc.CN_180426.txt"
outVname = "KSP-N3923-2-2018ku.V.lc.CN_180426.txt"
outIname = "KSP-N3923-2-2018ku.I.lc.CN_180426.txt"

#Fill out the following before making binned light curves
#output binned light curve filenames
#180423 saturated star photometry with satpix=25000
#180425 saturated star photometry with satpix=40000, relax centroid for SNR>2
#180426 don't bin saturated images. Use aperture photometry for very early.
binBname = "KSP-N3923-2-2018ku.B.lcbin.CN_180426.txt"
binVname = "KSP-N3923-2-2018ku.V.lcbin.CN_180426.txt"
binIname = "KSP-N3923-2-2018ku.I.lcbin.CN_180426.txt"

