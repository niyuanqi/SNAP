#################################################################
# Name:     ObjData.py                                          #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Important configuration data.                       #
#################################################################

#essential modules
import numpy as np

#window of interest
t1 = 85
t2 = 120
t1_early = 85
t2_early = 94.5

#supernova basic data
name = "KSP-N3923-2_2018ku"
RA = 177.757506
DEC = -28.744022

#Kasen MC modelling output files
logfile = "log.CN_180427.txt"
kasfile = "kasen.CN_180427.txt"

#early binned light curves for kasen monte carlo comparison
Bfile = "KSP-N3923-2-2018ku.B.lcbin.CN_180426.txt"
Vfile = "KSP-N3923-2-2018ku.V.lcbin.CN_180426.txt"
Ifile = "KSP-N3923-2-2018ku.I.lcbin.CN_180426.txt"
files = [Bfile, Vfile, Ifile]
#observed bands
band = ['B','V','i']
#band labels
Band = ['B','V','I']
#light curve files with limiting magnitudes at SN=1.0, 2.0, 3.0, 4.0, 5.0
limSNs = [1.0, 2.0, 3.0, 4.0, 5.0]
#[0.84134474606854293, 0.97724986805182079, 0.9986501019683699, 0.99996832875816688, 0.99999971334842808]
Bconffiles = ['KSP-N3923-2-2018ku.B.lcbin.CN_180426.S1.txt',
         'KSP-N3923-2-2018ku.B.lcbin.CN_180426.S2.txt',
         'KSP-N3923-2-2018ku.B.lcbin.CN_180426.S3.txt',
         'KSP-N3923-2-2018ku.B.lcbin.CN_180426.S4.txt',
         'KSP-N3923-2-2018ku.B.lcbin.CN_180426.S5.txt']
Vconffiles = ['KSP-N3923-2-2018ku.V.lcbin.CN_180426.S1.txt',
         'KSP-N3923-2-2018ku.V.lcbin.CN_180426.S2.txt',
         'KSP-N3923-2-2018ku.V.lcbin.CN_180426.S3.txt',
         'KSP-N3923-2-2018ku.V.lcbin.CN_180426.S4.txt',
         'KSP-N3923-2-2018ku.V.lcbin.CN_180426.S5.txt']
Iconffiles = ['KSP-N3923-2-2018ku.I.lcbin.CN_180426.S1.txt',
         'KSP-N3923-2-2018ku.I.lcbin.CN_180426.S2.txt',
         'KSP-N3923-2-2018ku.I.lcbin.CN_180426.S3.txt',
         'KSP-N3923-2-2018ku.I.lcbin.CN_180426.S4.txt',
         'KSP-N3923-2-2018ku.I.lcbin.CN_180426.S5.txt']
conffiles = np.array([Bconffiles, Vconffiles, Iconffiles])
#which sigs to plot
plot = { 1.0 : 0, 2.0 : 1, 3.0 : 2}
#style for each
style = ['k-', 'k--', 'k-.']
#monte carlo trials to use for each confidence level
#for max 3sig confidence, only ~1000 trials needed for MC.
nmc = 10000
#admittedly, this may not hold for 4sig confidence and higher, but speed...

#parameters for template fitting
#SNooPy format light curve
sn_file = "KSP-N3923-2_2018ku.txt"
#Phillips relation dataset (including sBV) from Burns
ph_file = "Phillips_dm15+st_results.dat"
#redshift
z = 0.005767
zerr = 0.000667
#Extinction coefficient (galactic) in each band S & F (2011)
#177h45m27.02s -28d44m38.48s Equ J2000
EBVgal = 0.0936
#CTIO B, CTIO V, SDSS i
Coefs = np.array([3.641, 2.682, 1.698])

#parameters from template fitting (for early light curve fitting)
#maximum epoch
Tmax = 103.23
Tmaxerr = 0.06
#max absolute magnitudes
maxes = np.array([-19.3285, -19.5142, -18.8032])

#parameters from early light curve fitting
t0 = -16.29
t0err = 0.11
t0obs = -16.38
t0obserr = 0.11

#parameters from Arnett fitting
#explosion parameters
m_ni = 1.04 #solar mass
e_51 = 0.76 #x10^51 ergs
m_ni_err = 0.10
e_51_err = 0.15
