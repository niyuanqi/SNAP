#################################################################
# Name:     ObjData.py                                          #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Important configuration data.                       #
#################################################################

#essential modules
import numpy as np

#window of interest
t1 = 278
t2 = 375
t1_early = 278
t2_early = 313
limcuts = [22,22,22]

#supernova basic data
name = "KSP-N59-1_2016bu"
RA = 3.2336166
DEC = -20.73819534

#light curve files
Bfile = "KSP-N59-1_2016bu.B.lc.CN_180712.txt"
Vfile = "KSP-N59-1_2016bu.V.lc.CN_180712.txt"
Ifile = "KSP-N59-1_2016bu.I.lc.CN_180712.txt"
files = [Bfile, Vfile, Ifile]
#observed bands
band = ['B','V','i']
#band labels
Band = ['B','V','I']

#parameters for template fitting
#SNooPy format light curve
sn_file = "KSP-N59-1_2016bu.txt"
#Phillips relation dataset (including sBV) from Burns
ph_file = "Phillips_dm15+st_results.dat"
#redshift
z = 0.115
zerr = 0.020
#Extinction coefficient (galactic) in each band S & F (2011)
#00h12m56.07s -20d44m17.50s Equ J2000
EBVgal = 0.0172
#CTIO B, CTIO V, SDSS i
Coefs = np.array([3.641, 2.682, 1.698])

#parameters from template fitting (for early light curve fitting)
#maximum epoch
Tmax = 315.857
Tmaxerr = 0.062
#max absolute magnitudes
maxes = np.array([-19.5610, -19.5104, -18.7906])

#early binned time series data files
binBfile = "KSP-N59-1_2016bu.B.lcbin.CN_180712.txt"
binVfile = "KSP-N59-1_2016bu.V.lcbin.CN_180712.txt"
binIfile = "KSP-N59-1_2016bu.I.lcbin.CN_180712.txt"
binfiles = [binBfile, binVfile, binIfile]

#parameters from early light curve fitting (for Arnett fitting)
t0 = -16.71
t0err = 0.57
t0obs = -18.63
t0obserr = 0.72

#parameters from Arnett fitting (for Kasen modelling)
#explosion parameters
m_ej = 1.25 #solar mass
e_51 = 0.90 #x10^51 ergs
m_ej_err = 0.15
e_51_err = 0.20
