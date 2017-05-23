#################################################################
# Name:     BinIm.py                                            #
# Author:   Yuan Qi Ni                                          #
# Version:  May 23, 2017                                        #
# Function: Program contains routines used for fits image       #
#           addition for improvement of SNR.                    #
#################################################################

#essential modules
import numpy as np

#essential imports
from SNAP.Astrometry import *

#function: bin images between two times in day of year float
def binTimes(t1, t2, out_name, delete_temp=True):
    
    #essential modules
    from glob import glob
    import subprocess
    import os

    #read all files
    files = glob('*.fits')
    #find files between t1 and t2
    times = np.zeros(len(files))
    for i in range(len(files)):
        ksp_time = filename.split('.')[3]
        times[i] = isot_day(ksp_isot(ksp_time))
    mask = np.logical_and(times<t2, times>t1)
    binfiles = files[mask]

    #get base output string
    out_base = out_name[:-4]
    wt_name = out_base+'weight.fits'
    xml_name = out_base+'xml'
    #swarp files between t1 and t2
    subprocess.call(['swarp','-COMBINE_TYPE','SUM','-IMAGEOUT_NAME',out_name,
                     '-WEIGHTOUT_NAME',wt_name,'-XML_NAME',xml_name]+binfiles)

#main function
if __name__ == "__main__":
    
    import argparse
    
    #receive arguments
    parser = argparse.ArgumentParser(description='make binned images.')
    parser.add_argument('t1', type=float, help='start time, day of year float')
    parser.add_argument('t2', type=float, help='end time, day of year float')
    parser.add_argument('out_name', type=str, help='output binned file name')
    args = parser.parse_args()
    
    #create binned image
    binTimes(args.t1, args.t2, args.out_name)
