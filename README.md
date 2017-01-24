# SNAP
SuperNova Analysis Package

Performs data analysis and processing on images containing transient sources.
Implements python routines for performing photometry, astrometry, etc.

# Processing
Contains routines for extracting light curves from raw images.
The flagship program is MagCalc.py

MagCalc.py:
Basic usage

python MagCalc.py -c phot -o N300-1.Q0.SN -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 N300-1.Q0.B.151010_1604.A.033278.005604N3646.0060.nh.crop.fits N300_1_Q0_SN.csv

# Analysis
Contains routines for analysing light curves
