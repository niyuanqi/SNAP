# SNAP
SuperNova Analysis Package

Performs data analysis and processing on images containing transient sources.
Implements python routines for performing photometry, astrometry, etc.

The flagship program is MagCalc.py.

LCgen.py is a sample script using which one may generate a MagCalc light curve.

The subpackage Analysis contains tools for dealing with such light curves.

To use any of these modules, add SNAP to your PYTHONPATH

*% PYTHONPATH=$PYTHONPATH:\<path containing SNAP directory\>*

Also requires latest astropy to function properly.


**MagCalc.py :**

Automatically performs differential photometry on given fits image files. It includes functions for automatic PSF fitting, planar background fitting, Kron aperture selection, photometric calibration using reference stars, and monte carlo limiting magnitude calculation. Agrees very well with Bertin and Arnout's SExtractor routine on uncrowded field point source photometry, but is superior at crowded field point source photometry. PSFs can be manipulated by a highly customizable set of flags. Can use provided reference star catalogs or can automatically query AAVSO (example given below). Can use a science image to perform reference star photometry while performing source photometry on a difference image with the same wcs and gain preferably constructed using DiffIm.py (example given below). 

Basic usage in command line (some samples)

*% python -m SNAP.MagCalc -c phot -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 example_image.fits catalog.csv*

*% python -m SNAP.MagCalc -c aavso -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 -f 16.0 example_image.fits -d example_diffIm.fits my_aavso_catalog_name.cat*

*% python -m SNAP.MagCalc -c phot -o N300-1.Q0.SN -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 N300-1.Q0.B.151010_1604.A.033278.005604N3646.0060.nh.crop.fits N300_1_Q0_SN.csv*

*% python -m SNAP.MagCalc -c diff -o KSP-N300-Nova -b 'B' -p 13.789218:-37.704572 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv N300-1.Q2.B.151009_0015.S.000859.005606N3754.0060.nh.fits N300-1.Q2.diff.cat*

*% python -m SNAP.MagCalc -c dprs -o KSP-OT-1 -b 'B' -p 140.92247:-21.969278 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv N2784-7.Q1.B.150402_2125.S.015081.092205N2208.0060.nh.fits N2784-7.Q1.DPRS.cat*

Try in terminal

*% python -m SNAP.MagCalc -h*
for explanation of flags and inputs

Basic usage in python routine (sample)

*% from SNAP.MagCalc import magnitude*

*% RAo, DECo, Io, SNo, Mo, Mo_err, Mlim = magnitude(image_fromfits, wcs_fromfits, 'dprs', '../N2784-7.Q1.DPRS.cat', (140.92247,-21.969278), radius=1000.0, name='KSP-OT-1', band='B', fwhm=5.0, limsnr=3.0, satmag=14.0, verbosity=0)*

Try in python shell on any imported modules from SNAP

*% from SNAP.MagCalc import \<module\>*

*% help(\<module\>)*
for explanation of functions and inputs

**DiffIm.py :**

Uses WCSremap and HOTPANTS routines (Andrew Becker) to subtract fits files and create image difference files. WCSremap matches images astrometrically, while HOTPANTS matches images photometrically (using convolution) for subtraction.

Basic usage (sample)

*% python -m SNAP.DiffIm srcfile_name reffile_name outfile_name*

Try in terminal

*% python -m SNAP.DiffIm -h*

for explanation of flags and inputs

**BinIm.py :**

Uses SWarp routine (Emmanuel Bertin) to create binned files with matched wcs.

Basic usage (sample)

*% python -m SNAP.BinIm t1 t2 outfile_name*

Try in terminal

*% python -m SNAP.BinIm -h*

for explanation of flags and inputs

**AutoSEx.py :**

Uses SExtractor routine (Bertin and Arnout) to detect objects in fits images and create catalog files.

Basic usage (sample)

*% python -m SNAP.AutoSEX srcfile_name*

*% python -m SNAP.AutoSEX srcfile_name --config_name myconfig.sex --xml_out --check_out OBJECTS*

Try in terminal

*% python -m SNAP.AutoSEx -h*

for explanation of flags and inputs

**MatchPhot.py :**

Operates on the output of SExtractor (Emmanuel Bertin) which takes in a fits image and compiles catalog file with fluxes of objects detected in the image. MatchPhot uses a list reference stars in the same field of view from AAVSO's APASS-DR9 catalog and extracts a photometric solution mapping SExtractor fluxes to their apparent magnitudes by matching bright reference stars from SExtractor catalog to their counterpart in APASS-DR9 catalog.

Basic usage (sample)

*% python -m SNAP.MatchPhot SExtractor_catname band -ref AAVSO_catname

Try in terminal

*% python -m SNAP.MatchPhot -h*

for explanation of flags and inputs

**Vizier.py :**

Contains python functions for querying/parsing vizier catalogs and currently supports USNO-B1, AAVSO-APASS-DR9.

**Astrometry.py :**

Contains python functions for computing astrometric quantities, like angles and lunar position.

**Catalog.py :**

Contains python functions for parsing various differential photometric reference star catalogs.

**Photometry.py :**

Contains python functions for PSF fitting, extraction, integration, etc.

## Analysis
Contains routines for analysing light curves, and miscellaneous tools.

**Cosmology.py :**

Contains python functions for computing cosmological distances and other quantities.

**FitsSandbox.py :**

Contains python functions for making fits test images.

**LCFitting.py :**

Contains python functions for fitting polynomials, supernova templates, etc to light curves of transients.

**LCRoutines.py :**

Contains python functions for reading and writing light curve files.

## Examples

Contains programs that use various SNAP routines to do things such as generate light curves from raw images, generate binned light curves, plot SN1a light curve, calibrate SN1a using light curve, etc.