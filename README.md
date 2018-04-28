# SNAP
SuperNova Analysis Package

Performs data analysis and processing on images containing transient sources.
Implements python routines for performing photometry, astrometry, etc. Also can analyse supernova type 1a light curves.

The flagship program is MagCalc.py.

cockpit-lc contains a compact set of routines which one may use to generate light curves from MagCalc. Keep the format of ObjData.py and replace values therein with your own. Then, run routines as outlined in README to generate a quick light curve. DataSetup.py synchronizes files from remote server, and generates file structure needed for cockpit-lc to work. CropFits.py crops raw files. LCgen.py generates light curve from cropped files using MagCalc.py. Can update light curve dynamically (picks up analysis where you left off, or when new data is available).

The subpackage Analysis contains tools for dealing with such light curves.

cockpit-sn1a uses the Analysis subpackage to analyse SN1a data. SNphillip.py can use SNooPy to fit for Phillips parameters (for fit for redshift using SNMCcalib.py), EarlyFit.py can fit power law to early light curve, ArnettFit.py can perform Arnett modelling, KasenCompare.py can compare early light curve to Kasen interaction models. The subdirectory kasen-mc also contains routines for extensive early light curve analysis using Kasen models.

To use any of these modules, add SNAP to your PYTHONPATH

*% PYTHONPATH=$PYTHONPATH:\<path containing SNAP directory\>*

Requires latest astropy to function properly.

Requires dill to use certain monte carlo intensive routines.

Requires SNooPy (Burns 2011) to analyse type 1a supernovae.


**MagCalc.py :**

Automatically performs differential photometry on given fits image files. It uses functions for automatic PSF fitting, planar background fitting, Kron aperture selection, photometric calibration using reference stars, and monte carlo limiting magnitude calculation. Agrees very well with Bertin and Arnout's SExtractor routine on uncrowded field point source photometry, but is superior at crowded field point source photometry. Operating conditions can be manipulated by a highly customizable set of flags. Can use provided reference star catalogs or can automatically query AAVSO (example given below). Can use a science image to perform reference star photometry while performing source photometry on a difference image with the same wcs and gain preferably constructed using DiffIm.py (example given below). Can perform either PSF photometry, automatic aperture photometry, or fixed aperture photometry. Can select how many degrees of freedom with which to fit source PSF. Can fit for background sky.

Basic usage in command line (some samples)

*% python -m SNAP.MagCalc -c phot -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 -f 16.0 --fit_sky conv_image.fits catalog.csv -d diff_image.fits*

*% python -m SNAP.MagCalc -c phot -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -a 0 -fwhm 5 -vvv -n 3.0 -s 14.0 -f 16.0 --fit_sky example_image.fits catalog.csv*

*% python -m SNAP.MagCalc -c aavso -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -a 10 -fwhm 5 -vvv -n 3.0 -s 14.0 -f 16.0 --fit_sky example_image.fits -d example_diffIm.fits my_aavso_catalog_name.cat*

*% python -m SNAP.MagCalc -c phot -o N300-1.Q0.SN -b 'B' -p 14.263303:-37.039900 -r 1000 -psf 1 -fwhm 5 -vvv -n 3.0 -s 14.0 -f 16.0 --fit_sky N300-1.Q0.B.151010_1604.A.033278.005604N3646.0060.nh.crop.fits N300_1_Q0_SN.csv*

*% python -m SNAP.MagCalc -c diff -o KSP-N300-Nova -b 'B' -p 13.789218:-37.704572 -r 1000 -psf 1 -fwhm 5 -n 3.0 -s 14.0 -f 16.0 --fit_sky -vv N300-1.Q2.B.151009_0015.S.000859.005606N3754.0060.nh.fits N300-1.Q2.diff.cat*

*% python -m SNAP.MagCalc -c dprs -o KSP-OT-1 -b 'B' -p 140.92247:-21.969278 -r 1000 -fwhm -psf 2 -fwhm 5 -n 3.0 -s 14.0 -f 16.0 --fit_sky -vv N2784-7.Q1.B.150402_2125.S.015081.092205N2208.0060.nh.fits N2784-7.Q1.DPRS.cat*

Try in terminal

*% python -m SNAP.MagCalc -h*
for explanation of flags and inputs

Basic usage in python routine (sample)

*% from SNAP.MagCalc import magnitude*

*% RAo, DECo, Io, SNo, Mo, Mo_err, Mlim = magnitude(image_fromfits, wcs_fromfits, 'dprs', '../N2784-7.Q1.DPRS.cat', (140.92247,-21.969278), radius=1000.0, name='KSP-OT-1', band='B', fwhm=5.0, limsnr=3.0, satmag=14.0, refMag=16.0, verbosity=0)*

Try in python shell on any imported modules from SNAP

*% from SNAP.MagCalc import \<module\>*

*% help(\<module\>)*
for explanation of functions and inputs

**DiffIm.py :**

Uses WCSremap and HOTPANTS routines (Andrew Becker) to subtract fits files and create image difference files. WCSremap matches images astrometrically, while HOTPANTS matches images photometrically (using convolution) for subtraction. Outputs a difference image, and a convolved image which is the science image photometrically matched to the difference image. When performing photometry, use convolved image for reference stars measurements and difference image for source measurements.

Basic usage (sample)

*% python -m SNAP.DiffIm srcfile_name reffile_name difffile_name convfile_name*

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

**StampIm.py :**

Contains functions for creating stamp png images and stamp collages from fits files.

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

Contains python functions for parsing various differential photometric reference star catalogs. Can automatically query aavso for a catalog, given catalog name will be name of saved file. Can also use a custom catalog given in "diff" format with columns ID,RA,DEC,B,Berr,V,Verr,i,ierr with 'NA' string denoting missing values.

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