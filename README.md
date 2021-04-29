# SuperNova Analysis Package (SNAP)

Contains routines for image manipulation (cropping, stacking, subtracting), photometry (PSF photometry, aperture photometry, estimating detection limits), and calibration (querying reference stars, differential photometry).
The main photometry program is MagCalc.py.
There are also many supernova-specific programs for analysing their light curves and fitting them with models as well as routines for transforming light curves between the KMTNet BVI filter system and standard Johnson/SDSS filters in the Analysis folder.

The latest manual is here -> ([KSPPhotManual_210429.pdf](KSPPhotManual_210429.pdf)).

It explains the standard workflow for how to use SNAP to process KSP images and extract light curves from them. It also provides detailed descriptions of how individual routines work, tests them on KMTNet data to establish behavioral benchmarks, and instructions on how you can setup or modify these routines/workflows for your purposes. Unfortunately, you need to download the pdf for the embedded hyperlinks to work.

You can also read below for some usage examples with a few programs.

---

Requires installing latest astropy for most things.

Requires installing synphot for synthetic photometry.

Requires installing dill to use certain monte carlo intensive routines.

Requires installing SNooPy (Burns 2011) to analyse type 1a supernovae.

To use any SNAP modules from command line, add SNAP to your PYTHONPATH as follows.

*% PYTHONPATH=$PYTHONPATH:\<path containing SNAP directory\>*

---

**MagCalc.py :**

Automatically performs differential photometry on given fits image files. It executes functions for automatic PSF fitting, planar background fitting, Kron aperture selection, photometric calibration using reference stars, and monte carlo limiting magnitude calculation. Agrees very well with Bertin and Arnout's SExtractor routine on uncrowded field point source photometry, but is superior at crowded field point source photometry. Operating conditions can be manipulated by a very customizable set of flags. Can use provided reference star catalogs or can automatically query AAVSO (example given below). Can use a science image to perform reference star photometry while performing source photometry on a difference image with the same wcs and gain, preferably constructed using DiffIm.py (example given below). Can perform either PSF photometry, automatic aperture photometry, or fixed aperture photometry. Can select how many degrees of freedom with which to fit source PSF. Can fit for background sky. Can handle multiple PSF fitting, when input name, psf, RAo, DECo are lists (aperture must be None, multiple aperture photometry is still under construction).

As of July 2018, MagCalc is able to perform multi-object PSF photometry.
The command line application of MagCalc has been preserved, and any old usage of MagCalc has been preserved (MagCalc will revert to single object photometry).
To use the new multi-object photometry, one need only replace the input parameters source name, ra, dec, psf, fitsky with python lists (not applicable to command line).
Each item in the list corresponds to each object's name, ra, dec, and psf to be used for source. Each object's fitsky parameter (boolean) indicates whether its annulus with participate in sky background fitting. Example: to have an annulus taken around every source, simply give a list of ones.
At sufficient verbosity, MagCalc will provide new plots. These are image plots of residuals of multi fit, and residuals of sky fit (including plot of annulus used). This is in addition to the old plots of fit cross-section (for each object).

Try the following line in terminal for an explanation of flags and inputs.

*% python -m SNAP.MagCalc -h*

Some examples of command line usage:

*% python -m SNAP.MagCalc -c phot -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 15.0 -f 16.0 --fit_sky conv_image.fits catalog.csv -d diff_image.fits*

*% python -m SNAP.MagCalc -c phot -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -a 0 -fwhm 5 -vvv -n 3.0 -s 15.0 -f 16.0 --fit_sky example_image.fits catalog.csv*

*% python -m SNAP.MagCalc -c aavso -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -a 10 -fwhm 5 -vvv -n 3.0 -s 15.0 -f 16.0 --fit_sky example_image.fits -d example_diffIm.fits my_aavso_catalog_name.cat*

*% python -m SNAP.MagCalc -c phot -o N300-1.Q0.SN -b 'B' -p 14.263303:-37.039900 -r 1000 -psf 1 -fwhm 5 -vvv -n 3.0 -s 15.0 -f 16.0 --fit_sky N300-1.Q0.B.151010_1604.A.033278.005604N3646.0060.nh.crop.fits N300_1_Q0_SN.csv*

*% python -m SNAP.MagCalc -c diff -o KSP-N300-Nova -b 'B' -p 13.789218:-37.704572 -r 1000 -psf 1 -fwhm 5 -n 3.0 -s 15.0 -f 16.0 --fit_sky -vv N300-1.Q2.B.151009_0015.S.000859.005606N3754.0060.nh.fits N300-1.Q2.diff.cat*

*% python -m SNAP.MagCalc -c dprs -o KSP-OT-1 -b 'B' -p 140.92247:-21.969278 -r 1000 -fwhm -psf 2 -fwhm 5 -n 3.0 -s 15.0 -f 16.0 --fit_sky -vv N2784-7.Q1.B.150402_2125.S.015081.092205N2208.0060.nh.fits N2784-7.Q1.DPRS.cat*

---

In python shell, you can try the following on any imported modules from SNAP for an explanation of functions and inputs.

*% from SNAP.MagCalc import \<module\>*

*% help(\<module\>)*

Some examples of MagCalc usage in python:

*% from SNAP.MagCalc import magnitude*

*% RAo, DECo, Io, SNo, Mo, Mo_err, Mlim = magnitude(image_fromfits, image_fromfits, wcs_fromfits, 'aavso', '../N2784-7.Q1.AAVSO.cat', (140.92247,-21.969278), radius=1000.0, aperture=15, name='KSP-OT-1', band='B', fwhm=5.0, limsnr=3.0, satmag=15.0, refMag=16.0, verbosity=0)*

---

**DiffIm.py :**

Uses WCSremap and HOTPANTS routines (Andrew Becker) to subtract fits files and create image difference files. WCSremap matches images astrometrically, while HOTPANTS matches images photometrically (using convolution) for subtraction. Outputs a difference image, and a convolved image which is the science image photometrically matched to the difference image. When performing photometry, use convolved image for reference stars measurements and difference image for source measurements.

Try the following line in terminal for an explanation of flags and inputs:

*% python -m SNAP.DiffIm -h*

Example usage:

*% python -m SNAP.DiffIm srcfile_name reffile_name difffile_name convfile_name*

---

**BinIm.py :**

Uses SWarp routine (Emmanuel Bertin) to create binned files with matched wcs.

Try the following line in terminal for an explanation of flags and inputs.

*% python -m SNAP.BinIm -h*

Example usage:

*% python -m SNAP.BinIm t1 t2 outfile_name*

---

**StampIm.py :**

Contains functions for creating stamp png images and stamp collages from fits files.

---

**AutoSEx.py :**

Uses SExtractor routine (Bertin and Arnout) to detect objects in fits images and create a photometry catalog.

Try the following line in terminal for an explanation of flags and inputs.

*% python -m SNAP.AutoSEx -h*

Example usage:

*% python -m SNAP.AutoSEX srcfile_name*

*% python -m SNAP.AutoSEX srcfile_name --config_name myconfig.sex --xml_out --check_out OBJECTS*

---

**MatchPhot.py :**

Operates on the output of SExtractor (Emmanuel Bertin) which takes in a fits image and compiles catalog file with fluxes of objects detected in the image. MatchPhot uses a list reference stars in the same field of view from AAVSO's APASS-DR9 catalog and extracts a photometric solution mapping SExtractor fluxes to their apparent magnitudes by matching bright reference stars from SExtractor catalog to their counterpart in APASS-DR9 catalog.

Try the following line in terminal for an explanation of flags and inputs.

*% python -m SNAP.MatchPhot -h*

Example usage:

*% python -m SNAP.MatchPhot SExtractor_catname band -ref AAVSO_catname

---

**Vizier.py :**

Contains functions for querying/parsing vizier catalogs and currently supports USNO-B1, AAVSO-APASS-DR9.

---

**Astrometry.py :**

Contains functions for computing astrometric quantities, like angles and lunar position.

---

**Catalog.py :**

Contains functions for parsing various differential photometric reference star catalogs. Can automatically query aavso for a catalog, given catalog name will be name of saved file. Can also use a custom catalog given in "diff" format with columns ID,RA,DEC,B,Berr,V,Verr,i,ierr with 'NA' string denoting missing values.

---

**Photometry.py :**

Contains functions for PSF fitting, extraction, integration, etc.

---

## Analysis
Code for analysing light curves, and miscellaneous tools.

**Cosmology.py :**

Contains functions for computing cosmological distances and other quantities.

**FitsSandbox.py :**

Contains functions for making fits test images.

**LCFitting.py :**

Contains functions for fitting polynomials, supernova templates, etc to light curves of transients.

**LCRoutines.py :**

Contains functions for reading and writing light curve files.

## Examples

Miscellaneous programs that do useful things such as make image histogram (for finding saturation limit), making image mask using Astroscrappy, testing multi-object photometry on single files, etc.

## cockpit-lc

Set of routines which one may use to generate light curves using MagCalc. Keep the format of ObjData.py and replace values therein with your own. Then, run routines as outlined in README to generate a quick light curve. DataSetup.py synchronizes files from remote server, and generates file structure needed for cockpit-lc to work. CropFits.py crops raw files. LCgen.py generates light curve from cropped files using MagCalc.py. Can update light curve dynamically (picks up analysis where you left off, or when new data is available). Copy the whole thing into an empty directory, and everything should work. As of July 2018, you can also perform multi-object psf photometry in cockpit-lc. Simply replace the corresponding parameters in ObjData.py with lists, and it will work out of the box if you use LCmgen.py. In fact, LCmgen.py also works on your old ObjData.py with single object inputs, hence it is supposed to be a replacement for LCgen.py. Advanced image subtraction is also possible in cockpit-lc using Diffgen.py. It uses MagCalc to retrieve PSF from image for calibration, uses Astroscrappy package (Curtis McCully) to create artifact masks, and then uses DiffIm.py to subtract reference image from image. Note, this is quite slow, and you can speed this up substatially by cropping both the science image and the reference image using CropIm.py to the same subsection of sky.

## cockpit-sn1a

Set of code that uses the Analysis subpackage to analyse SN1a data. SNphillip.py can use SNooPy to fit for Phillips parameters (for fit for redshift using SNMCcalib.py), EarlyFit.py can fit power law to early light curve, ArnettFit.py can perform Arnett modelling, KasenCompare.py can compare early light curve to Kasen interaction models. The subdirectory kasen-mc also contains routines for extensive early light curve analysis using Kasen models.