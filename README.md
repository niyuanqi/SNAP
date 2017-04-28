# SNAP
SuperNova Analysis Package

Performs data analysis and processing on images containing transient sources.
Implements python routines for performing photometry, astrometry, etc.

The flagship program is MagCalc.py.

LCgen.py is a samply script using which one may generate a MagCalc light curve.

The subpackage Analysis contains tools for dealing with such light curves.

**MagCalc.py :**

Automatically performs differential photometry on given fits image files. It includes functions for automatic PSF fitting, planar background fitting, Kron aperture selection, photometric calibration using reference stars, and monte carlo limiting magnitude calculation. Agrees very well with Bertin and Arnout's SExtractor routine on uncrowded field point source photometry, but is superior at crowded field point source photometry. PSFs can be manipulated by a highly customizable set of flags.

Basic usage in command line (some samples)

*% python -m SNAP.MagCalc -c phot -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 example.fits catalog.csv*

*% python -m SNAP.MagCalc -c diff -o KSP-N300-Nova -b 'B' -p 13.789218:-37.704572 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv N300-1.Q2.B.151009_0015.S.000859.005606N3754.0060.nh.fits N300-1.Q2.diff.cat*

*% python -m SNAP.MagCalc -c dprs -o KSP-OT-1 -b 'B' -p 140.92247:-21.969278 -r 1000 -fwhm 5 -n 3.0 -s 14.0 -vv N2784-7.Q1.B.150402_2125.S.015081.092205N2208.0060.nh.fits N2784-7.Q1.DPRS.cat*

Try 

*% python MagCalc.py -h*
for explanation of flags and inputs

Basic usage in python routine (sample)

*% RAo, DECo, Io, SNo, Mo, Mo_err, Mlim = magnitude(image, wcs, 'dprs', catname, (RA,DEC), radius=1000.0, name=name, band=band, fwhm=5.0, limsnr=SNRnoise, satmag=14.0, verbosity=0)*

Try in python shell on imported modules from SNAP

*% help(<module>)*
for explanation of functions and inputs

**DiffIm.py :**

Uses WCSremap and HOTPANTS routines (Andrew Becker) to subtract fits files and create image difference files. WCSremap matches images astrometrically, while HOTPANTS matches images photometrically for subtraction.

Basic usage (sample)

*% python DiffIm.py srcfile_name reffile_name outfile_name*

Try

*% python DiffIm.py -h*

for explanation of flags and inputs

**MatchPhot.py :**

Operates on the output of SExtractor (Emmanuel Bertin) which takes in a fits image and compiles catalog file with fluxes of objects detected in the image. MatchPhot uses a list reference stars in the same field of view from AAVSO's APASS-DR9 catalog and extracts a photometric solution mapping SExtractor fluxes to their apparent magnitudes by matching bright reference stars from SExtractor catalog to their counterpart in APASS-DR9 catalog.

Basic usage (sample)

*% python MatchPhot.py SExtractor_catname band -ref AAVSO_catname

Try

*% python MatchPhot.py -h*

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
