# SNAP
SuperNova Analysis Package

Performs data analysis and processing on images containing transient sources.
Implements python routines for performing photometry, astrometry, etc.

## Processing
Contains routines for extracting light curves from raw images.
The flagship program is MagCalc.py

**MagCalc.py :**

Automatically performs differential photometry on given fits image files. It includes functions for automatic PSF fitting, planar background fitting, Kron aperture selection, photometric calibration using reference stars, and monte carlo limiting magnitude calculation. Agrees very well with Bertin and Arnout's SExtractor routine on uncrowded field point source photometry, but is superior at crowded field point source photometry.

Basic usage (sample)

*% python MagCalc.py -c phot -o SOURCE_NAME -b 'B' -p 14.263303:-37.039900 -r 1000 -fwhm 5 -vvv -n 3.0 -s 14.0 example.fits catalog.csv*

Try 

*% python MagCalc.py -h*
for explanation of flags and inputs

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
