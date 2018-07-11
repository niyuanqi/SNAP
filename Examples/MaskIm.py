from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from astroscrappy import detect_cosmics

from SNAP.MagCalc import *
from SNAP.PSFlib import *
from ObjData import *

i = 2 #ref file number
#filename = "N2292-1.Q1.B.161023_1631-161025_0722.XCXA.064318N2550.00005.00005.FM48.BS0512.coadd.REF.fits"
#filename = "N2292-1.Q1.V.161024_1637-161025_1638.XCXA.064318N2550.00005.00005.FM42.BS0512.coadd.REF.fits"
filename = "N2292-1.Q1.I.161023_1745-161025_0835.XCXA.064317N2550.00005.00005.FM37.BS0512.coadd.REF.fits"
reflims = [40000, 40000, 40000]

#mask files
maskname='.'.join(filename.split('.')[:-1])+".mask.fits"
cleanname='.'.join(filename.split('.')[:-1])+".clean.fits"

#load image
hdu = fits.open(filename)
image = hdu[0].data
hdr = hdu[0].header
hdu.close()
#load wcs information
wcs = WCS(filename)

#info for MagCalc
fo = '.'.join(filename.split('.')[2:5])
band = fo[0]
#extract data from image
PSF, PSFerr, Med, Noise = magnitude(image, image, wcs, cattype, catname, (ra,dec), radius=size, psf=1, name=name, band=band, fwhm=5.0, limsnr=SNRnoise, satmag=satlvl, refmag=rellvl, fitsky=True, satpix=satpix, verbosity=2, diagnosis=True)

#convert to fwhm
fwhm = np.mean(E2moff_toFWHM(*PSF[:-1]))
#image size
imsize = np.mean(image.shape)
#lower limit:
llim = Med - 10*Noise

print ""
print "Image FWHM =",fwhm
print "Lower limit of valid counts =",llim
print ""

crmask, cleanarr = detect_cosmics(image, readnoise=Noise, satlevel=reflims[i], psffwhm=fwhm, psfsize=2.5*fwhm)

hdu = fits.PrimaryHDU(crmask.astype(uint8))
hdu.header = hdr
hdu.writeto(maskname)

hdu = fits.PrimaryHDU(cleanarr)
hdu.header = hdr
hdu.writeto(cleanname)
